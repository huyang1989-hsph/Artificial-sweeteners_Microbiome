# Conditional leave-one-participant-out stability of the first TWO sCCA components.
# Requires PMA. Sourcing this file does not fit a model.
# Candidate species, supplied transformations, and penalties remain fixed.
# PMA re-estimates centering/scaling within each training set when standardize=TRUE.
# This is a selection/pattern stability assessment, not held-out validation or
# a significance test. No elastic-net fitting or penalty tuning is performed.

SCCA_LOO_HELPER_VERSION <- "1.0"

.scca_loo_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

.scca_loo_cor <- function(a, b) {
  if (any(!is.finite(a)) || any(!is.finite(b)) ||
      stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
  as.numeric(stats::cor(a, b))
}

.scca_loo_cosine <- function(A, B) {
  a_norm <- sqrt(colSums(A^2))
  b_norm <- sqrt(colSums(B^2))
  .scca_loo_assert(all(is.finite(c(a_norm, b_norm))) &&
                     all(c(a_norm, b_norm) > 0),
                   "A requested component has a zero or invalid weight vector.")
  crossprod(sweep(A, 2L, a_norm, "/"), sweep(B, 2L, b_norm, "/"))
}

.scca_loo_align <- function(ref_u, ref_v, new_u, new_v) {
  # One-to-one matching of ONLY the first two components. The same sign flip
  # is applied to both weight vectors in a canonical pair. No rotations,
  # independent sign flips, or searching later components to improve agreement.
  similarity <- (.scca_loo_cosine(ref_u, new_u) +
                   .scca_loo_cosine(ref_v, new_v)) / 2
  direct <- abs(similarity[1L, 1L]) + abs(similarity[2L, 2L])
  swapped <- abs(similarity[1L, 2L]) + abs(similarity[2L, 1L])
  index <- if (swapped > direct) c(2L, 1L) else c(1L, 2L)
  matched <- similarity[cbind(1:2, index)]
  signs <- ifelse(matched < 0, -1, 1)
  list(u = sweep(new_u[, index, drop = FALSE], 2L, signs, "*"),
       v = sweep(new_v[, index, drop = FALSE], 2L, signs, "*"),
       index = index, signs = signs, matching_similarity = abs(matched),
       assignment_margin = abs(direct - swapped) / 2)
}

.scca_loo_measure <- function(X, Y, u, v, standardize) {
  if (standardize) {
    X_fit <- scale(X)
    Y_fit <- scale(Y)
  } else {
    X_fit <- X
    Y_fit <- Y
  }
  sx <- X_fit %*% u
  sy <- Y_fit %*% v
  .scca_loo_assert(all(is.finite(c(sx, sy))) &&
                     all(apply(sx, 2L, stats::sd) > 0) &&
                     all(apply(sy, 2L, stats::sd) > 0),
                   "A requested component has an undefined or constant score.")
  list(score_x = sx, score_y = sy,
       canonical_correlation = diag(stats::cor(sx, sy)),
       species_structure_loadings = stats::cor(X, sx),
       sweetener_profile = stats::cor(Y, sx))
}

.scca_loo_feature_names <- function(M, expected, label) {
  if (is.null(expected)) return(M)
  expected <- as.character(expected)
  .scca_loo_assert(length(expected) == ncol(M) && !anyNA(expected) &&
                     !anyDuplicated(expected) && setequal(colnames(M), expected),
                   paste0(label, " columns differ from the saved model's feature names."))
  M[, expected, drop = FALSE]
}

.scca_loo_quantiles <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, 5L))
  c(min(x), stats::quantile(x, c(.25, .5, .75), names = FALSE), max(x))
}

scca_loo_stability <- function(
    X, Y, reference_fit,
    component_names = c("Component_1", "Component_2"),
    analysis = "analysis", standardize = TRUE,
    output_dir = NULL, check_only = FALSE, omit_ids = NULL,
    progress_every = 25L, tolerance = 1e-6) {

  .scca_loo_assert(requireNamespace("PMA", quietly = TRUE),
                   "Install/load the PMA version used for your original sCCA analysis.")
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  .scca_loo_assert(is.numeric(X) && is.numeric(Y) &&
                     all(is.finite(X)) && all(is.finite(Y)),
                   "X and Y must be finite numeric matrices; no samples are silently removed.")
  .scca_loo_assert(nrow(X) >= 5L && nrow(X) == nrow(Y) &&
                     ncol(X) >= 2L && ncol(Y) >= 2L,
                   "Supply paired matrices with at least five participants and two columns each.")
  for (nm in c("X", "Y")) {
    M <- if (nm == "X") X else Y
    .scca_loo_assert(!is.null(colnames(M)) && !anyNA(colnames(M)) &&
                       !anyDuplicated(colnames(M)) && all(nzchar(colnames(M))),
                     paste(nm, "needs unique, nonmissing feature column names."))
    .scca_loo_assert(!is.null(rownames(M)) && !anyNA(rownames(M)) &&
                       !anyDuplicated(rownames(M)) && all(nzchar(rownames(M))),
                     paste(nm, "needs unique participant IDs as row names."))
  }
  .scca_loo_assert(setequal(rownames(X), rownames(Y)),
                   "X and Y must contain the same participant IDs. Prepare the original analytic sample first.")
  Y <- Y[match(rownames(X), rownames(Y)), , drop = FALSE]
  .scca_loo_assert(is.logical(standardize) && length(standardize) == 1L &&
                     !is.na(standardize) && is.logical(check_only) &&
                     length(check_only) == 1L && !is.na(check_only),
                   "standardize and check_only must each be TRUE or FALSE.")
  .scca_loo_assert(is.numeric(tolerance) && length(tolerance) == 1L &&
                     is.finite(tolerance) && tolerance > 0,
                   "tolerance must be a positive finite number.")
  .scca_loo_assert(is.numeric(progress_every) && length(progress_every) == 1L &&
                     is.finite(progress_every) && progress_every >= 0,
                   "progress_every must be a nonnegative number.")
  .scca_loo_assert(length(analysis) == 1L && !is.na(analysis) && nzchar(analysis) &&
                     length(component_names) == 2L && !anyNA(component_names) &&
                     !anyDuplicated(component_names) && all(nzchar(component_names)),
                   "Supply one analysis name and two distinct component names in the ORIGINAL order.")
  .scca_loo_assert(is.list(reference_fit) && is.matrix(reference_fit$u) &&
                     is.matrix(reference_fit$v) && ncol(reference_fit$u) >= 2L &&
                     ncol(reference_fit$v) >= 2L,
                   "reference_fit must be the full original PMA::CCA result, including u and v.")
  .scca_loo_assert(is.null(reference_fit$outcome),
                   "This helper supports the original unsupervised PMA::CCA analysis only.")
  for (nm in c("typex", "typez")) {
    .scca_loo_assert(is.null(reference_fit[[nm]]) ||
                       identical(reference_fit[[nm]], "standard"),
                     "This helper supports typex='standard' and typez='standard' only.")
  }
  for (nm in c("upos", "uneg", "vpos", "vneg")) {
    .scca_loo_assert(is.null(reference_fit[[nm]]) ||
                       identical(reference_fit[[nm]], FALSE),
                     "This helper expects the original unconstrained weight signs.")
  }
  X <- .scca_loo_feature_names(X, reference_fit$xnames, "X")
  Y <- .scca_loo_feature_names(Y, reference_fit$znames, "Y")
  .scca_loo_assert(nrow(reference_fit$u) == ncol(X) &&
                     nrow(reference_fit$v) == ncol(Y),
                   "Matrix columns must match ALL original sCCA candidates, including zero-weight features.")
  .scca_loo_assert(all(apply(X, 2L, stats::sd) > 0) &&
                     all(apply(Y, 2L, stats::sd) > 0),
                   "X or Y contains a constant column.")
  u_ref <- reference_fit$u[, 1:2, drop = FALSE]
  v_ref <- reference_fit$v[, 1:2, drop = FALSE]
  .scca_loo_assert(all(is.finite(c(u_ref, v_ref))),
                   "The original first two component weight vectors must be finite.")
  # Check nonzero norms before any fitting.
  .scca_loo_cosine(u_ref, u_ref)
  .scca_loo_cosine(v_ref, v_ref)
  penaltyx <- reference_fit$penaltyx
  penaltyz <- reference_fit$penaltyz
  .scca_loo_assert(is.numeric(penaltyx) && length(penaltyx) == 1L &&
                     is.finite(penaltyx) && penaltyx >= 0 && penaltyx <= 1 &&
                     is.numeric(penaltyz) && length(penaltyz) == 1L &&
                     is.finite(penaltyz) && penaltyz >= 0 && penaltyz <= 1,
                   "The original model must contain valid penaltyx and penaltyz settings.")
  K <- reference_fit$K
  if (is.null(K)) K <- ncol(reference_fit$u)
  niter <- reference_fit$niter
  if (is.null(niter)) niter <- 15L
  .scca_loo_assert(is.numeric(K) && length(K) == 1L && is.finite(K) &&
                     K >= 2L && K == as.integer(K) &&
                     K <= min(ncol(X), ncol(Y)) &&
                     is.numeric(niter) && length(niter) == 1L &&
                     is.finite(niter) && niter >= 1L && niter == as.integer(niter),
                   "Invalid original K or niter setting.")

  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  warning_rows <- list()
  fit_model <- function(A, B, stage) {
    tryCatch(withCallingHandlers(
      PMA::CCA(x = A, z = B, typex = "standard", typez = "standard",
               penaltyx = penaltyx, penaltyz = penaltyz, K = K,
               niter = niter, standardize = standardize, trace = FALSE),
      warning = function(w) {
        warning_rows[[length(warning_rows) + 1L]] <<-
          data.frame(stage = stage, message = conditionMessage(w))
      }), error = function(e) {
        stop(sprintf("sCCA failed at %s: %s. No fold has been silently discarded.",
                     stage, conditionMessage(e)), call. = FALSE)
      })
  }
  validate_fit <- function(fit, stage) {
    .scca_loo_assert(is.matrix(fit$u) && is.matrix(fit$v) &&
                       nrow(fit$u) == ncol(X) && nrow(fit$v) == ncol(Y) &&
                       ncol(fit$u) >= 2L && ncol(fit$v) >= 2L &&
                       all(is.finite(c(fit$u[, 1:2], fit$v[, 1:2]))),
                     paste("Invalid component weights at", stage))
  }
  reference <- .scca_loo_measure(X, Y, u_ref, v_ref, standardize)
  full_fit <- fit_model(X, Y, "full-data reference check")
  validate_fit(full_fit, "full-data reference check")
  aligned_full <- .scca_loo_align(u_ref, v_ref,
                                 full_fit$u[, 1:2, drop = FALSE],
                                 full_fit$v[, 1:2, drop = FALSE])
  stored_r <- reference_fit$cors
  has_stored_r <- !is.null(stored_r) && length(stored_r) >= 2L &&
    all(is.finite(stored_r[1:2]))
  baseline_check <- do.call(rbind, lapply(1:2, function(j) {
    u_diff <- abs(u_ref[, j] - aligned_full$u[, j])
    v_diff <- abs(v_ref[, j] - aligned_full$v[, j])
    set_match <- all((u_ref[, j] != 0) == (aligned_full$u[, j] != 0)) &&
      all((v_ref[, j] != 0) == (aligned_full$v[, j] != 0))
    coef_match <- all(u_diff <= tolerance * pmax(1, abs(u_ref[, j]))) &&
      all(v_diff <= tolerance * pmax(1, abs(v_ref[, j])))
    r_match <- if (has_stored_r) {
      abs(reference$canonical_correlation[j] - stored_r[j]) <= tolerance
    } else TRUE
    data.frame(analysis = analysis, component = component_names[j],
               matched_refit_component = aligned_full$index[j],
               joint_sign = aligned_full$signs[j],
               n_saved_species = sum(u_ref[, j] != 0),
               n_refit_species = sum(aligned_full$u[, j] != 0),
               max_species_weight_difference = max(u_diff),
               max_sweetener_weight_difference = max(v_diff),
               canonical_correlation = reference$canonical_correlation[j],
               stored_correlation_available = has_stored_r,
               stored_correlation_match = r_match,
               selected_set_match = set_match, coefficient_match = coef_match,
               baseline_passed = set_match && coef_match && r_match)
  }))
  rownames(baseline_check) <- NULL
  baseline_passed <- all(baseline_check$baseline_passed)
  if (!is.null(output_dir)) {
    utils::write.csv(baseline_check,
                      file.path(output_dir, "SCCA_full_data_reference_check.csv"), row.names = FALSE)
  }
  message(sprintf("%s: full-data reference check %d/2 components passed.",
                  analysis, sum(baseline_check$baseline_passed)))
  settings <- list(analysis = analysis, n_participants = nrow(X),
                   n_candidate_species = ncol(X), n_sweeteners = ncol(Y),
                   penaltyx = penaltyx, penaltyz = penaltyz, K = K, niter = niter,
                   standardize = standardize, PMA_version = as.character(utils::packageVersion("PMA")),
                   R_version = R.version.string, helper_version = SCCA_LOO_HELPER_VERSION)
  if (check_only) {
    return(list(baseline_check = baseline_check, baseline_passed = baseline_passed,
                reference = reference, settings = settings))
  }
  .scca_loo_assert(baseline_passed,
                   paste("Full-data weights or correlations do not reproduce the saved model.",
                         "Run check_only=TRUE and inspect baseline_check. No LOO fits were started."))

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  ids <- rownames(X)
  if (is.null(omit_ids)) omit_ids <- ids
  omit_ids <- as.character(omit_ids)
  .scca_loo_assert(length(omit_ids) > 0L && !anyNA(omit_ids) &&
                     !anyDuplicated(omit_ids) && all(omit_ids %in% ids),
                   "omit_ids must contain unique participant IDs found in X; omit for the full analysis.")
  b <- length(omit_ids)
  settings$n_refits <- b
  settings$is_complete_leave_one_out <- setequal(omit_ids, ids)
  if (!settings$is_complete_leave_one_out) {
    message(sprintf("Pilot/subset run: %d of %d participants will be omitted in turn.", b, n))
  }
  weights_x <- loadings_x <- array(NA_real_, c(p, 2L, b),
                                   dimnames = list(colnames(X), component_names, omit_ids))
  weights_y <- profiles_y <- array(NA_real_, c(q, 2L, b),
                                   dimnames = list(colnames(Y), component_names, omit_ids))
  fold_rows <- vector("list", 2L * b)

  for (i in seq_len(b)) {
    keep <- ids != omit_ids[i]
    A <- X[keep, , drop = FALSE]
    B <- Y[keep, , drop = FALSE]
    .scca_loo_assert(all(apply(A, 2L, stats::sd) > 0) &&
                       all(apply(B, 2L, stats::sd) > 0),
                     paste("A candidate becomes constant after omitting", omit_ids[i],
                           "- no columns or folds have been silently removed."))
    fit <- fit_model(A, B, paste("omission of", omit_ids[i]))
    validate_fit(fit, paste("omission of", omit_ids[i]))
    aligned <- .scca_loo_align(u_ref, v_ref, fit$u[, 1:2, drop = FALSE],
                               fit$v[, 1:2, drop = FALSE])
    measure <- .scca_loo_measure(A, B, aligned$u, aligned$v, standardize)
    weights_x[, , i] <- aligned$u
    weights_y[, , i] <- aligned$v
    loadings_x[, , i] <- measure$species_structure_loadings
    profiles_y[, , i] <- measure$sweetener_profile
    u_cos <- diag(.scca_loo_cosine(u_ref, aligned$u))
    v_cos <- diag(.scca_loo_cosine(v_ref, aligned$v))
    for (j in 1:2) {
      selected_ref <- u_ref[, j] != 0
      selected_new <- aligned$u[, j] != 0
      n_retained <- sum(selected_ref & selected_new)
      reference_loading <- reference$species_structure_loadings[, j]
      comparable <- reference_loading != 0
      fold_rows[[2L * (i - 1L) + j]] <- data.frame(
        analysis = analysis, omitted_participant = omit_ids[i],
        component = component_names[j], original_component = j,
        matched_refit_component = aligned$index[j],
        component_order_changed = aligned$index[j] != j,
        joint_sign = aligned$signs[j],
        matching_similarity = aligned$matching_similarity[j],
        assignment_margin = aligned$assignment_margin,
        canonical_correlation = measure$canonical_correlation[j],
        species_loading_similarity = .scca_loo_cor(
          reference_loading, measure$species_structure_loadings[, j]),
        sweetener_profile_similarity = .scca_loo_cor(
          reference$sweetener_profile[, j], measure$sweetener_profile[, j]),
        microbial_score_similarity = .scca_loo_cor(
          reference$score_x[keep, j], measure$score_x[, j]),
        sweetener_score_similarity = .scca_loo_cor(
          reference$score_y[keep, j], measure$score_y[, j]),
        species_weight_cosine = u_cos[j], sweetener_weight_cosine = v_cos[j],
        species_loading_direction_concordance = if (any(comparable)) {
          mean(sign(reference_loading[comparable]) ==
                 sign(measure$species_structure_loadings[comparable, j]))
        } else NA_real_,
        n_selected_species = sum(selected_new),
        n_reference_species_retained = n_retained,
        reference_species_retention = n_retained / sum(selected_ref),
        selection_jaccard = n_retained / sum(selected_ref | selected_new))
    }
    if (progress_every > 0 && (i %% progress_every == 0 || i == b)) {
      message(sprintf("%s: completed %d/%d sCCA refits.", analysis, i, b))
    }
  }
  fold_results <- do.call(rbind, fold_rows)
  rownames(fold_results) <- NULL
  metrics <- c("canonical_correlation", "species_loading_similarity",
               "sweetener_profile_similarity", "microbial_score_similarity",
               "sweetener_score_similarity", "species_weight_cosine",
               "sweetener_weight_cosine", "species_loading_direction_concordance",
               "reference_species_retention", "selection_jaccard", "n_selected_species",
               "assignment_margin")
  summary_rows <- list()
  for (j in 1:2) {
    ref_values <- c(canonical_correlation = unname(reference$canonical_correlation[j]),
                    species_loading_similarity = 1, sweetener_profile_similarity = 1,
                    microbial_score_similarity = 1, sweetener_score_similarity = 1,
                    species_weight_cosine = 1, sweetener_weight_cosine = 1,
                    species_loading_direction_concordance = 1,
                    reference_species_retention = 1, selection_jaccard = 1,
                    n_selected_species = sum(u_ref[, j] != 0), assignment_margin = NA_real_)
    d <- fold_results[fold_results$component == component_names[j], , drop = FALSE]
    for (metric in metrics) {
      values <- d[[metric]]
      quantiles <- .scca_loo_quantiles(values)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        analysis = analysis, component = component_names[j], metric = metric,
        n_refits = b, n_valid = sum(is.finite(values)),
        n_undefined = sum(!is.finite(values)), reference_value = unname(ref_values[metric]),
        median = quantiles[3L], q25 = quantiles[2L], q75 = quantiles[4L],
        minimum = quantiles[1L], maximum = quantiles[5L])
    }
  }
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL

  species_rows <- profile_rows <- overview_rows <- vector("list", 2L)
  for (j in 1:2) {
    W <- matrix(weights_x[, j, ], nrow = p, ncol = b)
    L <- matrix(loadings_x[, j, ], nrow = p, ncol = b)
    P <- matrix(profiles_y[, j, ], nrow = q, ncol = b)
    selected <- W != 0
    ref_selected <- u_ref[, j] != 0
    selection_frequency <- rowMeans(selected)
    loading_quantiles <- t(apply(L, 1L, .scca_loo_quantiles))
    loading_sign <- vapply(seq_len(p), function(k) {
      ref <- reference$species_structure_loadings[k, j]
      if (ref == 0) NA_real_ else mean(sign(L[k, ]) == sign(ref))
    }, numeric(1))
    weight_sign <- vapply(seq_len(p), function(k) {
      retained <- selected[k, ]
      if (!ref_selected[k] || !any(retained)) NA_real_ else {
        mean(sign(W[k, retained]) == sign(u_ref[k, j]))
      }
    }, numeric(1))
    species_rows[[j]] <- data.frame(
      analysis = analysis, component = component_names[j], species = colnames(X),
      reference_selected = ref_selected, reference_weight = u_ref[, j],
      reference_structure_loading = reference$species_structure_loadings[, j],
      n_refits = b, n_selected = rowSums(selected), selection_frequency = selection_frequency,
      weight_sign_agreement_when_selected = weight_sign,
      structure_loading_sign_agreement = loading_sign,
      median_structure_loading = loading_quantiles[, 3L],
      q25_structure_loading = loading_quantiles[, 2L],
      q75_structure_loading = loading_quantiles[, 4L])
    profile_quantiles <- t(apply(P, 1L, .scca_loo_quantiles))
    profile_rows[[j]] <- data.frame(
      analysis = analysis, component = component_names[j], sweetener = colnames(Y),
      reference_correlation = reference$sweetener_profile[, j],
      median_correlation = profile_quantiles[, 3L],
      q25_correlation = profile_quantiles[, 2L], q75_correlation = profile_quantiles[, 4L],
      minimum_correlation = profile_quantiles[, 1L], maximum_correlation = profile_quantiles[, 5L])
    d <- fold_results[fold_results$component == component_names[j], , drop = FALSE]
    overview_rows[[j]] <- data.frame(
      analysis = analysis, component = component_names[j], n_participants = n,
      n_refits = b, n_candidate_species = p,
      n_reference_selected_species = sum(ref_selected),
      n_reference_selected_every_fold = sum(ref_selected & selection_frequency == 1),
      proportion_reference_selected_every_fold = mean(selection_frequency[ref_selected] == 1),
      median_reference_retention_per_fold = stats::median(d$reference_species_retention),
      proportion_component_order_changed = mean(d$component_order_changed))
  }
  warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else {
    data.frame(stage = character(0), message = character(0))
  }
  result <- list(summary = summary, component_overview = do.call(rbind, overview_rows),
                 fold_results = fold_results, species_stability = do.call(rbind, species_rows),
                 sweetener_profiles = do.call(rbind, profile_rows), baseline_check = baseline_check,
                 baseline_passed = baseline_passed, warnings = warnings,
                 reference = list(u = u_ref, v = v_ref,
                                  canonical_correlation = reference$canonical_correlation,
                                  species_structure_loadings = reference$species_structure_loadings,
                                  sweetener_profile = reference$sweetener_profile),
                 aligned_weights_x = weights_x, aligned_weights_y = weights_y,
                 species_structure_loadings_by_fold = loadings_x,
                 sweetener_profiles_by_fold = profiles_y, settings = settings,
                 assessment = "Leave-one-out sCCA stability conditional on supplied preprocessing, candidate set, and fixed penalties")
  if (!is.null(output_dir)) {
    for (nm in c("summary", "component_overview", "fold_results", "species_stability",
                 "sweetener_profiles", "warnings")) {
      utils::write.csv(result[[nm]], file.path(output_dir, paste0("SCCA_LOO_", nm, ".csv")),
                        row.names = FALSE)
    }
    saveRDS(result, file.path(output_dir, "SCCA_LOO_results.rds"))
  }
  result
}
