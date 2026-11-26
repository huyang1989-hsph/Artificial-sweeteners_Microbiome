# Leave-one-participant-out ELN feature-selection stability
# Requires the already-installed glmnet package. Sourcing does not fit models.
#
# Scope: alpha, previously tuned lambda, candidate features, and supplied
# transformations remain fixed. This is conditional selection stability, not
# nested cross-validation or validation of the entire preprocessing/tuning pipeline.
# X must contain ALL original ELN candidate species, not only signature species.

ELN_LOOCV_HELPER_VERSION <- "3.0"

# Use the installed caret glmnet fitting callback, followed by the same
# default coefficient extraction as the original ELN.R. Never invoke an
# implicit exact=TRUE update or force a different lambda sequence.
.eln_fit_coefficients <- function(X, y, alpha_value, lambda_value,
                                  standardize, engine_fit, last = FALSE) {
  fit <- engine_fit(
    x = X, y = y, wts = NULL,
    param = data.frame(alpha = alpha_value, lambda = lambda_value),
    lev = NULL, last = last, classProbs = FALSE,
    family = "gaussian", standardize = standardize, intercept = TRUE
  )
  if (!is.null(fit$jerr) && fit$jerr != 0) {
    stop("glmnet reported a convergence issue; the fold has not been discarded.")
  }
  b <- as.matrix(stats::coef(fit, s = lambda_value, exact = FALSE))
  index <- match(colnames(X), rownames(b))
  if (anyNA(index)) stop("Unmatched species coefficients.")
  beta <- as.numeric(b[index, 1L])
  names(beta) <- colnames(X)
  if (any(!is.finite(beta))) stop("Nonfinite species coefficients.")
  beta
}

.eln_saved_coefficients <- function(entry) {
  if (!length(entry)) return(stats::setNames(numeric(0), character(0)))
  if (is.data.frame(entry) || is.matrix(entry)) {
    if (nrow(entry) != 1L) stop("Each saved sp_id entry must contain one coefficient row.")
    beta <- as.numeric(unlist(entry[1L, , drop = FALSE], use.names = FALSE))
    names(beta) <- colnames(entry)
  } else if (is.numeric(entry) && !is.null(names(entry))) {
    beta <- entry
  } else {
    stop("Unsupported saved coefficient format in sp_id.")
  }
  if (is.null(names(beta)) || anyNA(names(beta)) || anyDuplicated(names(beta)) ||
      any(!is.finite(beta))) stop("Saved coefficients need unique names and finite values.")
  beta[names(beta) != "(Intercept)" & beta != 0]
}

eln_loo_stability <- function(X, y, tuning, reference_species,
                              outcome = "outcome", standardize = TRUE,
                              progress_every = 50L) {
  if (!requireNamespace("glmnet", quietly = TRUE) ||
      !requireNamespace("caret", quietly = TRUE)) {
    stop("The glmnet and caret packages used in the original analysis are required.")
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("X must be numeric and finite; do not silently omit participants.")
  }
  if (!is.numeric(y) || length(y) != nrow(X) || any(!is.finite(y))) {
    stop("y must be a finite numeric vector aligned with the rows of X.")
  }
  if (is.null(colnames(X)) || anyDuplicated(colnames(X)) ||
      anyNA(colnames(X)) || any(!nzchar(colnames(X)))) {
    stop("X needs unique, nonmissing species column names.")
  }
  if (nrow(X) < 4L || ncol(X) < 2L || stats::sd(y) == 0) {
    stop("At least four participants, two candidates, and a variable outcome are required.")
  }
  tuning <- as.data.frame(tuning)
  if (!all(c("alpha", "lambda") %in% names(tuning))) {
    stop("tuning needs alpha and lambda columns.")
  }
  if (!nrow(tuning) || !is.numeric(tuning$alpha) || !is.numeric(tuning$lambda) ||
      any(!is.finite(as.matrix(tuning[, c("alpha", "lambda")])) ) ||
      any(tuning$alpha < 0 | tuning$alpha > 1) || any(tuning$lambda <= 0) ||
      anyDuplicated(tuning$alpha)) {
    stop("Supply one valid positive lambda per distinct alpha setting.")
  }
  reference_species <- unique(as.character(reference_species))
  if (anyNA(reference_species) || !all(reference_species %in% colnames(X))) {
    stop("Every reference species must occur among the candidate columns in X.")
  }
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE.")
  }

  n <- nrow(X)
  p <- ncol(X)
  k <- nrow(tuning)
  ids <- if (is.null(rownames(X))) as.character(seq_len(n)) else rownames(X)
  alpha_names <- paste0("alpha_", format(tuning$alpha, trim = TRUE))
  selected <- array(FALSE, dim = c(p, n, k),
                    dimnames = list(colnames(X), ids, alpha_names))
  ref_index <- match(reference_species, colnames(X))
  engine_fit <- caret::getModelInfo("glmnet", regex = FALSE)[[1L]]$fit

  for (a in seq_len(k)) {
    alpha_value <- as.numeric(tuning$alpha[a])
    lambda_value <- as.numeric(tuning$lambda[a])
    message(sprintf("%s: alpha=%g, fixed lambda=%g (%d leave-one-out fits)",
                    outcome, alpha_value, lambda_value, n))

    for (i in seq_len(n)) {
      X_train <- X[-i, , drop = FALSE]
      y_train <- y[-i]
      if (stats::sd(y_train) == 0) {
        stop(sprintf("%s becomes constant when participant %s is omitted.", outcome, ids[i]))
      }
      beta <- .eln_fit_coefficients(
        X_train, y_train, alpha_value, lambda_value,
        standardize = standardize, engine_fit = engine_fit, last = FALSE
      )
      selected[, i, a] <- beta != 0
      if (progress_every > 0 && (i %% progress_every == 0 || i == n)) {
        message(sprintf("  completed %d/%d", i, n))
      }
    }
  }

  # Per-species selection frequency separately for each alpha setting.
  by_alpha <- do.call(rbind, lapply(seq_len(k), function(a) {
    m <- matrix(selected[, , a], nrow = p, ncol = n)
    data.frame(outcome = outcome, alpha = tuning$alpha[a],
               lambda = tuning$lambda[a], species = colnames(X),
               reference_signature = colnames(X) %in% reference_species,
               n_folds = n, n_selected = rowSums(m),
               selection_frequency = rowMeans(m), row.names = NULL)
  }))

  # Reapply the original consensus rule inside each leave-one-out training set:
  # a species must have nonzero coefficients at every alpha setting in that fold.
  consensus <- apply(selected, c(1L, 2L), all)
  frequency <- data.frame(outcome = outcome, species = colnames(X),
                          reference_signature = colnames(X) %in% reference_species,
                          n_folds = n, n_selected = rowSums(consensus),
                          selection_frequency = rowMeans(consensus), row.names = NULL)
  frequency <- frequency[order(!frequency$reference_signature,
                               -frequency$selection_frequency, frequency$species), ]
  rownames(frequency) <- NULL

  retained <- if (length(ref_index)) colSums(consensus[ref_index, , drop = FALSE]) else rep(NA_real_, n)
  fold_summary <- data.frame(outcome = outcome, omitted_participant = ids,
                             n_signature_species = colSums(consensus),
                             n_reference_species_retained = retained,
                             proportion_reference_retained = if (length(ref_index)) {
                               retained / length(ref_index)
                             } else rep(NA_real_, n))
  ref_frequency <- rowMeans(consensus)[ref_index]
  n_union <- sum(rowSums(consensus) > 0)
  n_intersection <- sum(rowSums(consensus) == n)
  summary <- data.frame(
    outcome = outcome, n_participants = n, n_candidate_species = p, n_alpha = k,
    n_reference_signature_species = length(ref_index),
    n_loo_signature_union = n_union, n_loo_signature_intersection = n_intersection,
    proportion_intersection_of_union = if (n_union) n_intersection/n_union else NA_real_,
    median_reference_selection_frequency = if (length(ref_index)) stats::median(ref_frequency) else NA_real_,
    minimum_reference_selection_frequency = if (length(ref_index)) min(ref_frequency) else NA_real_,
    proportion_reference_selected_every_fold = if (length(ref_index)) mean(ref_frequency == 1) else NA_real_,
    median_reference_retention_per_fold = if (length(ref_index)) {
      stats::median(retained / length(ref_index))
    } else NA_real_
  )
  list(summary = summary, species_frequency = frequency,
       frequency_by_alpha = by_alpha, fold_summary = fold_summary,
       selected_by_fold_and_alpha = selected, consensus_by_fold = consensus,
       tuning = tuning, reference_species = reference_species,
       assessment = "Leave-one-out feature selection with preprocessing and tuned parameters fixed",
       helper_version = "3.0",
       coefficient_extraction = "caret glmnet callback; coef(s=lambda, exact=FALSE)",
       glmnet_version = as.character(utils::packageVersion("glmnet")),
       caret_version = as.character(utils::packageVersion("caret")))
}


eln_loo_from_saved <- function(X, data, eln_results, outcomes = NULL,
                               output_dir = NULL, standardize = TRUE,
                               progress_every = 50L, check_only = FALSE,
                               coefficient_tolerance = 1e-6) {
  if (!requireNamespace("glmnet", quietly = TRUE) ||
      !requireNamespace("caret", quietly = TRUE)) {
    stop("The glmnet and caret packages used in the original analysis are required.")
  }
  X <- as.matrix(X)
  if (!length(eln_results)) stop("Supply at least one saved ELN result object.")
  if (!is.numeric(X) || any(!is.finite(X)) || nrow(X) < 4L || ncol(X) < 2L ||
      is.null(colnames(X)) || anyDuplicated(colnames(X))) {
    stop("Supply a finite numeric candidate matrix with unique species column names.")
  }
  if (!is.logical(check_only) || length(check_only) != 1L || is.na(check_only)) {
    stop("check_only must be TRUE or FALSE.")
  }
  if (!is.numeric(coefficient_tolerance) || length(coefficient_tolerance) != 1L ||
      !is.finite(coefficient_tolerance) || coefficient_tolerance < 0) {
    stop("coefficient_tolerance must be one finite nonnegative number.")
  }
  if (is.null(rownames(X)) || anyDuplicated(rownames(X)) ||
      is.null(rownames(data)) || anyDuplicated(rownames(data)) ||
      !all(rownames(X) %in% rownames(data))) {
    stop("X and data need matching unique participant row names. Align by participant ID first.")
  }
  if (is.null(outcomes)) outcomes <- names(eln_results[[1L]]$sp_id)
  if (!length(outcomes) || !all(outcomes %in% names(data))) {
    stop(paste("Missing sweetener outcome columns:", paste(setdiff(outcomes, names(data)), collapse = ", ")))
  }
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  aligned <- data[match(rownames(X), rownames(data)), , drop = FALSE]
  results <- setNames(vector("list", length(outcomes)), outcomes)
  prepared <- setNames(vector("list", length(outcomes)), outcomes)
  engine_fit <- caret::getModelInfo("glmnet", regex = FALSE)[[1L]]$fit
  checks <- list()

  # Verify ALL requested outcomes before starting any leave-one-out fitting.
  for (outcome in outcomes) {
    valid <- vapply(eln_results, function(r) {
      !is.null(r$tuning_para) && outcome %in% rownames(r$tuning_para) &&
        outcome %in% names(r$sp_id)
    }, logical(1))
    if (!all(valid)) stop(sprintf("Saved results lack tuning/species information for %s.", outcome))
    tuning <- do.call(rbind, lapply(eln_results, function(r) {
      r$tuning_para[outcome, c("alpha", "lambda"), drop = FALSE]
    }))
    rownames(tuning) <- NULL
    original_sets <- lapply(eln_results, function(r) {
      names(.eln_saved_coefficients(r$sp_id[[outcome]]))
    })
    if (!all(unlist(original_sets, use.names = FALSE) %in% colnames(X))) {
      stop(sprintf("Some saved %s species are missing from X; supply the complete original candidate matrix.", outcome))
    }
    reference <- Reduce(intersect, original_sets)
    y <- aligned[[outcome]]
    if (!is.numeric(y) || any(!is.finite(y)) || stats::sd(y) == 0) {
      stop(sprintf("%s must be finite, numeric, and variable.", outcome))
    }
    if (any(!is.finite(as.matrix(tuning))) || anyDuplicated(tuning$alpha) ||
        any(tuning$alpha < 0 | tuning$alpha > 1) || any(tuning$lambda <= 0)) {
      stop(sprintf("Invalid saved tuning parameters for %s.", outcome))
    }
    for (a in seq_len(nrow(tuning))) {
      beta <- .eln_fit_coefficients(
        X, y, alpha_value = as.numeric(tuning$alpha[a]),
        lambda_value = as.numeric(tuning$lambda[a]),
        standardize = standardize, engine_fit = engine_fit, last = TRUE
      )
      saved <- .eln_saved_coefficients(eln_results[[a]]$sp_id[[outcome]])
      saved_full <- stats::setNames(rep(0, ncol(X)), colnames(X))
      saved_full[names(saved)] <- saved
      refit_species <- names(beta)[beta != 0]
      missing <- setdiff(names(saved), refit_species)
      additional <- setdiff(refit_species, names(saved))
      difference <- abs(beta - saved_full)
      coefficient_match <- all(difference <= coefficient_tolerance *
                                 pmax(1, abs(beta), abs(saved_full)))
      set_match <- !length(missing) && !length(additional)
      checks[[length(checks) + 1L]] <- data.frame(
        outcome = outcome, alpha = tuning$alpha[a], lambda = tuning$lambda[a],
        n_saved = length(saved), n_refit = length(refit_species),
        n_missing_saved = length(missing), n_additional_refit = length(additional),
        missing_saved = paste(missing, collapse = ";"),
        additional_refit = paste(additional, collapse = ";"),
        max_absolute_coefficient_difference = max(difference),
        selected_set_match = set_match, coefficient_match = coefficient_match,
        baseline_passed = set_match && coefficient_match,
        row.names = NULL
      )
    }
    prepared[[outcome]] <- list(tuning = tuning, reference = reference)
  }
  baseline_check <- do.call(rbind, checks)
  rownames(baseline_check) <- NULL
  baseline_passed <- all(baseline_check$baseline_passed)
  if (!is.null(output_dir)) {
    utils::write.csv(baseline_check, file.path(output_dir, "full_data_reference_check.csv"),
                      row.names = FALSE)
  }
  message(sprintf("Full-data reference check: %d/%d models match saved selections and coefficients.",
                  sum(baseline_check$baseline_passed), nrow(baseline_check)))
  if (check_only) {
    return(list(baseline_check = baseline_check, baseline_passed = baseline_passed,
                helper_version = "3.0"))
  }
  if (!baseline_passed) {
    stop("Full-data refits do not reproduce the saved models. Inspect full_data_reference_check.csv or rerun with check_only=TRUE. No LOO fits were started.")
  }

  for (outcome in outcomes) {
    tuning <- prepared[[outcome]]$tuning
    reference <- prepared[[outcome]]$reference
    ans <- eln_loo_stability(X = X, y = aligned[[outcome]], tuning = tuning,
                              reference_species = reference, outcome = outcome,
                              standardize = standardize,
                              progress_every = progress_every)
    ans$baseline_check <- baseline_check[baseline_check$outcome == outcome, , drop = FALSE]
    results[[outcome]] <- ans
    if (!is.null(output_dir)) {
      stem <- gsub("[^A-Za-z0-9_.-]", "_", outcome)
      saveRDS(ans, file.path(output_dir, paste0(stem, "_LOO.rds")))
      utils::write.csv(ans$species_frequency,
                        file.path(output_dir, paste0(stem, "_species_frequency.csv")), row.names = FALSE)
      utils::write.csv(ans$frequency_by_alpha,
                        file.path(output_dir, paste0(stem, "_frequency_by_alpha.csv")), row.names = FALSE)
      utils::write.csv(ans$fold_summary,
                        file.path(output_dir, paste0(stem, "_fold_summary.csv")), row.names = FALSE)
    }
  }
  combined <- do.call(rbind, lapply(results, function(z) z$summary))
  rownames(combined) <- NULL
  if (!is.null(output_dir)) {
    utils::write.csv(combined, file.path(output_dir, "LOO_selection_summary.csv"), row.names = FALSE)
  }
  list(summary = combined, results = results, baseline_check = baseline_check,
       baseline_passed = baseline_passed, helper_version = "3.0")
}
