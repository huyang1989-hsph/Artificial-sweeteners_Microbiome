# Short verification with simulated data; no study data are used.
# Run in the directory containing SCCA_LOO_stability.R.
source("SCCA_LOO_stability.R")
stopifnot(requireNamespace("PMA", quietly = TRUE))
set.seed(2026)
n <- 36L
latent <- matrix(rnorm(n * 2L), nrow = n)
X_test <- cbind(
  latent[, 1L] + rnorm(n, sd = .4),
  latent[, 1L] + rnorm(n, sd = .5),
  latent[, 1L] + rnorm(n, sd = .6),
  latent[, 2L] + rnorm(n, sd = .4),
  latent[, 2L] + rnorm(n, sd = .5),
  latent[, 2L] + rnorm(n, sd = .6),
  matrix(rnorm(n * 6L), nrow = n))
Y_test <- cbind(
  latent[, 1L] + rnorm(n, sd = .5),
  latent[, 1L] + rnorm(n, sd = .7),
  latent[, 2L] + rnorm(n, sd = .5),
  latent[, 2L] + rnorm(n, sd = .7))
rownames(X_test) <- rownames(Y_test) <- paste0("participant_", seq_len(n))
colnames(X_test) <- paste0("species_", seq_len(ncol(X_test)))
colnames(Y_test) <- paste0("sweetener_", seq_len(ncol(Y_test)))
fit_test <- PMA::CCA(X_test, Y_test, typex = "standard", typez = "standard",
                      penaltyx = .65, penaltyz = .85, K = 2,
                      standardize = TRUE, trace = FALSE)

# The original orientation/order is identifiable only up to a joint sign and
# component permutation. Verify those equivalences are handled correctly.
fit_relabelled <- fit_test
fit_relabelled$u <- fit_test$u[, c(2L, 1L), drop = FALSE]
fit_relabelled$v <- fit_test$v[, c(2L, 1L), drop = FALSE]
fit_relabelled$u[, 1L] <- -fit_relabelled$u[, 1L]
fit_relabelled$v[, 1L] <- -fit_relabelled$v[, 1L]
fit_relabelled$cors <- fit_test$cors[c(2L, 1L)]
relabelled_check <- scca_loo_stability(
  X_test, Y_test, fit_relabelled, analysis = "verification", check_only = TRUE)
stopifnot(relabelled_check$baseline_passed)

pilot <- scca_loo_stability(
  X_test, Y_test[n:1L, , drop = FALSE], fit_test,
  analysis = "verification", omit_ids = rownames(X_test)[1:3], progress_every = 0)
stopifnot(pilot$baseline_passed, nrow(pilot$fold_results) == 6L,
          all(pilot$component_overview$n_refits == 3L),
          !pilot$settings$is_complete_leave_one_out)

# Independently refit the first omission and compare its actual coefficient
# vectors to those saved by the helper, allowing the same order/sign freedoms.
manual <- PMA::CCA(X_test[-1L, , drop = FALSE], Y_test[-1L, , drop = FALSE],
                    typex = "standard", typez = "standard", penaltyx = .65,
                    penaltyz = .85, K = 2, standardize = TRUE, trace = FALSE)
manual_aligned <- .scca_loo_align(fit_test$u, fit_test$v, manual$u, manual$v)
stopifnot(max(abs(pilot$aligned_weights_x[, , 1L] - manual_aligned$u)) < 1e-10,
          max(abs(pilot$aligned_weights_y[, , 1L] - manual_aligned$v)) < 1e-10)
manual_sx <- scale(X_test[-1L, , drop = FALSE]) %*% manual_aligned$u
manual_sy <- scale(Y_test[-1L, , drop = FALSE]) %*% manual_aligned$v
stopifnot(max(abs(pilot$fold_results$canonical_correlation[1:2] -
                    diag(cor(manual_sx, manual_sy)))) < 1e-10)

# A changed data matrix must not pass as the saved full-data model.
bad_X <- -X_test
bad_check <- scca_loo_stability(
  bad_X, Y_test, fit_test, analysis = "verification", check_only = TRUE)
stopifnot(!bad_check$baseline_passed)

cat("Verification passed: component alignment, participant alignment, direct refit agreement, and mismatch detection.\n")
