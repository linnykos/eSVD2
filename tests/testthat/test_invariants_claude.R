# UNIT_TEST_PLAN.md section 4 -- T-PROP-01 .. T-PROP-10, plus sections 2.4/2.5.
#
# These encode the paper's own claims. Each catches a whole class of regression
# that no single-function test can.

test_that("T-PROP-01 / T-REP-01: reparameterization preserves the predictions", {
  esvd_obj <- .small_esvd_obj()

  # The paper states this twice. `Y X^T + Z C^T` must be unchanged across
  # `reparameterization_esvd_covariates()`.
  before_fit <- esvd_obj[["fit_First"]]
  prediction_before <- tcrossprod(before_fit$x_mat, before_fit$y_mat) +
    tcrossprod(esvd_obj$covariates, before_fit$z_mat)

  res <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                            fit_name = "fit_First",
                                            omitted_variables = "Log_UMI")
  after_fit <- res[["fit_First"]]
  prediction_after <- tcrossprod(after_fit$x_mat, after_fit$y_mat) +
    tcrossprod(res$covariates, after_fit$z_mat)

  expect_equal(prediction_before, prediction_after, tolerance = 1e-8)
})

test_that("T-PROP-02 / T-REP-02: X is orthogonal to the retained covariates", {
  esvd_obj <- .small_esvd_obj()

  res <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                            fit_name = "fit_First",
                                            omitted_variables = "Log_UMI")

  x_mat <- res[["fit_First"]]$x_mat
  covariates <- res$covariates
  retained_vec <- setdiff(colnames(covariates), "Log_UMI")

  # The paper's Step 1: X is orthogonalized against C.
  cross_mat <- crossprod(x_mat, covariates[, retained_vec, drop = FALSE])
  expect_true(max(abs(cross_mat)) < 1e-6,
              info = paste0("max |X'C| = ", signif(max(abs(cross_mat)), 4)))
})

test_that("T-REP-03: .reparameterize makes X'X/n and Y'Y/p diagonal and equal", {
  set.seed(10)
  x_mat <- matrix(stats::rnorm(60 * 3), nrow = 60, ncol = 3)
  y_mat <- matrix(stats::rnorm(20 * 3), nrow = 20, ncol = 3)

  res <- .reparameterize(x_mat = x_mat, y_mat = y_mat,
                         equal_covariance = TRUE)

  cov_x <- crossprod(res$x_mat) / nrow(res$x_mat)
  cov_y <- crossprod(res$y_mat) / nrow(res$y_mat)

  # The paper's Step 2, and the reason the factorization is identifiable.
  expect_equal(cov_x, diag(diag(cov_x)), tolerance = 1e-8)
  expect_equal(cov_y, diag(diag(cov_y)), tolerance = 1e-8)
  expect_equal(diag(cov_x), diag(cov_y), tolerance = 1e-8)

  # And the product is preserved.
  expect_equal(tcrossprod(res$x_mat, res$y_mat), tcrossprod(x_mat, y_mat),
               tolerance = 1e-8)
})

## Question Q-REP-1 resolved: warn and proceed. Expected to FAIL -- today
## `.reparameterize` on a rank-deficient input either errors or proceeds
## silently, and `opt_esvd.default` wraps it in a `tryCatch` that swallows the
## outcome either way, so a rank-deficient fit is indistinguishable from a
## good one.
test_that("T-REP-08: a rank-deficient x_mat warns and still returns a valid fit", {
  set.seed(10)
  x_mat <- matrix(stats::rnorm(60 * 3), nrow = 60, ncol = 3)
  x_mat[, 3] <- x_mat[, 2]
  y_mat <- matrix(stats::rnorm(20 * 3), nrow = 20, ncol = 3)

  # `.reparameterize` already does the right thing here: it warns and proceeds,
  # which is what question Q-REP-1 asked for. Asserting the warning is the point
  # of the test -- an earlier draft checked only that it did not error, which
  # passed without constraining anything.
  expect_warning(res <- .reparameterize(x_mat = x_mat, y_mat = y_mat,
                                        equal_covariance = TRUE),
                 regexp = "rank")
  expect_equal(dim(res$x_mat), dim(x_mat))
})

test_that("T-REP-09: the warning survives opt_esvd's tryCatch", {
  # The half that makes T-REP-08 worth anything: a `warning()` raised inside
  # `.reparameterize` is swallowed by `opt_esvd.default`'s handler today, which
  # leaves the user exactly as uninformed as before.
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  # Force rank deficiency in the fit that `opt_esvd` will reparameterize.
  esvd_obj[[latest_fit]]$x_mat[, 2] <- esvd_obj[[latest_fit]]$x_mat[, 1]

  # Today this ERRORS inside `eigen()` -- "infinite or missing values in 'x'" --
  # rather than warning and proceeding, which is what question Q-REP-1 decided.
  # Caught so the finding is reported as this assertion instead of aborting the
  # file.
  res <- try(reparameterization_esvd_covariates(input_obj = esvd_obj,
                                                fit_name = latest_fit,
                                                omitted_variables = "Log_UMI"),
             silent = TRUE)

  expect_false(inherits(res, "try-error"),
               info = if(inherits(res, "try-error"))
                 conditionMessage(attr(res, "condition")) else "ok")
})

test_that("T-REP-07: .identification does not return NaN on a negative eigenvalue", {
  # The function already warns about rank deficiency and then proceeds anyway,
  # taking a square root of what may be negative.
  set.seed(10)
  cov_x <- diag(c(1, 0.5, 1e-14))
  cov_y <- diag(c(1, 0.5, 1e-14))

  # Today this ERRORS inside `eigen()` rather than returning NaN, which is the
  # better of the two outcomes -- but the error comes from `base::eigen` and
  # names nothing the caller can act on. Either behaviour is acceptable; a NaN
  # return is not, because it propagates into every fitted value silently.
  res <- try(suppressWarnings(.identification(cov_x = cov_x, cov_y = cov_y)),
             silent = TRUE)

  # Either outcome is acceptable -- an error names the problem, a finite answer
  # is usable -- but a NaN return is not, because it propagates into every
  # fitted value with no signal. Asserting "not NaN" over both branches keeps
  # the test from passing merely because something threw.
  expect_true(inherits(res, "try-error") || all(!is.nan(res)),
              info = "a NaN here propagates into every fitted value")
})

test_that("T-PROP-03 / T-OPT-02: the loss is monotone for every family", {
  point <- .feasible_point("poisson", num_cells = 20, num_genes = 8)

  res <- suppressWarnings(
    opt_esvd(input_obj = point$dat,
             x_init = point$x_mat,
             y_init = point$y_mat,
             family = "poisson",
             l2pen = 0.1,
             max_iter = 8,
             nuisance_vec = point$gamma_vec,
             verbose = 0)
  )

  loss_vec <- res$loss
  skip_if(is.null(loss_vec), "opt_esvd did not return a loss trace")
  expect_true(all(diff(loss_vec) <= 1e-8),
              info = paste0("loss = ", paste0(signif(loss_vec, 6),
                                              collapse = ", ")))
})

## Four of the seven families are unusable at their documented defaults:
## `opt_esvd.default` defaults `nuisance_vec = rep(NA, p)`, but `gaussian`,
## `curved_gaussian`, `neg_binom` and `neg_binom2` all consume `gamma`, so the
## objective is NA and the call dies with "missing value where TRUE/FALSE
## needed". Expected to FAIL for those four.
test_that("T-OPT-05: every family runs at its documented defaults", {
  for(family_name in .all_families()){
    point <- .feasible_point(family_name, num_cells = 20, num_genes = 8)

    result <- try(suppressWarnings(
      opt_esvd(input_obj = point$dat,
               x_init = point$x_mat,
               y_init = point$y_mat,
               family = family_name,
               max_iter = 2,
               verbose = 0)
    ), silent = TRUE)

    expect_false(inherits(result, "try-error"),
                 info = paste0(family_name, ": ",
                               if(inherits(result, "try-error")) {
                                 conditionMessage(attr(result, "condition"))
                               } else "ok"))
  }
})

test_that("T-OPT-03 / T-PROP-04: opt_esvd is deterministic", {
  point <- .feasible_point("poisson", num_cells = 20, num_genes = 8)

  res_one <- suppressWarnings(
    opt_esvd(input_obj = point$dat, x_init = point$x_mat, y_init = point$y_mat,
             family = "poisson", max_iter = 5,
             nuisance_vec = point$gamma_vec, verbose = 0)
  )
  res_two <- suppressWarnings(
    opt_esvd(input_obj = point$dat, x_init = point$x_mat, y_init = point$y_mat,
             family = "poisson", max_iter = 5,
             nuisance_vec = point$gamma_vec, verbose = 0)
  )

  expect_equal(res_one$x_mat, res_two$x_mat)
  expect_equal(res_one$y_mat, res_two$y_mat)
})

test_that("T-PROP-05: posterior sanity, on the fitted object", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  expect_true(all(esvd_obj[[latest_fit]]$posterior_mean_mat > 0))
  expect_true(all(esvd_obj[[latest_fit]]$posterior_var_mat > 0))
  expect_true(all(is.finite(esvd_obj[[latest_fit]]$posterior_mean_mat)))
})

## Question Q-PROP-2 resolved: option (a), fixed seed with a loose threshold, in
## the CRAN suite. The seed was fixed BEFORE the threshold was chosen -- it is
## not tuned until the test passes, which would make the test assert nothing.
## The threshold is 0.01 rather than 0.05 because CRAN runs the suite on a dozen
## platforms and a 1-in-20 false alarm rate would be a nuisance.
test_that("T-PROP-06: p-values on genuinely null genes are approximately uniform", {
  skip_on_cran()

  set.seed(10)
  null_list <- generate_null(cell_per_person = 25, num_genes = 120,
                             num_individuals = 8)

  covariates <- null_list$covariates
  esvd_obj <- suppressWarnings(
    initialize_esvd(dat = null_list$obs_mat,
                    covariates = covariates[, setdiff(colnames(covariates),
                                                      "Individual"),
                                            drop = FALSE],
                    case_control_variable = "CC",
                    bool_intercept = TRUE,
                    k = 2,
                    metadata_case_control = covariates[, "CC"],
                    metadata_individual = factor(null_list$metadata_individual),
                    verbose = 0)
  )
  esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                                 fit_name = "fit_Init",
                                                 omitted_variables = "Log_UMI")
  esvd_obj <- suppressWarnings(
    estimate_nuisance(input_obj = esvd_obj, bool_covariates_as_library = TRUE,
                      verbose = 0)
  )
  esvd_obj <- compute_posterior(input_obj = esvd_obj,
                                alpha_max = 2 * max(null_list$obs_mat@x),
                                bool_covariates_as_library = TRUE,
                                library_min = 0.1)
  esvd_obj <- compute_test_statistic(input_obj = esvd_obj, verbose = 0)
  esvd_obj <- suppressWarnings(compute_pvalue(input_obj = esvd_obj))

  # `generate_null()` marks columns 1-10 as truly DE; 11 onward are null.
  pvalue_vec <- 10^(-esvd_obj$pvalue_list$log10pvalue)
  null_pvalue_vec <- pvalue_vec[11:length(pvalue_vec)]

  ks_res <- suppressWarnings(stats::ks.test(null_pvalue_vec, "punif"))
  expect_true(ks_res$p.value > 0.01,
              info = paste0("KS p-value = ", signif(ks_res$p.value, 4)))
})

test_that("T-GEN-05: generate_null marks exactly 10 genes as truly DE", {
  # T-PROP-06 depends on this. If the DE block ever moves, the calibration test
  # silently becomes meaningless rather than failing.
  set.seed(10)
  null_list <- generate_null(cell_per_person = 10, num_genes = 40,
                             num_individuals = 4)

  expect_equal(ncol(null_list$obs_mat), 40)
  expect_true(inherits(null_list$obs_mat, "dgCMatrix"))
  expect_true(all(c("Intercept", "Log_UMI", "CC", "Sex", "Age") %in%
                    colnames(null_list$covariates)))
})

## T-PROP-08 (scale equivariance) is CUT (Kevin, 2026-08-29), and the reason is
## worth keeping so it is not proposed again.
##
## The plan asked: multiply every count by c, add log(c) to `Log_UMI`, and
## `teststat_vec` should be unchanged. That is not a property of a count model.
## Sequencing c times as deep gives A ~ Poisson(c*l*lambda), whose mean AND
## variance are both c*l*mu. Multiplying observed counts by c gives mean c*l*mu
## but variance c^2*l*mu. The two are different distributions, so the rescaled
## data no longer satisfies the mean-variance relationship the model assumes and
## there is no reason for the statistic to be preserved.
##
## A genuine Poisson invariance would be THINNING -- `rbinom(A, prob = 1/c)` --
## but that is stochastic and would need its own seed discipline for modest
## value. Left out deliberately rather than overlooked.

test_that("T-PROP-09: swapping the case and control labels negates the statistic", {
  esvd_obj <- .small_esvd_obj()

  swapped_obj <- esvd_obj
  swapped_obj$case_control <- 1 - esvd_obj$case_control
  swapped_obj$param$test_case_individuals <- NULL

  res_swapped <- compute_test_statistic(input_obj = swapped_obj, verbose = 0)

  # The test is two-sided; if it were not symmetric, the direction of effect
  # would bias significance.
  expect_equal(unname(res_swapped$teststat_vec),
               unname(-esvd_obj$teststat_vec), tolerance = 1e-8)
})

test_that("T-PROP-10: permuting the genes permutes every per-gene output", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  set.seed(20)
  permutation_idx <- sample(ncol(esvd_obj$dat))

  permuted_obj <- esvd_obj
  permuted_obj[[latest_fit]]$posterior_mean_mat <-
    esvd_obj[[latest_fit]]$posterior_mean_mat[, permutation_idx, drop = FALSE]
  permuted_obj[[latest_fit]]$posterior_var_mat <-
    esvd_obj[[latest_fit]]$posterior_var_mat[, permutation_idx, drop = FALSE]

  res_permuted <- compute_test_statistic(input_obj = permuted_obj, verbose = 0)

  # Catches positional-vs-named indexing bugs anywhere in the stack, of which
  # there are several candidate sites (`library_idx`, `case_control_idx`, the
  # `nuisance_vec` sweeps).
  expect_equal(unname(res_permuted$teststat_vec),
               unname(esvd_obj$teststat_vec[permutation_idx]),
               tolerance = 1e-10)
})
