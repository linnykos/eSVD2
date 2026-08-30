# UNIT_TEST_PLAN.md section 2.7 -- T-POST-01 .. T-POST-11.
#
# `compute_posterior()` currently has three assertions, all about `dim()`.
# Everything below is new.

test_that("T-POST-01: posterior mean and variance are strictly positive and finite", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  mean_mat <- esvd_obj[[latest_fit]]$posterior_mean_mat
  var_mat <- esvd_obj[[latest_fit]]$posterior_var_mat

  expect_true(all(is.finite(mean_mat)))
  expect_true(all(is.finite(var_mat)))
  expect_true(all(mean_mat > 0))
  expect_true(all(var_mat > 0))
})

test_that("T-POST-02: var equals mean / SplusBeta, exactly", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  res <- compute_posterior(input_obj = esvd_obj,
                           alpha_max = 2 * max(esvd_obj$dat),
                           bool_covariates_as_library = TRUE,
                           bool_return_components = TRUE,
                           library_min = 0.1)

  mean_mat <- res[[latest_fit]]$posterior_mean_mat
  var_mat <- res[[latest_fit]]$posterior_var_mat
  numerator_mat <- res[[latest_fit]]$numerator_mat
  denominator_mat <- res[[latest_fit]]$denominator_mat

  # Equation 14: mean = num/denom, var = num/denom^2, so var = mean/denom.
  skip_if(is.null(numerator_mat), "bool_return_components did not store components")
  expect_equal(as.numeric(mean_mat),
               as.numeric(numerator_mat / denominator_mat), tolerance = 1e-10)
  expect_equal(as.numeric(var_mat),
               as.numeric(numerator_mat / denominator_mat^2), tolerance = 1e-10)
  expect_equal(as.numeric(var_mat),
               as.numeric(mean_mat / denominator_mat), tolerance = 1e-10)
})

test_that("T-POST-05: bool_return_components returns the documented debugging hook", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  res <- compute_posterior(input_obj = esvd_obj,
                           alpha_max = 2 * max(esvd_obj$dat),
                           bool_covariates_as_library = TRUE,
                           bool_return_components = TRUE,
                           library_min = 0.1)

  # An untested branch, and the components are the documented debugging hook.
  expect_true("numerator_mat" %in% names(res[[latest_fit]]))
  expect_true("denominator_mat" %in% names(res[[latest_fit]]))
})

test_that("T-POST-06: the two experimental booleans are mutually exclusive", {
  esvd_obj <- .small_esvd_obj()

  expect_error(compute_posterior(input_obj = esvd_obj,
                                 bool_adjust_covariates = TRUE,
                                 bool_covariates_as_library = TRUE))

  # And the experimental branch at least runs.
  expect_no_error(
    compute_posterior(input_obj = esvd_obj,
                      alpha_max = 2 * max(esvd_obj$dat),
                      bool_adjust_covariates = TRUE,
                      bool_covariates_as_library = FALSE,
                      library_min = 0.1)
  )
})

## Section 1.4, [verified]: `scale()` returns an n-by-1 MATRIX and moves
## `names()` to `rownames()`. So on the branch where the stabilization fires,
## `nuisance_vec` silently changes type and loses its names. Latent today,
## a bug on the next edit.
test_that("T-POST-07: nuisance_vec keeps its type and names on both stabilization branches", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]
  gene_vec <- colnames(esvd_obj$dat)

  for(bool_stabilize in c(TRUE, FALSE)){
    label <- paste0("bool_stabilize_underdispersion = ", bool_stabilize)

    res <- compute_posterior(input_obj = esvd_obj,
                             alpha_max = 2 * max(esvd_obj$dat),
                             bool_covariates_as_library = TRUE,
                             bool_stabilize_underdispersion = bool_stabilize,
                             library_min = 0.1)

    nuisance_vec <- res[[latest_fit]]$nuisance_vec
    expect_true(is.numeric(nuisance_vec), info = label)
    expect_false(is.matrix(nuisance_vec), info = label)
    expect_equal(names(nuisance_vec), gene_vec, info = label)
  }
})

## Question Q-POST-1 resolved: the CODE is right and the prose is wrong. The
## stabilization fires when `mean(log10(nuisance_vec)) > 0`, i.e. when the
## cohort is on average UNDER-dispersed -- `nuisance_vec` being the Gamma rate,
## a large value means little over-dispersion. The fix is a roxygen rewrite.
test_that("T-POST-08: stabilization fires when mean(log10(nuisance_vec)) > 0", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  # Force the under-dispersed regime: large rate = small over-dispersion.
  esvd_obj[[latest_fit]]$nuisance_vec[] <- 100
  res_under <- compute_posterior(input_obj = esvd_obj,
                                 alpha_max = 2 * max(esvd_obj$dat),
                                 bool_covariates_as_library = TRUE,
                                 bool_stabilize_underdispersion = TRUE,
                                 library_min = 0.1)
  res_under_off <- compute_posterior(input_obj = esvd_obj,
                                     alpha_max = 2 * max(esvd_obj$dat),
                                     bool_covariates_as_library = TRUE,
                                     bool_stabilize_underdispersion = FALSE,
                                     library_min = 0.1)
  expect_false(isTRUE(all.equal(
    as.numeric(res_under[[latest_fit]]$posterior_var_mat),
    as.numeric(res_under_off[[latest_fit]]$posterior_var_mat)
  )))

  # And the over-dispersed regime leaves it alone.
  esvd_obj[[latest_fit]]$nuisance_vec[] <- 0.01
  res_over <- compute_posterior(input_obj = esvd_obj,
                                alpha_max = 2 * max(esvd_obj$dat),
                                bool_covariates_as_library = TRUE,
                                bool_stabilize_underdispersion = TRUE,
                                library_min = 0.1)
  res_over_off <- compute_posterior(input_obj = esvd_obj,
                                    alpha_max = 2 * max(esvd_obj$dat),
                                    bool_covariates_as_library = TRUE,
                                    bool_stabilize_underdispersion = FALSE,
                                    library_min = 0.1)
  expect_equal(as.numeric(res_over[[latest_fit]]$posterior_var_mat),
               as.numeric(res_over_off[[latest_fit]]$posterior_var_mat))
})

test_that("T-POST-09: library_min actually binds", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  res_low <- compute_posterior(input_obj = esvd_obj,
                               alpha_max = 2 * max(esvd_obj$dat),
                               bool_covariates_as_library = TRUE,
                               bool_return_components = TRUE,
                               library_min = 0.1)
  res_high <- compute_posterior(input_obj = esvd_obj,
                                alpha_max = 2 * max(esvd_obj$dat),
                                bool_covariates_as_library = TRUE,
                                bool_return_components = TRUE,
                                library_min = 1e6)

  # `library_min` is the parameter whose default DISAGREES between the two
  # pipelines (section 1.5), so pinning that it binds at all is the
  # precondition for T-PG-03 meaning anything.
  expect_false(isTRUE(all.equal(
    as.numeric(res_low[[latest_fit]]$posterior_mean_mat),
    as.numeric(res_high[[latest_fit]]$posterior_mean_mat)
  )))
})

test_that("T-POST-10: alpha_max binds", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  res_loose <- compute_posterior(input_obj = esvd_obj,
                                 alpha_max = 1e6,
                                 bool_covariates_as_library = TRUE,
                                 library_min = 0.1)
  res_tight <- compute_posterior(input_obj = esvd_obj,
                                 alpha_max = 1,
                                 bool_covariates_as_library = TRUE,
                                 library_min = 0.1)

  expect_false(isTRUE(all.equal(
    as.numeric(res_loose[[latest_fit]]$posterior_mean_mat),
    as.numeric(res_tight[[latest_fit]]$posterior_mean_mat)
  )))
})

test_that("T-POST-11: nuisance_lower_quantile floors the low genes", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  res_default <- compute_posterior(input_obj = esvd_obj,
                                   alpha_max = 2 * max(esvd_obj$dat),
                                   bool_covariates_as_library = TRUE,
                                   library_min = 0.1,
                                   nuisance_lower_quantile = 0.01)
  res_median <- compute_posterior(input_obj = esvd_obj,
                                  alpha_max = 2 * max(esvd_obj$dat),
                                  bool_covariates_as_library = TRUE,
                                  library_min = 0.1,
                                  nuisance_lower_quantile = 0.5)

  # An untested parameter. Flooring half the genes at the median must change
  # the answer for at least those genes.
  expect_false(isTRUE(all.equal(
    as.numeric(res_default[[latest_fit]]$posterior_mean_mat),
    as.numeric(res_median[[latest_fit]]$posterior_mean_mat)
  )))
})

test_that("T-POST-12: posterior matrices carry the dimnames of dat", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  # `compute_test_statistic` names `teststat_vec` from these, and
  # `report_results` builds its `genes` column from that.
  expect_equal(colnames(esvd_obj[[latest_fit]]$posterior_mean_mat),
               colnames(esvd_obj$dat))
  expect_equal(rownames(esvd_obj[[latest_fit]]$posterior_mean_mat),
               rownames(esvd_obj$dat))
})
