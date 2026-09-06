context("Test compute_test_per_gene")

## compute_test_statistic is correct

test_that("compute_test_per_gene works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  covariates <- .get_object(eSVD_obj = eSVD_obj, what_obj = "covariates", which_fit = NULL)
  cc_vec <- covariates[,"case_control_1"]
  cc_levels <- sort(unique(cc_vec), decreasing = F)
  control_idx <- which(cc_vec == cc_levels[1])
  case_idx <- which(cc_vec == cc_levels[2])

  individual_vec <- metadata[,"individual"]
  control_individuals <- as.character(unique(individual_vec[control_idx]))
  case_individuals <- as.character(unique(individual_vec[case_idx]))

  res1 <- compute_test_statistic(input_obj = eSVD_obj$fit_First$posterior_mean_mat,
                                 posterior_var_mat = eSVD_obj$fit_First$posterior_var_mat,
                                 case_individuals = case_individuals,
                                 control_individuals = control_individuals,
                                 covariate_individual = "individual",
                                 individual_vec = individual_vec,
                                 metadata = metadata)

  # prep for compute_test_statistic
  eSVD_obj2 <- eSVD_obj
  eSVD_obj2$teststat_vec <- NULL
  eSVD_obj2$fit_First$nuisance_vec <- NULL
  eSVD_obj2$fit_First$posterior_mean_mat <- NULL
  eSVD_obj2$fit_First$posterior_var_mat <- NULL
  eSVD_obj2$fit_First$nuisance_vec <- NULL

  res2 <- compute_test_per_gene(input_obj = eSVD_obj,
                                alpha_max = 2 * max(dat),
                                bool_adjust_covariates = FALSE,
                                bool_covariates_as_library = TRUE,
                                bool_stabilize_underdispersion = TRUE,
                                library_min = 0.1,
                                nuisance_lower_quantile = 0.01,
                                pseudocount = 0,
                                verbose = 0 )

  # The fixture's posterior was computed with `alpha_max = 2 * max(dat)`
  # and `library_min = 0.1` (tests/assets/data_generation.R); with the
  # same settings here the two routes must agree gene by gene, not just in
  # the sum of their differences.
  expect_equal(unname(res2$teststat_vec), unname(res1$teststat_vec),
               tolerance = 1e-6)
  expect_equal(unname(res2$case_mean), unname(res1$case_mean),
               tolerance = 1e-6)
  expect_equal(unname(res2$control_mean), unname(res1$control_mean),
               tolerance = 1e-6)
})
