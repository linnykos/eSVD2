context("Test posterior")

## compute_posterior is correct

test_that("compute_posterior works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  eSVD_obj$fit_First$posterior_mean_mat <- NULL
  eSVD_obj$fit_First$posterior_var_mat <- NULL

  res <- compute_posterior(input_obj = eSVD_obj)

  expect_true(all(c("posterior_mean_mat", "posterior_var_mat") %in% names(res$fit_First)))
  expect_true(all(dim(res$fit_First$posterior_mean_mat) == dim(dat)))
  expect_true(all(dim(res$fit_First$posterior_var_mat) == dim(dat)))
})
