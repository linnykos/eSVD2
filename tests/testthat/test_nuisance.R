context("Test nuisance")

test_that("estimate_nuisance works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  eSVD_obj$teststat_vec <- NULL
  eSVD_obj$fit_First$nuisance_vec <- NULL

  res <- estimate_nuisance(input_obj = eSVD_obj,
                           verbose = 0)
  # plot(res$fit_First$nuisance_vec)

  expect_true("nuisance_vec" %in% names(res$fit_First))
  expect_true(all(res$fit_First$nuisance_vec > 0))
})
