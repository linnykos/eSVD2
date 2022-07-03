context("Test optimization")

test_that("opt_esvd works for eSVD_obj", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  eSVD_obj$fit_First <- NULL
  eSVD_obj$teststat_vec <- NULL
  eSVD_obj$latest_Fit <- "fit_Init"

  expect_true(inherits(eSVD_obj$dat, "dgCMatrix"))
  res <- opt_esvd(input_obj = eSVD_obj,
                  max_iter = 5)

  expect_true(is.list(res))
  expect_true(inherits(res, "eSVD"))
  expect_true(all(sort(names(res)) == sort(c("dat", "covariates", "param", "fit_Init",
                                             "fit_First", "latest_Fit"))))
  expect_true(inherits(res$fit_First, "eSVD_Fit"))
  expect_true(all(sort(names(res$fit_First)) == sort(c("x_mat", "y_mat",
                                                       "z_mat", "loss"))))
  expect_true(all(diff(res$fit_First$loss) < 0))
  expect_true(res$latest_Fit == "fit_First")


  ## ensures it works for dense matrices
  eSVD_obj$dat <- as.matrix(eSVD_obj$dat)
  res <- opt_esvd(input_obj = eSVD_obj,
                  max_iter = 5)

  expect_true(is.list(res))
  expect_true(inherits(res, "eSVD"))
})
