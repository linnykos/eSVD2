context("Test initialization")

## initialize_esvd is correct

test_that("initialize_esvd works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  res <- initialize_esvd(dat = dat,
                         covariates = covariates,
                         metadata_individual = metadata[,"individual"],
                         case_control_variable = "case_control_1",
                         k = 5,
                         lambda = 0.1,
                         offset_variables = NULL,
                         verbose = 0)

  expect_true(inherits(res, "eSVD"))
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("dat", "covariates",
                                             "fit_Init", "latest_Fit", "param",
                                             "case_control", "individual"))))
  expect_true(sum(abs(res$dat - dat)) <= 1e-4)
  expect_true(sum(abs(res$covariates - covariates)) <= 1e-4)
})

#########

## .initialize_coefficient is correct

test_that(".initialize_coefficient works for sparse matrices", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  res <- .initialize_coefficient(bool_intercept = F,
                                 covariates = covariates,
                                 dat = dat,
                                 lambda = 0.1,
                                 offset_variables = "Log_UMI")

  expect_true(all(dim(res) == c(ncol(dat), ncol(covariates))))
  expect_true(all(res[,"Intercept"] == 0))
  expect_true(all(res[,"Log_UMI"] == 1))
})
