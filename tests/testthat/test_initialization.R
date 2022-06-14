context("Test initialization")

## initialize_esvd is correct

test_that("initialize_esvd works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  res <- initialize_esvd(dat = dat,
                         covariates = covariates,
                         case_control_variable = "case_control",
                         k = 5,
                         lambda = 0.1,
                         mixed_effect_variables = colnames(covariates)[grep("individual", colnames(covariates))],
                         offset_variables = NULL,
                         verbose = 0)

  expect_true(inherits(res, "eSVD"))
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("dat", "covariates",
                                             "initial_Reg", "param"))))
  expect_true(sum(abs(res$dat - dat)) <= 1e-4)
  expect_true(sum(abs(res$covariates - cbind(1, covariates))) <= 1e-4)
  expect_true(all(sort(names(res$initial_Reg)) == sort(c("log_pval", "z_mat1", "z_mat2"))))
  expect_true(all(res$initial_Reg$log_pval <= 0))
  expect_true(length(names(res$initial_Reg$log_pval)) > 0)
  expect_true(all(names(res$initial_Reg$log_pval) == colnames(dat)))
  expect_true(mean(res$initial_Reg$log_pval[true_cc_status == 1]) >= mean(res$initial_Reg$log_pval[true_cc_status == 2]))
})

###########################

## apply_initial_threshold is correct

test_that("apply_initial_threshold works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  eSVD_obj$fit_Init <- NULL
  eSVD_obj$latest_Fit <- NULL
  eSVD_obj$fit_First <- NULL
  eSVD_obj$teststat_vec <- NULL
  eSVD_obj$param$init_pval_thres <- NULL

  n <- nrow(dat); p <- ncol(dat); k <- eSVD_obj$param$init_k
  res <- apply_initial_threshold(eSVD_obj = eSVD_obj,
                                 pval_thres = 0.1)

  expect_true(inherits(res, "eSVD"))
  expect_true(is.list(res))
  expect_true("init_pval_thres" %in% names(res$param))
  expect_true(res$param$init_pval_thres == 0.1)
  expect_true(all(sort(names(res)) == sort(c("dat", "covariates",
                                             "initial_Reg", "param", "fit_Init",
                                             "latest_Fit"))))
  expect_true(all(sort(names(res$fit_Init)) == sort(c("x_mat", "y_mat",
                                                      "z_mat"))))
  expect_true(all(dim(res$fit_Init$x_mat) == c(n,k)))
  expect_true(all(dim(res$fit_Init$y_mat) == c(p,k)))
  expect_true(all(dim(res$fit_Init$z_mat) == c(p,ncol(eSVD_obj$covariates))))
  expect_true(length(colnames(res$fit_Init$z_mat)) > 0)
  expect_true(all(colnames(res$fit_Init$z_mat) == c("Intercept", colnames(covariates))))
})

########################

## .initialize_coefficient is correct

test_that(".initialize_coefficient works for sparse matrices", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  res <- .initialize_coefficient(bool_intercept = T,
                                 case_control_variable = "case_control",
                                 covariates = covariates,
                                 dat = dat,
                                 lambda = 0.1,
                                 mixed_effect_variables = colnames(covariates)[grep("individual", colnames(covariates))],
                                 offset_variables = NULL)

  expect_true(inherits(res, "initial_Reg"))
  expect_true(all(res$log_pval <= 0))
  expect_true(all(sort(names(res)) == sort(c("log_pval", "z_mat1", "z_mat2"))))
  expect_true(mean(res$log_pval[true_cc_status == 1]) >= mean(res$log_pval[true_cc_status == 2]))
})
