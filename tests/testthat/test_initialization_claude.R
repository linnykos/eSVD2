# Tests for initialize_esvd() beyond the legacy test_initialization.R.

## The single-covariate branch of `.initialize_coefficient()` (reached when
## only zero or one covariate is left to estimate after removing the intercept
## and the offsets) cannot use glmnet, which needs two predictors. It used to
## fall back to `stats::glm()` WITHOUT the offset, so a design with only the
## case-control indicator initialized every gene's intercept without adjusting
## for sequencing depth. The oracle is `stats::glm()` with the offset.
test_that("T-INIT-09: the GLM fallback keeps the library-size offset", {
  set.seed(10)
  n <- 80
  p <- 3
  dat <- matrix(stats::rpois(n * p, lambda = 4), nrow = n, ncol = p)
  dimnames(dat) <- list(paste0("c", seq_len(n)), paste0("g", seq_len(p)))
  covariate_df <- data.frame(CC = factor(rep(c(0, 1), each = n / 2)))
  covariates <- format_covariates(dat = dat, covariate_df = covariate_df)
  expect_equal(colnames(covariates), c("Intercept", "Log_UMI", "CC_1"))
  offset_vec <- covariates[, "Log_UMI"]

  # One covariate to estimate, with an intercept.
  z_mat <- .initialize_coefficient(bool_intercept = TRUE,
                                   covariates = covariates,
                                   dat = dat,
                                   lambda = 0.1,
                                   offset_variables = "Log_UMI",
                                   verbose = 0)
  expect_true(all(z_mat[, "Log_UMI"] == 1))
  for(j in seq_len(p)){
    glm_fit <- stats::glm(dat[, j] ~ covariates[, "CC_1"] + offset(offset_vec),
                          family = stats::poisson)
    expect_equal(unname(z_mat[j, c("Intercept", "CC_1")]),
                 unname(stats::coef(glm_fit)),
                 tolerance = 1e-6, info = paste0("gene ", j))
  }
  # With the offset in place the intercept is on the per-unit-depth scale,
  # far below log(mean count); without it the two would nearly coincide.
  expect_true(all(z_mat[, "Intercept"] < log(mean(dat)) - 1))

  # One covariate, no intercept.
  z_mat_noint <- .initialize_coefficient(bool_intercept = FALSE,
                                         covariates = covariates,
                                         dat = dat,
                                         lambda = 0.1,
                                         offset_variables = "Log_UMI",
                                         verbose = 0)
  expect_true(all(z_mat_noint[, "Intercept"] == 0))
  glm_fit <- stats::glm(dat[, 1] ~ 0 + covariates[, "CC_1"] + offset(offset_vec),
                        family = stats::poisson)
  expect_equal(unname(z_mat_noint[1, "CC_1"]), unname(stats::coef(glm_fit)),
               tolerance = 1e-6)

  # Zero covariates to estimate: intercept only.
  covariates0 <- covariates[, c("Intercept", "Log_UMI")]
  z_mat0 <- .initialize_coefficient(bool_intercept = TRUE,
                                    covariates = covariates0,
                                    dat = dat,
                                    lambda = 0.1,
                                    offset_variables = "Log_UMI",
                                    verbose = 0)
  expect_equal(dim(z_mat0), c(p, 2))
  for(j in seq_len(p)){
    glm_fit <- stats::glm(dat[, j] ~ 1 + offset(offset_vec),
                          family = stats::poisson)
    expect_equal(unname(z_mat0[j, "Intercept"]), unname(stats::coef(glm_fit)),
                 tolerance = 1e-6, info = paste0("gene ", j))
  }

  # Zero covariates and no intercept: nothing to estimate, nothing to fit.
  z_mat00 <- .initialize_coefficient(bool_intercept = FALSE,
                                     covariates = covariates0,
                                     dat = dat,
                                     lambda = 0.1,
                                     offset_variables = "Log_UMI",
                                     verbose = 0)
  expect_true(all(z_mat00[, "Intercept"] == 0))
  expect_true(all(z_mat00[, "Log_UMI"] == 1))
})

## The full initializer on the same minimal design runs end to end, and the
## sparse and dense inputs agree.
test_that("T-INIT-10: initialize_esvd runs with only the case-control covariate", {
  set.seed(10)
  n <- 80
  p <- 6
  dat <- matrix(stats::rpois(n * p, lambda = 4), nrow = n, ncol = p)
  dimnames(dat) <- list(paste0("c", seq_len(n)), paste0("g", seq_len(p)))
  individual_vec <- factor(rep(paste0("i", 1:8), each = n / 8))
  covariate_df <- data.frame(CC = factor(rep(c(0, 1), each = n / 2)))
  covariates <- format_covariates(dat = dat, covariate_df = covariate_df)

  res_dense <- initialize_esvd(dat = dat,
                               covariates = covariates,
                               metadata_individual = individual_vec,
                               bool_intercept = TRUE,
                               case_control_variable = "CC_1",
                               k = 2,
                               lambda = 0.1)
  dat_sparse <- methods::as(methods::as(dat * 1.0, "dMatrix"), "CsparseMatrix")
  res_sparse <- initialize_esvd(dat = dat_sparse,
                                covariates = covariates,
                                metadata_individual = individual_vec,
                                bool_intercept = TRUE,
                                case_control_variable = "CC_1",
                                k = 2,
                                lambda = 0.1)

  expect_true(inherits(res_dense, "eSVD"))
  expect_equal(res_dense$fit_Init$z_mat, res_sparse$fit_Init$z_mat,
               tolerance = 1e-8)
  expect_true(all(is.finite(res_dense$fit_Init$x_mat)))
})
