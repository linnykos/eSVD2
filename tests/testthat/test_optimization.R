context("Test optimization")

test_that("opt_esvd works for poisson", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- log(tcrossprod(x_mat, y_mat))

  # Simulate data
  dat <- generate_data(nat_mat, family = "poisson", nuisance_param_vec = NA,
                       library_size_vec = 1)

  # Simulate a covariate whose coefficients are fixed, i.e., an offset
  covariate1 <- rnorm(n)
  coef1 <- rep(1.23, p)
  # Simulate another covariate whose coefficients are to be estimated
  covariate2 <- rnorm(n)
  coef2 <- rep(0.1, p)
  c_mat <- cbind(covariate1, covariate2)
  z_mat <- cbind(coef1, coef2)

  res <- opt_esvd(x_init = x_mat,
                  y_init = y_mat,
                  dat = dat,
                  z_init = z_mat,
                  covariates = c_mat,
                  max_iter = 5)

  expect_true(is.list(res))
  expect_true(class(res) == "eSVD")
  expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat", "covariates", "loss",
                                             "z_mat","param"))))
})
