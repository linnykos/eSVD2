context("Test optimization helper functions")

test_that("opt_x and opt_yz works", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- log(tcrossprod(x_mat, y_mat))

  # Simulate data
  dat <- generate_data(nat_mat,
                       family = "poisson",
                       nuisance_param_vec = NA,
                       library_size_vec = 1)

  # Simulate a covariate whose coefficients are fixed, i.e., an offset
  covariate1 <- rnorm(n)
  coef1 <- rep(1.23, p)
  # Simulate another covariate whose coefficients are to be estimated
  covariate2 <- rnorm(n)
  coef2 <- rep(0.1, p)
  # Combine X and covariates
  c_mat <- cbind(covariate1, covariate2)
  xc_mat <- cbind(x_mat, c_mat)
  # Combine Y and Z
  z_mat <- cbind(coef1, coef2)
  yz_mat <- cbind(y_mat, z_mat)
  # fixed_cols stands for the columns in YZ that need to be fixed
  # In this case the first column of Z, (k+1)-th column of YZ
  fixed_cols <- k + 1

  loader <- data_loader(dat)
  family <- esvd_family("poisson")
  library_size_vec <- rep(1, n)
  nuisance_param_vec <- rep(NA, p)
  l2penx <- 0.001
  l2peny <- 0.01
  l2penz <- 0.1

  # Compute initial loss
  loss0 <- objfn_all_r(
    XC = x_mat, YZ = y_mat, k = k, loader = loader, family = family,
    s = library_size_vec, gamma = nuisance_param_vec,
    l2penx = l2penx, l2peny = l2peny, l2penz = l2penz
  )

  # Optimize X
  cat("\n")
  newxc <- opt_x(xc_mat, yz_mat, k, loader = loader, family = family,
                 s = library_size_vec, gamma = nuisance_param_vec,
                 l2penx = l2penx, verbose = 0)

  # Test whether C has remained the same
  expect_equal(max(abs(xc_mat[, -(1:k)] - newxc[, -(1:k)])), 0)

  # Compute new loss
  loss1 <- objfn_all_r(
    XC = newxc, YZ = yz_mat, k = k, loader = loader, family = family,
    s = library_size_vec, gamma = nuisance_param_vec,
    l2penx = l2penx, l2peny = l2peny, l2penz = l2penz
  )

  # Optimize Y
  newyz <- opt_yz(yz_mat, newxc, k, fixed_cols = fixed_cols,
                  loader = loader, family = family,
                  s = library_size_vec, gamma = nuisance_param_vec,
                  l2peny = l2peny, l2penz = l2penz,
                  verbose = 0)

  # Test whether the first column of Z has remained the same
  expect_equal(max(abs(yz_mat[, k + 1] - newyz[, k + 1])), 0)

  # Compute new loss
  loss2 <- objfn_all_r(
    XC = newxc, YZ = newyz, k = k, loader = loader, family = family,
    s = library_size_vec, gamma = nuisance_param_vec,
    l2penx = l2penx, l2peny = l2peny, l2penz = l2penz
  )

  # cat("loss0 = ", loss0, ", loss1 = ", loss1, ", loss2 = ", loss2, "\n", sep = "")
  expect_lt(loss1, loss0)
  expect_lt(loss2, loss1)
})
