# context("Test Poisson")
#
# ## .evaluate_objective.poisson is correct
#
# test_that(".evaluate_objective.poisson works and is correct", {
#   set.seed(123)
#   n <- 100
#   p <- 150
#   k <- 5
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- log(tcrossprod(x_mat, y_mat))
#   library_size_vec <- 1:n
#
#   ## Simulate data
#   dat <- generate_data(nat_mat, family = "poisson", library_size_vec = library_size_vec)
#
#   res1 <- .evaluate_objective.poisson(dat, x_mat, y_mat, library_size_vec = rep(1, n))
#   res2 <- .evaluate_objective.poisson(dat, x_mat, y_mat, library_size_vec = library_size_vec)
#
#   expect_true(is.numeric(res1))
#   expect_true(is.numeric(res2))
#
#   res1b <- 0
#   for(i in 1:n) {
#     for(j in 1:p) {
#       res1b <- res1b + (exp(c(x_mat[i, ] %*% y_mat[j, ])) -
#                         dat[i, j] * c(x_mat[i, ] %*% y_mat[j, ])) / (n * p)
#     }
#   }
#
#   res2b <- 0
#   for(i in 1:n) {
#     for(j in 1:p) {
#       res2b <- res2b + (exp(c(x_mat[i, ] %*% y_mat[j, ]) + log(library_size_vec[i])) -
#                         dat[i, j] * c(x_mat[i, ] %*% y_mat[j, ])) / (n * p)
#     }
#   }
#
#   expect_true(abs(res1 - res1b) <= 1e-4)
#   expect_true(abs(res2 - res2b) <= 1e-4)
# })
#
# #################
#
# ## .evaluate_objective_single.poisson is correct
#
# test_that(".evaluate_objective_single.poisson works and is correct in the x direction", {
#   set.seed(123)
#   n <- 100
#   p <- 150
#   k <- 5
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- log(tcrossprod(x_mat, y_mat))
#   library_size_vec <- 1:n
#
#   ## Simulate data
#   dat <- generate_data(nat_mat, family = "poisson", library_size_vec = library_size_vec)
#
#   i <- 1
#   res1 <- .evaluate_objective_single.poisson(x_mat[i, ], y_mat, dat[i, ], library_size_vec = rep(1, n)[i])
#   res2 <- .evaluate_objective_single.poisson(x_mat[i, ], y_mat, dat[i, ], library_size_vec = library_size_vec[i])
#
#   expect_true(is.numeric(res1))
#   expect_true(is.numeric(res2))
#
#   res1b <- 0
#   for(j in 1:p) {
#     res1b <- res1b + (exp(c(x_mat[i, ] %*% y_mat[j,])) -
#                       dat[i, j] * c(x_mat[i, ] %*% y_mat[j, ])) / p
#   }
#   res2b <- 0
#   for(j in 1:p) {
#     res2b <- res2b + (exp(c(x_mat[i, ] %*% y_mat[j, ]) + log(library_size_vec[i])) -
#                       dat[i, j] * c(x_mat[i, ] %*% y_mat[j, ])) / p
#   }
#
#   expect_true(abs(res1 - res1b) <= 1e-4)
#   expect_true(abs(res2 - res2b) <= 1e-4)
# })
#
# test_that(".evaluate_objective_single.poisson works and is correct in the y direction", {
#   set.seed(123)
#   n <- 100
#   p <- 150
#   k <- 5
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- log(tcrossprod(x_mat, y_mat))
#   library_size_vec <- 1:n
#
#   ## Simulate data
#   dat <- generate_data(nat_mat, family = "poisson", library_size_vec = library_size_vec)
#
#   j <- 1
#   res1 <- .evaluate_objective_single.poisson(y_mat[j, ], x_mat, dat[, j], library_size_vec = rep(1, n))
#   res2 <- .evaluate_objective_single.poisson(y_mat[j, ], x_mat, dat[, j], library_size_vec = library_size_vec)
#
#   expect_true(is.numeric(res1))
#   expect_true(is.numeric(res2))
#
#   res1b <- 0
#   for(i in 1:n) {
#     res1b <- res1b + (exp(c(x_mat[i, ] %*% y_mat[j, ])) -
#                       dat[i, j] * c(x_mat[i, ] %*% y_mat[j, ])) / n
#   }
#   res2b <- 0
#   for(i in 1:n) {
#     res2b <- res2b + (exp(c(x_mat[i, ] %*% y_mat[j, ]) + log(library_size_vec[i])) -
#                       dat[i, j] * c(x_mat[i, ] %*% y_mat[j, ])) / n
#   }
#
#   expect_true(abs(res1 - res1b) <= 1e-4)
#   expect_true(abs(res2 - res2b) <= 1e-4)
# })
#
# ###################################
#
# ## .gradient_vec.poisson is correct
#
# test_that(".gradient_vec.poisson works and is correct in the x direction", {
#   set.seed(123)
#   n <- 100
#   p <- 150
#   k <- 5
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- log(tcrossprod(x_mat, y_mat))
#   library_size_vec <- 1:n
#
#   ## Simulate data
#   dat <- generate_data(nat_mat, family = "poisson", library_size_vec = library_size_vec)
#
#   i <- 1
#   res1 <- .gradient_vec.poisson(x_mat[i, ], y_mat, dat[i, ], library_size_vec = rep(1, n)[i])
#   res2 <- .gradient_vec.poisson(x_mat[i, ], y_mat, dat[i, ], library_size_vec = library_size_vec[i])
#
#   expect_true(is.numeric(res1))
#   expect_true(is.numeric(res2))
#
#   res1b <- rep(0, k)
#   for(j in 1:p) {
#     res1b <- res1b + y_mat[j, ] * (exp(c(x_mat[i, ] %*% y_mat[j, ])) - dat[i, j]) / p
#   }
#   res2b <- rep(0, k)
#   for(j in 1:p) {
#     res2b <- res2b + y_mat[j, ] * (exp(c(x_mat[i, ] %*% y_mat[j, ]) + log(library_size_vec[i])) - dat[i, j]) / p
#   }
#
#   expect_true(sum(abs(res1 - res1b)) <= 1e-4)
#   expect_true(sum(abs(res2 - res2b)) <= 1e-4)
# })
#
# test_that(".gradient_vec.poisson works and is correct in the y direction", {
#   set.seed(123)
#   n <- 100
#   p <- 150
#   k <- 5
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- log(tcrossprod(x_mat, y_mat))
#   library_size_vec <- 1:n
#
#   ## Simulate data
#   dat <- generate_data(nat_mat, family = "poisson", library_size_vec = library_size_vec)
#
#   j <- 1
#   res1 <- .gradient_vec.poisson(y_mat[j, ], x_mat, dat[, j], library_size_vec = rep(1, n))
#   res2 <- .gradient_vec.poisson(y_mat[j, ], x_mat, dat[, j], library_size_vec = library_size_vec)
#
#   expect_true(is.numeric(res1))
#   expect_true(is.numeric(res2))
#
#   res1b <- rep(0, k)
#   for(i in 1:n) {
#     res1b <- res1b + x_mat[i, ] * (exp(c(x_mat[i, ] %*% y_mat[j, ])) - dat[i, j]) / n
#   }
#   res2b <- rep(0, k)
#   for(i in 1:n) {
#     res2b <- res2b + x_mat[i, ] * (exp(c(x_mat[i, ] %*% y_mat[j, ]) + log(library_size_vec[i])) - dat[i, j]) / n
#   }
#
#   expect_true(sum(abs(res1 - res1b)) <= 1e-4)
#   expect_true(sum(abs(res2 - res2b)) <= 1e-4)
# })
#
# ########################################
#
# ## .hessian_vec.poisson is correct
#
# test_that(".hessian_vec.poisson works and is correct in the x direction", {
#   set.seed(123)
#   n <- 100
#   p <- 150
#   k <- 5
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- log(tcrossprod(x_mat, y_mat))
#   library_size_vec <- 1:n
#
#   ## Simulate data
#   dat <- generate_data(nat_mat, family = "poisson", library_size_vec = library_size_vec)
#
#   i <- 1
#   res1 <- .hessian_vec.poisson(x_mat[i, ], y_mat, dat[i, ], library_size_vec = rep(1, n)[i])
#   res2 <- .hessian_vec.poisson(x_mat[i, ], y_mat, dat[i, ], library_size_vec = library_size_vec[i])
#
#   expect_true(is.numeric(res1))
#   expect_true(is.numeric(res2))
#
#   res1b <- matrix(0, k, k)
#   for(j in 1:p) {
#     res1b <- res1b + (exp(c(x_mat[i, ] %*% y_mat[j, ]))) * tcrossprod(y_mat[j, ]) / p
#   }
#   res2b <- matrix(0, k, k)
#   for(j in 1:p) {
#     res2b <- res2b + (exp(c(x_mat[i, ] %*% y_mat[j, ]) + log(library_size_vec[i]))) * tcrossprod(y_mat[j, ]) / p
#   }
#
#   expect_true(sum(abs(res1 - res1b)) <= 1e-4)
#   expect_true(sum(abs(res2 - res2b)) <= 1e-4)
# })
#
# test_that(".hessian_vec.poisson works and is correct in the y direction", {
#   set.seed(123)
#   n <- 100
#   p <- 150
#   k <- 5
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- log(tcrossprod(x_mat, y_mat))
#   library_size_vec <- 1:n
#
#   ## Simulate data
#   dat <- generate_data(nat_mat, family = "poisson", library_size_vec = library_size_vec)
#
#   j <- 1
#   res1 <- .hessian_vec.poisson(y_mat[j, ], x_mat, dat[, j], library_size_vec = rep(1, n))
#   res2 <- .hessian_vec.poisson(y_mat[j, ], x_mat, dat[, j], library_size_vec = library_size_vec)
#
#   expect_true(is.numeric(res1))
#   expect_true(is.numeric(res2))
#
#   res1b <- matrix(0, k, k)
#   for(i in 1:n) {
#     res1b <- res1b + (exp(c(x_mat[i, ] %*% y_mat[j, ]))) * tcrossprod(x_mat[i, ]) / n
#   }
#   res2b <- matrix(0, k, k)
#   for(i in 1:n) {
#     res2b <- res2b + (exp(c(x_mat[i, ] %*% y_mat[j, ]) + log(library_size_vec[i]))) * tcrossprod(x_mat[i, ]) / n
#   }
#
#   expect_true(sum(abs(res1 - res1b)) <= 1e-4)
#   expect_true(sum(abs(res2 - res2b)) <= 1e-4)
# })
