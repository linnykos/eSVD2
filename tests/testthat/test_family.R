context("Test distribution family functions")

run_test <- function(dat, u_mat, v_mat, family, nuisance_param_vec = rep(1, nrow(v_mat)), library_size_vec = rep(1, nrow(u_mat)))
{
    n <- nrow(dat)
    p <- ncol(dat)
    library_size_vec <- .parse_library_size(dat, library_size_vec)

    # Test objective functions
    loss1 <- family$objfn_all(dat, u_mat, v_mat, nuisance_param_vec, library_size_vec)
    loss2 <- numeric(n)
    for(i in 1:n)
    {
        loss2[i] <- family$objfn(u_mat[i, ], v_mat, dat[i, ], nuisance_param_vec,
                                 library_size = library_size_vec[i])
    }
    loss3 <- numeric(p)
    for(j in 1:p)
    {
        loss3[j] <- family$objfn(v_mat[j, ], u_mat, dat[, j], nuisance_param_vec[j],
                                 library_size = library_size_vec)
    }
    # The three loss values should be equal
    expect_lt(loss1 - mean(loss2), 1e-8)
    expect_lt(loss1 - mean(loss3), 1e-8)

    # Test gradients
    for(i in 1:n)
    {
        grad1 <- family$grad(u_mat[i, ], v_mat, dat[i, ],
                             nuisance_param_vec = nuisance_param_vec,
                             library_size = library_size_vec[i])
        grad2 <- numDeriv::grad(family$objfn, u_mat[i, ], other_mat = v_mat, dat_vec = dat[i, ],
                                nuisance_param_vec = nuisance_param_vec,
                                library_size = library_size_vec[i])
        expect_lt(max(abs(grad1 - grad2)), 1e-6)
    }
    for(j in 1:p)
    {
        grad1 <- family$grad(v_mat[j, ], u_mat, dat[, j],
                             nuisance_param_vec = nuisance_param_vec[j],
                             library_size = library_size_vec)
        grad2 <- numDeriv::grad(family$objfn, v_mat[j, ], other_mat = u_mat, dat_vec = dat[, j],
                                nuisance_param_vec = nuisance_param_vec[j],
                                library_size = library_size_vec)
        expect_lt(max(abs(grad1 - grad2)), 1e-6)
    }

    # Test Hessians
    for(i in 1:n)
    {
        hess1 <- family$hessian(u_mat[i, ], v_mat, dat[i, ],
                                nuisance_param_vec = nuisance_param_vec,
                                library_size = library_size_vec[i])
        hess2 <- numDeriv::hessian(family$objfn, u_mat[i, ], other_mat = v_mat, dat_vec = dat[i, ],
                                   nuisance_param_vec = nuisance_param_vec,
                                   library_size = library_size_vec[i])
        expect_lt(max(abs(hess1 - hess2)), 1e-6)
    }
    for(j in 1:p)
    {
        hess1 <- family$hessian(v_mat[j, ], u_mat, dat[, j],
                                nuisance_param_vec = nuisance_param_vec[j],
                                library_size = library_size_vec)
        hess2 <- numDeriv::hessian(family$objfn, v_mat[j, ], other_mat = u_mat, dat_vec = dat[, j],
                                   nuisance_param_vec = nuisance_param_vec[j],
                                   library_size = library_size_vec)
        expect_lt(max(abs(hess1 - hess2)), 1e-6)
    }
}


######################## Gaussian ########################

test_that("Functions for Gaussian distribution", {
    # Simulate natural parameter matrix
    set.seed(123)
    n <- 10
    p <- 15
    k <- 2
    nuisance_param_vec <- rep(2, p); library_size_vec = rep(1, n)
    u_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
    v_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
    nat_mat <- tcrossprod(u_mat, v_mat)

    # Simulate data with default library size (all one)
    dat <- eSVD2::generate_data(
        nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
        library_size_vec = library_size_vec
    )

    # Test
    run_test(dat, u_mat, v_mat, .gaussian, nuisance_param_vec = nuisance_param_vec,
             library_size_vec = library_size_vec)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .gaussian, nuisance_param_vec = nuisance_param_vec,
             library_size_vec = library_size_vec)

    # Simulate data with a library size vector
    library_size_vec <- sample(10:20, n, replace = TRUE)
    dat <- eSVD2::generate_data(
        nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
        library_size_vec = library_size_vec
    )

    # Test
    run_test(dat, u_mat, v_mat, .gaussian, nuisance_param_vec = nuisance_param_vec,
             library_size_vec = library_size_vec)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .gaussian, nuisance_param_vec = nuisance_param_vec,
             library_size_vec = library_size_vec)
})

######################## Curved Gaussian ########################

test_that("Functions for curved-Gaussian distribution", {
    # Simulate natural parameter matrix
    set.seed(123)
    n <- 10
    p <- 15
    k <- 2
    scalar <- 2
    u_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
    v_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
    nat_mat <- tcrossprod(u_mat, v_mat)

    # Simulate data with default library size (all one)
    dat <- eSVD2::generate_data(
        nat_mat, family = "curved_gaussian", nuisance_param_vec = scalar,
        library_size_vec = 1, tol = 1e-3
    )

    # Test
    run_test(dat, u_mat, v_mat, .curved_gaussian, nuisance_param_vec = scalar,
             library_size_vec = 1)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .curved_gaussian, nuisance_param_vec = scalar,
             library_size_vec = 1)

    # Simulate data with a library size vector
    library_size_vec <- sample(10:20, n, replace = TRUE)
    dat <- eSVD2::generate_data(
        nat_mat, family = "curved_gaussian", nuisance_param_vec = scalar,
        library_size_vec = library_size_vec, tol = 1e-3
    )

    # Test
    run_test(dat, u_mat, v_mat, .curved_gaussian, nuisance_param_vec = scalar,
             library_size_vec = library_size_vec)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .curved_gaussian, nuisance_param_vec = scalar,
             library_size_vec = library_size_vec)
})

######################## Exponential ########################

test_that("Functions for exponential distribution", {
    # Simulate natural parameter matrix
    set.seed(123)
    n <- 10
    p <- 15
    k <- 2
    u_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
    v_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
    nat_mat <- tcrossprod(u_mat, v_mat)

    # Simulate data with default library size (all one)
    dat <- eSVD2::generate_data(
        nat_mat, family = "exponential", nuisance_param_vec = NA,
        library_size_vec = 1
    )

    # Test
    run_test(dat, u_mat, v_mat, .exponential, nuisance_param_vec = NA,
             library_size_vec = 1)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .exponential, nuisance_param_vec = NA,
             library_size_vec = 1)

    # Simulate data with a library size vector
    library_size_vec <- sample(10:20, n, replace = TRUE)
    dat <- eSVD2::generate_data(
        nat_mat, family = "exponential", nuisance_param_vec = NA,
        library_size_vec = library_size_vec
    )

    # Test
    run_test(dat, u_mat, v_mat, .exponential, nuisance_param_vec = NA,
             library_size_vec = library_size_vec)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .exponential, nuisance_param_vec = NA,
             library_size_vec = library_size_vec)
})

######################## Poisson ########################

test_that("Functions for Poisson distribution", {
    # Simulate natural parameter matrix
    set.seed(123)
    n <- 10
    p <- 15
    k <- 2
    u_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
    v_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
    nat_mat <- tcrossprod(u_mat, v_mat)

    # Simulate data with default library size (all one)
    dat <- eSVD2::generate_data(
        nat_mat, family = "poisson", nuisance_param_vec = NA,
        library_size_vec = 1
    )

    # Test
    run_test(dat, u_mat, v_mat, .poisson, nuisance_param_vec = NA,
             library_size_vec = 1)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .poisson, nuisance_param_vec = NA,
             library_size_vec = 1)

    # Simulate data with a library size vector
    library_size_vec <- sample(10:20, n, replace = TRUE)
    dat <- eSVD2::generate_data(
        nat_mat, family = "poisson", nuisance_param_vec = NA,
        library_size_vec = library_size_vec
    )

    # Test
    run_test(dat, u_mat, v_mat, .poisson, nuisance_param_vec = NA,
             library_size_vec = library_size_vec)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .poisson, nuisance_param_vec = NA,
             library_size_vec = library_size_vec)
})

######################## Negative binomial ########################

test_that("Functions for negative binomial distribution", {
    # Simulate data
    set.seed(123)
    n <- 10
    p <- 15
    k <- 2
    scalar <- 10
    u_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
    v_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
    nat_mat <- tcrossprod(u_mat, v_mat)
    dat <- eSVD2::generate_data(
        nat_mat, family = "neg_binom", nuisance_param_vec = scalar,
        library_size_vec = 1
    )

    # Test
    run_test(dat, u_mat, v_mat, .neg_binom, nuisance_param_vec = scalar,
             library_size_vec = 1)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .neg_binom, nuisance_param_vec = scalar,
             library_size_vec = 1)
})

######################## Bernoulli ########################

test_that("Functions for Bernoulli distribution", {
    # Simulate data
    set.seed(123)
    n <- 10
    p <- 15
    k <- 2
    u_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
    v_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
    nat_mat <- tcrossprod(u_mat, v_mat)
    dat <- eSVD2::generate_data(
        nat_mat, family = "bernoulli", nuisance_param_vec = NA,
        library_size_vec = 1
    )

    # Test
    run_test(dat, u_mat, v_mat, .bernoulli, nuisance_param_vec = NA,
             library_size_vec = 1)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .bernoulli, nuisance_param_vec = NA,
             library_size_vec = 1)
})
