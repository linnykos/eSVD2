context("Test distribution family functions")

run_test <- function(dat, u_mat, v_mat, family, ...)
{
    n <- nrow(dat)
    p <- ncol(dat)

    # Test objective functions
    loss1 <- family$objfn_all(dat, u_mat, v_mat, ...)
    loss2 <- numeric(n)
    for(i in 1:n)
    {
        loss2[i] <- family$objfn(u_mat[i, ], v_mat, dat[i, ], ...)
    }
    loss3 <- numeric(p)
    for(j in 1:p)
    {
        loss3[j] <- family$objfn(v_mat[j, ], u_mat, dat[, j], ...)
    }
    # The three loss values should be equal
    expect_lt(loss1 - mean(loss2), 1e-8)
    expect_lt(loss1 - mean(loss3), 1e-8)

    # Test gradients
    for(i in 1:n)
    {
        grad1 <- family$grad(u_mat[i, ], v_mat, dat[i, ], ...)
        grad2 <- numDeriv::grad(family$objfn, u_mat[i, ], other_mat = v_mat, dat_vec = dat[i, ], ...)
        expect_lt(max(abs(grad1 - grad2)), 1e-6)
    }
    for(j in 1:p)
    {
        grad1 <- family$grad(v_mat[j, ], u_mat, dat[, j], ...)
        grad2 <- numDeriv::grad(family$objfn, v_mat[j, ], other_mat = u_mat, dat_vec = dat[, j], ...)
        expect_lt(max(abs(grad1 - grad2)), 1e-6)
    }

    # Test Hessians
    for(i in 1:n)
    {
        hess1 <- family$hessian(u_mat[i, ], v_mat, dat[i, ], ...)
        hess2 <- numDeriv::hessian(family$objfn, u_mat[i, ], other_mat = v_mat, dat_vec = dat[i, ], ...)
        expect_lt(max(abs(hess1 - hess2)), 1e-6)
    }
    for(j in 1:p)
    {
        hess1 <- family$hessian(v_mat[j, ], u_mat, dat[, j], ...)
        hess2 <- numDeriv::hessian(family$objfn, v_mat[j, ], other_mat = u_mat, dat_vec = dat[, j], ...)
        expect_lt(max(abs(hess1 - hess2)), 1e-6)
    }
}


######################## Curved Gaussian ########################

test_that("Functions for curved-Gaussian distribution", {
    # Simulate data
    set.seed(123)
    n <- 10
    p <- 15
    k <- 2
    scalar <- 2
    u_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
    v_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
    nat_mat <- tcrossprod(u_mat, v_mat)
    dat <- eSVD2::generate_data(
        nat_mat, family = "curved_gaussian", nuisance_param_vec = scalar,
        library_size_vec = NA, tol = 1e-3
    )

    # Test
    run_test(dat, u_mat, v_mat, .curved_gaussian, scalar = 2)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .curved_gaussian, scalar = 2)
})

######################## Exponential ########################

test_that("Functions for exponential distribution", {
    # Simulate data
    set.seed(123)
    n <- 10
    p <- 15
    k <- 2
    u_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
    v_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
    nat_mat <- tcrossprod(u_mat, v_mat)
    dat <- eSVD2::generate_data(
        nat_mat, family = "exponential", nuisance_param_vec = NA, library_size_vec = NA
    )

    # Test
    run_test(dat, u_mat, v_mat, .exponential)

    # Test missing values
    dat[sample(length(dat), n * p * 0.1)] <- NA
    run_test(dat, u_mat, v_mat, .exponential)
})
