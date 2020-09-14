# Distribution: exponential
# Natural parameter: m_{ij} = u_i^Tv_j
# Relation to canonical parameters: m_{ij} = -lambda_{ij}, where E(a_{ij}) = 1/lambda_{ij}
# optimization problem: -log(-m_{ij}) - a_{ij}*m_{ij}

.evaluate_objective.exponential <- function(dat, u_mat, v_mat, ...) {
    # Check dimensions
    n <- nrow(dat)
    p <- ncol(dat)
    stopifnot(
        ncol(u_mat) == ncol(v_mat),
        nrow(u_mat) == n,
        nrow(v_mat) == p
    )

    # Compute natural parameters
    nat_mat <- tcrossprod(u_mat, v_mat)
    stopifnot(all(nat_mat < 0))

    # Only compute likelihood on non-missing data
    idx <- which(!is.na(dat))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_mat[idx]
    dat_vals <- dat[idx]
    negloglik <- -log(-nat_vals) - nat_vals * dat_vals
    sum(negloglik) / n / p
}

.evaluate_objective_single.exponential <- function(dat_vec, current_vec, other_mat, ...) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- dat_vec[idx]
    negloglik <- -log(-nat_vals) - nat_vals * dat_vals
    sum(negloglik) / length(dat_vec)
}

.gradient_vec.exponential <- function(dat_vec, current_vec, other_mat, ...) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- dat_vec[idx]
    other_vals <- other_mat[idx, , drop = FALSE]
    grad <- other_vals * (-1 / nat_vals - dat_vals)

    colSums(grad) / length(dat_vec)
}

.hessian_vec.exponential <- function(dat_vec, current_vec, other_mat, ...) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- dat_vec[idx]
    other_vals <- other_mat[idx, , drop = FALSE]

    term1 <- t(other_vals) %*% diag(1 / nat_vals^2) %*% other_vals

    term1 / length(dat_vec)
}

.feasibility.exponential <- function(current_vec, other_mat, ...) {
    nat_vec <- c(other_mat %*% current_vec)
    all(nat_vec < 0)
}



### Test correctness ###
# tests = function(dat, u_mat, v_mat)
# {
#     n = nrow(dat)
#     p = ncol(dat)
#
#     # Test objective functions
#     loss1 = .evaluate_objective.exponential(dat, u_mat, v_mat)
#     loss2 = numeric(n)
#     for(i in 1:n)
#     {
#         loss2[i] = .evaluate_objective_single.exponential(
#             dat[i, ], u_mat[i, ], v_mat
#         )
#     }
#     loss3 = numeric(p)
#     for(j in 1:p)
#     {
#         loss3[j] = .evaluate_objective_single.exponential(
#             dat[, j], v_mat[j, ], u_mat
#         )
#     }
#     stopifnot(
#         abs(loss1 - mean(loss2)) < 1e-8,
#         abs(loss1 - mean(loss3)) < 1e-8
#     )
#
#     # Test gradients
#     for(i in 1:n)
#     {
#         grad1 = .gradient_vec.exponential(
#             dat[i, ], u_mat[i, ], v_mat
#         )
#         grad2 = numDeriv::grad(function(x) {
#             .evaluate_objective_single.exponential(dat[i, ], x, v_mat)
#         }, u_mat[i, ])
#         stopifnot(max(abs(grad1 - grad2)) < 1e-6)
#     }
#     for(j in 1:p)
#     {
#         grad1 = .gradient_vec.exponential(
#             dat[, j], v_mat[j, ], u_mat
#         )
#         grad2 = numDeriv::grad(function(x) {
#             .evaluate_objective_single.exponential(dat[, j], x, u_mat)
#         }, v_mat[j, ])
#         stopifnot(max(abs(grad1 - grad2)) < 1e-6)
#     }
#
#     # Test Hessians
#     for(i in 1:n)
#     {
#         hess1 = .hessian_vec.exponential(
#             dat[i, ], u_mat[i, ], v_mat
#         )
#         hess2 = numDeriv::hessian(function(x) {
#             .evaluate_objective_single.exponential(dat[i, ], x, v_mat)
#         }, u_mat[i, ])
#         stopifnot(max(abs(hess1 - hess2)) < 1e-6)
#     }
#     for(j in 1:p)
#     {
#         hess1 = .hessian_vec.exponential(
#             dat[, j], v_mat[j, ], u_mat
#         )
#         hess2 = numDeriv::hessian(function(x) {
#             .evaluate_objective_single.exponential(dat[, j], x, u_mat)
#         }, v_mat[j, ])
#         stopifnot(max(abs(hess1 - hess2)) < 1e-6)
#     }
# }
# # Simulate data
# set.seed(123)
# n = 10
# p = 15
# k = 2
# u_mat = matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
# v_mat = -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
# nat_mat = tcrossprod(u_mat, v_mat)
# dat = eSVD2::generate_data(
#     nat_mat, family = "exponential", nuisance_param_vec = NA, library_size_vec = NA
# )
#
# # Test
# tests(dat, u_mat, v_mat)
#
# # Test missing values
# dat[sample(length(dat), n * p * 0.1)] = NA
# tests(dat, u_mat, v_mat)



.evaluate_objective_mat.exponential <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat < 0))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum(-log(-nat_mat[idx]) - nat_mat[idx]*dat[idx])
}

.gradient_mat.exponential <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat < 0))

  n <- nrow(dat); p <- ncol(dat)

  (-1/(nat_mat) - dat)/(n*p)
}



