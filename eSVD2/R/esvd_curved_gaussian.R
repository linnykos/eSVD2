# Distribution: one-parameter Gaussian where sd = mean/scalar
# Natural parameter: m_{ij} = u_i^Tv_j
# Relation to canonical parameters: m_{ij} = 1/mu_{ij}
# Optimization problem: -log(m_{ij}) - scalar^2*a_{ij}^2*(-m_{ij}^2)/2 - scalar^2*a_{ij}*m_{ij}

.evaluate_objective.curved_gaussian <- function(dat, u_mat, v_mat, scalar = 2, ...) {
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
    stopifnot(all(nat_mat > 0))

    # Only compute likelihood on non-missing data
    idx <- which(!is.na(dat))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_mat[idx]
    dat_vals <- dat[idx]
    negloglik <- -log(nat_vals) +
        scalar^2 * dat_vals^2 * nat_vals^2 / 2 -
        scalar^2 * dat_vals * nat_vals
    sum(negloglik) / n / p
}

.evaluate_objective_single.curved_gaussian <- function(dat_vec, current_vec, other_mat, scalar = 2, ...) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec > 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- dat_vec[idx]
    negloglik <- -log(nat_vals) +
        scalar^2 * dat_vals^2 * nat_vals^2 / 2 -
        scalar^2 * dat_vals * nat_vals
    sum(negloglik) / length(dat_vec)
}

.gradient_vec.curved_gaussian <- function(dat_vec, current_vec, other_mat, scalar = 2, ...) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec > 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- dat_vec[idx]
    other_vals <- other_mat[idx, , drop = FALSE]
    grad <- other_vals * (-1 / nat_vals + scalar^2 * dat_vals^2 * nat_vals -
                              scalar^2 * dat_vals)
    colSums(grad) / length(dat_vec)
}



### Test correctness ###
# tests = function(dat, u_mat, v_mat, scalar)
# {
#     n = nrow(dat)
#     p = ncol(dat)
#
#     # Test objective functions
#     loss1 = .evaluate_objective.curved_gaussian(dat, u_mat, v_mat, scalar = scalar)
#     loss2 = numeric(n)
#     for(i in 1:n)
#     {
#         loss2[i] = .evaluate_objective_single.curved_gaussian(
#             dat[i, ], u_mat[i, ], v_mat, scalar = scalar
#         )
#     }
#     loss3 = numeric(p)
#     for(j in 1:p)
#     {
#         loss3[j] = .evaluate_objective_single.curved_gaussian(
#             dat[, j], v_mat[j, ], u_mat, scalar = scalar
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
#         grad1 = .gradient_vec.curved_gaussian(
#             dat[i, ], u_mat[i, ], v_mat, scalar = scalar
#         )
#         grad2 = numDeriv::grad(function(x) {
#             .evaluate_objective_single.curved_gaussian(dat[i, ], x, v_mat, scalar = scalar)
#         }, u_mat[i, ])
#         stopifnot(max(abs(grad1 - grad2)) < 1e-6)
#     }
# }
# # Simulate data
# set.seed(123)
# n = 10
# p = 15
# k = 2
# scalar = 2
# u_mat = matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
# v_mat = matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
# nat_mat = tcrossprod(u_mat, v_mat)
# dat = eSVD2::generate_data(
#     nat_mat, family = "curved_gaussian", nuisance_param_vec = scalar,
#     library_size_vec = NA, tol = 1e-3
# )
#
# # Test
# tests(dat, u_mat, v_mat, scalar = scalar)
#
# # Test missing values
# dat[sample(length(dat), n * p * 0.1)] = NA
# tests(dat, u_mat, v_mat, scalar = scalar)



.evaluate_objective_mat.curved_gaussian <- function(dat, nat_mat, scalar = 2, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat > 0))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum(-log(nat_mat[idx]) -
        nat_mat[idx]*dat[idx]*scalar^2 +
        nat_mat[idx]^2*dat[idx]^2*scalar^2/2)
}

.gradient_mat.curved_gaussian <- function(dat, nat_mat, scalar = 2, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat > 0))

  n <- nrow(dat); p <- ncol(dat)

  (-1/(nat_mat) - scalar^2*dat + scalar^2*dat^2*nat_mat)/(n*p)
}
