# Distribution: exponential
# Natural parameter: m_{ij} = u_i^Tv_j
# Relation to canonical parameters: m_{ij} = -lambda_{ij}, where E(a_{ij}) = 1/lambda_{ij}
# Optimization problem: -log(-m_{ij}) - a_{ij}*m_{ij}

.evaluate_objective.exponential <- function(
    dat, u_mat, v_mat, nuisance_param_vec, library_size_vec, ...
) {
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
    dat_vals <- (dat / library_size_vec)[idx]
    negloglik <- -log(-nat_vals) - nat_vals * dat_vals
    sum(negloglik) / n / p
}

# length(library_size) == 1 if current vector is u
# length(library_size) == n if current vector is v
.evaluate_objective_single.exponential <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- (dat_vec / library_size)[idx]
    negloglik <- -log(-nat_vals) - nat_vals * dat_vals
    sum(negloglik) / length(dat_vec)
}

.gradient_vec.exponential <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- (dat_vec / library_size)[idx]
    other_vals <- other_mat[idx, , drop = FALSE]
    grad <- other_vals * (-1 / nat_vals - dat_vals)

    colSums(grad) / length(dat_vec)
}

.hessian_vec.exponential <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- (dat_vec / library_size)[idx]
    other_vals <- other_mat[idx, , drop = FALSE]

    term1 <- t(other_vals) %*% diag(1 / nat_vals^2) %*% other_vals

    term1 / length(dat_vec)
}

.feasibility.exponential <- function(current_vec, other_mat, ...) {
    nat_vec <- c(other_mat %*% current_vec)
    all(nat_vec < 0)
}

.exponential <- structure(
    list(
        objfn_all = .evaluate_objective.exponential,
        objfn     = .evaluate_objective_single.exponential,
        grad      = .gradient_vec.exponential,
        hessian   = .hessian_vec.exponential,
        feas      = .feasibility.exponential
    ),
    class = "esvd_family"
)



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
