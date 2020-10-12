# Distribution: Gaussian
# Natural parameter: m_{ij} = u_i^Tv_j
# Relation to canonical parameters: m_{ij} = mu_{ij}
# Optimization problem: (m_{ij} - a_{ij})^2

.evaluate_objective.gaussian <- function(dat, u_mat, v_mat, ...) {
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

    # Only compute likelihood on non-missing data
    idx <- which(!is.na(dat))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_mat[idx]
    dat_vals <- dat[idx]
    negloglik <- (nat_vals - dat_vals)^2
    sum(negloglik) /  n / p
}

.evaluate_objective_single.gaussian <- function(current_vec, other_mat, dat_vec, ...) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- dat_vec[idx]
    negloglik <- sum((nat_vals - dat_vals)^2)
    sum(negloglik) / length(dat_vec)
}

.gradient_vec.gaussian <- function(current_vec, other_mat, dat_vec, ...) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- dat_vec[idx]
    other_vals <- other_mat[idx, , drop = FALSE]
    # Broadcast (nat_vals - dat_vals) to each column of other_vals
    grad <- 2 * other_vals * (nat_vals - dat_vals)
    colSums(grad) / length(dat_vec)
}

.hessian_vec.gaussian <- function(current_vec, other_mat, dat_vec, ...) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- nat_vec[idx]
    dat_vals <- dat_vec[idx]
    other_vals <- other_mat[idx, , drop = FALSE]

    term1 <- 2 * crossprod(other_vals)

    term1 / length(dat_vec)
}

.feasibility.gaussian <- function(current_vec, other_mat, ...) {
    rep(TRUE, nrow(other_mat))
}

.gaussian <- structure(
    list(
        objfn_all = .evaluate_objective.gaussian,
        objfn     = .evaluate_objective_single.gaussian,
        grad      = .gradient_vec.gaussian,
        hessian   = .hessian_vec.gaussian,
        feas      = .feasibility.gaussian
    ),
    class = "esvd_family"
)



.evaluate_objective_mat.gaussian <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum((nat_mat[idx] - dat[idx])^2)
}

.gradient_mat.gaussian <- function(dat, nat_mat, scalar = 2, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)))

  n <- nrow(dat); p <- ncol(dat)

  1/(n*p) * 2 * (nat_mat - dat)
}
