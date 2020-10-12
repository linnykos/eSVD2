# Distribution: Poisson
# Natural parameter: m_{ij} = u_i^Tv_j
# Relation to canonical parameters: m_{ij} = log(lambda_{ij}), where E(a_{ij}) = lambda_{ij}
# Optimization problem: exp(m_{ij}) - a_{ij}*mu_{ij}

.evaluate_objective.poisson <- function(
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

    # Only compute likelihood on non-missing data
    idx <- which(!is.na(dat))
    stopifnot(length(idx) > 0)

    log_library_size_vals <- rep(log(library_size_vec), p)[idx]
    nat_vals <- nat_mat[idx] + log_library_size_vals
    dat_vals <- dat[idx]
    negloglik <- exp(nat_vals) - nat_vals * dat_vals
    sum(negloglik) / n / p
}

# length(library_size) == 1 if current vector is u
# length(library_size) == n if current vector is v
.evaluate_objective_single.poisson <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- (nat_vec + log(library_size))[idx]
    dat_vals <- dat_vec[idx]
    negloglik <- exp(nat_vals) - nat_vals * dat_vals
    sum(negloglik) / length(dat_vec)
}

.gradient_vec.poisson <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- (nat_vec + log(library_size))[idx]
    dat_vals <- dat_vec[idx]
    other_vals <- other_mat[idx, , drop = FALSE]
    grad <- other_vals * (exp(nat_vals) - dat_vals)

    colSums(grad) / length(dat_vec)
}

.hessian_vec.poisson <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    nat_vals <- (nat_vec + log(library_size))[idx]
    dat_vals <- dat_vec[idx]
    other_vals <- other_mat[idx, , drop = FALSE]

    term1 <- t(other_vals) %*% diag(exp(nat_vals)) %*% other_vals

    term1 / length(dat_vec)
}

.feasibility.poisson <- function(current_vec, other_mat, ...) {
    rep(TRUE, nrow(other_mat))
}

.poisson <- structure(
    list(
        objfn_all = .evaluate_objective.poisson,
        objfn     = .evaluate_objective_single.poisson,
        grad      = .gradient_vec.poisson,
        hessian   = .hessian_vec.poisson,
        feas      = .feasibility.poisson
    ),
    class = "esvd_family"
)



.evaluate_objective_mat.poisson <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat > 0))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum(exp(nat_mat[idx]) - nat_mat[idx]*dat[idx])
}

.gradient_mat.poisson <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat > 0))

  n <- nrow(dat); p <- ncol(dat)

  (exp(nat_mat) - dat)/(n*p)
}
