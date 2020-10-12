# Distribution: negative binomial
# Natural parameter: m_{ij} = u_i^Tv_j
# Relation to canonical parameters: m_{ij} = log(p_{ij})
# Optimization problem: -r*log(1-exp(m_{ij})) - a_{ij}*mu_{ij}

.evaluate_objective.neg_binom <- function(dat, u_mat, v_mat, scalar, ...) {
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
    negloglik <- -scalar * log(1 - exp(nat_vals)) - nat_vals * dat_vals
    sum(negloglik) / n / p
}

.evaluate_objective_single.neg_binom <- function(current_vec, other_mat, dat_vec, scalar, ...) {
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
    negloglik <- -scalar * log(1 - exp(nat_vals)) - nat_vals * dat_vals
    sum(negloglik) / length(dat_vec)
}


# f(x) = log(1 - exp(x))
# f'(x) = -exp(x) / (1 - exp(x)) = -1 / (exp(-x) - 1) = 1 / (1 - exp(-x))
# f''(x) = -exp(-x) / (1 - exp(-x))^2 = [1 / (1 - exp(-x))] * [-exp(-x) / (1 - exp(-x))]
#        = f'(x) * [1 - f'(x)]

.gradient_vec.neg_binom <- function(current_vec, other_mat, dat_vec, scalar, ...) {
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
    grad <- other_vals * (scalar / (exp(-nat_vals) - 1) - dat_vals)
    colSums(grad) / length(dat_vec)
}

.hessian_vec.neg_binom <- function(current_vec, other_mat, dat_vec, scalar, ...) {
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

    fpx = 1 / (1 - exp(nat_vals))
    fppx = fpx * (1 - fpx)
    term1 <- t(other_vals) %*% diag(-scalar * fppx) %*% other_vals

    term1 / length(dat_vec)
}

.feasibility.neg_binom <- function(current_vec, other_mat, ...) {
    nat_vec <- c(other_mat %*% current_vec)
    all(nat_vec < 0)
}

.neg_binom <- structure(
    list(
        objfn_all = .evaluate_objective.neg_binom,
        objfn     = .evaluate_objective_single.neg_binom,
        grad      = .gradient_vec.neg_binom,
        hessian   = .hessian_vec.neg_binom,
        feas      = .feasibility.neg_binom
    ),
    class = "esvd_family"
)



.evaluate_objective_mat.neg_binom <- function(dat, nat_mat, scalar, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat < 0))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum(-scalar * log(1-exp(nat_mat[idx])) - nat_mat[idx]*dat[idx])
}

.gradient_mat.neg_binom <- function(dat, nat_mat, scalar, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat < 0))

  n <- nrow(dat); p <- ncol(dat)
  (scalar * exp(nat_mat)/(1-exp(nat_mat)) - dat)/(n*p)
}
