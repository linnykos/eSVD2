# Distribution: Bernoulli

.evaluate_objective.bernoulli <- function(
  dat, x_mat, y_mat, nuisance_param_vec, library_size_vec, ...
) {
  # Check dimensions
  n <- nrow(dat)
  p <- ncol(dat)
  stopifnot(
    ncol(x_mat) == ncol(y_mat),
    nrow(x_mat) == n,
    nrow(y_mat) == p
  )

  # Compute natural parameters
  nat_mat <- tcrossprod(x_mat, y_mat)

  # Only compute likelihood on non-missing data
  idx <- which(!is.na(dat))
  stopifnot(length(idx) > 0)

  nat_vals <- nat_mat[idx]
  dat_vals <- dat[idx]
  negloglik <- log(1 + exp(nat_vals)) - nat_vals * dat_vals
  mean(negloglik)
}

.evaluate_objective_single.bernoulli <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat)
  )

  nat_vec <- c(other_mat %*% current_vec)
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  nat_vals <- nat_vec[idx]
  dat_vals <- dat_vec[idx]
  negloglik <- log(1 + exp(nat_vals)) - nat_vals * dat_vals
  mean(negloglik)
}

# f(x) = log(1 + exp(x))
# f'(x) = exp(x) / (1 + exp(x)) = 1 / (1 + exp(-x)) = plogis(x)
# f''(x) = exp(-x) / (1 + exp(-x))^2 = [1 / (1 + exp(-x))] * [exp(-x) / (1 + exp(-x))]
#        = f'(x) * [1 - f'(x)]

.gradient_vec.bernoulli <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
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
  grad <- other_vals * (stats::plogis(nat_vals) - dat_vals)

  colMeans(grad)
}

.hessian_vec.bernoulli <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size, ...
) {
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

  fpx = stats::plogis(nat_vals)
  fppx = fpx * (1 - fpx)
  hess <- crossprod(.mult_vec_mat(fppx, other_vals), other_vals) / length(idx)
  hess
}

.feasibility.bernoulli <- function(current_vec, other_mat, ...) {
  TRUE
}

.bernoulli <- structure(
  list(
    objfn_all = .evaluate_objective.bernoulli,
    objfn     = .evaluate_objective_single.bernoulli,
    grad      = .gradient_vec.bernoulli,
    hessian   = .hessian_vec.bernoulli,
    feas      = .feasibility.bernoulli
  ),
  class = "esvd_family"
)
