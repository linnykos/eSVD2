# Distribution: Poisson

.evaluate_objective.poisson <- function(
  dat, x_mat, y_mat, nuisance_param_vec, library_size_vec, ...
) {
  # Check dimensions
  n <- nrow(dat)
  p <- ncol(dat)
  stopifnot(
    length(library_size_vec) == n,
    ncol(x_mat) == ncol(y_mat),
    nrow(x_mat) == n,
    nrow(y_mat) == p
  )

  # Compute natural parameters
  nat_mat <- tcrossprod(x_mat, y_mat)

  # Only compute likelihood on non-missing data
  idx <- which(!is.na(dat))
  stopifnot(length(idx) > 0)
  # `negloglik` contains NAs with the same locations as those in `dat`
  negloglik <- exp(nat_mat + log(library_size_vec)) - nat_mat * dat
  mean(negloglik[idx])
}

# If current vector is x: length(library_size_vec) == 1
# If current vector is y: length(library_size_vec) == n
.evaluate_objective_single.poisson <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat)
  )

  nat_vec <- c(other_mat %*% current_vec)
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  negloglik <- exp(nat_vec + log(library_size_vec)) - nat_vec * dat_vec
  mean(negloglik[idx])
}

.gradient_vec.poisson <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat)
  )

  nat_vec <- c(other_mat %*% current_vec)
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  grad <- other_mat * (exp(nat_vec + log(library_size_vec)) - dat_vec)
  # `grad` contains NA rows corresponding to the NAs in `dat_vec`
  grad <- grad[idx, , drop = FALSE]
  colMeans(grad)
}

.hessian_vec.poisson <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat)
  )

  nat_vec <- c(other_mat %*% current_vec)
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  exp_term <- exp(nat_vec + log(library_size_vec))
  exp_other_term <- .mult_vec_mat(exp_term, other_mat)
  exp_other_vals <- exp_other_term[idx, , drop = FALSE]
  other_vals <- other_mat[idx, , drop = FALSE]
  hess <- crossprod(exp_other_vals, other_vals) / length(idx)
  hess
}

.feasibility.poisson <- function(current_vec, other_mat, ...) {
  TRUE
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
