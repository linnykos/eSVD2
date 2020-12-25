# Distribution: one-parameter Gaussian where sd = mean/scalar

.evaluate_objective.curved_gaussian <- function(
  dat, x_mat, y_mat, nuisance_param_vec, library_size_vec, ...
) {
  # Check dimensions
  n <- nrow(dat)
  p <- ncol(dat)
  stopifnot(
    length(nuisance_param_vec) == p,
    length(library_size_vec) == n,
    ncol(x_mat) == ncol(y_mat),
    nrow(x_mat) == n,
    nrow(y_mat) == p
  )

  # Compute natural parameters
  nat_mat <- tcrossprod(x_mat, y_mat)
  stopifnot(all(nat_mat > 0))

  # Only compute likelihood on non-missing data
  idx <- which(!is.na(dat))
  stopifnot(length(idx) > 0)
  # `negloglik` contains NAs with the same locations as those in `dat`
  negloglik <- -log(nat_mat) +
    .mult_vec_mat(1 / (2 * library_size_vec), .mult_mat_vec(dat^2, nuisance_param_vec^2)) * nat_mat^2 -
    .mult_mat_vec(dat, nuisance_param_vec^2) * nat_mat
  mean(negloglik[idx])
}

# If current vector is x: length(library_size_vec) == 1
# If current vector is y: length(library_size_vec) == n
.evaluate_objective_single.curved_gaussian <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  check_dim_x <- length(library_size_vec) == 1 &&
                 length(nuisance_param_vec) == nrow(other_mat)
  check_dim_y <- length(library_size_vec) == nrow(other_mat) &&
                 length(nuisance_param_vec) == 1
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat),
    check_dim_x || check_dim_y
  )

  nat_vec <- c(other_mat %*% current_vec)
  stopifnot(all(nat_vec > 0))
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  negloglik <- -log(nat_vec) +
    nuisance_param_vec^2 * dat_vec^2 * nat_vec^2 / (2 * library_size_vec) -
    nuisance_param_vec^2 * dat_vec * nat_vec
  mean(negloglik[idx])
}

.gradient_vec.curved_gaussian <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  check_dim_x <- length(library_size_vec) == 1 &&
                 length(nuisance_param_vec) == nrow(other_mat)
  check_dim_y <- length(library_size_vec) == nrow(other_mat) &&
                 length(nuisance_param_vec) == 1
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat),
    check_dim_x || check_dim_y
  )

  nat_vec <- c(other_mat %*% current_vec)
  stopifnot(all(nat_vec > 0))
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  grad <- other_mat * (-1 / nat_vec + nuisance_param_vec^2 * dat_vec^2 * nat_vec / library_size_vec -
                       nuisance_param_vec^2 * dat_vec)
  # `grad` contains NA rows corresponding to the NAs in `dat_vec`
  grad <- grad[idx, ]
  colMeans(grad)
}

.hessian_vec.curved_gaussian <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  check_dim_x <- length(library_size_vec) == 1 &&
                 length(nuisance_param_vec) == nrow(other_mat)
  check_dim_y <- length(library_size_vec) == nrow(other_mat) &&
                 length(nuisance_param_vec) == 1
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat),
    check_dim_x || check_dim_y
  )

  nat_vec <- c(other_mat %*% current_vec)
  stopifnot(all(nat_vec > 0))
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  other_vals <- other_mat[idx, , drop = FALSE]
  term1_vals <- .mult_vec_mat(1 / nat_vec^2, other_mat)[idx, , drop = FALSE]
  term2_vals <- .mult_vec_mat(nuisance_param_vec^2 * dat_vec^2 / library_size_vec, other_mat)[idx, , drop = FALSE]
  hess <- crossprod(term1_vals + term2_vals, other_vals) / length(idx)
  hess
}

.feasibility.curved_gaussian <- function(current_vec, other_mat, ...) {
  nat_vec <- c(other_mat %*% current_vec)
  all(nat_vec > 0)
}

.curved_gaussian <- structure(
  list(
    objfn_all = .evaluate_objective.curved_gaussian,
    objfn     = .evaluate_objective_single.curved_gaussian,
    grad      = .gradient_vec.curved_gaussian,
    hessian   = .hessian_vec.curved_gaussian,
    feas      = .feasibility.curved_gaussian
  ),
  class = "esvd_family"
)
