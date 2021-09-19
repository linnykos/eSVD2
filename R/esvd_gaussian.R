# Distribution: Gaussian

.evaluate_objective.gaussian <- function(
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

  # Compute natural parameters, Theta = XY'
  nat_mat <- tcrossprod(x_mat, y_mat)

  # Only compute likelihood on non-missing data
  idx <- which(!is.na(dat))
  stopifnot(length(idx) > 0)
  # `dat_vals` and `negloglik` contain NAs with the same locations as those in `dat`
  dat_vals <- (dat / library_size_vec)
  negloglik <- .mult_mat_vec((nat_mat - dat_vals)^2, 1 / nuisance_param_vec^2)
  mean(negloglik[idx])
}

# If current vector is x:
#   length(library_size_vec) == 1, length(nuisance_param_vec) == p
# If current vector is y:
#   length(library_size_vec) == n, length(nuisance_param_vec) == 1
.evaluate_objective_single.gaussian <- function(
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
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  dat_vals <- (dat_vec / library_size_vec)
  negloglik <- (nat_vec - dat_vals)^2 / nuisance_param_vec^2
  mean(negloglik[idx])
}

.gradient_vec.gaussian <- function(
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
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  dat_vals <- (dat_vec / library_size_vec)
  other_vals <- other_mat
  # Broadcast `nat_vals - dat_vals` to each column of `other_vals`
  grad <- 2 * other_vals * (nat_vec - dat_vals) / nuisance_param_vec^2
  # `grad` contains NA rows corresponding to the NAs in `dat_vec`
  grad <- grad[idx, , drop = FALSE]
  colMeans(grad)
}

.hessian_vec.gaussian <- function(
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

  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  if(length(nuisance_param_vec) == 1) {
    other_vals <- other_mat / nuisance_param_vec
  } else {
    other_vals <- .mult_vec_mat(1 / nuisance_param_vec, other_mat)
  }

  other_vals <- other_vals[idx, , drop = FALSE]
  hess <- 2 * crossprod(other_vals) / length(idx)
  hess
}

.feasibility.gaussian <- function(current_vec, other_mat, ...) {
  TRUE
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



# See eSVD2_writing/writeup/Writeup4
#
# Log-density for the whole data matrix [n x p]
.log_prob.gaussian <- function(A, theta, s, gamma)
{
  res <- -0.5 * (theta - A / s)^2 * s
  sweep(res, 2, gamma^2, "/")
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.gaussian <- function(Ai, thetai, si, gamma)
{
  -0.5 * (thetai - Ai / si)^2 * si / gamma^2
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.gaussian <- function(Aj, thetaj, s, gammaj)
{
  -0.5 * (thetaj - Aj / s)^2 * s / gammaj^2
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.gaussian <- function(Ai, thetai, si, gamma)
{
  -(si * thetai - Ai) / gamma^2
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.gaussian <- function(Aj, thetaj, s, gammaj)
{
  -(s * thetaj - Aj) / gammaj^2
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.gaussian <- function(Ai, thetai, si, gamma)
{
  -si / gamma^2
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.gaussian <- function(Aj, thetaj, s, gammaj)
{
  -s / gammaj^2
}

# Feasibility of the natural parameter
.feasibility.gaussian <- function(theta)
{
  TRUE
}

.esvd.gaussian <- structure(
  list(
    name           = "gaussian",
    log_prob       = .log_prob.gaussian,
    log_prob_row   = .log_prob_row.gaussian,
    log_prob_col   = .log_prob_col.gaussian,
    dlog_prob_row  = .dlog_prob_row.gaussian,
    dlog_prob_col  = .dlog_prob_col.gaussian,
    d2log_prob_row = .d2log_prob_row.gaussian,
    d2log_prob_col = .d2log_prob_col.gaussian,
    feasibility    = .feasibility.gaussian,
    feas_always    = TRUE
  ),
  class = "esvd_family"
)
