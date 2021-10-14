# Distribution: exponential

.evaluate_objective.exponential <- function(
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
  stopifnot(all(nat_mat < 0))

  # Only compute likelihood on non-missing data
  idx <- which(!is.na(dat))
  stopifnot(length(idx) > 0)
  negloglik <- -library_size_vec * log(-nat_mat) - dat * nat_mat
  mean(negloglik[idx])
}

# If current vector is x: length(library_size_vec) == 1
# If current vector is y: length(library_size_vec) == n
.evaluate_objective_single.exponential <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat)
  )

  nat_vec <- c(other_mat %*% current_vec)
  stopifnot(all(nat_vec < 0))
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  negloglik <- -library_size_vec * log(-nat_vec) - dat_vec * nat_vec
  mean(negloglik[idx])
}

.gradient_vec.exponential <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat)
  )

  nat_vec <- c(other_mat %*% current_vec)
  stopifnot(all(nat_vec < 0))
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  grad <- other_mat * (-library_size_vec / nat_vec - dat_vec)
  # `grad` contains NA rows corresponding to the NAs in `dat_vec`
  grad <- grad[idx, , drop = FALSE]
  colMeans(grad)
}

.hessian_vec.exponential <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat)
  )

  nat_vec <- c(other_mat %*% current_vec)
  stopifnot(all(nat_vec < 0))
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  other_vals <- other_mat[idx, , drop = FALSE]
  term1_vals <- .mult_vec_mat(library_size_vec / nat_vec^2, other_mat)[idx, , drop = FALSE]
  hess <- crossprod(term1_vals, other_vals) / length(idx)
  hess
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



# See eSVD2_writing/writeup/Writeup4
#
# Log-density for the whole data matrix [n x p]
.log_prob.exponential <- function(A, theta, s, gamma)
{
  A * theta + s * log(-theta)
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.exponential <- function(Ai, thetai, si, gamma)
{
  Ai * thetai + si * log(-thetai)
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.exponential <- function(Aj, thetaj, s, gammaj)
{
  Aj * thetaj + s * log(-thetaj)
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.exponential <- function(Ai, thetai, si, gamma)
{
  Ai + si / thetai
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.exponential <- function(Aj, thetaj, s, gammaj)
{
  Aj + s / thetaj
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.exponential <- function(Ai, thetai, si, gamma)
{
  -si / thetai^2
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.exponential <- function(Aj, thetaj, s, gammaj)
{
  -s / thetaj^2
}

# Feasibility of the natural parameter
.feasibility.exponential <- function(theta)
{
  all(theta < 0)
}

.esvd.exponential <- structure(
  list(
    log_prob       = .log_prob.exponential,
    log_prob_row   = .log_prob_row.exponential,
    log_prob_col   = .log_prob_col.exponential,
    dlog_prob_row  = .dlog_prob_row.exponential,
    dlog_prob_col  = .dlog_prob_col.exponential,
    d2log_prob_row = .d2log_prob_row.exponential,
    d2log_prob_col = .d2log_prob_col.exponential,
    feasibility    = .feasibility.exponential,
    feas_always    = FALSE,
    domain         = c(-Inf, 0)
  ),
  class = "esvd_family"
)
.esvd.exponential <- list2env(.esvd.exponential)
