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
  grad <- grad[idx, , drop = FALSE]
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



# See eSVD2_writing/writeup/Writeup4
#
# Log-density for the whole data matrix [n x p]
.log_prob.curved_gaussian <- function(A, theta, s, gamma)
{
  gamma2A <- sweep(A, 2, gamma^2, "*")
  gamma2A * theta - 0.5 * gamma2A * A * theta^2 / s + log(theta)
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.curved_gaussian <- function(Ai, thetai, si, gamma)
{
  gamma2A <- gamma^2 * Ai
  gamma2A * thetai - 0.5 * gamma2A * Ai * thetai^2 / si + log(thetai)
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.curved_gaussian <- function(Aj, thetaj, s, gammaj)
{
  gamma2A <- gammaj^2 * Aj
  gamma2A * thetaj - 0.5 * gamma2A * Aj * thetaj^2 / s + log(thetaj)
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.curved_gaussian <- function(Ai, thetai, si, gamma)
{
  gamma2A <- gamma^2 * Ai
  gamma2A - gamma2A * Ai * thetai / si + 1 / thetai
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.curved_gaussian <- function(Aj, thetaj, s, gammaj)
{
  gamma2A <- gammaj^2 * Aj
  gamma2A - gamma2A * Aj * thetaj / s + 1 / thetaj
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.curved_gaussian <- function(Ai, thetai, si, gamma)
{
  -(gamma * Ai)^2 / si - 1 / thetai^2
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.curved_gaussian <- function(Aj, thetaj, s, gammaj)
{
  -(gammaj * Aj)^2 / s - 1 / thetaj^2
}

# Feasibility of the natural parameter
.feasibility.curved_gaussian <- function(theta)
{
  all(theta > 0)
}

.esvd.curved_gaussian <- structure(
  list(
    log_prob       = .log_prob.curved_gaussian,
    log_prob_row   = .log_prob_row.curved_gaussian,
    log_prob_col   = .log_prob_col.curved_gaussian,
    dlog_prob_row  = .dlog_prob_row.curved_gaussian,
    dlog_prob_col  = .dlog_prob_col.curved_gaussian,
    d2log_prob_row = .d2log_prob_row.curved_gaussian,
    d2log_prob_col = .d2log_prob_col.curved_gaussian,
    feasibility    = .feasibility.curved_gaussian,
    feas_always    = FALSE
  ),
  class = "esvd_family"
)
