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



# Softplus function
# log(1 + exp(x)) = log(1 + exp(-|x|)) + max(x, 0)
softplus <- function(x)
{
  pmax(x, 0) + log(1 + exp(-abs(x)))
}

# See eSVD2_writing/writeup/Writeup4
#
# Log-density for the whole data matrix [n x p]
.log_prob.bernoulli <- function(A, theta, s, gamma)
{
  A * theta - softplus(theta)
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.bernoulli <- function(Ai, thetai, si, gamma)
{
  Ai * thetai - softplus(thetai)
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.bernoulli <- function(Aj, thetaj, s, gammaj)
{
  Aj * thetaj - softplus(thetaj)
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.bernoulli <- function(Ai, thetai, si, gamma)
{
  Ai - stats::plogis(thetai)
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.bernoulli <- function(Aj, thetaj, s, gammaj)
{
  Aj - stats::plogis(thetaj)
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.bernoulli <- function(Ai, thetai, si, gamma)
{
  p <- stats::plogis(thetai)
  -p * (1 - p)
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.bernoulli <- function(Aj, thetaj, s, gammaj)
{
  p <- stats::plogis(thetaj)
  -p * (1 - p)
}

# Feasibility of the natural parameter
.feasibility.gaussian <- function(theta)
{
  TRUE
}

.esvd.bernoulli <- structure(
  list(
    name           = "bernoulli",
    log_prob       = .log_prob.bernoulli,
    log_prob_row   = .log_prob_row.bernoulli,
    log_prob_col   = .log_prob_col.bernoulli,
    dlog_prob_row  = .dlog_prob_row.bernoulli,
    dlog_prob_col  = .dlog_prob_col.bernoulli,
    d2log_prob_row = .d2log_prob_row.bernoulli,
    d2log_prob_col = .d2log_prob_col.bernoulli,
    feasibility    = .feasibility.bernoulli,
    feas_always    = TRUE
  ),
  class = "esvd_family"
)
