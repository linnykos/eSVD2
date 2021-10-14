# Distribution: Gaussian

# See eSVD2_writing/writeup/2021-05-20-covariates.pdf
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

# Initialize the natural parameter from data
# theta = mu = mean
.dat_to_nat.gaussian <- function(A, gamma, tol = 1e-3) {
  A
}

# Convert natural parameter to mean
.nat_to_canon.gaussian <- function(theta) {
  theta
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
    feas_always    = TRUE,
    domain         = c(-Inf, Inf),
    dat_to_nat     = .dat_to_nat.gaussian,
    nat_to_canon   = .nat_to_canon.gaussian
  ),
  class = "esvd_family"
)
.esvd.gaussian <- list2env(.esvd.gaussian)
