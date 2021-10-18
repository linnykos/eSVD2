# Distribution: exponential
# Canonical parameter: lambda (mean = lambda)
# Natural parameter: theta = -1/lambda

# See eSVD2_writing/writeup/2021-05-20-covariates.pdf
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

# Initialize the natural parameter from data
# theta = -1/lambda = -1/mean
.dat_to_nat.exponential <- function(A, gamma, tol = 1e-3) {
  res <- -1 / (A + tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Convert natural parameter to canonical parameter
.nat_to_canon.exponential <- function(theta) {
   -1 / theta
}

.esvd.exponential <- structure(
  list(
    name           = "exponential",
    log_prob       = .log_prob.exponential,
    log_prob_row   = .log_prob_row.exponential,
    log_prob_col   = .log_prob_col.exponential,
    dlog_prob_row  = .dlog_prob_row.exponential,
    dlog_prob_col  = .dlog_prob_col.exponential,
    d2log_prob_row = .d2log_prob_row.exponential,
    d2log_prob_col = .d2log_prob_col.exponential,
    feasibility    = .feasibility.exponential,
    feas_always    = FALSE,
    domain         = c(-Inf, 0),
    dat_to_nat     = .dat_to_nat.exponential,
    nat_to_canon   = .nat_to_canon.exponential
  ),
  class = "esvd_family"
)
.esvd.exponential <- list2env(.esvd.exponential)
