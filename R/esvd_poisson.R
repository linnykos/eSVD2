# Distribution: Poisson
# Canonical parameter: lambda
# Natural parameter: theta = log(lambda)

# See eSVD2_writing/writeup/2021-05-20-covariates.pdf
#
# Log-density for the whole data matrix [n x p]
.log_prob.poisson <- function(A, theta, s, gamma)
{
  A * theta - exp(log(s) + theta)
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.poisson <- function(Ai, thetai, si, gamma)
{
  Ai * thetai - exp(log(si) + thetai)
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.poisson <- function(Aj, thetaj, s, gammaj)
{
  Aj * thetaj - exp(log(s) + thetaj)
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.poisson <- function(Ai, thetai, si, gamma)
{
  Ai - exp(log(si) + thetai)
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.poisson <- function(Aj, thetaj, s, gammaj)
{
  Aj - exp(log(s) + thetaj)
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.poisson <- function(Ai, thetai, si, gamma)
{
  -exp(log(si) + thetai)
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.poisson <- function(Aj, thetaj, s, gammaj)
{
  -exp(log(s) + thetaj)
}

# Feasibility of the natural parameter
.feasibility.poisson <- function(theta)
{
  TRUE
}

# Initialize the natural parameter from data
# theta = log(lambda) = log(mean)
.dat_to_nat.poisson <- function(A, gamma, tol = 1e-3) {
  res <- log(A + tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Convert natural parameter to canonical parameter
.nat_to_canon.poisson <- function(theta) {
  exp(theta)
}

.esvd.poisson <- structure(
  list(
    name           = "poisson",
    log_prob       = .log_prob.poisson,
    log_prob_row   = .log_prob_row.poisson,
    log_prob_col   = .log_prob_col.poisson,
    dlog_prob_row  = .dlog_prob_row.poisson,
    dlog_prob_col  = .dlog_prob_col.poisson,
    d2log_prob_row = .d2log_prob_row.poisson,
    d2log_prob_col = .d2log_prob_col.poisson,
    feasibility    = .feasibility.poisson,
    feas_always    = TRUE,
    domain         = c(-Inf, Inf),
    dat_to_nat     = .dat_to_nat.poisson,
    nat_to_canon   = .nat_to_canon.poisson
  ),
  class = "esvd_family"
)
.esvd.poisson <- list2env(.esvd.poisson)
