# Distribution: Bernoulli
# Canonical parameter: p
# Natural parameter: theta = log(p/(1-p))

# Softplus function
# log(1 + exp(x)) = log(1 + exp(-|x|)) + max(x, 0)
softplus <- function(x)
{
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# See eSVD2_writing/writeup/2021-05-20-covariates.pdf
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
.feasibility.bernoulli <- function(theta)
{
  TRUE
}

# Initialize the natural parameter from data
# Do a simple mapping: A=0 ==> theta=+1 (p ~= 0.73)
#                      A=1 ==> theta=-1 (p ~= 0.27)
.dat_to_nat.bernoulli <- function(A, gamma, tol = 1e-3) {
  res <- ifelse(A > 0.5, 1, -1)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Convert natural parameter to canonical parameter
.nat_to_canon.bernoulli <- function(theta) {
  stats::plogis(theta)
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
    feas_always    = TRUE,
    domain         = c(-Inf, Inf),
    dat_to_nat     = .dat_to_nat.bernoulli,
    nat_to_canon   = .nat_to_canon.bernoulli
  ),
  class = "esvd_family"
)
.esvd.bernoulli <- list2env(.esvd.bernoulli)
