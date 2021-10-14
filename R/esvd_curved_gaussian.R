# Distribution: one-parameter Gaussian where sd = mean/scalar

# See eSVD2_writing/writeup/2021-05-20-covariates.pdf
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

# Initialize the natural parameter from data
# theta = 1 / mu = 1 / mean
.dat_to_nat.curved_gaussian <- function(A, gamma, tol = 1e-3) {
  res <- 1 / pmax(A, tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Convert natural parameter to mean
.nat_to_canon.curved_gaussian <- function(theta) {
  1 / theta
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
    feas_always    = FALSE,
    domain         = c(0, Inf),
    dat_to_nat     = .dat_to_nat.curved_gaussian,
    nat_to_canon   = .nat_to_canon.curved_gaussian
  ),
  class = "esvd_family"
)
.esvd.curved_gaussian <- list2env(.esvd.curved_gaussian)
