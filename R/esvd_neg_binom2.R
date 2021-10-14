# Distribution: negative binomial

# See eSVD2_writing/writeup/2021-05-20-covariates.pdf
#
# Log-density for the whole data matrix [n x p]
.log_prob.neg_binom2 <- function(A, theta, s, gamma)
{
  exptheta <- exp(theta)
  logexpthetag <- log(sweep(exptheta, 2, gamma, "+"))
  Ag <- sweep(A, 2, gamma, "+")
  A*theta - Ag*logexpthetag
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.neg_binom2 <- function(Ai, thetai, si, gamma)
{
  Ai * thetai - (Ai + gamma) * log(exp(thetai) + gamma)
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.neg_binom2 <- function(Aj, thetaj, s, gammaj)
{
  Aj * thetaj - (Aj + gammaj) * log(exp(thetaj) + gammaj)
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.neg_binom2 <- function(Ai, thetai, si, gamma)
{
  exptheta <- exp(thetai)
  Ai - (Ai + gamma)*(exptheta / (exptheta + gamma))
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.neg_binom2 <- function(Aj, thetaj, s, gammaj)
{
  exptheta <- exp(thetaj)
  Aj - (Aj + gammaj) * (exptheta / (exptheta + gammaj))
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.neg_binom2 <- function(Ai, thetai, si, gamma)
{
  exptheta <- exp(thetai)
  -(Ai + gamma) * exptheta * gamma / (exptheta + gamma)^2
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.neg_binom2 <- function(Aj, thetaj, s, gammaj)
{
  exptheta <- exp(thetaj)
  -(Aj + gammaj) * exptheta * gammaj / (exptheta + gammaj)^2
}

# Feasibility of the natural parameter
.feasibility.neg_binom2 <- function(theta)
{
  TRUE
}

# Initialize the natural parameter from data
# theta = log(p)
.dat_to_nat.neg_binom2 <- function(A, gamma, tol = 1e-3) {
  res <- log(A + tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Convert natural parameter to mean
.nat_to_canon.neg_binom2 <- function(theta) {
  exp(theta)
}

.esvd.neg_binom2 <- structure(
  list(
    name           = "neg_binom2",
    log_prob       = .log_prob.neg_binom2,
    log_prob_row   = .log_prob_row.neg_binom2,
    log_prob_col   = .log_prob_col.neg_binom2,
    dlog_prob_row  = .dlog_prob_row.neg_binom2,
    dlog_prob_col  = .dlog_prob_col.neg_binom2,
    d2log_prob_row = .d2log_prob_row.neg_binom2,
    d2log_prob_col = .d2log_prob_col.neg_binom2,
    feasibility    = .feasibility.neg_binom2,
    feas_always    = TRUE,
    domain         = c(-Inf, Inf),
    dat_to_nat     = .dat_to_nat.neg_binom2,
    nat_to_canon   = .nat_to_canon.neg_binom2
  ),
  class = "esvd_family"
)
.esvd.neg_binom2 <- list2env(.esvd.neg_binom2)
