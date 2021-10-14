# Distribution: negative binomial
# Log-density for the whole data matrix [n x p]
.log_prob.neg_binom <- function(A, theta, s, gamma)
{
  slog1theta <- log(1 - exp(theta)) * s
  sweep(slog1theta, 2, gamma, "*") + A * theta
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.neg_binom <- function(Ai, thetai, si, gamma)
{
  Ai * thetai + si * gamma * log(1 - exp(thetai))
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.neg_binom <- function(Aj, thetaj, s, gammaj)
{
  Aj * thetaj + s * gammaj * log(1 - exp(thetaj))
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.neg_binom <- function(Ai, thetai, si, gamma)
{
  exptheta <- exp(thetai)
  Ai - si * gamma * exptheta / (1 - exptheta)
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.neg_binom <- function(Aj, thetaj, s, gammaj)
{
  exptheta <- exp(thetaj)
  Aj - s * gammaj * exptheta / (1 - exptheta)
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.neg_binom <- function(Ai, thetai, si, gamma)
{
  exptheta <- exp(thetai)
  -si * gamma * exptheta / (1 - exptheta)^2
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.neg_binom <- function(Aj, thetaj, s, gammaj)
{
  exptheta <- exp(thetaj)
  -s * gammaj * exptheta / (1 - exptheta)^2
}

# Feasibility of the natural parameter
.feasibility.neg_binom <- function(theta)
{
  all(theta < 0)
}

.dat_to_nat.neg_binom <- function(A, gamma, tol = 1e-3){
  res <- sapply(1:ncol(A), function(j){
    (A[,j] + tol)/gamma[j]
  })
  res <- log(res/(1+res))

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

.nat_to_canon.neg_binom <- function(theta){
  exp(theta)
}

.esvd.neg_binom <- structure(
  list(
    name           = "neg_binom",
    log_prob       = .log_prob.neg_binom,
    log_prob_row   = .log_prob_row.neg_binom,
    log_prob_col   = .log_prob_col.neg_binom,
    dlog_prob_row  = .dlog_prob_row.neg_binom,
    dlog_prob_col  = .dlog_prob_col.neg_binom,
    d2log_prob_row = .d2log_prob_row.neg_binom,
    d2log_prob_col = .d2log_prob_col.neg_binom,
    feasibility    = .feasibility.neg_binom,
    feas_always    = FALSE,
    domain         = c(-Inf, 0),
    dat_to_nat     = .dat_to_nat.neg_binom,
    nat_to_canon   = .nat_to_canon.neg_binom
  ),
  class = "esvd_family"
)
.esvd.neg_binom <- list2env(.esvd.neg_binom)
