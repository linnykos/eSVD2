# Gaussian
# theta = mu = mean
.dat_to_nat.gaussian <- function(A, gamma, tol = 1e-3) {
  A
}

# Curved Gaussian
# theta = 1 / mu = 1 / mean
.dat_to_nat.curved_gaussian <- function(A, gamma, tol = 1e-3) {
  res <- 1 / pmax(A, tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Exponential
# theta = -1/lambda = -1/mean
.dat_to_nat.exponential <- function(A, gamma, tol = 1e-3) {
  res <- -1 / (A + tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Poisson
# theta = log(lambda) = log(mean)
.dat_to_nat.poisson <- function(A, gamma, tol = 1e-3) {
  res <- log(A + tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Negative binomial
# theta = log(p), mean = r * p / (1 - p)
.dat_to_nat.neg_binom <- function(A, gamma, tol = 1e-3) {
  res <- sweep(A + tol, 2, gamma, "/")
  res <- log(res / (1 + res))

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Negative binomial 2
# theta = log(p)
.dat_to_nat.neg_binom2 <- function(A, gamma, tol = 1e-3) {
  res <- log(A + tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

# Bernoulli
# Do a simple mapping: A=0 ==> theta=+1 (p ~= 0.73)
#                      A=1 ==> theta=-1 (p ~= 0.27)
.dat_to_nat.bernoulli <- function(A, gamma, tol = 1e-3) {
  res <- ifelse(A > 0.5, 1, -1)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}



## Convert natural parameter to canonical parameter

# This function returns a function for the given family
.nat_to_canon <- function(family)
{
  expr <- switch(family,
                 gaussian        = quote(theta),
                 curved_gaussian = quote(1 / theta),
                 exponential     = quote(-1 / theta),
                 poisson         = quote(exp(theta)),
                 neg_binom       = quote(exp(theta)),
                 neg_binom2      = quote(exp(theta)),
                 bernoulli       = quote(stats::plogis(theta)))
  fn <- function(theta) { theta }
  body(fn)[[2]] <- expr
  fn
}



#' Internal Constructor for distribution family
#'
#' @param family string
#'
#' @return an object with conversion functions
esvd_family <- function(family)
{
  family <- as.character(family)
  obj <- .esvd_family(family)

  # Additional functions
  obj$feasibility <- function(nat)
  {
    obj$feas_always ||
      (all(nat > obj$domain[1]) && all(nat < obj$domain[2]))
  }
  obj$dat_to_nat <- get(sprintf(".dat_to_nat.%s", family))
  obj$nat_to_canon <- .nat_to_canon(family)

  structure(obj, class = "esvd_family")
}
