# Compute W'DW, W [m x n], D [m x m] is diagonal, d [m x 1] the diagonal elements of D
wtdw <- function(w, d)
{
  # Broadcast d to each column of w
  dw <- w * d
  crossprod(dw, w)
}

# Objective function value for all the parameters (X, Y, B)
#
# A:      data matrix [n x p]
# X:      latent factor for cells [n x k]
# Y:      latent factor for genes [p x k]
# Z:      covariates [n x r] or NULL
# B:      regression coefficients for Z [p x r] or NULL
# s:      library size vector [n x 1]
# gamma:  nuisance parameter vector [p x 1]
# offset: offset vector [n x 1]
objfn_all <- function(X, Y, B, Z, A, family, s, gamma, offset, ...)
{
  # theta = XY'+ZB' + offset
  #
  # We allow Z and B to be NULL, meaning no covariates are in the model
  theta <- tcrossprod(cbind(X, Z), cbind(Y, B))  # [n x p]
  theta <- sweep(theta, 1, offset, "+")

  # Call the family function log_prob() to compute entry-wise log-density value
  log_prob <- family$log_prob(A, theta, s, gamma)  # [n x p]
  # Some entries in A may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(A))
  stopifnot(length(obs_ind) > 0)
  # Compute negative log-likelihood function value
  -mean(log_prob[obs_ind])
}

# Objective function value for the i-th row of X
objfn_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, ...)
{
  if(!is.null(family$cpp_functions))
  {
    return(objfn_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma, offseti))
  }

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offseti     # [p x 1]
  log_prob <- family$log_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # Some entries in Ai may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Ai))
  stopifnot(length(obs_ind) > 0)
  # Compute negative log-likelihood function value
  -mean(log_prob[obs_ind])
}

# Objective function value for the j-th row of Y
objfn_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, ...)
{
  if(!is.null(family$cpp_functions))
  {
    return(objfn_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj, offset))
  }

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset      # [n x 1]
  log_prob <- family$log_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # Some entries in Aj may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Aj))
  stopifnot(length(obs_ind) > 0)
  # Compute negative log-likelihood function value
  -mean(log_prob[obs_ind])
}

# Gradient for the i-th row of X
grad_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, ...)
{
  if(!is.null(family$cpp_functions))
  {
    return(grad_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma, offseti))
  }

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offseti       # [p x 1]
  dlog_prob <- family$dlog_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # Some entries in Ai may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Ai))
  stopifnot(length(obs_ind) > 0)
  # Compute the gradient
  Ysub <- Y[obs_ind, , drop = FALSE]
  prod <- crossprod(Ysub, dlog_prob[obs_ind])
  -drop(prod) / length(obs_ind)
}

# Gradient for the j-th row of Y
grad_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, ...)
{
  if(!is.null(family$cpp_functions))
  {
    return(grad_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj, offset))
  }

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset        # [n x 1]
  dlog_prob <- family$dlog_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # Some entries in Aj may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Aj))
  stopifnot(length(obs_ind) > 0)
  # Compute the gradient
  Xsub <- X[obs_ind, , drop = FALSE]
  prod <- crossprod(Xsub, dlog_prob[obs_ind])
  -drop(prod) / length(obs_ind)
}

# Hessian for the i-th row of X
hessian_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, ...)
{
  if(!is.null(family$cpp_functions))
  {
    return(hessian_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma, offseti))
  }

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offseti         # [p x 1]
  d2log_prob <- family$d2log_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # Some entries in Ai may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Ai))
  stopifnot(length(obs_ind) > 0)
  # Compute the Hessian
  Ysub <- Y[obs_ind, , drop = FALSE]
  diag <- d2log_prob[obs_ind]
  -wtdw(Ysub, diag) / length(obs_ind)
}

# Hessian for the j-th row of Y
hessian_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, ...)
{
  if(!is.null(family$cpp_functions))
  {
    return(hessian_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj, offset))
  }

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset          # [n x 1]
  d2log_prob <- family$d2log_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # Some entries in Aj may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Aj))
  stopifnot(length(obs_ind) > 0)
  # Compute the Hessian
  Xsub <- X[obs_ind, , drop = FALSE]
  diag <- d2log_prob[obs_ind]
  -wtdw(Xsub, diag) / length(obs_ind)
}

# Move direction for the i-th row of X, d = -inv(H) * g
direction_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, ...)
{
  if(!is.null(family$cpp_functions))
  {
    return(direction_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma, offseti))
  }

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offseti         # [p x 1]
  dlog_prob <- family$dlog_prob_row(Ai, thetai, si, gamma)    # [p x 1]
  d2log_prob <- family$d2log_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # Some entries in Ai may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Ai))
  stopifnot(length(obs_ind) > 0)
  # Compute the Hessian and gradient
  Ysub <- Y[obs_ind, , drop = FALSE]
  g <- crossprod(Ysub, -dlog_prob[obs_ind])
  direc <- -solve(wtdw(Ysub, -d2log_prob[obs_ind]), g)

  # # H = VY' * VY / N
  # # g = Y' * w / N
  # # d = -inv(H) * g = -inv(VY' * VY) * (Y' * w)
  # # Let VY = QR, then VY' * VY = R'Q'QR = R'R
  # # d = -inv(R) * inv(R') * (Y' * w)
  # weights <- -d2log_prob[obs_ind]
  # VY <- Ysub * sqrt(weights)
  # g <- crossprod(Ysub, -dlog_prob[obs_ind])
  # R <- qr.R(qr(VY))
  # direc <- -backsolve(R, forwardsolve(t(R), g))

  list(grad = g / length(obs_ind), direction = direc)
}

# Move direction for the j-th row of Y
direction_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, ...)
{
  if(!is.null(family$cpp_functions))
  {
    return(direction_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj, offset))
  }

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset          # [n x 1]
  dlog_prob <- family$dlog_prob_col(Aj, thetaj, s, gammaj)    # [n x 1]
  d2log_prob <- family$d2log_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # Some entries in Aj may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Aj))
  stopifnot(length(obs_ind) > 0)
  # Compute the Hessian and gradient
  Xsub <- X[obs_ind, , drop = FALSE]
  g <- crossprod(Xsub, -dlog_prob[obs_ind])
  direc <- -solve(wtdw(Xsub, -d2log_prob[obs_ind]), g)

  # # H = VX' * VX / N
  # # g = X' * w / N
  # # d = -inv(H) * g = -inv(VX' * VX) * (X' * w)
  # # Let VX = QR, then VX' * VX = R'Q'QR = R'R
  # # d = -inv(R) * inv(R') * (X' * w)
  # weights <- -d2log_prob[obs_ind]
  # VX <- Xsub * sqrt(weights)
  # g <- crossprod(Xsub, -dlog_prob[obs_ind])
  # R <- qr.R(qr(VX))
  # direc <- -backsolve(R, forwardsolve(t(R), g))

  list(grad = g / length(obs_ind), direction = direc)
}

# Feasibility of the i-th row of X
# thetai = Xi * Y' + Zi * B' + offseti
feas_Xi <- function(Xi, Y, B, Zi, family, offseti, ...)
{
  if(family$feas_always)
    return(TRUE)

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offseti
  family$feasibility(thetai)
}

# Feasibility of the j-th row of Y
# thetaj = X * Yj' + Z * Bj' + offset
feas_Yj <- function(Yj, X, Bj, Z, family, offset, ...)
{
  if(family$feas_always)
    return(TRUE)

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset
  family$feasibility(thetaj)
}

#################
.parse_library_size <- function(dat, library_size_vec) {
  stopifnot(length(library_size_vec) %in% c(1, nrow(dat)))
  n <- nrow(dat)

  if(any(is.na(library_size_vec))){
    library_size_vec <- rowSums(dat)
  } else if(length(library_size_vec) == 1) {
    library_size_vec <- rep(library_size_vec[1], n)
  }

  library_size_vec/min(library_size_vec)
}

.string_to_distr_funcs <- function(family)
{
  switch(family,
         bernoulli = .esvd.bernoulli,
         curved_gaussian = .esvd.curved_gaussian,
         exponential = .esvd.exponential,
         gaussian = .esvd.gaussian,
         neg_binom = .esvd.neg_binom,
         neg_binom2 = .esvd.neg_binom2,
         poisson = .esvd.poisson,
         stop("unsupported distribution family"))
}
