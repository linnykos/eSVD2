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
# l2pen:  ridge penalty [1]
objfn_all <- function(X, Y, B, Z, A, family, s, gamma, offset_mat, l2pen, ...)
{
  # theta = XY'+ZB' + offset
  #
  # We allow Z and B to be NULL, meaning no covariates are in the model
  XZ <- cbind(X, Z)
  YB <- cbind(Y, B)
  theta <- tcrossprod(XZ, YB)  # [n x p]
  theta <- theta + offset_mat

  # Call the family function log_prob() to compute entry-wise log-density value
  log_prob <- family$log_prob(A, theta, s, gamma)  # [n x p]
  # Some entries in A may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(A))
  stopifnot(length(obs_ind) > 0)
  # Compute negative log-likelihood function value
  nll <- -mean(log_prob[obs_ind])
  # Ridge penalty
  ridge <- if(l2pen > 0) l2pen / length(obs_ind) * (sum(X^2) + sum(YB^2)) else 0
  nll + ridge
}

# Objective function value for the i-th row of X
objfn_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, offset, l2pen, ...)
{
  # Call the C++ version if the given family has C++ implementation
  if(!is.null(family$cpp_functions))
  {
    storage.mode(Xi) <- storage.mode(Y) <- storage.mode(Ai) <- storage.mode(gamma) <- "double"
    return(objfn_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, l2pen))
  }

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offset     # [p x 1]
  log_prob <- family$log_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # Some entries in Ai may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Ai))
  stopifnot(length(obs_ind) > 0)
  # Compute negative log-likelihood function value
  nll <- -mean(log_prob[obs_ind])
  # Ridge penalty
  ridge <- if(l2pen > 0) l2pen * sum(Xi^2) / length(obs_ind) else 0
  nll + ridge
}

# Objective function value for the j-th row of Y
objfn_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, l2pen, ...)
{
  # Call the C++ version if the given family has C++ implementation
  if(!is.null(family$cpp_functions))
  {
    storage.mode(Yj) <- storage.mode(X) <- storage.mode(Aj) <- storage.mode(s) <- storage.mode(offset) <- "double"
    return(objfn_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, l2pen))
  }

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset      # [n x 1]
  log_prob <- family$log_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # Some entries in Aj may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Aj))
  stopifnot(length(obs_ind) > 0)
  # Compute negative log-likelihood function value
  nll <- -mean(log_prob[obs_ind])
  # Ridge penalty
  ridge <- if(l2pen > 0) l2pen * sum(Yj^2) / length(obs_ind) else 0
  nll + ridge
}

# Gradient for the i-th row of X
grad_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, offset, l2pen, ...)
{
  # Call the C++ version if the given family has C++ implementation
  if(!is.null(family$cpp_functions))
  {
    storage.mode(Xi) <- storage.mode(Y) <- storage.mode(Ai) <- storage.mode(gamma) <- "double"
    return(grad_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, l2pen))
  }

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offset       # [p x 1]
  dlog_prob <- family$dlog_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # Some entries in Ai may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Ai))
  stopifnot(length(obs_ind) > 0)
  # Compute the gradient
  Ysub <- Y[obs_ind, , drop = FALSE]
  prod <- crossprod(Ysub, dlog_prob[obs_ind])
  (-drop(prod) + 2 * l2pen * Xi) / length(obs_ind)
}

# Gradient for the j-th row of Y
grad_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, l2pen, ...)
{
  # Call the C++ version if the given family has C++ implementation
  if(!is.null(family$cpp_functions))
  {
    storage.mode(Yj) <- storage.mode(X) <- storage.mode(Aj) <- storage.mode(s) <- storage.mode(offset) <- "double"
    return(grad_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, l2pen))
  }

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset        # [n x 1]
  dlog_prob <- family$dlog_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # Some entries in Aj may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Aj))
  stopifnot(length(obs_ind) > 0)
  # Compute the gradient
  Xsub <- X[obs_ind, , drop = FALSE]
  prod <- crossprod(Xsub, dlog_prob[obs_ind])
  (-drop(prod) + 2 * l2pen * Yj) / length(obs_ind)
}

# Hessian for the i-th row of X
hessian_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, offset, l2pen, ...)
{
  # Call the C++ version if the given family has C++ implementation
  if(!is.null(family$cpp_functions))
  {
    storage.mode(Xi) <- storage.mode(Y) <- storage.mode(Ai) <- storage.mode(gamma) <- "double"
    return(hessian_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, l2pen))
  }

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offset         # [p x 1]
  d2log_prob <- family$d2log_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # Some entries in Ai may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Ai))
  stopifnot(length(obs_ind) > 0)
  # Compute the Hessian
  Ysub <- Y[obs_ind, , drop = FALSE]
  d2 <- d2log_prob[obs_ind]
  hess <- -wtdw(Ysub, d2)
  diag(hess) <- diag(hess) + 2 * l2pen
  hess / length(obs_ind)
}

# Hessian for the j-th row of Y
hessian_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, l2pen, ...)
{
  # Call the C++ version if the given family has C++ implementation
  if(!is.null(family$cpp_functions))
  {
    storage.mode(Yj) <- storage.mode(X) <- storage.mode(Aj) <- storage.mode(s) <- storage.mode(offset) <- "double"
    return(hessian_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, l2pen))
  }

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset          # [n x 1]
  d2log_prob <- family$d2log_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # Some entries in Aj may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Aj))
  stopifnot(length(obs_ind) > 0)
  # Compute the Hessian
  Xsub <- X[obs_ind, , drop = FALSE]
  d2 <- d2log_prob[obs_ind]
  hess <- -wtdw(Xsub, d2)
  diag(hess) <- diag(hess) + 2 * l2pen
  hess / length(obs_ind)
}

# Move direction for the i-th row of X, d = -inv(H) * g
direction_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, offset, l2pen, ...)
{
  # Call the C++ version if the given family has C++ implementation
  if(!is.null(family$cpp_functions))
  {
    storage.mode(Xi) <- storage.mode(Y) <- storage.mode(Ai) <- storage.mode(gamma) <- "double"
    return(direction_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma, offseti, l2pen))
  }

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offset          # [p x 1]
  dlog_prob <- family$dlog_prob_row(Ai, thetai, si, gamma)    # [p x 1]
  d2log_prob <- family$d2log_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # Some entries in Ai may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Ai))
  stopifnot(length(obs_ind) > 0)
  # Compute the Hessian and gradient
  Ysub <- Y[obs_ind, , drop = FALSE]
  d2 <- d2log_prob[obs_ind]
  hess <- -wtdw(Ysub, d2)
  diag(hess) <- diag(hess) + 2 * l2pen
  g <- -crossprod(Ysub, dlog_prob[obs_ind])
  g <- drop(g) + 2 * l2pen * Xi
  # Fall back to gradient direction if Hessian is singular
  direc <- tryCatch(-solve(hess, g), error = function(e) -g / length(obs_ind))

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
direction_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, l2pen, ...)
{
  # Call the C++ version if the given family has C++ implementation
  if(!is.null(family$cpp_functions))
  {
    storage.mode(Yj) <- storage.mode(X) <- storage.mode(Aj) <- storage.mode(s) <- storage.mode(offset) <- "double"
    return(direction_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj, offset, l2pen))
  }

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj)) + offset          # [n x 1]
  dlog_prob <- family$dlog_prob_col(Aj, thetaj, s, gammaj)    # [n x 1]
  d2log_prob <- family$d2log_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # Some entries in Aj may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(Aj))
  stopifnot(length(obs_ind) > 0)
  # Compute the Hessian and gradient
  Xsub <- X[obs_ind, , drop = FALSE]
  d2 <- d2log_prob[obs_ind]
  hess <- -wtdw(Xsub, d2)
  diag(hess) <- diag(hess) + 2 * l2pen
  g <- -crossprod(Xsub, dlog_prob[obs_ind])
  g <- drop(g) + 2 * l2pen * Yj
  # Fall back to gradient direction if Hessian is singular
  direc <- tryCatch(-solve(hess, g), error = function(e) -g / length(obs_ind))

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
feas_Xi <- function(Xi, Y, B, Zi, family, offset, ...)
{
  if(family$feas_always)
    return(TRUE)

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi)) + offset
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
