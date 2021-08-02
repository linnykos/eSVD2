
# Compute W'DW, W [m x n], D [m x m] is diagonal, d [m x 1] the diagonal elements of D
wtdw <- function(w, d)
{
  # Broadcast d to each column of w
  dw <- w * d
  crossprod(dw, w)
}

# Objective function value for all the parameters (X, Y, B)
#
# A:     data matrix [n x p]
# X:     latent factor for cells [n x k]
# Y:     latent factor for genes [p x k]
# Z:     covariates [n x r] or NULL
# B:     regression coefficients for Z [p x r] or NULL
# s:     library size vector [n x 1]
# gamma: nuisance parameter vector [p x 1]
objfn_all <- function(X, Y, B, Z, A, family, s, gamma, ...)
{
  # theta = XY'+ZB'
  #
  # We allow Z and B to be NULL, meaning no covariates are in the model
  theta <- tcrossprod(cbind(X, Z), cbind(Y, B))  # [n x p]

  # Call the family function log_prob() to compute entry-wise log-density value
  log_prob <- family$log_prob(A, theta, s, gamma)  # [n x p]
  # Some entries in A may be NAs, and we only aggregate non-missing values
  obs_ind <- which(!is.na(A))
  stopifnot(length(obs_ind) > 0)
  # Compute negative log-likelihood function value
  -mean(log_prob[obs_ind])
}

# Objective function value for the i-th row of X
objfn_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, ...)
{
  # thetai <- drop(cbind(Y, B) %*% c(Xi, Zi))               # [p x 1]
  # log_prob <- family$log_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # # Some entries in Ai may be NAs, and we only aggregate non-missing values
  # obs_ind <- which(!is.na(Ai))
  # stopifnot(length(obs_ind) > 0)
  # # Compute negative log-likelihood function value
  # -mean(log_prob[obs_ind])
  objfn_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma)
}

# Objective function value for the j-th row of Y
objfn_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, ...)
{
  # thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj))               # [n x 1]
  # log_prob <- family$log_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # # Some entries in Aj may be NAs, and we only aggregate non-missing values
  # obs_ind <- which(!is.na(Aj))
  # stopifnot(length(obs_ind) > 0)
  # # Compute negative log-likelihood function value
  # -mean(log_prob[obs_ind])
  objfn_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj)
}

# Gradient for the i-th row of X
grad_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, ...)
{
  # thetai <- drop(cbind(Y, B) %*% c(Xi, Zi))                 # [p x 1]
  # dlog_prob <- family$dlog_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # # Some entries in Ai may be NAs, and we only aggregate non-missing values
  # obs_ind <- which(!is.na(Ai))
  # stopifnot(length(obs_ind) > 0)
  # # Compute the gradient
  # Ysub <- Y[obs_ind, , drop = FALSE]
  # prod <- crossprod(Ysub, dlog_prob[obs_ind])
  # -drop(prod) / length(obs_ind)
  grad_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma)
}

# Gradient for the j-th row of Y
grad_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, ...)
{
  # thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj))                 # [n x 1]
  # dlog_prob <- family$dlog_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # # Some entries in Aj may be NAs, and we only aggregate non-missing values
  # obs_ind <- which(!is.na(Aj))
  # stopifnot(length(obs_ind) > 0)
  # # Compute the gradient
  # Xsub <- X[obs_ind, , drop = FALSE]
  # prod <- crossprod(Xsub, dlog_prob[obs_ind])
  # -drop(prod) / length(obs_ind)
  grad_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj)
}

# Hessian for the i-th row of X
hessian_Xi <- function(Xi, Y, B, Zi, Ai, family, si, gamma, ...)
{
  # thetai <- drop(cbind(Y, B) %*% c(Xi, Zi))                   # [p x 1]
  # d2log_prob <- family$d2log_prob_row(Ai, thetai, si, gamma)  # [p x 1]
  # # Some entries in Ai may be NAs, and we only aggregate non-missing values
  # obs_ind <- which(!is.na(Ai))
  # stopifnot(length(obs_ind) > 0)
  # # Compute the Hessian
  # Ysub <- Y[obs_ind, , drop = FALSE]
  # diag <- d2log_prob[obs_ind]
  # -wtdw(Ysub, diag) / length(obs_ind)
  hessian_Xi_impl(Xi, Y, B, Zi, Ai, family, si, gamma)
}

# Hessian for the j-th row of Y
hessian_Yj <- function(Yj, X, Bj, Z, Aj, family, s, gammaj, ...)
{
  # thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj))                   # [n x 1]
  # d2log_prob <- family$d2log_prob_col(Aj, thetaj, s, gammaj)  # [n x 1]
  # # Some entries in Aj may be NAs, and we only aggregate non-missing values
  # obs_ind <- which(!is.na(Aj))
  # stopifnot(length(obs_ind) > 0)
  # # Compute the Hessian
  # Xsub <- X[obs_ind, , drop = FALSE]
  # diag <- d2log_prob[obs_ind]
  # -wtdw(Xsub, diag) / length(obs_ind)
  hessian_Yj_impl(Yj, X, Bj, Z, Aj, family, s, gammaj)
}

# Feasibility of the i-th row of X
# thetai = Xi * Y' + Zi * B'
feas_Xi <- function(Xi, Y, B, Zi, family, ...)
{
  if(family$feas_always)
    return(TRUE)

  thetai <- drop(cbind(Y, B) %*% c(Xi, Zi))
  family$feasibility(thetai)
}

# Feasibility of the j-th row of Y
# thetaj = X * Yj' + Z * Bj'
feas_Yj <- function(Yj, X, Bj, Z, family, ...)
{
  if(family$feas_always)
    return(TRUE)

  thetaj <- drop(cbind(X, Z) %*% c(Yj, Bj))
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
         poisson = .esvd.poisson,
         stop("unsupported distribution family"))
}
