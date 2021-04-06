# Optimization for eSVD

# Optimize u given v
opt_u_given_v <- function(x0, y_mat, dat, opt_fun,
                          family_funcs, nuisance_param_vec,
                          library_size_vec, verbose = 0, ...)
{
  n <- nrow(dat)
  p <- ncol(dat)
  x_mat <- x0
  for(i in 1:n)
  {
    if(verbose >= 2)
      cat("===== Optimizing Row ", i, " of U =====\n", sep = "")
    opt <- opt_fun(
      x0[i, ], family_funcs$objfn, family_funcs$grad, family_funcs$hessian, family_funcs$feas,
      eps_rel = 1e-3, verbose = (verbose >= 3),
      other_mat = y_mat, dat_vec = dat[i, ],
      nuisance_param_vec = nuisance_param_vec,
      library_size = library_size_vec[i], ...
    )
    x_mat[i, ] <- opt$x
    if(verbose >= 3)
      cat("==========\n\n")
  }
  x_mat
}

# Optimize v given u
opt_v_given_u <- function(y0, x_mat, dat, opt_fun,
                          family_funcs, nuisance_param_vec,
                          library_size_vec, verbose = 0, ...)
{
  n <- nrow(dat)
  p <- ncol(dat)
  y_mat <- y0
  for(j in 1:p)
  {
    if(verbose >= 2)
      cat("===== Optimizing Row ", j, " of V =====\n", sep = "")
    opt <- opt_fun(
      y0[j, ], family_funcs$objfn, family_funcs$grad, family_funcs$hessian, family_funcs$feas,
      eps_rel = 1e-3, verbose = (verbose >= 3),
      other_mat = x_mat, dat_vec = dat[, j],
      nuisance_param_vec = nuisance_param_vec[j],
      library_size = library_size_vec, ...
    )
    y_mat[j, ] <- opt$x
    if(verbose >= 3)
      cat("==========\n\n")
  }
  y_mat
}

#######################

#' Main optimization function for eSVD
#'
#' @param x_init initial estimate of a \eqn{n \times k} matrix for the \eqn{k}-dimensional
#' embedding for the \eqn{n} cells
#' @param y_init initial estimate of a \eqn{p \times k} matrix for the \eqn{k}-dimensional
#' embedding for the \eqn{p} genes
#' @param dat dataset where the \eqn{n} rows represent cells and \eqn{p} columns represent genes
#' @param family A character string, one of \code{"gaussian"}, \code{"exponential"},
#'               \code{"poisson"}, \code{"neg_binom"}, \code{"curved_gaussian"},
#'               and \code{"bernoulli"}.
#' @param method A character string indicating the optimization method,
#'               either \code{"newton"} or \code{"lbfgs"}.
#' @param nuisance_param_vec either \code{NA} or a single numeric or a length-\eqn{p}
#' vector of numerics representing nuisance parameters (for \code{family = "neg_binom"} or
#' \code{family = "curved_gausian"}).
#' @param library_size_vec either \code{NA} or a single numeric (default is \code{1}) or
#' a length-\eqn{n} vector of numerics.
#' If \code{NA}, the library size will be estimated.
#' @param max_iter positive integer
#' @param verbose no-negtaive integer
#' @param ... additional parameters, not used
#'
#' @return a list with elements \code{x_mat} and \code{y_mat} (and others), representing the two
#' latent matrices
#' @export
opt_esvd <- function(x_init, y_init, dat, family, method = c("newton", "lbfgs"),
                     nuisance_param_vec = NA, library_size_vec = 1,
                     max_iter = 100, tol = 1e-6,
                     verbose = 0, ...)
{
  stopifnot(
    nrow(x_init) == nrow(dat),
    nrow(y_init) == ncol(dat),
    ncol(x_init) == ncol(y_init),
    is.character(family),
    length(which(!is.na(dat))) > 0
  )

  n <- nrow(dat); p <- nrow(dat); k <- ncol(x_init)
  family_funcs <- .string_to_distr_funcs(family)
  library_size_vec <- .parse_library_size(dat, library_size_vec)

  if(all(!is.na(nuisance_param_vec)) & length(nuisance_param_vec) == 1) {
    nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
  }

  method <- match.arg(method)
  opt_fun <- if(method == "newton") constr_newton else constr_lbfgs

  x_mat <- x_init
  y_mat <- y_init
  losses <- c()
  for(i in 1:max_iter)
  {
    if(verbose >= 1)
      cat("========== eSVD Iter ", i, " ==========\n\n", sep = "")
    # Optimize u given v
    x_mat <- opt_u_given_v(x_mat, y_mat, dat, opt_fun,
                           family_funcs, nuisance_param_vec,
                           library_size_vec, verbose, ...)
    # Orthogonalize u
    # x_mat = sqrt(n) * svd(x_mat)$u
    # Optimize v given u
    y_mat <- opt_v_given_u(y_mat, x_mat, dat, opt_fun,
                           family_funcs, nuisance_param_vec,
                           library_size_vec, verbose, ...)
    # Loss function
    loss <- family_funcs$objfn_all(dat, x_mat, y_mat, nuisance_param_vec,
                                   library_size_vec, ...)
    losses <- c(losses, loss)
    if(verbose >= 1)
      cat("========== eSVD Iter ", i, ", loss = ", loss, " ==========\n\n", sep = "")
    # Orthogonalize v
    # y_mat = sqrt(p) * svd(y_mat)$u

    # Convergence test
    resid <- abs(losses[i] - losses[i - 1])
    thresh <- tol * max(1, abs(losses[i - 1]))
    if(i >=2 && resid <= thresh)
      break
  }

  list(x_mat = x_mat, y_mat = y_mat, loss = losses,
       nuisance_param_vec = nuisance_param_vec, library_size_vec = library_size_vec)
}
