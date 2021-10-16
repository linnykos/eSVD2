# Optimization for eSVD

# Optimize X given Y and B
opt_x <- function(X0, Y, B, Z, A,
                  family, s, gamma, offset_vec, l2pen,
                  opt_fun,
                  verbose = 0, ...)
{
  if(!is.null(family$cpp_functions) && identical(opt_fun, constr_newton))
  {
    return(opt_x_cpp(X0, Y, B, Z, A, family, s, gamma, offset_vec, l2pen, verbose))
  }

  n <- nrow(A)
  X <- X0
  # Optimize each row of X
  for(i in 1:n)
  {
    if(verbose >= 2)
      cat("===== Optimizing Row ", i, " of X =====\n", sep = "")
    Zi <- if(is.null(Z)) NULL else Z[i, ]

    opt <- opt_fun(
      x0 = X0[i, ],
      f = objfn_Xi,
      gr = grad_Xi,
      hn = hessian_Xi,
      direc = direction_Xi,
      feas = feas_Xi,
      eps_rel = 1e-3,
      verbose = (verbose >= 3),
      Y = Y, B = B, Zi = Zi, Ai = A[i, ],
      family = family, si = s[i], gamma = gamma, offseti = offset_vec[i], l2pen = l2pen, ...
    )

    X[i, ] <- opt$x
    if(verbose >= 3)
      cat("==========\n\n")
  }
  X
}

# Optimize Y and B given X
opt_yb <- function(YB0, XZ, A, family, s, gamma, offset_vec, l2pen, opt_fun, verbose = 0, ...)
{
  if(!is.null(family$cpp_functions) && identical(opt_fun, constr_newton))
  {
    return(opt_yb_cpp(YB0, XZ, A, family, s, gamma, offset_vec, l2pen, verbose))
  }

  p <- ncol(A)
  YB <- YB0
  # Optimize each row of Y and B
  for(j in 1:p)
  {
    if(verbose >= 2)
      cat("===== Optimizing Row ", j, " of Y =====\n", sep = "")

    opt <- opt_fun(
      x0 = YB0[j, ],
      f = objfn_Yj,
      gr = grad_Yj,
      hn = hessian_Yj,
      direc = direction_Yj,
      feas = feas_Yj,
      eps_rel = 1e-3,
      verbose = (verbose >= 3),
      X = XZ, Bj = NULL, Z = NULL, Aj = A[, j],
      family = family, s = s, gammaj = gamma[j], offset = offset_vec, l2pen = l2pen, ...
    )

    YB[j, ] <- opt$x
    if(verbose >= 3)
      cat("==========\n\n")
  }
  YB
}

#######################

#' Main optimization function for eSVD
#'
#' @param x_init             initial estimate of a \eqn{n \times k}{n × k} matrix for the \eqn{k}-dimensional
#'                           embedding for the \eqn{n} cells
#' @param y_init             initial estimate of a \eqn{p \times k}{p × k} matrix for the \eqn{k}-dimensional
#'                           embedding for the \eqn{p} genes
#' @param dat                dataset where the \eqn{n} rows represent cells and \eqn{p} columns represent genes
#' @param family             a character string, one of \code{"gaussian"}, \code{"exponential"},
#'                           \code{"poisson"}, \code{"neg_binom"}, \code{"curved_gaussian"},
#'                           and \code{"bernoulli"}
#' @param method             a character string indicating the optimization method,
#'                           either \code{"newton"} or \code{"lbfgs"}
#' @param b_init             initial estimate of a \eqn{n \times d}{n × d} matrix for the \eqn{d}-dimensional
#'                           embedding for the \eqn{n} cells
#' @param covariates         an \eqn{n \times d}{n × d} matrix representing the additional \eqn{d} covariates,
#'                           or \code{NULL} if no covariate is given
#' @param nuisance_param_vec either \code{NA} or a single numeric or a length-\eqn{p}
#'                           vector of numerics representing nuisance parameters (for \code{family = "neg_binom"}
#'                           or \code{family = "curved_gausian"})
#' @param library_size_vec   either \code{NA} or a single numeric (default is \code{1}) or
#'                           a length-\eqn{n} vector of numerics. If \code{NA}, the library size will be estimated
#' @param offset_vec         a vector of length-\eqn{n} that represents a constant amount added to each row of the
#'                           natural parameter matrix
#' @param l2pen              the ridge penalty parameter for (X, Y, B)
#' @param reparameterize     reparameterize \code{x_mat} and \code{y_mat} after every iteration
#' @param reestimate_nuisance a boolean for whether the nuisance parameter is reestimate after every iteration
#' @param global_estimate    a boolean for whether or not the same nuisance parameter is used for all genes
#' @param min_nuisance       a small positive number for the minimum nuisance parameter value
#' @param max_nuisance       a large positive number for the maximum nuisance parameter value
#' @param max_iter           a positive integer giving the maximum number of iterations
#' @param tol                a small positive number for the tolerance of optimization error
#' @param verbose            a non-negative integer to indicate the verbosity of messages
#' @param ...                additional parameters, currently not used
#'
#' @return A list with elements \code{x_mat} and \code{y_mat} (and others), representing the two
#'         latent matrices.
#' @export
opt_esvd <- function(x_init,
                     y_init,
                     dat,
                     family = "gaussian",
                     method = c("newton", "lbfgs"),
                     b_init = NULL,
                     covariates = NULL,
                     nuisance_param_vec = NA,
                     library_size_vec = NA,
                     offset_vec = rep(0, nrow(x_init)),
                     l2pen = 0,
                     reparameterize = F,
                     reestimate_nuisance = F,
                     global_estimate = F,
                     min_nuisance = 0.01,
                     max_nuisance = 1e5,
                     max_iter = 100,
                     tol = 1e-6,
                     verbose = 0,
                     ...)
{
  n <- nrow(dat)
  p <- ncol(dat)
  k <- ncol(x_init)
  stopifnot(
    nrow(x_init) == n,
    nrow(y_init) == p,
    ncol(y_init) == k,
    is.character(family),
    sum(!is.na(dat)) > 0
  )

  # C++ code needs double type
  storage.mode(x_init) <- "double"
  storage.mode(y_init) <- "double"
  storage.mode(dat) <- "double"

  family <- .string_to_distr_funcs(family)
  library_size_vec <- .parse_library_size(dat, library_size_vec)
  if(all(!is.na(nuisance_param_vec)) && length(nuisance_param_vec) == 1) {
    nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
  }
  library_size_vec <- as.numeric(library_size_vec)
  nuisance_param_vec <- as.numeric(nuisance_param_vec)

  method <- match.arg(method)
  opt_fun <- if(method == "newton") constr_newton else constr_lbfgs

  if(is.null(covariates))
  {
    b_mat <- NULL
  } else {
    if(is.null(b_init))
    {
      r <- ncol(covariates)
      b_mat <- matrix(0, p, r)
    } else {
      b_mat <- b_init
    }
  }
  x_mat <- x_init
  y_mat <- y_init

  losses <- c()
  for(i in 1:max_iter)
  {
    if(verbose >= 1)
      cat("========== eSVD Iter ", i, " ==========\n\n", sep = "")
    # Optimize X given Y and B
    x_mat <- opt_x(X0 = x_mat,
                   Y = y_mat,
                   B = b_mat,
                   Z = covariates,
                   A = dat,
                   family = family,
                   s = library_size_vec,
                   gamma = nuisance_param_vec,
                   offset_vec = offset_vec,
                   l2pen = l2pen,
                   opt_fun = opt_fun,
                   verbose = verbose, ...)

    # Optimize Y and B given X
    yb_mat <- cbind(y_mat, b_mat)
    xz_mat <- cbind(x_mat, covariates)
    yb_mat <- opt_yb(yb_mat,
                     XZ = xz_mat,
                     A = dat,
                     family = family,
                     s = library_size_vec,
                     gamma = nuisance_param_vec,
                     offset_vec = offset_vec,
                     l2pen = l2pen,
                     opt_fun = opt_fun,
                     verbose = verbose, ...)

    # Split Y and B
    if(is.null(covariates))
    {
      y_mat <- yb_mat
      b_mat <- NULL
    } else {
      y_mat <- yb_mat[, 1:k, drop = FALSE]
      b_mat <- yb_mat[, -(1:k), drop = FALSE]
    }

    if(reestimate_nuisance & family$name == "neg_binom2"){
      theta_mat <- tcrossprod(yb_mat, cbind(x_mat, covariates))
      theta_mat <- sweep(theta_mat, 2, offset_vec, "+")

      glmgampoi_res <- glmGamPoi::overdispersion_mle(y = Matrix::t(dat),
                                                     mean = exp(theta_mat),
                                                     global_estimate = global_estimate)

      if(global_estimate){
        nuisance_param_vec <- rep(glmgampoi_res$estimate, p)
        if(verbose >= 1)
          cat("Updated nuisance value: ", round(glmgampoi_res$estimate, 2), "\n", sep = "")
      } else {
        nuisance_param_vec <- glmgampoi_res$estimate
        nuisance_param_vec <- pmin(pmax(nuisance_param_vec, min_nuisance), max_nuisance)
        if(verbose >= 1){
          cat("Updated nuisance values:\n", sep = "")
          print(round(stats::quantile(nuisance_param_vec),2))
        }
      }
    }

    if(reparameterize){
      tmp <- .reparameterize(x_mat, y_mat, equal_covariance = T)
      x_mat <- tmp$x_mat; y_mat <- tmp$y_mat
    }

    # Loss function
    loss <- objfn_all(X = x_mat,
                      Y = y_mat,
                      B = b_mat,
                      Z = covariates,
                      A = dat,
                      family = family,
                      s = library_size_vec,
                      gamma = nuisance_param_vec,
                      offset = offset_vec,
                      l2pen = l2pen)
    losses <- c(losses, loss)
    if(verbose >= 1)
      cat("========== eSVD Iter ", i, ", loss = ", loss, " ==========\n\n", sep = "")

    # Convergence test
    resid <- abs(losses[i] - losses[i - 1])
    thresh <- tol * max(1, abs(losses[i - 1]))
    if(i >=2 && resid <= thresh)
      break
  }

  tmp <- .reparameterize(x_mat, y_mat, equal_covariance = T)
  x_mat <- tmp$x_mat; y_mat <- tmp$y_mat

  list(x_mat = x_mat,
       y_mat = y_mat,
       covariates = covariates,
       b_mat = b_mat,
       loss = losses,
       nuisance_param_vec = nuisance_param_vec,
       library_size_vec = library_size_vec,
       offset_vec = offset_vec)
}
