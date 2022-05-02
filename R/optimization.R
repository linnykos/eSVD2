# Optimization for eSVD

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
#' @param run_cpp            whether to use C++ code
#' @param verbose            a non-negative integer to indicate the verbosity of messages
#' @param ...                additional parameters, currently not used
#'
#' @return A list with elements \code{x_mat} and \code{y_mat} (and others), representing the two
#'         latent matrices.
#' @export
opt_esvd <- function(input_obj, ...) UseMethod("opt_esvd")

#' @export
opt_esvd.eSVD <- function(input_obj,
                          l2pen = 0.1,
                          max_iter = 100,
                          method = c("newton", "lbfgs"),
                          offset_variables = NULL,
                          tol = 1e-6,
                          verbose = 0,
                          fit_name = "fit_First",
                          fit_previous = "fit_Init",
                          ...){
  dat <- .get_object(eSVD_obj = input_obj, what_obj = "dat", which_fit = NULL)
  covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
  x_mat <- .get_object(eSVD_obj = input_obj, what_obj = "x_mat", which_fit = fit_previous)
  y_mat <- .get_object(eSVD_obj = input_obj, what_obj = "y_mat", which_fit = fit_previous)
  z_mat <- .get_object(eSVD_obj = input_obj, what_obj = "z_mat", which_fit = fit_previous)

  param <- .opt_esvd_format_param(family = "poisson",
                                  l2pen = l2pen,
                                  max_iter = max_iter,
                                  method = method,
                                  offset_variables = offset_variables,
                                  tol = tol,
                                  prefix = paste0(fit_name, "_"))
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)

  res <- opt_esvd.default(input_obj = dat,
                          x_init = x_mat,
                          y_init = y_mat,
                          z_init = z_mat,
                          covariates = covariates,
                          family = "poisson",
                          l2pen = l2pen,
                          max_iter = max_iter,
                          method = method[1],
                          offset_variables = offset_variables,
                          tol = tol,
                          verbose = verbose, ...)

  input_obj[[fit_name]] <- .form_esvd_fit(
    x_mat = res$x_mat,
    y_mat = res$y_mat,
    z_mat = res$z_mat,
    loss = res$loss
  )
  input_obj[["latest_Fit"]] <- fit_name
  input_obj
}

#' @export
opt_esvd.default <- function(input_obj,
                             x_init,
                             y_init,
                             z_init = NULL,
                             covariates = NULL,
                             family = "poisson",
                             l2pen = 0.1,
                             max_iter = 100,
                             method = c("newton", "lbfgs"),
                             offset_variables = NULL,
                             tol = 1e-6,
                             verbose = 0,
                             ...)
{
  n <- nrow(input_obj)
  p <- ncol(input_obj)
  k <- ncol(x_init)
  stopifnot(
    nrow(x_init) == n, nrow(y_init) == p, ncol(y_init) == k,
    is.character(family), sum(!is.na(input_obj)) > 0,
    method %in% c("newton", "lbfgs")
  )
  if(!all(is.null(offset_variables))){
    stopifnot(is.character(offset_variables),
              all(offset_variables %in% colnames(covariates)))
  }

  # Convert family string to internal family object, e.g. `.esvd.poisson` and `.esvd.neg_binom2`
  family_str <- family
  family <- esvd_family(family_str)
  param <- .opt_esvd_format_param(family = family_str,
                                  l2pen = l2pen,
                                  max_iter = max_iter,
                                  method = method,
                                  offset_variables = offset_variables,
                                  tol = tol)

  # Parse optimization method
  method <- match.arg(method, choice = c("newton", "lbfgs"))
  opt_fun <- if(method == "newton") constr_newton else constr_lbfgs
  loader <- data_loader(input_obj)

  # Initialize embedding matrices
  z_mat <- .opt_esvd_setup_z_mat(covariates = covariates,
                                 p = ncol(input_obj),
                                 z_init = z_init)

  xc_mat <- cbind(x_init, covariates)
  yz_mat <- cbind(y_init, z_mat)
  fixed_cols <- which(colnames(yz_mat) %in% offset_variables)

  losses <- c()
  for(i in seq_len(max_iter))
  {
    if(verbose >= 1) cat("========== eSVD Iter ", i, " ==========\n\n", sep = "")
    # Optimize X given Y and B
    # [[note to self: When all the coding is done, fix the notation -- Z should be the coefficients (currently B), while C should be the covariates (currently)]]
    xc_mat <- opt_x(
      XC_init = xc_mat,
      YZ = yz_mat,
      k = k,
      loader = loader,
      family = family,
      s = rep(1, n),
      gamma = rep(NA, p),
      l2penx = l2pen,
      verbose = verbose)

    # Optimize Y and B given X
    yz_mat <- opt_yz(
      YZ_init = yz_mat,
      XC = xc_mat,
      k = k,
      fixed_cols = fixed_cols,
      loader = loader,
      family = family,
      s = rep(1, n),
      gamma = rep(NA, p),
      l2peny = l2pen,
      l2penz = l2pen,
      verbose = verbose)

    # Loss function
    loss <- objfn_all_r(
      XC = xc_mat,
      YZ = yz_mat,
      k = k,
      loader = loader,
      family = family,
      s = rep(1, n),
      gamma = rep(NA, p),
      l2penx = l2pen,
      l2peny = l2pen,
      l2penz = l2pen
    )

    losses <- c(losses, loss)
    if(verbose >= 1) cat("========== eSVD Iter ", i, ", loss = ", loss, " ==========\n\n", sep = "")

    # Convergence test
    if(i >= 2){
      resid <- abs(losses[i] - losses[i - 1])
      thresh <- tol * max(1, abs(losses[i - 1]))
      if(resid <= thresh) break()
    }
  }

  x_mat <- xc_mat[,1:k, drop = F]
  y_mat <- yz_mat[,1:k, drop = F]
  if(k <= ncol(yz_mat)){
    z_mat <- yz_mat[,(k+1):ncol(yz_mat), drop = F]
  }
  tmp <- tryCatch(.reparameterize(x_mat, y_mat, equal_covariance = T),
                  error = function(e){list(x_mat = x_mat, y_mat = y_mat)})
  x_mat <- tmp$x_mat
  y_mat <- tmp$y_mat

  tmp <- .opt_esvd_format_matrices(covariates = covariates,
                                   dat = input_obj,
                                   x_mat = x_mat,
                                   y_mat = y_mat,
                                   z_mat = z_mat)

  list(x_mat = x_mat,
       y_mat = y_mat,
       covariates = covariates,
       z_mat = z_mat,
       loss = losses,
       param = param)
}
