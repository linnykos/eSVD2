#' Optimize eSVD
#'
#' Generic function interface
#'
#' @param input_obj Main object
#' @param ...       Additional parameters
#'
#' @return Output dependent on class of \code{input_obj}
#' @export
opt_esvd <- function(input_obj, ...) {UseMethod("opt_esvd")}

#' Optimize eSVD for eSVD objects
#'
#' @param input_obj         \code{eSVD} object outputed from \code{apply_initial_threshold}.
#' @param fit_name          String for the name of that will become the current fit when
#'                          storing the results in \code{input_obj}.
#' @param fit_previous      String for the name of the previous fit that this function will
#'                          grab the initialization values from.
#' @param l2pen             Small positive number for the amount of penalization for both the cells'
#'                          and the genes' latent vectors as well as the coefficients.
#' @param max_iter          Positive integer for number of iterations.
#' @param offset_variables  A vector of strings depicting which column names in \code{input_obj$covariate}
#'                          be treated as an offset during the optimization (i.e., their coefficients will not change
#'                          throughout the optimization).
#' @param tol               Small positive number to differentiate between zero and non-zero.
#' @param verbose           Integer.
#' @param ...               Additional parameters.
#'
#' @return \code{eSVD} object with added elements with name to whatever
#' \code{fit_name} was set to.
#' @export
opt_esvd.eSVD <- function(input_obj,
                          fit_name = "fit_First",
                          fit_previous = "fit_Init",
                          l2pen = 0.1,
                          max_iter = 100,
                          offset_variables = NULL,
                          tol = 1e-6,
                          verbose = 0,
                          ...){
  dat <- .get_object(eSVD_obj = input_obj, what_obj = "dat", which_fit = NULL)
  covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
  x_mat <- .get_object(eSVD_obj = input_obj, what_obj = "x_mat", which_fit = fit_previous)
  y_mat <- .get_object(eSVD_obj = input_obj, what_obj = "y_mat", which_fit = fit_previous)
  z_mat <- .get_object(eSVD_obj = input_obj, what_obj = "z_mat", which_fit = fit_previous)

  param <- .opt_esvd_format_param(family = "poisson",
                                  l2pen = l2pen,
                                  max_iter = max_iter,
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


#' Optimize eSVD for matrices or sparse matrices.
#'
#' @param input_obj          Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
#'                           and \eqn{p} columns represent genes.
#'                           The rows and columns of the matrix should be named.
#' @param x_init             Initial matrix of the cells' latent vectors that is \eqn{n} rows and \eqn{k}
#'                           columns. The row names should be the same as \code{input_obj}.
#' @param y_init             Initial matrix of the genes' latent vectors that is \eqn{p} rows and \eqn{k}
#'                           columns. The row names should be the same as the column names of \code{input_obj}.
#' @param z_init             Initial matrix of the genes' coefficient vectors that is \eqn{p} rows and \code{ncol(covariates)}
#'                           columns. The row names should be the same as the column names of \code{input_obj},
#'                           and the column names should be the same as \code{covariates}.
#' @param covariates         \code{matrix} object with \eqn{n} rows with the same rownames as \code{input_obj} where the columns
#'                           represent the different covariates.
#'                           Notably, this should contain only numerical columns (i.e., all categorical
#'                           variables should have already been split into numerous indicator variables).
#' @param family             String among \code{"gaussian"}, \code{"curved_gaussian"},
#'                           \code{"exponential"}, \code{"poisson"}, \code{"neg_binom"},
#'                           \code{"neg_binom2"}, or \code{"bernoulli"}. Notably, with exception of
#'                           \code{"neg_binom2"}, all the other families are parameterized such that
#'                           eSVD is fitting the dot product to be the canonical parameter of these
#'                           expoential-family distributions. For \code{"neg_binom2"}, the dot
#'                           product is the log-mean of the distribution (i.e., similar to the canonical
#'                           parameterization of the Poisson family).
#' @param l2pen              Small positive number for the amount of penalization for both the cells'
#'                           and the genes' latent vectors as well as the coefficients.
#' @param library_multipler  Vector of positive numerics of length \eqn{n}. It is the multiplier
#'                           such that the variance of cell \code{i}'s entries is the mean of
#'                           cell \code{i}'s entries times the square-root of cell \code{i}'s
#'                           value in \code{library_multipler} (entry-wise). This is used as
#'                           an alternative interpretation of how library-size affects a cell's
#'                           gene expression (instead of using the library size as a covariate to be
#'                           regressed out).
#' @param max_iter           Positive integer for number of iterations.
#' @param nuisance_vec       Vector of non-negative numerics (or \code{NA}'s) of length \eqn{p},
#'                           representing each gene's nuisance parameter when using an exponential-family
#'                           distribution that requires one.
#'                           It is used only when \code{family} is  \code{"curved_gaussian"} or
#'                           \code{"neg_binom"} or \code{"neg_binom2"}.
#' @param offset_variables   A vector of strings depicting which column names in \code{input_obj$covariate}
#'                           be treated as an offset during the optimization (i.e., their coefficients will not change
#'                           throughout the optimization).
#' @param tol                Small positive number to differentiate between zero and non-zero.
#' @param verbose            Integer.
#' @param ...                Additional parameters.
#'
#' @return a \code{list} with elements \code{x_mat}, \code{y_mat},
#' \code{z_mat}, \code{library_multiplier}, \code{loss}, \code{nuisance_vec}
#' and \code{param}.
#' @export
opt_esvd.default <- function(input_obj,
                             x_init,
                             y_init,
                             z_init = NULL,
                             covariates = NULL,
                             family = "poisson",
                             l2pen = 0.1,
                             library_multipler = rep(1, nrow(input_obj)),
                             max_iter = 100,
                             nuisance_vec = rep(NA, ncol(input_obj)),
                             offset_variables = NULL,
                             tol = 1e-6,
                             verbose = 0,
                             ...)
{
  n <- nrow(input_obj)
  p <- ncol(input_obj)
  k <- ncol(x_init)
  stopifnot(
    inherits(input_obj, c("matrix", "dgCMatrix")),
    nrow(x_init) == n, nrow(y_init) == p, ncol(y_init) == k,
    is.character(family), sum(!is.na(input_obj)) > 0
  )
  if(!all(is.null(offset_variables))){
    stopifnot(is.character(offset_variables),
              all(offset_variables %in% colnames(covariates)))
  }

  # Convert family string to internal family object
  family_str <- as.character(family)
  family <- esvd_family(family_str)
  param <- .opt_esvd_format_param(family = family_str,
                                  l2pen = l2pen,
                                  max_iter = max_iter,
                                  offset_variables = offset_variables,
                                  tol = tol)

  # Load the data
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
    # Optimize X given C, Y, and Z
    xc_mat <- opt_x(
      XC_init = xc_mat,
      YZ = yz_mat,
      k = k,
      loader = loader,
      family = family,
      s = library_multipler,
      gamma = nuisance_vec,
      l2penx = l2pen,
      verbose = verbose)

    # Optimize Y and Z given X
    yz_mat <- opt_yz(
      YZ_init = yz_mat,
      XC = xc_mat,
      k = k,
      fixed_cols = fixed_cols,
      loader = loader,
      family = family,
      s = library_multipler,
      gamma = nuisance_vec,
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
      s = library_multipler,
      gamma = nuisance_vec,
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

  list(x_mat = tmp$x_mat,
       y_mat = tmp$y_mat,
       covariates = covariates,
       z_mat = tmp$z_mat,
       library_multipler = library_multipler,
       loss = losses,
       nuisance_vec = nuisance_vec,
       param = param)
}
