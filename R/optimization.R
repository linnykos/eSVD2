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
#' @param input_obj         \code{eSVD} object output from \code{initialize_esvd}
#'                          (or an earlier \code{opt_esvd}).
#' @param fit_name          String for the name of that will become the current fit when
#'                          storing the results in \code{input_obj}.
#' @param fit_previous      String for the name of the previous fit that this function will
#'                          grab the initialization values from.
#' @param l2pen             Small positive number for the amount of penalization for both the cells'
#'                          and the genes' latent vectors as well as the coefficients.
#' @param max_iter          Positive integer for number of iterations.
#' @param offset_variables  A vector of strings depicting which column names in \code{input_obj$covariates}
#'                          be treated as an offset during the optimization (i.e., their coefficients will not change
#'                          throughout the optimization).
#' @param tol               Relative tolerance of the stopping rule: iteration stops once the objective
#'                          changes by at most \code{tol} times \code{max(1, |previous objective|)}.
#' @param verbose           Integer.
#' @param ...               Additional parameters.
#'
#' @return \code{eSVD} object with added elements with name to whatever
#' \code{fit_name} was set to.
#' @examples
#' set.seed(10)
#' sim <- generate_null(cell_per_person = 15, num_genes = 40,
#'                      num_individuals = 8)
#' esvd_obj <- initialize_esvd(dat = sim$obs_mat,
#'                             covariates = sim$covariates,
#'                             metadata_individual = sim$metadata_individual,
#'                             case_control_variable = "CC",
#'                             bool_intercept = TRUE,
#'                             k = 2,
#'                             lambda = 0.1)
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_Init",
#'                                                omitted_variables = "Log_UMI")
#' esvd_obj <- opt_esvd(input_obj = esvd_obj,
#'                      max_iter = 5,
#'                      offset_variables = setdiff(colnames(esvd_obj$covariates), "CC"),
#'                      fit_name = "fit_First",
#'                      fit_previous = "fit_Init")
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_First",
#'                                                omitted_variables = "Log_UMI")
#' esvd_obj$latest_Fit
#' esvd_obj$fit_First$loss
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
#'                           eSVD is fitting the dot product to be the natural parameter of these
#'                           exponential-family distributions (see the "Natural parameter" line at the
#'                           top of each \code{src/family_*.cpp}). For \code{"neg_binom2"}, the dot
#'                           product is the log-mean of the distribution (i.e., similar to the canonical
#'                           parameterization of the Poisson family).
#' @param l2pen              Small positive number for the amount of penalization for both the cells'
#'                           and the genes' latent vectors as well as the coefficients.
#' @param library_multipler  Vector of positive numerics of length \eqn{n}: the per-cell multiplier
#'                           \eqn{s_i} in the likelihood, entering by family (for example, the
#'                           Poisson mean is \eqn{s_i e^{\theta_{ij}}}, and the Gaussian has mean
#'                           \eqn{s_i\theta_{ij}} and variance \eqn{s_i\gamma_j^2});
#'                           \code{"neg_binom2"} ignores it. This is used as
#'                           an alternative interpretation of how library-size affects a cell's
#'                           gene expression (instead of using the library size as a covariate to be
#'                           regressed out).
#' @param max_iter           Positive integer for number of iterations.
#' @param nuisance_vec       Vector of positive numerics of length \eqn{p},
#'                           representing each gene's nuisance parameter when using an exponential-family
#'                           distribution that requires one: the standard deviation for
#'                           \code{"gaussian"}, the mean divided by the standard deviation (the
#'                           inverse coefficient of variation) for \code{"curved_gaussian"}, and the size (number of failures) for
#'                           \code{"neg_binom"} and \code{"neg_binom2"}. It is ignored by
#'                           \code{"poisson"}, \code{"exponential"} and \code{"bernoulli"}.
#'                           The default \code{NULL} uses \code{1} for every gene, which is a
#'                           placeholder rather than an estimate; supply your own values for
#'                           the families that use it.
#' @param offset_variables   A vector of strings depicting which column names in \code{input_obj$covariates}
#'                           be treated as an offset during the optimization (i.e., their coefficients will not change
#'                           throughout the optimization).
#' @param tol                Relative tolerance of the stopping rule: iteration stops once the objective
#'                           changes by at most \code{tol} times \code{max(1, |previous objective|)}.
#' @param verbose            Integer
#' @param ...                Additional parameters
#'
#' @return a \code{list} with elements \code{x_mat}, \code{y_mat},
#' \code{z_mat}, \code{covariates}, \code{library_multipler}, \code{loss}
#' (the objective after every iteration), \code{nuisance_vec} and
#' \code{param}. A warning is raised if any row or column update ended in a
#' failed line search.
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
                             nuisance_vec = NULL,
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
  # Four of the seven families consume `nuisance_vec` inside the objective,
  # and an NA there used to surface as "missing value where TRUE/FALSE
  # needed" from the line search. A default of 1 lets every family run at
  # its defaults; it is documented as a placeholder, not an estimate.
  if(anyNA(x_init) || anyNA(y_init) || (!is.null(z_init) && anyNA(z_init)) ||
     (!is.null(covariates) && anyNA(covariates))){
    stop("`x_init`, `y_init`, `z_init` and `covariates` must not contain NA")
  }
  if(is.null(nuisance_vec)) nuisance_vec <- rep(1, p)
  if(length(nuisance_vec) != p || !is.numeric(nuisance_vec) ||
     any(!is.finite(nuisance_vec)) || any(nuisance_vec <= 0)){
    stop("`nuisance_vec` must be a vector of ", p,
         " finite positive numerics (one per gene), or NULL")
  }
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
  # Line-search failures are counted in C++ and handed back as an attribute
  # (a C++ frame must not raise an R warning; see constrained_newton.cpp).
  # They are summed over all iterations and reported once, below.
  num_linesearch_failed <- 0
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
    num_linesearch_failed <- num_linesearch_failed +
      .attr_or_zero(xc_mat, "num_linesearch_failed")

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
    num_linesearch_failed <- num_linesearch_failed +
      .attr_or_zero(yz_mat, "num_linesearch_failed")

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

    if(!is.finite(loss)){
      stop("the eSVD objective became non-finite at iteration ", i,
           " (loss = ", loss, "). The fit has diverged; a smaller step ",
           "(larger `l2pen`) or a better initialization is needed, and ",
           "collinear covariates are a common cause")
    }
    losses <- c(losses, loss)
    if(verbose >= 1) cat("========== eSVD Iter ", i, ", loss = ", loss, " ==========\n\n", sep = "")

    # Convergence test
    if(i >= 2){
      resid <- abs(losses[i] - losses[i - 1])
      thresh <- tol * max(1, abs(losses[i - 1]))
      if(verbose >= 2) {
        print(paste0("Residual of loss: ", resid))
        print(paste0("Threshold for termination: ", thresh))
      }
      if(resid <= thresh) break()
    }
  }

  if(num_linesearch_failed > 0){
    warning("the Newton line search failed for ", num_linesearch_failed,
            " row/column update(s) over ", length(losses), " iteration(s), ",
            "leaving those rows or columns at their last accepted step; the ",
            "fit may not have converged")
  }

  # Subsetting drops the C++ attribute, so it does not leak into the fit.
  x_mat <- xc_mat[,1:k, drop = FALSE]
  y_mat <- yz_mat[,1:k, drop = FALSE]
  if(k < ncol(yz_mat)){
    z_mat <- yz_mat[,(k+1):ncol(yz_mat), drop = FALSE]
  }
  tmp <- tryCatch(.reparameterize(x_mat, y_mat, equal_covariance = TRUE),
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
