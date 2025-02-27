#' Optimize X given C, Y and Z
#'
#' @param XC_init Initial value for the `[X C]` matrix, X [n x k], C [n x r]
#' @param YZ      The `[Y Z]` matrix, Y [p x k], Z [p x r]
#' @param k       Number of columns in X and Y
#' @param loader  The data loader, typically returned by data_loader()
#' @param family  A family object, typically returned by esvd_family()
#' @param s       The library size vector, [n x 1]
#' @param gamma   The nuisance parameter vector, [p x 1]
#' @param l2penx  The l2 penalty parameter for X, a scalar
#' @param verbose Verbosity parameter
#' @param inplace Whether the input XC_init will be modified and returned
#' @param ...     Additional parameters
#'
opt_x <- function(XC_init, YZ, k, loader, family, s, gamma, l2penx,
                  verbose = 0, inplace = FALSE, ...)
{
  storage.mode(YZ) <- "double"
  .opt_x(XC_init, YZ, k, loader, family, s, gamma, l2penx, verbose, inplace)
}

#' Optimize Y and Z given X and C
#'
#' @param YZ_init    Initial value for the `[Y Z]` matrix, Y [p x k], Z [p x r]
#' @param XC         The `[X C]` matrix, X [n x k], C [n x r]
#' @param k          Number of columns in X and Y
#' @param fixed_cols Which columns in YZ need to be fixed
#' @param loader     The data loader, typically returned by data_loader()
#' @param family     A family object, typically returned by esvd_family()
#' @param s          The library size vector, [n x 1]
#' @param gamma      The nuisance parameter vector, [p x 1]
#' @param l2peny     The l2 penalty parameter for Y, a scalar
#' @param l2penz     The l2 penalty parameter for Z, a scalar
#' @param verbose    Verbosity parameter
#' @param inplace    Whether the input XC_init will be modified and returned
#' @param ...     Additional parameters
#'
opt_yz <- function(YZ_init, XC, k, fixed_cols, loader, family, s, gamma,
                   l2peny, l2penz, verbose = 0, inplace = FALSE, ...)
{
  # YZind will be passed to C++, should be zero-based
  YZind <- setdiff(1:ncol(YZ_init), fixed_cols) - 1
  storage.mode(XC) <- "double"
  .opt_yz(YZ_init, XC, k, YZind, loader, family, s, gamma, l2peny, l2penz,
          verbose, inplace)
}

##########################

.opt_esvd_format_param <- function(family,
                                   l2pen,
                                   max_iter,
                                   offset_variables,
                                   tol,
                                   prefix = "") {
  res <- list(family = family,
              l2pen = l2pen,
              max_iter = max_iter,
              offset_variables = offset_variables,
              tol = tol)
  names(res) <- paste0(prefix, names(res))
  res
}

# Initialize the Z matrix according to covariates
.opt_esvd_setup_z_mat <- function(covariates, p, z_init) {
  if(is.null(covariates))
  {
    z_mat <- NULL
  } else {
    if(is.null(z_init))
    {
      z_mat <- matrix(0, nrow = p, ncol = ncol(covariates))
    } else {
      z_mat <- z_init
    }
  }
}

# Set row and column names of output matrices
.opt_esvd_format_matrices <- function(covariates, dat, x_mat, y_mat, z_mat) {
  stopifnot(nrow(dat) == nrow(x_mat),
            ncol(dat) == nrow(y_mat),
            ncol(x_mat) == ncol(y_mat))
  if(length(rownames(dat)) > 0) rownames(x_mat) <- rownames(dat)
  if(length(colnames(dat)) > 0) rownames(y_mat) <- colnames(dat)

  if(!all(is.null(covariates))){
    stopifnot(nrow(covariates) == nrow(dat),
              nrow(z_mat) == ncol(dat),
              ncol(covariates) == ncol(z_mat))

    if(length(colnames(covariates)) > 0) colnames(z_mat) <- colnames(covariates)
    if(length(colnames(dat)) > 0) rownames(z_mat) <- colnames(dat)
  }

  list(x_mat = x_mat, y_mat = y_mat, z_mat = z_mat)
}
