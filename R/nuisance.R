#' Estimate nuisance values
#'
#' Generic function interface
#'
#' @param input_obj Main object
#' @param ...       Additional parameters
#'
#' @return Output dependent on class of \code{input_obj}
#' @export
estimate_nuisance <- function(input_obj, ...) {UseMethod("estimate_nuisance")}

#' Estimate nuisance values for eSVD objects
#'
#' Assumes a Gamma-Poisson model where the mean and variance are proportionally
#' related.
#'
#' @param input_obj \code{eSVD} object outputed from \code{opt_esvd.eSVD}.
#'                  Specifically, the nuisance parameters will be estimated
#'                  based on the fit in \code{input_obj[[input_obj[["latest_Fit"]]]]}.
#' @param verbose   Integer.
#' @param ...       Additional parameters.
#'
#' @return \code{eSVD} object with \code{nuisance_vec} appended to the list in
#' \code{input_obj[[input_obj[["latest_Fit"]]]]}.
#' @export
estimate_nuisance.eSVD <- function(input_obj,
                                   verbose = 0, ...){
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"))

  dat <- .get_object(eSVD_obj = input_obj, what_obj = "dat", which_fit = NULL)
  covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
  latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
  case_control_variable <- .get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
  x_mat <-.get_object(eSVD_obj = input_obj, what_obj = "x_mat", which_fit = latest_Fit)
  y_mat <-.get_object(eSVD_obj = input_obj, what_obj = "y_mat", which_fit = latest_Fit)
  z_mat <-.get_object(eSVD_obj = input_obj, what_obj = "z_mat", which_fit = latest_Fit)

  nat_mat1 <- tcrossprod(x_mat, y_mat)
  nat_mat2 <- tcrossprod(covariates[,case_control_variable,drop = F],
                         z_mat[,case_control_variable,drop = F])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)

  offset_variable <- setdiff(colnames(covariates), case_control_variable)
  library_mat <- exp(tcrossprod(
    covariates[,offset_variable], z_mat[,offset_variable]
  ))

  nuisance_vec <- estimate_nuisance.default(
    input_obj = dat,
    mean_mat = mean_mat_nolib,
    library_mat = library_mat,
    verbose = verbose
  )

  input_obj[[latest_Fit]]$nuisance_vec <- nuisance_vec
  input_obj
}

#' Estimate nuisance values for matrix or sparse matrices.
#'
#' @param input_obj    Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
#'                     and \eqn{p} columns represent genes.
#'                     The rows and columns of the matrix should be named.
#' @param mean_mat     A \code{matrix} of \eqn{n} rows and \eqn{p} columns that represents the
#'                     expected value of each entry.
#' @param library_mat  A \code{matrix} of \eqn{n} rows and \eqn{p} columns that represents the
#'                     library size of each entry.
#' @param verbose      Integer.
#' @param ...          Additional parameters.
#'
#' @return Vector of length \eqn{p}
#' @export
estimate_nuisance.default <- function(input_obj,
                                      mean_mat,
                                      library_mat,
                                      verbose = 0, ...){
  stopifnot(inherits(input_obj, c("matrix", "dgCMatrix")),
            is.matrix(mean_mat), is.matrix(library_mat),
            all(dim(mean_mat) == dim(input_obj)),
            all(dim(library_mat) == dim(input_obj)))

  p <- ncol(input_obj)
  nuisance_vec <- sapply(1:p, function(j){
    if(verbose ==1 && p > 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(j)

    val <- tryCatch(
      gamma_rate(x = as.numeric(input_obj[,j]),
                 mu = mean_mat[,j],
                 s = library_mat[,j]),
      # exp(log_gamma_rate(x = as.numeric(input_obj[,j]),
      #                mu = mean_mat[,j],
      #                s = library_mat[,j])),
      error = function(c) {
        if(verbose > 0) print(paste0("Nuisance estimation failed at variable ", j))
        0
      })

    val
  })
  if(length(colnames(input_obj)) > 0) names(nuisance_vec) <- colnames(input_obj)

  nuisance_vec
}

