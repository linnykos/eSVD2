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
                                   bool_covariates_as_library = F,
                                   bool_library_includes_interept = T,
                                   bool_use_log = F,
                                   min_val =  1e-4,
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

  library_size_variable <- .get_object(input_obj,
                                       which_fit = "param",
                                       what_obj = "init_library_size_variable")
  library_size_variables <- library_size_variable
  if(bool_covariates_as_library) library_size_variables <- c(library_size_variables, setdiff(colnames(covariates), c("Intercept", case_control_variable)))
  if(bool_library_includes_interept) library_size_variables <- c("Intercept", library_size_variables)

  library_idx <- which(colnames(covariates) %in% library_size_variables)

  nat_mat1 <- tcrossprod(x_mat, y_mat)
  nat_mat2 <- tcrossprod(covariates[,-library_idx], z_mat[,-library_idx])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)

  library_mat <- exp(tcrossprod(
    covariates[,library_idx], z_mat[,library_idx]
  ))

  nuisance_vec <- estimate_nuisance.default(
    input_obj = dat,
    mean_mat = mean_mat_nolib,
    library_mat = library_mat,
    bool_use_log = bool_use_log,
    min_val = min_val,
    verbose = verbose
  )

  input_obj[[latest_Fit]]$nuisance_vec <- nuisance_vec
  param <- .format_param_nuisance(bool_covariates_as_library = bool_covariates_as_library,
                                  bool_library_includes_interept = bool_library_includes_interept,
                                  bool_use_log = bool_use_log,
                                  min_val = min_val)
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)

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
                                      bool_use_log = F,
                                      min_val =  1e-4,
                                      verbose = 0, ...){
  stopifnot(inherits(input_obj, c("matrix", "dgCMatrix")),
            is.matrix(mean_mat), is.matrix(library_mat),
            all(dim(mean_mat) == dim(input_obj)),
            all(dim(library_mat) == dim(input_obj)))

  p <- ncol(input_obj)
  nuisance_vec <- sapply(1:p, function(j){
    if(verbose ==1 && p > 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(paste0(j, " of ", p))

    .nuisance_in_sequence(j = j,
                          mu = mean_mat[,j],
                          s = library_mat[,j],
                          x = as.numeric(input_obj[,j]),
                          bool_use_log = bool_use_log,
                          verbose = verbose)
  })
  if(length(colnames(input_obj)) > 0) names(nuisance_vec) <- colnames(input_obj)

  pmax(nuisance_vec, min_val)
}

.nuisance_in_sequence <- function(j, mu, s, x, bool_use_log, verbose){
  if(!bool_use_log){
    res <- tryCatch(
      gamma_rate(x = x,
                 mu = mu,
                 s = s),
      warning = function(e){NULL},
      error = function(e){NULL})
    if(!all(is.null(res))) {return(res)}
  }

  res <- tryCatch(
    exp(log_gamma_rate(x = x,
                       mu = mu,
                       s = s)),
    warning = function(e){NULL},
    error = function(e){NULL})
  if(!all(is.null(res))) {return(res)}

  if(verbose > 0) warning(paste0("Nuisance estimation failed at variable ", j))
  return(0)
}


.format_param_nuisance <- function(bool_covariates_as_library,
                                   bool_library_includes_interept,
                                   bool_use_log,
                                   min_val) {
  list(nuisance_bool_covariates_as_library = bool_covariates_as_library,
       nuisance_bool_library_includes_interept = bool_library_includes_interept,
       nuisance_bool_use_log = bool_use_log,
       nuisance_min_val = min_val)
}


