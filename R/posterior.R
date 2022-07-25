#' Compute posterior according to Gamma-Poisson model
#'
#' Generic function interface
#'
#' @param input_obj Main object
#' @param ...       Additional parameters
#'
#' @return Output dependent on class of \code{input_obj}
#' @export
compute_posterior <- function(input_obj, ...) {UseMethod("compute_posterior")}

#' Compute posterior according to Gamma-Poisson model for eSVD object
#'
#' The posterior is computed based on whatever \code{input_obj$latest_Fit} is set to.
#'
#' @param input_obj                 \code{eSVD} object outputed from \code{opt_esvd.eSVD}.
#' @param alpha_max                 Maximum value of numerator when computing posterior.
#' @param nuisance_lower_quantile   All the nuisance values that are smaller than this quantile
#'                                  are set to this quantile.
#' @param ...                       Additional parameters.
#'
#' @return \code{eSVD} object with \code{posterior_mean_mat}
#' and \code{posterior_var_mat} appended to the list in
#' \code{input_obj[[input_obj[["latest_Fit"]]]]}.
#' @export
compute_posterior.eSVD <- function(input_obj,
                                   alpha_max = 50,
                                   bool_adjust_covariates = F,
                                   bool_covariates_as_library = T,
                                   nuisance_lower_quantile = 0.01,
                                   ...){
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"))

  dat <- .get_object(eSVD_obj = input_obj, what_obj = "dat", which_fit = NULL)
  covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
  latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
  esvd_res <- .get_object(eSVD_obj = input_obj, what_obj = NULL, which_fit = latest_Fit)
  nuisance_vec <- .get_object(eSVD_obj = input_obj, what_obj = "nuisance", which_fit = latest_Fit)
  stopifnot(all(nuisance_vec >= 0))

  param <- .format_param_posterior(alpha_max = alpha_max,
                                   bool_adjust_covariates = bool_adjust_covariates,
                                   bool_covariates_as_library = bool_covariates_as_library,
                                   nuisance_lower_quantile = nuisance_lower_quantile)
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)
  case_control_variable <- .get_object(eSVD_obj = input_obj, which_fit = "param", what_obj = "init_case_control_variable")
  library_size_variable <- .get_object(eSVD_obj = input_obj, which_fit = "param", what_obj = "init_library_size_variable")
  bool_library_includes_interept <- .get_object(eSVD_obj = input_obj, which_fit = "param", what_obj = "nuisance_bool_library_includes_interept")

  res <- compute_posterior.default(
    input_obj = dat,
    case_control_variable = case_control_variable,
    covariates = covariates,
    esvd_res = esvd_res,
    library_size_variable = library_size_variable,
    nuisance_vec = nuisance_vec,
    alpha_max = alpha_max,
    bool_adjust_covariates = bool_adjust_covariates,
    bool_covariates_as_library = bool_covariates_as_library,
    bool_library_includes_interept = bool_library_includes_interept,
    nuisance_lower_quantile = nuisance_lower_quantile
  )

  input_obj[[latest_Fit]]$posterior_mean_mat <- res$posterior_mean_mat
  input_obj[[latest_Fit]]$posterior_var_mat <- res$posterior_var_mat

  input_obj
}

#' Compute posterior according to Gamma-Poisson model for matrices and sparse matrices.
#'
#' @param input_obj                Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
#'                                 and \eqn{p} columns represent genes.
#'                                 The rows and columns of the matrix should be named.
#' @param case_control_variable    A string of the column name of \code{covariates} which depicts the case-control
#'                                 status of each cell. Notably, this should be a binary variable where a \code{1}
#'                                 is hard-coded to describe case, and a \code{0} to describe control.
#' @param covariates               \code{matrix} object with \eqn{n} rows with the same rownames as \code{dat} where the columns
#'                                 represent the different covariates.
#'                                 Notably, this should contain only numerical columns (i.e., all categorical
#'                                 variables should have already been split into numerous indicator variables).
#' @param esvd_res                 Output of \code{opt_esvd.default}, notably with elements
#'                                 \code{x_mat}, \code{y_mat} and \code{z_mat}
#' @param nuisance_vec             Vector of non-negative numerics of length \code{ncol(input_obj)}, such as
#'                                 the output of \code{estimate_nuisance.default}.
#' @param alpha_max                Maximum value of numerator when computing posterior.
#' @param nuisance_lower_quantile  All the nuisance values that are smaller than this quantile
#'                                 are set to this quantile.
#' @param ...                      Additional parameters.
#'
#' @return A \code{list} of elements \code{posterior_mean_mat}
#' and \code{posterior_var_mat}
#' @export
compute_posterior.default <- function(input_obj,
                                      case_control_variable,
                                      covariates,
                                      esvd_res,
                                      library_size_variable,
                                      nuisance_vec,
                                      alpha_max = 50,
                                      bool_adjust_covariates = F,
                                      bool_covariates_as_library = T,
                                      bool_library_includes_interept = T,
                                      nuisance_lower_quantile = 0.01,
                                      ...){
  stopifnot(inherits(input_obj, c("matrix", "dgCMatrix")),
            is.list(esvd_res),
            all(c("x_mat", "y_mat", "z_mat") %in% names(esvd_res)),
            case_control_variable %in% colnames(covariates),
            library_size_variable %in% colnames(covariates),
            all(colnames(covariates) == colnames(esvd_res$z_mat)),
            nrow(esvd_res$x_mat) == nrow(covariates),
            nrow(esvd_res$y_mat) == nrow(esvd_res$z_mat),
            ncol(esvd_res$x_mat) == ncol(esvd_res$y_mat),
            !bool_adjust_covariates | !bool_covariates_as_library)

  if(inherits(input_obj, "dgCMatrix")) input_obj <- as.matrix(input_obj)
  case_control_idx <- which(colnames(covariates) == case_control_variable)

  library_size_variables <- library_size_variable
  if(bool_covariates_as_library) library_size_variables <- unique(c(library_size_variables, setdiff(colnames(covariates), c("Intercept", case_control_variable))))
  if(bool_library_includes_interept) library_size_variables <-  unique(c("Intercept", library_size_variables))

  library_idx <- which(colnames(covariates) %in% library_size_variables)
  idx_vec <- c(case_control_idx, library_idx)

  nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
  nat_mat2 <- tcrossprod(covariates[,-library_idx], esvd_res$z_mat[,-library_idx])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)
  library_mat <- exp(tcrossprod(
    covariates[,library_idx], esvd_res$z_mat[,library_idx]
  ))

  nuisance_vec <- pmax(nuisance_vec,
                       stats::quantile(nuisance_vec, probs = nuisance_lower_quantile))
  Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
                 STATS = nuisance_vec, FUN = "*")
  AplusAlpha <- input_obj + Alpha
  # adjust the Alpha's based on the confounding covariates
  if(bool_adjust_covariates){
    tmp <- log(AplusAlpha)
    nat_mat_confounder <- tcrossprod(covariates[,-idx_vec,drop = F],
                                     esvd_res$z_mat[,-idx_vec,drop = F])
    AplusAlpha <- exp(tmp - nat_mat_confounder)
  }

  if(!is.null(alpha_max)) AplusAlpha <- pmin(AplusAlpha, alpha_max)

  SplusBeta <- sweep(library_mat, MARGIN = 2,
                     STATS = nuisance_vec, FUN = "+")
  posterior_mean_mat <- AplusAlpha/SplusBeta
  posterior_var_mat <- AplusAlpha/SplusBeta^2

  rownames(posterior_mean_mat) <- rownames(input_obj)
  rownames(posterior_var_mat) <- rownames(input_obj)
  colnames(posterior_mean_mat) <- colnames(input_obj)
  colnames(posterior_var_mat) <- colnames(input_obj)

  list(posterior_mean_mat = posterior_mean_mat,
       posterior_var_mat = posterior_var_mat)
}


.format_param_posterior <- function(alpha_max,
                                    bool_adjust_covariates,
                                    bool_covariates_as_library,
                                    nuisance_lower_quantile) {
  list(posterior_alpha_max = alpha_max,
       posterior_bool_adjust_covariates = bool_adjust_covariates,
       posterior_bool_covariates_as_library = bool_covariates_as_library,
       posterior_nuisance_lower_quantile = nuisance_lower_quantile)
}

