#' @export
compute_posterior <- function(input_obj, ...) UseMethod("compute_posterior")

#' @export
compute_posterior.eSVD <- function(input_obj,
                                   alpha_max = 50,
                                   nuisance_lower_quantile = 0.01,
                                   ...){
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"))

  dat <- .get_object(eSVD_obj = input_obj, what_obj = "dat", which_fit = NULL)
  covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
  latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
  esvd_res <- .get_object(eSVD_obj = input_obj, what_obj = NULL, which_fit = latest_Fit)
  case_control_variable <- .get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
  nuisance_vec <-.get_object(eSVD_obj = input_obj, what_obj = "nuisance", which_fit = latest_Fit)
  stopifnot(all(nuisance_vec >= 0))

  param <- .format_param_posterior(alpha_max = alpha_max,
                                   nuisance_lower_quantile = nuisance_lower_quantile)
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)

  res <- compute_posterior.default(
    input_obj = dat,
    case_control_variable = case_control_variable,
    covariates = covariates,
    esvd_res = esvd_res,
    nuisance_vec = nuisance_vec,
    alpha_max = alpha_max,
    nuisance_lower_quantile = nuisance_lower_quantile
  )

  input_obj[[latest_Fit]]$posterior_mean_mat <- res$posterior_mean_mat
  input_obj[[latest_Fit]]$posterior_var_mat <- res$posterior_var_mat

  input_obj
}

#' @export
compute_posterior.default <- function(input_obj,
                                      case_control_variable,
                                      covariates,
                                      esvd_res,
                                      nuisance_vec,
                                      alpha_max = 50,
                                      nuisance_lower_quantile = 0.01,
                                      ...){
  stopifnot(inherits(input_obj, c("matrix", "dgCMatrix")),
            is.list(esvd_res),
            all(c("x_mat", "y_mat", "z_mat") %in% names(esvd_res)),
            case_control_variable %in% colnames(covariates),
            all(colnames(covariates) == colnames(esvd_res$z_mat)),
            nrow(esvd_res$x_mat) == nrow(covariates),
            nrow(esvd_res$y_mat) == nrow(esvd_res$z_mat),
            ncol(esvd_res$x_mat) == ncol(esvd_res$y_mat))

  if(inherits(input_obj, "dgCMatrix")) input_obj <- as.matrix(input_obj)
  offset_var <- setdiff(colnames(covariates), case_control_variable)

  nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
  nat_mat2 <- tcrossprod(covariates[,case_control_variable,drop = F],
                         esvd_res$z_mat[,case_control_variable,drop = F])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)
  library_mat <- exp(tcrossprod(
    covariates[,offset_var],
    esvd_res$z_mat[,offset_var]
  ))

  nuisance_vec <- pmax(nuisance_vec,
                       stats::quantile(nuisance_vec, probs = nuisance_lower_quantile))
  Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
                 STATS = nuisance_vec, FUN = "*")
  Alpha <- pmin(Alpha, alpha_max)
  AplusAlpha <- input_obj + Alpha
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
                                    nuisance_lower_quantile) {
  list(posterior_alpha_max = alpha_max,
       posterior_nuisance_lower_quantile = nuisance_lower_quantile)
}

