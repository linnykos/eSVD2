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
#' @param input_obj                       \code{eSVD} object output from \code{opt_esvd.eSVD}.
#' @param alpha_max                       Upper bound on the fitted mean without the library, \eqn{\mu_{ji}}, before it forms the prior \eqn{\mu_{ji}\beta_j}, default is \code{1e3}.
#' @param bool_adjust_covariates          Boolean to adjust the numerator in the posterior by the donor covariates, default is \code{FALSE}.
#'                                        This parameter is experimental, and we have not yet encountered a scenario where it is useful to be set to be \code{TRUE}.
#' @param bool_covariates_as_library      Boolean to include the donor covariates effects in the adjusted library size, default is \code{TRUE}
#' @param bool_return_components          Boolean to return the numerator and denominator of the posterior terms as well (themselves cell-by-gene matrices), default is \code{FALSE}
#' @param bool_stabilize_underdispersion  Boolean, default \code{TRUE}. \code{nuisance_vec}
#'                                        holds Gamma \emph{rates}, which share the units of the library
#'                                        size, so a geometric mean above 1 means less over-dispersion
#'                                        than a unit rate in those units (not under-dispersion: a finite
#'                                        rate is always over-dispersed). When that happens
#'                                        (\code{mean(log10(nuisance_vec)) > 0}), every gene's rate,
#'                                        capped ones included, is divided by the one geometric mean so that
#'                                        \code{mean(log10(nuisance_vec))} becomes 0; otherwise the
#'                                        rates are left alone
#' @param library_min                     All covariate-adjusted library size smaller than this value are set to this value, default is 0.1.
#' @param nuisance_lower_quantile         All the nuisance values that are smaller than this quantile
#'                                        are set to this quantile, default is 0.01; one number in \eqn{[0, 1]}, or \code{NULL} for no floor.
#' @param pseudocount                     The additional count that is added to the count matrix, default is 0.
#' @param ...                             Additional parameters.
#'
#' @return \code{eSVD} object with \code{posterior_mean_mat}
#' and \code{posterior_var_mat} appended to the list in
#' \code{input_obj[[input_obj[["latest_Fit"]]]]}.
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
#' esvd_obj <- estimate_nuisance(input_obj = esvd_obj)
#' esvd_obj <- compute_posterior(input_obj = esvd_obj,
#'                               alpha_max = 2 * max(sim$obs_mat))
#' dim(esvd_obj$fit_First$posterior_mean_mat)
#' @export
compute_posterior.eSVD <- function(input_obj,
                                   alpha_max = 1e3,
                                   bool_adjust_covariates = FALSE,
                                   bool_covariates_as_library = TRUE,
                                   bool_return_components = FALSE,
                                   bool_stabilize_underdispersion = TRUE,
                                   library_min = 0.1,
                                   nuisance_lower_quantile = 0.01,
                                   pseudocount = 0,
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
                                   bool_return_components = bool_return_components,
                                   bool_stabilize_underdispersion = bool_stabilize_underdispersion,
                                   library_min = library_min,
                                   nuisance_lower_quantile = nuisance_lower_quantile,
                                   pseudocount = pseudocount)
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)
  # `.combine_two_named_lists` keeps an entry that is already there, so a
  # second call with other settings would leave the first call's settings
  # in `param`, and `recompute_pvalue` would repeat the posterior with them.
  input_obj$param[names(param)] <- param
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
    bool_return_components = bool_return_components,
    bool_stabilize_underdispersion = bool_stabilize_underdispersion,
    library_min = library_min,
    nuisance_lower_quantile = nuisance_lower_quantile,
    pseudocount = pseudocount
  )

  input_obj[[latest_Fit]]$posterior_mean_mat <- res$posterior_mean_mat
  input_obj[[latest_Fit]]$posterior_var_mat <- res$posterior_var_mat
  if(bool_return_components){
    input_obj[[latest_Fit]]$numerator_mat <- res$numerator_mat
    input_obj[[latest_Fit]]$denominator_mat <- res$denominator_mat
  }

  input_obj
}

#' Compute posterior according to Gamma-Poisson model for matrices and sparse matrices.
#'
#' @param input_obj                       Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
#'                                        and \eqn{p} columns represent genes.
#'                                        The rows and columns of the matrix should be named.
#' @param case_control_variable           A string of the column name of \code{covariates} which depicts the case-control
#'                                        status of each cell. Notably, this should be a binary variable where a \code{1}
#'                                        is hard-coded to describe case, and a \code{0} to describe control.
#' @param covariates                      \code{matrix} object with \eqn{n} rows with the same rownames as \code{dat} where the columns
#'                                        represent the different covariates.
#'                                        Notably, this should contain only numerical columns (i.e., all categorical
#'                                        variables should have already been split into numerous indicator variables).
#' @param esvd_res                        Output of \code{opt_esvd.default}, notably with elements
#'                                        \code{x_mat}, \code{y_mat} and \code{z_mat}
#' @param library_size_variable           A string of the variable name (which must be in \code{covariates}) of which variable denotes the sequenced (i.e., observed) library size.
#' @param nuisance_vec                    Vector of non-negative numerics of length \code{ncol(input_obj)}, such as
#'                                        the output of \code{estimate_nuisance.default}. These are Gamma
#'                                        \emph{rates} \eqn{\beta_j = 1/\gamma_j}: the posterior mean is
#'                                        \eqn{(A_{ji} + \mu_{ji}\beta_j)/(\ell_{ji} + \beta_j)}, so a larger
#'                                        value pulls the posterior harder towards the fitted mean.
#' @param alpha_max                       Upper bound on the fitted mean without the library, \eqn{\mu_{ji}}, before it forms the prior \eqn{\mu_{ji}\beta_j}, default is \code{1e3}.
#' @param bool_adjust_covariates          Boolean to adjust the numerator in the posterior by the donor covariates, default is \code{FALSE}.
#'                                        This parameter is experimental, and we have not yet encountered a scenario where it is useful to be set to be \code{TRUE}.
#' @param bool_covariates_as_library      Boolean to include the donor covariates effects in the adjusted library size, default is \code{TRUE}.
#' @param bool_library_includes_interept  Boolean if the intercept term from the eSVD matrix factorization should be included in the calculation for the covariate-adjusted library size, default is \code{TRUE}.
#' @param bool_return_components          Boolean to return the numerator and denominator of the posterior terms as well (themselves cell-by-gene matrices), default is \code{FALSE}.
#' @param bool_stabilize_underdispersion  Boolean, default \code{TRUE}. \code{nuisance_vec}
#'                                        holds Gamma \emph{rates}, which share the units of the library
#'                                        size, so a geometric mean above 1 means less over-dispersion
#'                                        than a unit rate in those units (not under-dispersion: a finite
#'                                        rate is always over-dispersed). When that happens
#'                                        (\code{mean(log10(nuisance_vec)) > 0}), every gene's rate,
#'                                        capped ones included, is divided by the one geometric mean so that
#'                                        \code{mean(log10(nuisance_vec))} becomes 0; otherwise the
#'                                        rates are left alone.
#' @param library_min                     All covariate-adjusted library size smaller than this value are set to this value, default is 0.1.
#' @param nuisance_lower_quantile         All the nuisance values that are smaller than this quantile
#'                                        are set to this quantile; one number in \eqn{[0, 1]}, or \code{NULL} for no floor.
#' @param pseudocount                     The additional count that is added to the count matrix, default is 0.
#' @param ...                             Additional parameters.
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
                                      alpha_max = 1e3,
                                      bool_adjust_covariates = FALSE,
                                      bool_covariates_as_library = TRUE,
                                      bool_library_includes_interept = TRUE,
                                      bool_return_components = FALSE,
                                      bool_stabilize_underdispersion = TRUE,
                                      library_min = 0.1,
                                      nuisance_lower_quantile = 0.01,
                                      pseudocount = 0,
                                      ...){
  .check_posterior_args(alpha_max = alpha_max,
                        bool_adjust_covariates = bool_adjust_covariates,
                        bool_covariates_as_library = bool_covariates_as_library,
                        bool_stabilize_underdispersion = bool_stabilize_underdispersion,
                        library_min = library_min,
                        nuisance_lower_quantile = nuisance_lower_quantile,
                        pseudocount = pseudocount)
  stopifnot(inherits(input_obj, c("matrix", "dgCMatrix")),
            is.list(esvd_res),
            all(c("x_mat", "y_mat", "z_mat") %in% names(esvd_res)),
            case_control_variable %in% colnames(covariates),
            library_size_variable %in% colnames(covariates),
            all(colnames(covariates) == colnames(esvd_res$z_mat)),
            nrow(esvd_res$x_mat) == nrow(covariates),
            nrow(esvd_res$y_mat) == nrow(esvd_res$z_mat),
            ncol(esvd_res$x_mat) == ncol(esvd_res$y_mat))

  if(inherits(input_obj, "dgCMatrix")) input_obj <- as.matrix(input_obj)
  if(!is.null(pseudocount) && pseudocount > 0) input_obj <- input_obj + pseudocount
  if(is.null(case_control_variable)) case_control_variable <- numeric(0)
  case_control_idx <- which(colnames(covariates) == case_control_variable)

  # The one rule for which columns form the library size, shared with
  # `estimate_nuisance` and `compute_test_per_gene`. They agree when they are
  # passed the same library booleans, as they are by default.
  library_idx <- .nuisance_library_idx(
    covariates = covariates,
    case_control_variable = case_control_variable,
    library_size_variable = library_size_variable,
    bool_covariates_as_library = bool_covariates_as_library,
    bool_library_includes_interept = bool_library_includes_interept
  )
  idx_vec <- c(case_control_idx, library_idx)

  nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
  nat_mat2 <- tcrossprod(covariates[,-library_idx], esvd_res$z_mat[,-library_idx])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)
  library_mat <- exp(tcrossprod(
    covariates[,library_idx], esvd_res$z_mat[,library_idx]
  ))
  if(!is.null(library_min)) library_mat <- pmax(library_mat, library_min)

  # `stats::quantile(x, probs = NULL)` is `numeric(0)`, and `pmax()` with it
  # would empty `nuisance_vec`, so NULL has to skip the floor explicitly.
  if(!is.null(nuisance_lower_quantile)){
    nuisance_vec <- pmax(nuisance_vec,
                         stats::quantile(nuisance_vec,
                                         probs = nuisance_lower_quantile))
  }
  # Recentre log10(nuisance_vec) at 0 when the geometric mean of the rates is
  # above 1 (`nuisance_vec` is the Gamma rate, in the units of the library
  # size, so that is LESS over-dispersion than a unit rate; it is not
  # under-dispersion). Every rate, capped ones included, is divided by the
  # same factor. Written as a
  # plain subtraction rather than `scale()`, which returns a p x 1 matrix and
  # moves `names()` to `rownames()`.
  if(bool_stabilize_underdispersion && mean(log10(nuisance_vec)) > 0) {
    nuisance_vec <- 10^(log10(nuisance_vec) - mean(log10(nuisance_vec)))
  }

  if(!is.null(alpha_max)) mean_mat_nolib <- pmin(mean_mat_nolib, alpha_max)
  Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
                 STATS = nuisance_vec, FUN = "*")
  AplusAlpha <- input_obj + Alpha
  # adjust the Alpha's based on the confounding covariates
  if(bool_adjust_covariates){
    tmp <- log(AplusAlpha)
    nat_mat_confounder <- tcrossprod(covariates[,-idx_vec,drop = FALSE],
                                     esvd_res$z_mat[,-idx_vec,drop = FALSE])
    AplusAlpha <- exp(tmp - nat_mat_confounder)
  }


  SplusBeta <- sweep(library_mat, MARGIN = 2,
                     STATS = nuisance_vec, FUN = "+")
  posterior_mean_mat <- AplusAlpha/SplusBeta
  posterior_var_mat <- AplusAlpha/SplusBeta^2

  rownames(posterior_mean_mat) <- rownames(input_obj)
  rownames(posterior_var_mat) <- rownames(input_obj)
  colnames(posterior_mean_mat) <- colnames(input_obj)
  colnames(posterior_var_mat) <- colnames(input_obj)

  if(!bool_return_components){
    list(posterior_mean_mat = posterior_mean_mat,
         posterior_var_mat = posterior_var_mat)
  } else {
    list(posterior_mean_mat = posterior_mean_mat,
         posterior_var_mat = posterior_var_mat,
         numerator_mat = AplusAlpha,
         denominator_mat = SplusBeta)
  }
}

#' Check the settings the two posterior paths share
#'
#' Called first by \code{compute_posterior.default} and by
#' \code{compute_test_per_gene}, so that the matrix path and the per-gene
#' path refuse the same inputs with the same messages.
#'
#' @param alpha_max                       \code{NULL} (no bound) or one number above 0; \code{Inf} is allowed.
#' @param bool_adjust_covariates          One \code{TRUE} or \code{FALSE}.
#' @param bool_covariates_as_library      One \code{TRUE} or \code{FALSE}; not \code{TRUE} together with
#'                                        \code{bool_adjust_covariates}.
#' @param bool_stabilize_underdispersion  One \code{TRUE} or \code{FALSE}.
#' @param library_min                     \code{NULL} (no floor) or one positive finite number.
#' @param nuisance_lower_quantile         \code{NULL} (no floor) or one number in \eqn{[0, 1]}.
#' @param pseudocount                     \code{NULL} or one non-negative finite number.
#'
#' @return \code{NULL}, invisibly; an error names the first setting that fails.
#' @noRd
.check_posterior_args <- function(alpha_max,
                                  bool_adjust_covariates,
                                  bool_covariates_as_library,
                                  bool_stabilize_underdispersion,
                                  library_min,
                                  nuisance_lower_quantile,
                                  pseudocount){
  .describe <- function(x){
    if(is.null(x)) return("NULL")
    paste0(deparse(x), collapse = "")
  }
  .is_one_number <- function(x){
    is.numeric(x) && length(x) == 1 && !is.na(x)
  }

  bool_list <- list(bool_adjust_covariates = bool_adjust_covariates,
                    bool_covariates_as_library = bool_covariates_as_library,
                    bool_stabilize_underdispersion = bool_stabilize_underdispersion)
  for(arg_name in names(bool_list)){
    x <- bool_list[[arg_name]]
    if(!is.logical(x) || length(x) != 1 || is.na(x)){
      stop("`", arg_name, "` must be one TRUE or FALSE; received ", .describe(x))
    }
  }
  if(bool_adjust_covariates && bool_covariates_as_library){
    stop("`bool_adjust_covariates` and `bool_covariates_as_library` cannot ",
         "both be TRUE: the first divides the posterior numerator by the ",
         "covariate effects that the second puts in the library size")
  }

  if(!is.null(alpha_max) && (!.is_one_number(alpha_max) || alpha_max <= 0)){
    stop("`alpha_max` must be NULL or one number above 0; received ",
         .describe(alpha_max))
  }
  if(!is.null(library_min) && (!.is_one_number(library_min) ||
                               !is.finite(library_min) || library_min <= 0)){
    stop("`library_min` must be NULL or one positive finite number; received ",
         .describe(library_min))
  }
  if(!is.null(nuisance_lower_quantile) &&
     (!.is_one_number(nuisance_lower_quantile) ||
      nuisance_lower_quantile < 0 || nuisance_lower_quantile > 1)){
    stop("`nuisance_lower_quantile` must be NULL or one number between 0 and 1; ",
         "received ", .describe(nuisance_lower_quantile))
  }
  if(!is.null(pseudocount) && (!.is_one_number(pseudocount) ||
                               !is.finite(pseudocount) || pseudocount < 0)){
    stop("`pseudocount` must be NULL or one non-negative finite number; ",
         "received ", .describe(pseudocount))
  }

  invisible()
}

.format_param_posterior <- function(alpha_max,
                                    bool_adjust_covariates,
                                    bool_covariates_as_library,
                                    bool_return_components,
                                    bool_stabilize_underdispersion,
                                    library_min,
                                    nuisance_lower_quantile,
                                    pseudocount) {
  list(posterior_alpha_max = alpha_max,
       posterior_bool_adjust_covariates = bool_adjust_covariates,
       posterior_bool_covariates_as_library = bool_covariates_as_library,
       posterior_bool_return_components = bool_return_components,
       posterior_bool_stabilize_underdispersion = bool_stabilize_underdispersion,
       posterior_library_min = library_min,
       posterior_nuisance_lower_quantile = nuisance_lower_quantile,
       posterior_pseudocount = pseudocount)
}

