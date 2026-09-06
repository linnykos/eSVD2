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

#' Estimate nuisance values for eSVD objects (i.e., over-dispersion)
#'
#' Assumes a Gamma-Poisson model where the mean and variance are proportionally
#' related: \eqn{\lambda_{ji} \sim \mathrm{Gamma}(\mathrm{mean} = \mu_{ji},
#' \mathrm{var} = \gamma_j \mu_{ji})}. The estimated \code{nuisance_vec} is
#' the Gamma \emph{rate} \eqn{\beta_j = 1/\gamma_j}, the reciprocal of the
#' over-dispersion \eqn{\gamma_j} in Lin, Qiu and Roeder (2024), so a
#' \strong{larger} value means \strong{less} over-dispersion.
#'
#' @param input_obj                       \code{eSVD} object output from \code{opt_esvd.eSVD}.
#'                                        Specifically, the nuisance parameters will be estimated
#'                                        based on the fit in \code{input_obj[[input_obj[["latest_Fit"]]]]}.
#' @param bool_covariates_as_library      Boolean to adjust the numerator in the posterior by the donor covariates, default is \code{FALSE}.
#'                                        This parameter is experimental, and we have not yet encountered a scenario where it is useful to be set to be \code{TRUE}.
#' @param bool_library_includes_interept  Boolean if the intercept term from the eSVD matrix factorization should be included in the calculation for the covariate-adjusted library size, default is \code{TRUE}.
#' @param bool_use_log                    Boolean if the nuisance (i.e., over-dispersion) parameter should be estimated on the log scale, default is \code{FALSE}.
#' @param min_val                         Minimum value of the nuisance parameter.
#' @param verbose                         Integer.
#' @param ...                             Additional parameters.
#'
#' @return \code{eSVD} object with \code{nuisance_vec} appended to the list in
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
#' summary(esvd_obj$fit_First$nuisance_vec)
#' @export
estimate_nuisance.eSVD <- function(input_obj,
                                   bool_covariates_as_library = FALSE,
                                   bool_library_includes_interept = TRUE,
                                   bool_use_log = FALSE,
                                   min_val =  1e-4,
                                   verbose = 0, ...){
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"))

  dat <- .get_object(eSVD_obj = input_obj, what_obj = "dat", which_fit = NULL)
  covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
  latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
  case_control_variable <- .get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
  if(is.null(case_control_variable)) case_control_variable <- numeric(0)
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

  res <- .estimate_nuisance_matrix(
    dat = dat,
    mean_mat = mean_mat_nolib,
    library_mat = library_mat,
    bool_use_log = bool_use_log,
    min_val = min_val,
    verbose = verbose
  )

  input_obj[[latest_Fit]]$nuisance_vec <- res$nuisance_vec
  param <- .format_param_nuisance(bool_covariates_as_library = bool_covariates_as_library,
                                  bool_library_includes_interept = bool_library_includes_interept,
                                  bool_use_log = bool_use_log,
                                  min_val = min_val,
                                  num_failed = res$num_failed)
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)

  input_obj
}

#' Estimate nuisance values for matrix or sparse matrices.
#'
#' Each gene's value is the maximum-likelihood Gamma \emph{rate}
#' \eqn{\beta_j = 1/\gamma_j} under \eqn{A_{ji} \mid \lambda_{ji} \sim
#' \mathrm{Poisson}(\ell_{ji} \lambda_{ji})}, \eqn{\lambda_{ji} \sim
#' \mathrm{Gamma}(\mu_{ji} \beta_j, \beta_j)}, with \eqn{\mu} from
#' \code{mean_mat} and \eqn{\ell} from \code{library_mat}. Larger values
#' mean less over-dispersion.
#'
#' @param input_obj    Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
#'                     and \eqn{p} columns represent genes.
#'                     The rows and columns of the matrix should be named.
#' @param mean_mat     A \code{matrix} of \eqn{n} rows and \eqn{p} columns that represents the
#'                     expected value of each entry.
#' @param library_mat  A \code{matrix} of \eqn{n} rows and \eqn{p} columns that represents the
#'                     library size of each entry.
#' @param bool_use_log Boolean if the nuisance (i.e., over-dispersion) parameter should be estimated on the log scale, default is \code{FALSE}.
#' @param min_val      Minimum value of the nuisance parameter.
#' @param verbose      Integer.
#' @param ...          Additional parameters.
#'
#' @return Numeric vector of length \eqn{p} (named by \code{colnames(input_obj)}
#' when present) of Gamma rates. A gene whose estimation fails on both
#' routes gets \code{min_val}, and a warning reports how many genes did.
#' @export
estimate_nuisance.default <- function(input_obj,
                                      mean_mat,
                                      library_mat,
                                      bool_use_log = FALSE,
                                      min_val =  1e-4,
                                      verbose = 0, ...){
  res <- .estimate_nuisance_matrix(dat = input_obj,
                                   mean_mat = mean_mat,
                                   library_mat = library_mat,
                                   bool_use_log = bool_use_log,
                                   min_val = min_val,
                                   verbose = verbose)

  res$nuisance_vec
}

#' Estimate every gene's nuisance parameter, counting the failures
#'
#' Shared by both \code{estimate_nuisance} methods. A gene falls through to
#' \code{min_val} when both \code{gamma_rate} and \code{log_gamma_rate} fail;
#' that used to be indistinguishable from a successful fit, since the
#' \code{0} returned on failure was clamped to \code{min_val} and the warning
#' fired only under \code{verbose > 0}. The count is returned so the
#' \code{eSVD} method can record it in \code{param$nuisance_num_failed}, and
#' a warning is raised whenever it is positive.
#'
#' @inheritParams estimate_nuisance.default
#' @param dat  The count matrix (\code{input_obj} of the methods).
#'
#' @returns List with \code{nuisance_vec} (named by \code{colnames(dat)})
#' and \code{num_failed}.
#' @noRd
.estimate_nuisance_matrix <- function(dat,
                                      mean_mat,
                                      library_mat,
                                      bool_use_log,
                                      min_val,
                                      verbose){
  stopifnot(inherits(dat, c("matrix", "dgCMatrix")),
            is.matrix(mean_mat), is.matrix(library_mat))
  if(!all(dim(mean_mat) == dim(dat)) || !all(dim(library_mat) == dim(dat))){
    stop("`mean_mat` (", paste0(dim(mean_mat), collapse = " x "),
         ") and `library_mat` (", paste0(dim(library_mat), collapse = " x "),
         ") must both have the dimensions of the count matrix (",
         paste0(dim(dat), collapse = " x "), ")")
  }

  p <- ncol(dat)
  nuisance_vec <- sapply(1:p, function(j){
    if(verbose ==1 && p > 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(paste0(j, " of ", p))

    .nuisance_in_sequence(j = j,
                          mu = mean_mat[,j],
                          s = library_mat[,j],
                          x = as.numeric(dat[,j]),
                          bool_use_log = bool_use_log,
                          verbose = verbose)
  })

  failed_idx <- which(is.na(nuisance_vec))
  num_failed <- length(failed_idx)
  if(num_failed > 0){
    warning("nuisance estimation failed for ", num_failed, " of ", p,
            " gene(s); those genes are set to `min_val` = ", min_val)
    nuisance_vec[failed_idx] <- 0
  }

  if(length(colnames(dat)) > 0) names(nuisance_vec) <- colnames(dat)

  list(nuisance_vec = pmax(nuisance_vec, min_val),
       num_failed = num_failed)
}

# Returns NA when both routes fail; the caller counts and clamps.
.nuisance_in_sequence <- function(j, mu, s, x, bool_use_log, verbose){
  if(!bool_use_log){
    res <- tryCatch(
      gamma_rate(x = x,
                 mu = mu,
                 s = s),
      warning = function(e){NULL},
      error = function(e){NULL})
    if(!is.null(res) && length(res) == 1 && is.finite(res) && res > 0) {return(res)}
  }

  res <- tryCatch(
    exp(log_gamma_rate(x = x,
                       mu = mu,
                       s = s)),
    warning = function(e){NULL},
    error = function(e){NULL})
  if(!is.null(res) && length(res) == 1 && is.finite(res) && res > 0) {return(res)}

  if(verbose > 0) print(paste0("Nuisance estimation failed at variable ", j))
  return(NA_real_)
}


.format_param_nuisance <- function(bool_covariates_as_library,
                                   bool_library_includes_interept,
                                   bool_use_log,
                                   min_val,
                                   num_failed = 0) {
  list(nuisance_bool_covariates_as_library = bool_covariates_as_library,
       nuisance_bool_library_includes_interept = bool_library_includes_interept,
       nuisance_bool_use_log = bool_use_log,
       nuisance_min_val = min_val,
       nuisance_num_failed = num_failed)
}


