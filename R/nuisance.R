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
#' \strong{The cap.} Each gene's rate is the maximum-likelihood estimate,
#' bounded above: \eqn{\hat\beta_j = \min\{\hat\beta_j^{\mathrm{MLE}},\ c
#' \cdot \mathrm{median}_i\, \ell_{ji}\}}, with \eqn{c} =
#' \code{cap_multiplier} and \eqn{\ell_{ji}} the covariate-adjusted library
#' size of gene \eqn{j} in cell \eqn{i}. The rate has the units of the library
#' size, which is why the bound is a multiple of it. The bound is a
#' calibration device and not an estimate of the smallest over-dispersion in
#' the data: a gene whose counts are no more variable around the fit than a
#' Poisson has no finite maximum-likelihood rate, its posterior would collapse
#' onto the fit, and its test statistic would be inflated. The bound limits
#' how much more weight the posterior may give the fit for one gene than for
#' the typical gene. \code{cap_multiplier = Inf} removes it, which is the
#' behaviour of version 1.1.0. \code{cap_multiplier = 1} is close to, but not
#' the same as, the bound of the version that accompanied Lin, Qiu and Roeder
#' (2024), which was the largest library size of the gene and not a multiple
#' of its median. \code{recompute_pvalue} changes the multiplier on a
#' fitted object without estimating anything again.
#'
#' \strong{The floor} \code{min_val} is in the same units: no gene's rate is
#' below \code{min_val} times the median of its library size. On the
#' unit-free scale of \code{plot_nuisance} every rate therefore lies between
#' \code{min_val} and \code{cap_multiplier}, and no gene is above its cap.
#' \code{min_val} must be below \code{cap_multiplier}.
#'
#' \strong{The status} of each gene records what happened to it:
#' \describe{
#'   \item{\code{estimated}}{the gene keeps its maximum-likelihood rate.}
#'   \item{\code{capped}}{a finite maximum-likelihood rate exists and is above
#'     the bound.}
#'   \item{\code{boundary}}{no finite maximum-likelihood rate exists: the
#'     score statistic for over-dispersion, \eqn{D_j = \sum_i [(A_{ji} -
#'     \ell_{ji}\mu_{ji})^2 - A_{ji}] / \mu_{ji}}, is not positive. Such a
#'     gene is set to the bound, whatever value the iterations of the
#'     estimate stopped at (which is what \code{nuisance_mle_vec} holds for
#'     it). With \code{cap_multiplier = Inf} it keeps that value.}
#'   \item{\code{failed}}{both estimation routes failed and the rate is the
#'     floor, \code{min_val} times the median of its library size.}
#' }
#'
#' @param input_obj                       \code{eSVD} object output from \code{opt_esvd.eSVD}.
#'                                        Specifically, the nuisance parameters will be estimated
#'                                        based on the fit in \code{input_obj[[input_obj[["latest_Fit"]]]]}.
#' @param bool_covariates_as_library      Boolean to adjust the numerator in the posterior by the donor covariates, default is \code{FALSE}.
#'                                        This parameter is experimental, and we have not yet encountered a scenario where it is useful to be set to be \code{TRUE}.
#' @param bool_library_includes_interept  Boolean if the intercept term from the eSVD matrix factorization should be included in the calculation for the covariate-adjusted library size, default is \code{TRUE}.
#' @param bool_use_log                    Boolean if the nuisance (i.e., over-dispersion) parameter should be estimated on the log scale, default is \code{FALSE}.
#' @param cap_multiplier                  One positive number, default \code{10}: no gene's rate may exceed this
#'                                        many times the median over cells of its library size.
#'                                        \code{Inf} means no cap.
#' @param min_val                         Minimum value of the rate, as a multiple of the median over cells of the
#'                                        gene's library size (the units of \code{cap_multiplier}); must be below
#'                                        \code{cap_multiplier}. Default \code{1e-4}.
#' @param verbose                         Integer.
#' @param ...                             Additional parameters.
#'
#' @return \code{eSVD} object with the following named vectors, one entry per
#' gene, added to the list in \code{input_obj[[input_obj[["latest_Fit"]]]]}:
#' \code{nuisance_vec} (the rate after the cap, which is what every later
#' stage uses), \code{nuisance_mle_vec} (the rate before the cap),
#' \code{nuisance_library_median_vec} (the median over cells of the gene's
#' library size, so the cap is \code{cap_multiplier} times this),
#' \code{nuisance_status} (a factor with levels \code{estimated},
#' \code{capped}, \code{boundary}, \code{failed}; see above),
#' \code{gene_mean_count_vec} (the gene's mean count over cells) and
#' \code{gene_sparsity_vec} (the fraction of cells in which its count is
#' zero). \code{input_obj$param} records \code{nuisance_cap_multiplier},
#' \code{nuisance_num_capped} (the number of genes whose rate the cap
#' replaced, the \code{boundary} genes among them included),
#' \code{nuisance_num_boundary} and \code{nuisance_num_failed}; calling the
#' function again overwrites them. \code{plot_nuisance} draws these.
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
#' table(esvd_obj$fit_First$nuisance_status)
#' esvd_obj$param$nuisance_num_capped
#' @export
estimate_nuisance.eSVD <- function(input_obj,
                                   bool_covariates_as_library = FALSE,
                                   bool_library_includes_interept = TRUE,
                                   bool_use_log = FALSE,
                                   cap_multiplier = 10,
                                   min_val =  1e-4,
                                   verbose = 0, ...){
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"))
  .check_cap_multiplier(cap_multiplier)
  .check_min_val(min_val, cap_multiplier)

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
  library_idx <- .nuisance_library_idx(
    covariates = covariates,
    case_control_variable = case_control_variable,
    library_size_variable = library_size_variable,
    bool_covariates_as_library = bool_covariates_as_library,
    bool_library_includes_interept = bool_library_includes_interept
  )

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
    cap_multiplier = cap_multiplier,
    min_val = min_val,
    verbose = verbose
  )

  input_obj[[latest_Fit]]$nuisance_vec <- res$nuisance_vec
  input_obj[[latest_Fit]]$nuisance_mle_vec <- res$nuisance_mle_vec
  input_obj[[latest_Fit]]$nuisance_library_median_vec <- res$library_median_vec
  input_obj[[latest_Fit]]$nuisance_status <- res$nuisance_status

  # The two summaries of the counts that `plot_nuisance` draws the rate
  # against. Stored here because `eSVD(bool_diet = TRUE)` drops `dat`.
  gene_mean_count_vec <- Matrix::colMeans(dat)
  gene_sparsity_vec <- 1 - Matrix::colSums(dat != 0) / nrow(dat)
  names(gene_mean_count_vec) <- colnames(dat)
  names(gene_sparsity_vec) <- colnames(dat)
  input_obj[[latest_Fit]]$gene_mean_count_vec <- gene_mean_count_vec
  input_obj[[latest_Fit]]$gene_sparsity_vec <- gene_sparsity_vec

  param <- .format_param_nuisance(bool_covariates_as_library = bool_covariates_as_library,
                                  bool_library_includes_interept = bool_library_includes_interept,
                                  bool_use_log = bool_use_log,
                                  cap_multiplier = cap_multiplier,
                                  min_val = min_val,
                                  num_boundary = res$num_boundary,
                                  num_capped = res$num_capped,
                                  num_failed = res$num_failed)
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)
  # `.combine_two_named_lists` keeps an entry that is already there, so a
  # second call at another cap would leave the first call's multiplier and
  # counts beside rates they do not describe.
  input_obj$param[names(param)] <- param

  input_obj
}

#' Estimate nuisance values for matrix or sparse matrices.
#'
#' Each gene's value is the maximum-likelihood Gamma \emph{rate}
#' \eqn{\beta_j = 1/\gamma_j} under \eqn{A_{ji} \mid \lambda_{ji} \sim
#' \mathrm{Poisson}(\ell_{ji} \lambda_{ji})}, \eqn{\lambda_{ji} \sim
#' \mathrm{Gamma}(\mu_{ji} \beta_j, \beta_j)}, with \eqn{\mu} from
#' \code{mean_mat} and \eqn{\ell} from \code{library_mat}. Larger values
#' mean less over-dispersion. The estimate is then bounded above by
#' \code{cap_multiplier} times the median of the gene's column of
#' \code{library_mat}; see \code{estimate_nuisance.eSVD} for why.
#'
#' @param input_obj    Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
#'                     and \eqn{p} columns represent genes.
#'                     The rows and columns of the matrix should be named.
#' @param mean_mat     A \code{matrix} of \eqn{n} rows and \eqn{p} columns that represents the
#'                     expected value of each entry.
#' @param library_mat  A \code{matrix} of \eqn{n} rows and \eqn{p} columns that represents the
#'                     library size of each entry.
#' @param bool_use_log Boolean if the nuisance (i.e., over-dispersion) parameter should be estimated on the log scale, default is \code{FALSE}.
#' @param cap_multiplier One positive number, default \code{10}: no gene's rate may exceed this
#'                     many times the median of its column of \code{library_mat}.
#'                     \code{Inf} means no cap.
#' @param min_val      Minimum value of the rate, as a multiple of the median
#'                     over cells of the gene's library size (the units of
#'                     \code{cap_multiplier}); must be below
#'                     \code{cap_multiplier}. Default \code{1e-4}.
#' @param verbose      Integer.
#' @param ...          Additional parameters.
#'
#' @return Numeric vector of length \eqn{p} (named by \code{colnames(input_obj)}
#' when present) of Gamma rates, after the cap. A gene whose estimation fails
#' on both routes gets \code{min_val}, and a warning reports how many genes
#' did. The rates before the cap and each gene's status are returned by the
#' \code{eSVD} method only.
#' @export
estimate_nuisance.default <- function(input_obj,
                                      mean_mat,
                                      library_mat,
                                      bool_use_log = FALSE,
                                      cap_multiplier = 10,
                                      min_val =  1e-4,
                                      verbose = 0, ...){
  .check_cap_multiplier(cap_multiplier)
  .check_min_val(min_val, cap_multiplier)

  res <- .estimate_nuisance_matrix(dat = input_obj,
                                   mean_mat = mean_mat,
                                   library_mat = library_mat,
                                   bool_use_log = bool_use_log,
                                   cap_multiplier = cap_multiplier,
                                   min_val = min_val,
                                   verbose = verbose)

  res$nuisance_vec
}

#' Estimate every gene's nuisance parameter, counting the failures
#'
#' Shared by both \code{estimate_nuisance} methods. A gene falls through to
#' the floor, \code{min_val} times its median library size, when both
#' \code{gamma_rate} and \code{log_gamma_rate} fail; that used to be
#' indistinguishable from a successful fit, since the \code{0} returned on
#' failure was clamped to \code{min_val} and the warning
#' fired only under \code{verbose > 0}. The count is returned so the
#' \code{eSVD} method can record it in \code{param$nuisance_num_failed}, and
#' a warning is raised whenever it is positive.
#'
#' The cap is applied here, in R, and not in \code{gamma_rate}, which stays a
#' maximum-likelihood estimator.
#'
#' @inheritParams estimate_nuisance.default
#' @param dat  The count matrix (\code{input_obj} of the methods).
#'
#' @returns List with \code{library_median_vec}, \code{nuisance_mle_vec},
#' \code{nuisance_status}, \code{nuisance_vec} (each named by
#' \code{colnames(dat)}), \code{num_boundary}, \code{num_capped} and
#' \code{num_failed}.
#' @noRd
.estimate_nuisance_matrix <- function(dat,
                                      mean_mat,
                                      library_mat,
                                      bool_use_log,
                                      cap_multiplier,
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
  # One pass over the genes: the estimate and the boundary statistic share
  # the column, which on a `dgCMatrix` is the expensive part. `vapply` keeps
  # the 2 x p shape when p is 1.
  per_gene_mat <- vapply(seq_len(p), function(j){
    if(verbose ==1 && p > 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(paste0(j, " of ", p))

    x_vec <- as.numeric(dat[,j])
    mu_vec <- mean_mat[,j]
    s_vec <- library_mat[,j]
    c(.nuisance_in_sequence(j = j,
                            mu = mu_vec,
                            s = s_vec,
                            x = x_vec,
                            bool_use_log = bool_use_log,
                            verbose = verbose),
      .compute_boundary_statistic(x_vec = x_vec,
                                  mu_vec = mu_vec,
                                  s_vec = s_vec))
  }, FUN.VALUE = numeric(2))
  nuisance_vec <- per_gene_mat[1, ]
  boundary_statistic_vec <- per_gene_mat[2, ]

  # `as.numeric()`: whether `colMedians` carries the column names depends on
  # the version of matrixStats. The names are set from `dat` below.
  library_median_vec <- as.numeric(matrixStats::colMedians(library_mat))

  bool_failed_vec <- is.na(nuisance_vec)
  failed_idx <- which(bool_failed_vec)
  num_failed <- length(failed_idx)
  if(num_failed > 0){
    warning("nuisance estimation failed for ", num_failed, " of ", p,
            " gene(s); those genes are set to the floor, `min_val` = ",
            min_val, " times their median library size")
    nuisance_vec[failed_idx] <- 0
  }
  # The floor shares the units of the cap and is below it (`.check_min_val`),
  # so no gene can end above its cap.
  nuisance_mle_vec <- pmax(nuisance_vec, min_val * library_median_vec)
  # A statistic that is not finite (a zero or non-finite fitted mean) says
  # nothing about the boundary.
  bool_boundary_vec <- is.finite(boundary_statistic_vec) &
    boundary_statistic_vec <= 0

  if(length(colnames(dat)) > 0){
    names(nuisance_mle_vec) <- colnames(dat)
    names(library_median_vec) <- colnames(dat)
  }

  cap_res <- .apply_nuisance_cap(nuisance_mle_vec = nuisance_mle_vec,
                                 library_median_vec = library_median_vec,
                                 bool_boundary_vec = bool_boundary_vec,
                                 bool_failed_vec = bool_failed_vec,
                                 cap_multiplier = cap_multiplier,
                                 min_val = min_val)

  list(library_median_vec = library_median_vec,
       nuisance_mle_vec = nuisance_mle_vec,
       nuisance_status = cap_res$nuisance_status,
       nuisance_vec = cap_res$nuisance_vec,
       num_boundary = cap_res$num_boundary,
       num_capped = cap_res$num_capped,
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
                                   cap_multiplier,
                                   min_val,
                                   num_boundary = 0,
                                   num_capped = 0,
                                   num_failed = 0) {
  list(nuisance_bool_covariates_as_library = bool_covariates_as_library,
       nuisance_bool_library_includes_interept = bool_library_includes_interept,
       nuisance_bool_use_log = bool_use_log,
       nuisance_cap_multiplier = cap_multiplier,
       nuisance_min_val = min_val,
       nuisance_num_boundary = num_boundary,
       nuisance_num_capped = num_capped,
       nuisance_num_failed = num_failed)
}


