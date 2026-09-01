#' Compute the degree of freedom
#'
#' This is an intermediary function used in \code{compute_pvalue}
#'
#' @param input_obj \code{eSVD} object outputed from \code{compute_test_statistic}.
#' @param min_cells_per_individual Minimum number of cells an individual must
#'                                 contribute; see \code{compute_test_statistic}.
#'
#' @return a named vector of degree-of-freedom values, one for each gene
#' @noRd
.compute_df <- function(input_obj, min_cells_per_individual = 3){
  # The Welch degrees of freedom are recomputed from the posterior matrices,
  # which need `dat` for their dimensions. `eSVD(bool_diet = TRUE)` removes
  # `dat` and takes the fused `compute_test_per_gene` route instead, so
  # reaching here without it means the two were mixed.
  if(is.null(input_obj[["dat"]])){
    stop("`input_obj$dat` is missing, and `compute_pvalue` needs it to ",
         "recompute the per-individual statistics. An object built with ",
         "`eSVD(bool_diet = TRUE)` has already had its p-values computed by ",
         "`compute_test_per_gene`; otherwise keep `dat` on the object")
  }
  stopifnot(all(!is.null(input_obj[["case_control"]])) && all(input_obj[["case_control"]] %in% c(0,1)) && length(input_obj[["case_control"]]) == nrow(input_obj[["dat"]]),
            all(!is.null(input_obj[["individual"]])) && all(is.factor(input_obj[["individual"]])) && length(input_obj[["individual"]]) == nrow(input_obj[["dat"]]))

  cc_vec <- input_obj[["case_control"]]
  cc_levels <- sort(unique(cc_vec), decreasing = F)
  stopifnot(length(cc_levels) == 2)
  control_idx <- which(cc_vec == cc_levels[1])
  case_idx <- which(cc_vec == cc_levels[2])

  latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
  posterior_mean_mat <- .get_object(eSVD_obj = input_obj, what_obj = "posterior_mean_mat", which_fit = latest_Fit)
  posterior_var_mat <- .get_object(eSVD_obj = input_obj, what_obj = "posterior_var_mat", which_fit = latest_Fit)

  individual_vec <- input_obj[["individual"]]
  control_individuals <- unique(individual_vec[control_idx])
  case_individuals <- unique(individual_vec[case_idx])

  # The Welch-Satterthwaite denominator below contains (v/n)^2/(n-1), so an arm
  # with one individual makes `df_vec` zero and every downstream `stats::pt()`
  # NaN. Stopping here names the cause; without it the failure surfaces much
  # later as "missing values and NaN's not allowed".
  .check_cohort_is_testable(case_individuals = case_individuals,
                            control_individuals = control_individuals,
                            individual_vec = individual_vec,
                            min_cells_per_individual = min_cells_per_individual)

  tmp <- .determine_individual_indices(case_individuals = case_individuals,
                                               control_individuals = control_individuals,
                                               individual_vec = individual_vec)
  all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
  avg_mat <- .construct_averaging_matrix(idx_list = all_indiv_idx,
                                                 n = nrow(posterior_mean_mat))
  avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
  avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

  case_row_idx <- 1:length(case_individuals)
  control_row_idx <- (length(case_individuals)+1):nrow(avg_posterior_mean_mat)
  case_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[case_row_idx,,drop = F])
  control_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[control_row_idx,,drop = F])
  case_gaussian_var <- .compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[case_row_idx,,drop = F],
    avg_posterior_var_mat = avg_posterior_var_mat[case_row_idx,,drop = F]
  )
  control_gaussian_var <- .compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[control_row_idx,,drop = F],
    avg_posterior_var_mat = avg_posterior_var_mat[control_row_idx,,drop = F]
  )

  n1 <- length(case_individuals)
  n2 <- length(control_individuals)

  # see https://www.theopeneducator.com/doe/hypothesis-Testing-Inferential-Statistics-Analysis-of-Variance-ANOVA/Two-Sample-T-Test-Unequal-Variance
  numerator_vec <- (case_gaussian_var/n1 + control_gaussian_var/n2)^2
  denominator_vec <- (case_gaussian_var/n1)^2/(n1-1) + (control_gaussian_var/n2)^2/(n2-1)
  df_vec <- numerator_vec/denominator_vec
  names(df_vec) <- names(case_gaussian_var)

  df_vec
}

#' Map Welch t-statistics to Gaussian statistics
#'
#' Implements \eqn{\hat{Z}_j = \Phi^{-1}(F_{df_j}(\hat{T}_j))}. The composition
#' is evaluated on the log scale and mirrored through zero, because
#' \code{qnorm(pt(t, df))} saturates at exactly \code{+Inf} once
#' \code{pt(t, df)} rounds to \code{1} (\code{t = 40, df = 18} is enough),
#' whereas the lower tail stays finite far longer. One strongly up-regulated
#' gene used to knock the whole empirical null over to its crudest fallback
#' this way (CRAN_READINESS.md 1.1).
#'
#' @param teststat_vec  Numeric vector of t-statistics.
#' @param df_vec        Numeric vector of degrees of freedom, same length.
#'
#' @returns Numeric vector with \code{names(teststat_vec)}.
#' @noRd
.t_to_gaussian <- function(teststat_vec, df_vec){
  stopifnot(length(teststat_vec) == length(df_vec))

  lower_log_prob <- stats::pt(-abs(teststat_vec), df = df_vec, log.p = TRUE)
  gaussian_vec <- sign(teststat_vec) * -stats::qnorm(lower_log_prob, log.p = TRUE)
  names(gaussian_vec) <- names(teststat_vec)

  gaussian_vec
}

#' Compute p-values
#'
#' Converts the Welch statistics in \code{input_obj$teststat_vec} to Gaussian
#' statistics via their Welch-Satterthwaite degrees of freedom, fits an
#' empirical null with \code{multtest}, and returns two-sided p-values on the
#' \eqn{-\log_{10}} scale with Benjamini-Hochberg adjustment.
#'
#' @param input_obj   \code{eSVD} object outputed from \code{compute_test_statistic}.
#' @param min_cells_per_individual  Minimum number of cells an individual must
#'                    contribute; see \code{compute_test_statistic}.
#' @param verbose     Integer.
#' @param ...         Additional parameters.
#'
#' @return \code{eSVD} object with added element \code{"pvalue_list"}, a list
#' with \code{df_vec}, \code{fdr_vec}, \code{gaussian_teststat},
#' \code{log10pvalue} (\eqn{-\log_{10}} of the two-sided p-value),
#' \code{method} (which empirical-null estimator ran; see \code{multtest}),
#' \code{null_mean} and \code{null_sd}.
#' @export
compute_pvalue <- function(input_obj,
                           min_cells_per_individual = 3,
                           verbose = 0,
                           ...){
  if(is.null(input_obj[["teststat_vec"]])){
    stop("`input_obj` has no `teststat_vec`; run `compute_test_statistic` ",
         "before `compute_pvalue`")
  }

  df_vec <- .compute_df(input_obj = input_obj,
                        min_cells_per_individual = min_cells_per_individual)

  teststat_vec <- input_obj$teststat_vec
  stopifnot(length(teststat_vec) == length(df_vec))
  gaussian_teststat <- .t_to_gaussian(teststat_vec = teststat_vec,
                                      df_vec = df_vec)

  fdr_res <- multtest(gaussian_teststat)
  fdr_vec <- fdr_res$fdr_vec
  names(fdr_vec) <- names(gaussian_teststat)
  log10pvalue_vec <- fdr_res$logpvalue_vec
  names(log10pvalue_vec) <- names(teststat_vec)

  pvalue_list <- list(
    df_vec = df_vec,
    fdr_vec = fdr_vec,
    gaussian_teststat = gaussian_teststat,
    log10pvalue = log10pvalue_vec,
    method = fdr_res$method,
    null_mean = fdr_res$null_mean,
    null_sd = fdr_res$null_sd
  )

  input_obj[["pvalue_list"]] <- pvalue_list
  input_obj
}
