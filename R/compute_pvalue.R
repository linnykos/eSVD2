#' Compute the degree of freedom
#'
#' This is an intermediary function used in \code{compute_pvalue}
#'
#' @param input_obj \code{eSVD} object outputed from \code{compute_test_statistic}.
#'
#' @return a named vector of degree-of-freedom values, one for each gene
.compute_df <- function(input_obj){
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

#' Compute p-values
#'
#' @param input_obj   \code{eSVD} object outputed from \code{compute_test_statistic}.
#' @param verbose     Integer.
#' @param ...         Additional parameters.
#'
#' @return \code{eSVD} object with added element \code{"pvalue_list"}
#' @export
compute_pvalue <- function(input_obj,
                           verbose = 0,
                           ...){
  df_vec <- .compute_df(input_obj = input_obj)
  names(df_vec) <- names(df_vec)

  teststat_vec <- input_obj$teststat_vec
  p <- length(teststat_vec)
  gaussian_teststat <- sapply(1:p, function(j){
    stats::qnorm(stats::pt(teststat_vec[j], df = df_vec[j]))
  })
  names(gaussian_teststat) <- names(teststat_vec)

  fdr_res <- multtest(gaussian_teststat)
  fdr_vec <- fdr_res$fdr_vec
  names(fdr_vec) <- names(gaussian_teststat)
  null_mean <- fdr_res$null_mean
  null_sd <- fdr_res$null_sd
  log10pvalue_vec <- sapply(gaussian_teststat, function(x){
    if(x < null_mean) {
      Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
    } else {
      Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
    }
  })
  log10pvalue_vec <- -(log10pvalue_vec/log(10) + log10(2))
  names(log10pvalue_vec) <- names(teststat_vec)

  pvalue_list <- list(
    df_vec = df_vec,
    fdr_vec = fdr_vec,
    gaussian_teststat = gaussian_teststat,
    log10pvalue = log10pvalue_vec,
    null_mean = null_mean,
    null_sd = null_sd
  )

  input_obj[["pvalue_list"]] <- pvalue_list
  input_obj
}
