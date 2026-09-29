#' Report the DE results from eSVD
#'
#' @param input_obj \code{eSVD} object output from \code{compute_pvalue}.
#'
#' @returns a data frame with one row per gene (in the order of
#' \code{input_obj$teststat_vec}) and columns \code{genes}, \code{logFC}
#' (\eqn{\log_2} of the case-to-control ratio of posterior means),
#' \code{logFC_se} (the standard error of \code{logFC}, also on the
#' \eqn{\log_2} scale; see \code{compute_log_fold_change} for what it is a
#' standard error of and how it compares with those of 'DESeq2', 'dreamlet'
#' and 'NEBULA'),
#' \code{log10pvalue} (\eqn{-\log_{10}} of the two-sided p-value),
#' \code{pvalue} and \code{pvalue_adj} (Benjamini-Hochberg). \code{pvalue}
#' underflows to \code{0} for any gene with \code{log10pvalue} above about
#' 308, so genes should be ranked by \code{log10pvalue}, which keeps the
#' distinction. \code{logFC / logFC_se} is not the statistic behind
#' \code{pvalue}, which is computed on the linear scale and calibrated with
#' an empirical null. \code{logFC_se} is \code{NA}, with a warning, for an
#' object built before version 1.1.0, which does not store it.
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
#' esvd_obj <- compute_test_statistic(input_obj = esvd_obj)
#' esvd_obj <- compute_pvalue(input_obj = esvd_obj)
#' result_df <- report_results(esvd_obj)
#' utils::head(result_df[order(result_df$log10pvalue, decreasing = TRUE), ])
#' @export
report_results <- function(input_obj){
  if(all(c("pvalue_list", "case_mean", "control_mean", "teststat_vec") %in% names(input_obj))){
    log10pvalue <- input_obj$pvalue_list$log10pvalue
    pvalue <- 10^(-log10pvalue)
    pvalue_adj <- input_obj$pvalue_list$fdr_vec
    genes <- names(input_obj$teststat_vec)

    if(all(c("log2fc_vec", "log2fc_se_vec") %in% names(input_obj))){
      # Both are read from the object, so the two columns always describe
      # the same genes in the same way (NA together, never Inf beside NA).
      logFC <- input_obj[["log2fc_vec"]]
      logFC_se <- input_obj[["log2fc_se_vec"]]
      # The columns are joined to the others by position.
      stopifnot(length(logFC) == length(genes),
                length(logFC_se) == length(genes),
                identical(names(logFC), genes),
                identical(names(logFC_se), genes))
    } else {
      # An object saved by a version before 1.1.0 has the arm means and no
      # standard error. It still gets its data frame, with the column present
      # and empty, so code that reads the other columns keeps working.
      warning("`input_obj` has no `log2fc_vec` and `log2fc_se_vec`, so ",
              "`logFC_se` is NA. The object was built before eSVD2 1.1.0; ",
              "see `compute_log_fold_change`")
      logFC <- log2(input_obj$case_mean / input_obj$control_mean)
      logFC_se <- rep(NA_real_, length(logFC))
    }

    df <- data.frame(genes = genes,
                     logFC = unname(logFC),
                     logFC_se = unname(logFC_se),
                     log10pvalue = unname(log10pvalue),
                     pvalue = unname(pvalue),
                     pvalue_adj = unname(pvalue_adj))
    rownames(df) <- df$genes
    return(df)

  } else {
    message("input_obj does not have all the results computed yet")
    invisible()
  }
}
