#' Report the DE results from eSVD
#'
#' @param input_obj \code{eSVD} object outputed from \code{compute_pvalue}.
#'
#' @returns a data frame with one row per gene (in the order of
#' \code{input_obj$teststat_vec}) and columns \code{genes}, \code{logFC}
#' (\eqn{\log_2} of the case-to-control ratio of posterior means),
#' \code{log10pvalue} (\eqn{-\log_{10}} of the two-sided p-value),
#' \code{pvalue} and \code{pvalue_adj} (Benjamini-Hochberg). \code{pvalue}
#' underflows to \code{0} for any gene with \code{log10pvalue} above about
#' 308, so genes should be ranked by \code{log10pvalue}, which keeps the
#' distinction.
#' @export
report_results <- function(input_obj){
  if(all(c("pvalue_list", "case_mean", "control_mean", "teststat_vec") %in% names(input_obj))){
    log10pvalue <- input_obj$pvalue_list$log10pvalue
    pvalue <- 10^(-log10pvalue)
    pvalue_adj <- input_obj$pvalue_list$fdr_vec
    logFC <- log2(input_obj$case_mean / input_obj$control_mean)
    genes <- names(input_obj$teststat_vec)

    df <- data.frame(genes = genes,
                     logFC = unname(logFC),
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
