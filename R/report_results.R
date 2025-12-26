#' Report the DE results from eSVD
#'
#' @param input_obj \code{eSVD} object outputed from \code{compute_pvalue}.
#'
#' @returns a data frame
#' @export
report_results <- function(input_obj){
  if(all(c("pvalue_list", "case_mean", "control_mean", "teststat_vec") %in% names(input_obj))){
    pvalue <- 10^(-input_obj$pvalue_list$log10pvalue)
    pvalue_adj <- input_obj$pvalue_list$fdr_vec
    logFC <- log2(input_obj$case_mean / input_obj$control_mean)
    genes <- names(input_obj$teststat_vec)

    df <- data.frame(genes = genes,
                     logFC = logFC,
                     pvalue = pvalue,
                     pvalue_adj = pvalue_adj)
    rownames(df) <- df$genes
    return(df)

  } else {
    message("input_obj does not have all the results computed yet")
    invisible()
  }
}
