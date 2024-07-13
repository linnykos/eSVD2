#' Fisher's exact test
#' 
#' Perform a Fisher's exact test. This formula is from https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
#'
#' @param set1_genes A vector of characters of your resulting gene list
#' @param set2_genes A vector of characters of a gene list you are comparing against
#' @param all_genes A vector of characters of all the genes in your "universe". This should include all the genes in \code{set1_genes} and \code{set2_genes}.
#' @param verbose Numeric
#'
#' @return A list of: \code{pvalue} for the p-value,
#' \code{set_bg_len} for the number of genes in \code{all_genes} but not in \code{set2_genes},
#' \code{set_overlap_len} for the number of genes in both \code{set1_genes} and \code{set2_genes},
#' \code{set1_len} for number of genes in \code{set1_genes}, and
#' \code{set2_len} for number of genes in \code{set2_genes}.
#' @export
fisher_test <- function(set1_genes,
                        set2_genes,
                        all_genes,
                        verbose = 0){
  stopifnot(all(set1_genes %in% all_genes),
            all(set2_genes %in% all_genes),
            length(set1_genes) > 0,
            length(set2_genes) > 0,
            length(all_genes) > 0,
            is.character(set1_genes),
            is.character(set2_genes),
            is.character(all_genes))
  
  m <- length(set2_genes)
  n <- length(all_genes) - m
  k <- length(set1_genes)
  x <- length(intersect(set1_genes, set2_genes))
  
  if(verbose > 0){
    paste0("#Other: ", m, 
           ", #Bg: ", n, 
           ", #Select: ", k, 
           ", #Intersect: ", x,
           ", #Expected: ", round(m*(k/length(all_genes)),1) )
  }
  
  pvalue <- sum(sapply(x:k, function(i){
    stats::dhyper(x = i, m = m, n = n, k = k, log = F)
  }))
  
  list(pvalue = pvalue,
       set1_len = k,
       set2_len = m,
       set_bg_len = n,
       set_overlap_len = x)
}