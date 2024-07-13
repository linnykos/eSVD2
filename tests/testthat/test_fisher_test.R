context("Test for Fisher enrichment")

test_that("fisher_test works", {
  all_genes <- as.character(1:30)
  set1_genes <- as.character(1:15)
  set2_genes <- as.character(c(1:12, 16:18))

  res <- fisher_test(set1_genes = set1_genes,
                     set2_genes = set2_genes,
                     all_genes = all_genes)

  stopifnot(abs(res$pvalue - 0.0014) <= 1e-2)
})
