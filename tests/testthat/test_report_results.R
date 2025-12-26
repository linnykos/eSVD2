context("Test report_results")

## report_results is correct

test_that("report_results works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  covariates <- .get_object(eSVD_obj = eSVD_obj, what_obj = "covariates", which_fit = NULL)
  cc_vec <- covariates[,"case_control_1"]
  cc_levels <- sort(unique(cc_vec), decreasing = F)
  control_idx <- which(cc_vec == cc_levels[1])
  case_idx <- which(cc_vec == cc_levels[2])

  individual_vec <- metadata[,"individual"]
  control_individuals <- as.character(unique(individual_vec[control_idx]))
  case_individuals <- as.character(unique(individual_vec[case_idx]))

  eSVD_obj <- compute_test_statistic(input_obj = eSVD_obj)
  eSVD_obj <- compute_pvalue(input_obj = eSVD_obj)

  res <- report_results(eSVD_obj)

  expect_true(sum(abs(stats::p.adjust(res$pvalue, method = "BH") - res$pvalue_adj)) <= 1e-6)
  expect_true(is.data.frame(res))
  expect_true(all(sort(names(res)) == c("genes", "logFC", "pvalue", "pvalue_adj")))

})
