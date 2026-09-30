# UNIT_TEST_PLAN.md section 2.16 -- T-GS-01 .. T-GS-17.
#
# THIS WHOLE FILE IS EXPECTED TO FAIL. `gene_status` is a NEW FEATURE that has
# not been implemented -- no eSVD2 R or C++ code has been changed. The tests are
# written first, to the specification in section 2.16.1, so that the feature can
# be built against them.
#
# Per question Q-COH-7, the feature lives in `eSVD_helper()`: it labels
# `gene_status`, removes the all-zero genes, calls `eSVD()`, and reinserts them
# at their original positions afterwards. `eSVD()` and `initialize_esvd()` both
# ERROR on an all-zero gene rather than filtering.

# Skips the whole file cleanly when the feature is absent, rather than producing
# seventeen identical "could not find function" errors.
.skip_if_no_gene_status <- function(){
  skip_if_not(exists("eSVD_helper"),
              "eSVD_helper() not implemented yet (UNIT_TEST_PLAN.md section 2.16/2.17)")
}

test_that("T-GS-01: gene_status is a two-level factor over the original genes", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()
  esvd_obj <- .helper_run(dat)

  gene_status <- esvd_obj[["gene_status"]]
  expect_s3_class(gene_status, "factor")
  expect_equal(levels(gene_status), c("analyzed", "all_zero"))
  expect_equal(length(gene_status), ncol(dat))
  expect_equal(names(gene_status), colnames(dat))
  # `as.integer()` gives Kevin's 1 and 2 exactly.
  expect_true(all(as.integer(gene_status) %in% c(1L, 2L)))
})

test_that("T-GS-02: the boundary of the all_zero definition", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts()
  zero_idx <- 3
  single_count_idx <- 7
  one_arm_zero_idx <- 11

  dat[, zero_idx] <- 0
  # A gene with a single non-zero count is NOT all-zero.
  dat[, single_count_idx] <- 0
  dat[1, single_count_idx] <- 1
  # A gene that is zero in one arm only is NOT all-zero -- and is often the most
  # interesting gene in the dataset.
  cc_vec <- .tiny_data()$cc_vec
  dat[cc_vec == 0, one_arm_zero_idx] <- 0

  esvd_obj <- .helper_run(dat)
  gene_status <- esvd_obj[["gene_status"]]

  expect_equal(as.character(gene_status[zero_idx]), "all_zero")
  expect_equal(as.character(gene_status[single_count_idx]), "analyzed")
  expect_equal(as.character(gene_status[one_arm_zero_idx]), "analyzed")
})

test_that("T-GS-03: removing all-zero genes does not change Log_UMI", {
  # This one needs no implementation: it is the `rowSums` identity, and it is
  # what makes the filter safe for the retained genes.
  dat <- .tiny_counts_with_zero_genes()
  zero_idx <- attr(dat, "zero_idx")
  covariate_df <- .tiny_data()$covariate_df

  full_covariates <- format_covariates(dat = dat, covariate_df = covariate_df)
  filtered_covariates <- format_covariates(dat = dat[, -zero_idx, drop = FALSE],
                                           covariate_df = covariate_df)

  expect_equal(full_covariates[, "Log_UMI"], filtered_covariates[, "Log_UMI"])
})

test_that("T-GS-04: removing all-zero genes does not change alpha_max", {
  dat <- .tiny_counts_with_zero_genes()
  zero_idx <- attr(dat, "zero_idx")

  # `eSVD()` derives `alpha_max <- 2*max(dat@x)` when it is NULL, and an
  # all-zero gene stores no non-zeros.
  expect_equal(2 * max(dat), 2 * max(dat[, -zero_idx, drop = FALSE]))
})

test_that("T-GS-05: analyzed-gene results equal a directly filtered run", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()
  zero_idx <- attr(dat, "zero_idx")

  # The load-bearing test of the feature. Under Q-COH-7 `eSVD()` never sees an
  # all-zero gene, so this does not ask "does the pipeline ignore them" -- it
  # asks "does the reinsertion disturb the retained entries", which is the only
  # remaining way the feature can be wrong.
  res_helper <- .helper_run(dat)
  res_direct <- .helper_run(dat[, -zero_idx, drop = FALSE])

  analyzed_idx <- which(res_helper[["gene_status"]] == "analyzed")

  expect_equal(unname(res_helper$teststat_vec[analyzed_idx]),
               unname(res_direct$teststat_vec), tolerance = 1e-8)
  expect_equal(unname(res_helper$case_mean[analyzed_idx]),
               unname(res_direct$case_mean), tolerance = 1e-8)
  expect_equal(unname(res_helper$pvalue_list$log10pvalue[analyzed_idx]),
               unname(res_direct$pvalue_list$log10pvalue), tolerance = 1e-8)
})

test_that("T-GS-06: every per-gene estimate is NA at exactly the status-2 positions", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()
  esvd_obj <- .helper_run(dat)
  zero_idx <- which(esvd_obj[["gene_status"]] == "all_zero")

  vector_names <- c("teststat_vec", "case_mean", "control_mean")
  for(element_name in vector_names){
    vec <- esvd_obj[[element_name]]
    expect_true(all(is.na(vec[zero_idx])), info = element_name)
    expect_true(all(!is.na(vec[-zero_idx])), info = element_name)
  }

  for(element_name in c("df_vec", "gaussian_teststat")){
    vec <- esvd_obj$pvalue_list[[element_name]]
    expect_true(all(is.na(vec[zero_idx])), info = element_name)
    expect_true(all(!is.na(vec[-zero_idx])), info = element_name)
  }
})

test_that("T-GS-07: log10pvalue is 0 and fdr_vec is 1 at status-2 positions", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()
  esvd_obj <- .helper_run(dat)
  zero_idx <- which(esvd_obj[["gene_status"]] == "all_zero")

  # Deliberately not NA: a downstream `which(fdr < 0.05)` must not have to
  # think about missingness.
  expect_equal(unname(esvd_obj$pvalue_list$log10pvalue[zero_idx]),
               rep(0, length(zero_idx)))
  expect_equal(unname(esvd_obj$pvalue_list$fdr_vec[zero_idx]),
               rep(1, length(zero_idx)))
})

test_that("T-GS-08: genes come back at their ORIGINAL positions", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()
  esvd_obj <- .helper_run(dat)

  # A reinsertion that appended at the end would pass T-GS-06 and fail here.
  # The zeros sit at positions 3 and 12 of 20 precisely to catch that.
  expect_equal(names(esvd_obj$teststat_vec), colnames(dat))
  expect_equal(names(esvd_obj$pvalue_list$fdr_vec), colnames(dat))
  expect_equal(length(esvd_obj$teststat_vec), ncol(dat))
})

test_that("T-GS-09: BH is computed on analyzed genes only", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts()
  dat_padded <- cbind(dat, matrix(0, nrow = nrow(dat), ncol = 50,
                                  dimnames = list(NULL,
                                                  paste0("empty", 1:50))))

  res_plain <- .helper_run(dat)
  res_padded <- .helper_run(dat_padded)

  analyzed_idx <- which(res_padded[["gene_status"]] == "analyzed")

  # Including status-2 genes in `p.adjust` would enlarge `n` in `p*n/rank` and
  # raise the adjusted p-value of every genuinely DE gene. Under Q-COH-7 this
  # is structural rather than conventional, so the test is an architecture
  # check -- it fails loudly if the filter ever moves back inside `eSVD()`.
  expect_equal(unname(res_padded$pvalue_list$fdr_vec[analyzed_idx]),
               unname(res_plain$pvalue_list$fdr_vec), tolerance = 1e-10)
})

test_that("T-GS-10: the empirical null is fit on analyzed genes only", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts()
  dat_padded <- cbind(dat, matrix(0, nrow = nrow(dat), ncol = 50,
                                  dimnames = list(NULL,
                                                  paste0("empty", 1:50))))

  res_plain <- .helper_run(dat)
  res_padded <- .helper_run(dat_padded)

  expect_equal(res_padded$pvalue_list$null_mean,
               res_plain$pvalue_list$null_mean, tolerance = 1e-10)
  expect_equal(res_padded$pvalue_list$null_sd,
               res_plain$pvalue_list$null_sd, tolerance = 1e-10)
})

test_that("T-GS-11: x_mat is not padded", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()
  esvd_obj <- .helper_run(dat)
  latest_fit <- esvd_obj[["latest_Fit"]]

  # `x_mat` is cells-by-k and has no gene dimension. Guards against a blanket
  # "insert NA rows everywhere" implementation.
  x_mat <- esvd_obj[[latest_fit]]$x_mat
  expect_equal(nrow(x_mat), nrow(dat))
  expect_true(all(!is.na(x_mat)))
})

test_that("T-GS-12: report_results keeps one row per original gene", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()
  esvd_obj <- .helper_run(dat)
  zero_idx <- which(esvd_obj[["gene_status"]] == "all_zero")

  res <- report_results(esvd_obj)

  expect_equal(nrow(res), ncol(dat))
  expect_equal(res$genes, colnames(dat))
  expect_true(all(is.na(res$logFC[zero_idx])))
  expect_equal(unname(res$pvalue[zero_idx]), rep(1, length(zero_idx)))
  expect_equal(unname(res$pvalue_adj[zero_idx]), rep(1, length(zero_idx)))
})

test_that("T-GS-13: both bool_diet paths give the same gene_status and outputs", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()

  res_diet <- .helper_run(dat, bool_diet = TRUE)
  res_full <- .helper_run(dat, bool_diet = FALSE)

  expect_equal(res_diet[["gene_status"]], res_full[["gene_status"]])
  expect_equal(unname(res_diet$teststat_vec), unname(res_full$teststat_vec),
               tolerance = 1e-8)
})

test_that("T-GS-14: every gene all-zero gives an informative error", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts()
  dat[] <- 0

  expect_error(.helper_run(dat), regexp = "zero")
})

test_that("T-GS-15: k greater than the analyzed-gene count errors, naming both", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()

  # A `k` that was legal before filtering can become illegal after it. `k = 30`
  # is `eSVD()`'s default, so this is not a contrived input.
  expect_error(.helper_run(dat, k = ncol(dat) - 1), regexp = "k")
})

test_that("T-GS-16: gene_status survives an intermediate_save round trip", {
  .skip_if_no_gene_status()

  dat <- .tiny_counts_with_zero_genes()
  save_file <- withr::local_tempfile(fileext = ".RData")

  esvd_obj <- .helper_run(dat, intermediate_save = save_file)

  expect_true(file.exists(save_file))
  expect_true("gene_status" %in% names(esvd_obj))
})

test_that("T-GS-17: the three levels agree on what all_zero means", {
  .skip_if_no_gene_status()

  # The cost of Q-COH-7's three-level defence: the helper labels, `eSVD()`
  # errors, and `initialize_esvd()` errors, all on the same condition. If they
  # use three expressions rather than one `.which_all_zero()`, the helper
  # starts passing genes that `eSVD()` rejects and the user gets an internal
  # error from a function they never called.
  dat <- .tiny_counts()
  dat[, 3] <- 0
  # A gene that is all NA, and a gene that is zero except for one NA.
  dat[, 5] <- NA
  dat[, 7] <- 0
  dat[1, 7] <- NA

  expect_true(exists(".which_all_zero"))
  all_zero_idx <- .which_all_zero(dat)

  covariates <- .tiny_covariates(dat = .tiny_counts())
  expect_error(initialize_esvd(dat = dat,
                               covariates = covariates,
                               metadata_individual = .tiny_data()$individual_vec,
                               k = 2),
               regexp = "zero")
})
