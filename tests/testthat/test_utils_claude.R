# UNIT_TEST_PLAN.md section 2.14 -- T-UTIL-01 .. T-UTIL-12.
#
# `utils.R`, `data_management.R` and `report_results.R` have zero coverage.
# `.get_object()` in particular is the accessor every pipeline function routes
# through, with fourteen branches and no tests.

test_that("T-UTIL-01: .mult_mat_vec and .mult_vec_mat are the identities they claim", {
  set.seed(10)
  mat <- matrix(stats::rnorm(20), nrow = 5, ncol = 4)
  vec_col <- stats::rnorm(4)
  vec_row <- stats::rnorm(5)

  # The whole point of these helpers is to be a fast identity, so assert the
  # identity rather than a stored answer.
  expect_equal(.mult_mat_vec(mat, vec_col), mat %*% diag(vec_col))
  expect_equal(.mult_vec_mat(vec_row, mat), diag(vec_row) %*% mat)
})

test_that("T-UTIL-01b: .mult_mat_vec rejects a length mismatch", {
  mat <- matrix(0, nrow = 5, ncol = 4)
  expect_error(.mult_mat_vec(mat, rep(1, 3)))
  expect_error(.mult_vec_mat(rep(1, 3), mat))
})

test_that("T-UTIL-02: .nonzero_col matches a dense recomputation", {
  set.seed(10)
  dense_mat <- matrix(stats::rpois(60, lambda = 0.6), nrow = 12, ncol = 5)
  dense_mat[, 3] <- 0
  sparse_mat <- methods::as(methods::as(dense_mat * 1.0, "dMatrix"),
                            "CsparseMatrix")

  for(col_idx in seq_len(ncol(dense_mat))){
    label <- paste0("column ", col_idx)
    expected_idx <- which(dense_mat[, col_idx] != 0)

    expect_equal(as.numeric(.nonzero_col(sparse_mat, col_idx,
                                         bool_value = FALSE)),
                 as.numeric(expected_idx),
                 info = label)
    expect_equal(as.numeric(.nonzero_col(sparse_mat, col_idx,
                                         bool_value = TRUE)),
                 as.numeric(dense_mat[expected_idx, col_idx]),
                 info = label)
  }

  # The all-zero column is the case that exercises the `val1 == val2` early
  # return, and it is common in real single-cell data.
  expect_equal(.nonzero_col(sparse_mat, 3, bool_value = FALSE), numeric(0))
  expect_equal(.nonzero_col(sparse_mat, 3, bool_value = TRUE), numeric(0))
})

test_that("T-UTIL-03: .combine_two_named_lists preserves a NULL value under its name", {
  list1 <- list(a = 1, b = 2)
  list2 <- list(c = NULL, d = 4)

  res <- .combine_two_named_lists(list1, list2)

  # The `TEMP_NAME` trick exists so that a NULL *value* keeps its name. If the
  # name were dropped, `.get_object(which_fit = "param")` later fails with
  # `"what_obj is not found" is not TRUE`, which names nothing.
  expect_true("c" %in% names(res))
  expect_null(res[["c"]])
  expect_equal(res[["d"]], 4)

  # An existing name is not overwritten.
  res2 <- .combine_two_named_lists(list(a = 1), list(a = 99, b = 2))
  expect_equal(res2[["a"]], 1)
  expect_equal(res2[["b"]], 2)
})

test_that("T-UTIL-04: .get_object returns the right thing for every branch", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  expect_equal(.get_object(esvd_obj, "dat", NULL), esvd_obj$dat)
  expect_equal(.get_object(esvd_obj, "covariates", NULL), esvd_obj$covariates)
  expect_equal(.get_object(esvd_obj, "latest_Fit", NULL), latest_fit)
  expect_equal(.get_object(esvd_obj, "x_mat", latest_fit),
               esvd_obj[[latest_fit]]$x_mat)
  expect_equal(.get_object(esvd_obj, "y_mat", latest_fit),
               esvd_obj[[latest_fit]]$y_mat)
  expect_equal(.get_object(esvd_obj, "z_mat", latest_fit),
               esvd_obj[[latest_fit]]$z_mat)
  expect_equal(.get_object(esvd_obj, "nuisance", latest_fit),
               esvd_obj[[latest_fit]]$nuisance_vec)
  expect_equal(.get_object(esvd_obj, "posterior_mean_mat", latest_fit),
               esvd_obj[[latest_fit]]$posterior_mean_mat)
  expect_equal(.get_object(esvd_obj, "posterior_var_mat", latest_fit),
               esvd_obj[[latest_fit]]$posterior_var_mat)

  # The `NULL` what_obj branch returns the whole fit.
  expect_s3_class(.get_object(esvd_obj, NULL, latest_fit), "eSVD_Fit")

  # And the param branch.
  expect_equal(.get_object(esvd_obj, "init_k", "param"), 2)
})

## Section 1.8: `stopifnot("what_obj is not found")` passes a CHARACTER to
## `stopifnot`, which is not a logical, so the error reads
## `"what_obj is not found" is not TRUE` and never names the key that was asked
## for. Expected to FAIL until that becomes a `stop()`.
test_that("T-UTIL-05: an unrecognized what_obj errors with a message naming the key", {
  esvd_obj <- .small_esvd_obj()

  expect_error(.get_object(esvd_obj, "not_a_real_key", "fit_First"),
               regexp = "not_a_real_key")
})

test_that("T-UTIL-06: .get_object guards which_fit for the initial_Reg keys", {
  esvd_obj <- .small_esvd_obj()

  expect_error(.get_object(esvd_obj, "z_mat1", "fit_First"))
  expect_error(.get_object(esvd_obj, "z_mat2", "fit_First"))
  expect_error(.get_object(esvd_obj, "log_pval", "fit_First"))
})

test_that("T-UTIL-07: report_results on an incomplete object messages and returns NULL", {
  esvd_obj <- .small_esvd_obj()
  esvd_obj$pvalue_list <- NULL

  expect_message(res <- report_results(esvd_obj))
  expect_null(res)
})

test_that("T-UTIL-08: report_results logFC and pvalue are the documented transforms", {
  esvd_obj <- .small_esvd_obj()
  res <- report_results(esvd_obj)

  expect_s3_class(res, "data.frame")
  expect_equal(res$logFC,
               unname(log2(esvd_obj$case_mean / esvd_obj$control_mean)))
  expect_equal(res$pvalue, unname(10^(-esvd_obj$pvalue_list$log10pvalue)))
  expect_equal(res$pvalue_adj, unname(esvd_obj$pvalue_list$fdr_vec))
  expect_equal(res$genes, names(esvd_obj$teststat_vec))
})

test_that("T-UTIL-09: fisher_test matches stats::fisher.test", {
  all_genes <- paste0("gene_", 1:100)

  # Several configurations rather than one, and an external oracle rather than
  # the hardcoded 1e-2 constant the existing test uses.
  config_grid <- expand.grid(size1 = c(10, 20, 40),
                             size2 = c(10, 30),
                             overlap = c(0, 3, 8))

  for(row_idx in seq_len(nrow(config_grid))){
    size1 <- config_grid$size1[row_idx]
    size2 <- config_grid$size2[row_idx]
    overlap <- config_grid$overlap[row_idx]
    if(overlap > min(size1, size2)) next
    label <- paste0("size1=", size1, " size2=", size2, " overlap=", overlap)

    set1_genes <- all_genes[seq_len(size1)]
    set2_genes <- c(all_genes[seq_len(overlap)],
                    all_genes[size1 + seq_len(size2 - overlap)])

    res <- fisher_test(set1_genes = set1_genes,
                       set2_genes = set2_genes,
                       all_genes = all_genes)

    contingency_mat <- matrix(
      c(overlap,
        size1 - overlap,
        size2 - overlap,
        length(all_genes) - size1 - size2 + overlap),
      nrow = 2, ncol = 2
    )
    expected <- stats::fisher.test(contingency_mat,
                                   alternative = "greater")$p.value

    expect_equal(res$pvalue, as.numeric(expected), tolerance = 1e-8,
                 info = label)
    expect_equal(res$set_overlap_len, overlap, info = label)
    expect_equal(res$set1_len, size1, info = label)
  }
})

## `fisher_test(verbose = 1)` computes a `paste0(...)` and discards it, so the
## whole verbose branch does nothing. Expected to FAIL.
test_that("T-UTIL-10: fisher_test(verbose = 1) actually emits its diagnostic", {
  all_genes <- paste0("gene_", 1:100)

  expect_output(fisher_test(set1_genes = all_genes[1:10],
                            set2_genes = all_genes[5:20],
                            all_genes = all_genes,
                            verbose = 1))
})

test_that("T-UTIL-11: fisher_test errors when a set is not a subset of all_genes", {
  all_genes <- paste0("gene_", 1:100)

  expect_error(fisher_test(set1_genes = c("not_a_gene"),
                           set2_genes = all_genes[1:10],
                           all_genes = all_genes))
})

test_that("T-UTIL-12: print.esvd_data_loader prints and returns invisibly", {
  set.seed(10)
  dat <- matrix(stats::rpois(20, lambda = 3) * 1.0, nrow = 5, ncol = 4)
  loader <- data_loader(dat)

  expect_output(print(loader))
  expect_invisible(print(loader))
})
