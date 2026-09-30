# The fused per-gene path against the matrix path -- T-PGENE-01.
#
# `compute_test_per_gene()` repeats, gene by gene, what `compute_posterior()`,
# `compute_test_statistic()` and `compute_pvalue()` do on whole matrices. The
# two paths must agree on every setting that changes which covariate columns
# form the library size, because each path had its own copy of that rule.
#
# Oracles are tagged per test: [invariant] an identity the design implies.

context("Test compute_test_per_gene (claude)")

## [invariant] the per-gene path equals the matrix path, on every setting of
## the library. `bool_library_includes_interept` is read from `param`, so it
## is set by running `estimate_nuisance()` again. Pinned before the inline
## library rule of `compute_test_per_gene()` was routed through
## `.nuisance_library_idx()`; T-POST-13 pins the matrix path against the
## column sets written by hand, so the two together pin this path against
## them as well. `bool_adjust_covariates` is refused beside
## `bool_covariates_as_library`, so that cell of the grid is left out.
test_that("T-PGENE-01: compute_test_per_gene equals the matrix path on every library setting", {
  esvd_obj <- .small_esvd_obj()
  alpha_max <- 2 * max(esvd_obj$dat)

  grid <- expand.grid(bool_covariates_as_library = c(TRUE, FALSE),
                      bool_library_includes_interept = c(TRUE, FALSE),
                      bool_adjust_covariates = c(FALSE, TRUE))
  grid <- grid[!(grid$bool_adjust_covariates & grid$bool_covariates_as_library), ]
  expect_equal(nrow(grid), 6)

  for(i in seq_len(nrow(grid))){
    label <- paste0(names(grid), " = ", unlist(grid[i, ]), collapse = ", ")
    obj_nuisance <- suppressWarnings(estimate_nuisance(
      input_obj = esvd_obj,
      bool_covariates_as_library = TRUE,
      bool_library_includes_interept = grid$bool_library_includes_interept[i]
    ))

    matrix_obj <- compute_posterior(
      input_obj = obj_nuisance,
      alpha_max = alpha_max,
      bool_adjust_covariates = grid$bool_adjust_covariates[i],
      bool_covariates_as_library = grid$bool_covariates_as_library[i],
      library_min = 0.1
    )
    matrix_obj <- compute_test_statistic(input_obj = matrix_obj, verbose = 0)
    matrix_obj <- .muffle_locfdr_fallback(compute_pvalue(input_obj = matrix_obj))

    per_gene_obj <- .muffle_locfdr_fallback(compute_test_per_gene(
      input_obj = obj_nuisance,
      alpha_max = alpha_max,
      bool_adjust_covariates = grid$bool_adjust_covariates[i],
      bool_covariates_as_library = grid$bool_covariates_as_library[i],
      library_min = 0.1
    ))

    for(what in c("teststat_vec", "case_mean", "control_mean", "case_var",
                  "control_var", "log2fc_vec", "log2fc_se_vec")){
      expect_equal(per_gene_obj[[what]], matrix_obj[[what]],
                   tolerance = 1e-8, info = paste0(label, ": ", what))
    }
    expect_equal(per_gene_obj$pvalue_list$df_vec,
                 matrix_obj$pvalue_list$df_vec,
                 tolerance = 1e-8, info = label)
    expect_equal(per_gene_obj$pvalue_list$gaussian_teststat,
                 matrix_obj$pvalue_list$gaussian_teststat,
                 tolerance = 1e-8, info = label)
  }

  # The grid moves the statistic: without this, the six agreements above
  # could all be the same number.
  obj_a <- .muffle_locfdr_fallback(compute_test_per_gene(
    input_obj = esvd_obj, alpha_max = alpha_max,
    bool_covariates_as_library = TRUE, library_min = 0.1))
  obj_b <- .muffle_locfdr_fallback(compute_test_per_gene(
    input_obj = esvd_obj, alpha_max = alpha_max,
    bool_covariates_as_library = FALSE, library_min = 0.1))
  expect_false(isTRUE(all.equal(obj_a$teststat_vec, obj_b$teststat_vec)))
})
