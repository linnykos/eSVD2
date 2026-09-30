# The fused per-gene path against the matrix path -- T-PGENE-01 .. T-PGENE-04.
#
# `compute_test_per_gene()` repeats, gene by gene, what `compute_posterior()`,
# `compute_test_statistic()` and `compute_pvalue()` do on whole matrices. The
# two paths must agree on every setting that changes which covariate columns
# form the library size, because each path had its own copy of that rule.
#
# The two paths must also refuse the same inputs (T-PGENE-02, -03): before
# 1.2.0 the per-gene path ran on a combination the matrix path refuses, and
# the two treated `nuisance_lower_quantile = NULL` differently.
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

# ---- refusals ----------------------------------------------------------------

# The message of the error `expr` raises, or NA when it raises none.
.error_message <- function(expr){
  tryCatch({
    force(expr)
    NA_character_
  }, error = function(cnd){conditionMessage(cnd)})
}

## [invariant] the matrix path refuses `bool_adjust_covariates` beside
## `bool_covariates_as_library` (T-POST-06); the per-gene path ran on it.
## Both must refuse it, with the same message, naming both arguments.
test_that("T-PGENE-02: both paths refuse bool_adjust_covariates with bool_covariates_as_library", {
  esvd_obj <- .small_esvd_obj()

  per_gene_msg <- .error_message(compute_test_per_gene(
    input_obj = esvd_obj,
    bool_adjust_covariates = TRUE,
    bool_covariates_as_library = TRUE
  ))
  matrix_msg <- .error_message(compute_posterior(
    input_obj = esvd_obj,
    bool_adjust_covariates = TRUE,
    bool_covariates_as_library = TRUE
  ))

  expect_false(is.na(per_gene_msg))
  expect_match(per_gene_msg, "bool_adjust_covariates", fixed = TRUE)
  expect_match(per_gene_msg, "bool_covariates_as_library", fixed = TRUE)
  expect_identical(per_gene_msg, matrix_msg)
})

## [invariant] each nonsensical setting is refused by both paths, by name,
## one setting at a time; the documented values (`NULL` for "off", and
## `alpha_max = Inf`) are accepted by both.
test_that("T-PGENE-03: both paths refuse the same nonsensical settings, naming them", {
  esvd_obj <- .small_esvd_obj()

  bad_list <- list(
    list(arg = "bool_adjust_covariates", value = NA),
    list(arg = "bool_covariates_as_library", value = "TRUE"),
    list(arg = "bool_stabilize_underdispersion", value = c(TRUE, FALSE)),
    list(arg = "alpha_max", value = 0),
    list(arg = "alpha_max", value = -1),
    list(arg = "alpha_max", value = NA_real_),
    list(arg = "library_min", value = 0),
    list(arg = "library_min", value = -1),
    list(arg = "library_min", value = Inf),
    list(arg = "library_min", value = NA_real_),
    list(arg = "nuisance_lower_quantile", value = -0.1),
    list(arg = "nuisance_lower_quantile", value = 1.5),
    list(arg = "nuisance_lower_quantile", value = NA_real_),
    list(arg = "nuisance_lower_quantile", value = c(0.1, 0.2)),
    list(arg = "pseudocount", value = -1),
    list(arg = "pseudocount", value = NA_real_)
  )

  for(bad in bad_list){
    label <- paste0(bad$arg, " = ", paste0(deparse(bad$value), collapse = ""))
    arg_list <- stats::setNames(list(bad$value), bad$arg)

    per_gene_msg <- .error_message(do.call(
      compute_test_per_gene, c(list(input_obj = esvd_obj), arg_list)
    ))
    matrix_msg <- .error_message(do.call(
      compute_posterior, c(list(input_obj = esvd_obj), arg_list)
    ))

    expect_false(is.na(per_gene_msg), info = paste0("per gene: ", label))
    expect_false(is.na(matrix_msg), info = paste0("matrix: ", label))
    expect_match(per_gene_msg, paste0("`", bad$arg, "`"), fixed = TRUE,
                 info = paste0("per gene: ", label))
    expect_match(matrix_msg, paste0("`", bad$arg, "`"), fixed = TRUE,
                 info = paste0("matrix: ", label))
  }

  good_list <- list(list(alpha_max = NULL),
                    list(alpha_max = Inf),
                    list(library_min = NULL),
                    list(nuisance_lower_quantile = NULL),
                    list(pseudocount = NULL))
  for(good in good_list){
    label <- paste0(names(good), " = ", deparse(good[[1]]))
    expect_no_error(.muffle_locfdr_fallback(do.call(
      compute_test_per_gene, c(list(input_obj = esvd_obj), good)
    )), message = paste0("per gene: ", label))
    expect_no_error(do.call(
      compute_posterior, c(list(input_obj = esvd_obj), good)
    ), message = paste0("matrix: ", label))
  }
})

## [invariant] `nuisance_lower_quantile = NULL` means no floor, on both paths:
## the quantile at 0 is the smallest rate, so flooring there changes nothing.
## The matrix path used to compute `stats::quantile(x, probs = NULL)`, which
## is `numeric(0)`, so `pmax()` emptied `nuisance_vec`.
test_that("T-PGENE-04: nuisance_lower_quantile = NULL skips the floor on both paths", {
  esvd_obj <- .small_esvd_obj()
  alpha_max <- 2 * max(esvd_obj$dat)

  .run_both <- function(nuisance_lower_quantile){
    matrix_obj <- compute_posterior(input_obj = esvd_obj,
                                    alpha_max = alpha_max,
                                    nuisance_lower_quantile = nuisance_lower_quantile)
    matrix_obj <- compute_test_statistic(input_obj = matrix_obj, verbose = 0)
    per_gene_obj <- .muffle_locfdr_fallback(compute_test_per_gene(
      input_obj = esvd_obj,
      alpha_max = alpha_max,
      nuisance_lower_quantile = nuisance_lower_quantile
    ))
    list(matrix = matrix_obj$teststat_vec, per_gene = per_gene_obj$teststat_vec)
  }
  null_list <- .run_both(NULL)
  zero_list <- .run_both(0)

  expect_equal(length(null_list$matrix), ncol(esvd_obj$dat))
  expect_equal(null_list$per_gene, null_list$matrix, tolerance = 1e-8)
  expect_equal(null_list$matrix, zero_list$matrix, tolerance = 1e-12)
  expect_equal(null_list$per_gene, zero_list$per_gene, tolerance = 1e-12)
})

## Finding of the code review (2026-09-29). `eSVD()` reaches the posterior only
## after the initialization, both rounds of `opt_esvd` and the nuisance
## estimate, so a setting the posterior refuses must be refused on entry, or
## a long fit is lost to it. [invariant] nothing is fitted: the intermediate
## file `eSVD()` writes after the initialization does not exist.
test_that("T-PGENE-05: eSVD() refuses the posterior settings before fitting anything", {
  skip_if_not_installed("SeuratObject")
  dat <- .tiny_counts()

  bad_list <- list(
    list(bool_adjust_covariates = TRUE, bool_covariates_as_library = TRUE),
    list(library_min = 0),
    list(pseudocount = -1)
  )
  for(bad in bad_list){
    label <- paste0(names(bad), " = ", unlist(bad), collapse = ", ")
    save_path <- tempfile(fileext = ".RData")
    msg <- .error_message(do.call(.esvd_run,
                                  c(list(dat = dat,
                                         intermediate_save = save_path),
                                    bad)))

    expect_false(is.na(msg), info = label)
    expect_match(msg, paste0("`", names(bad)[1], "`"), fixed = TRUE,
                 info = label)
    expect_false(file.exists(save_path), info = label)
  }
})
