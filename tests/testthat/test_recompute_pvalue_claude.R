# `recompute_pvalue()` -- T-REDO-01 .. T-REDO-12.
#
# The factorization does not depend on the cap on the nuisance rate, and the
# maximum-likelihood rates are stored beside the capped ones. So the test can
# be redone at another cap from a fitted object: apply the cap to the stored
# rates, then repeat the posterior, the test statistic and the p-values with
# the settings recorded in `param`.
#
# The oracle throughout is the same analysis run FROM SCRATCH at the new cap.
# The two go through the same code path on the same inputs, so they are
# compared at a tolerance of 1e-10 and not at a statistical one: the pipeline
# amplifies a perturbation of 1e-9 in the fit to order 1 in a statistic, so a
# looser tolerance would hide a redo that rebuilt slightly different
# covariates.
#
# An object built with `bool_diet = TRUE` has neither counts nor covariates;
# they are rebuilt from the Seurat object and the arguments `eSVD()` recorded.
#
# Fitted fixtures are cached in `.fixture_cache` (helper-fixtures.R).

# ---- fixtures ---------------------------------------------------------------

# Everything the redo is documented to replace.
.redo_result_elements <- function(){
  c("teststat_vec", "case_mean", "control_mean", "case_var", "control_var",
    "log2fc_vec", "log2fc_se_vec", "pvalue_list")
}

.expect_same_results <- function(res, reference, label){
  for(element_name in .redo_result_elements()){
    expect_equal(res[[element_name]], reference[[element_name]],
                 tolerance = 1e-10,
                 info = paste0(label, ", ", element_name))
  }

  fit <- res[[res$latest_Fit]]
  fit_reference <- reference[[reference$latest_Fit]]
  for(element_name in c("nuisance_vec", "nuisance_mle_vec",
                        "nuisance_library_median_vec", "nuisance_status",
                        "posterior_mean_mat", "posterior_var_mat")){
    expect_equal(fit[[element_name]], fit_reference[[element_name]],
                 tolerance = 1e-10,
                 info = paste0(label, ", ", element_name))
  }

  for(element_name in c("nuisance_cap_multiplier", "nuisance_num_capped",
                        "nuisance_num_boundary")){
    expect_equal(res$param[[element_name]], reference$param[[element_name]],
                 info = paste0(label, ", ", element_name))
  }
}

# `eSVD()` on F-TINY, cached by its settings.
.redo_esvd_obj <- function(bool_diet, cap_multiplier){
  cache_key <- paste0("redo_esvd_", bool_diet, "_", cap_multiplier)
  if(is.null(.fixture_cache[[cache_key]])){
    .fixture_cache[[cache_key]] <- .esvd_run(dat = .tiny_counts(),
                                             bool_diet = bool_diet,
                                             cap_multiplier = cap_multiplier)
  }
  .fixture_cache[[cache_key]]
}

# A cohort that `eSVD_helper()` has to filter in both directions: two all-zero
# genes (positions 3 and 12), and one individual left with two cells, which
# the helper drops. The Seurat object returned is the UNFILTERED one, which is
# what a user holds and passes to the redo.
.redo_helper_seurat <- function(){
  dat <- .tiny_counts_with_zero_genes()
  attr(dat, "zero_idx") <- NULL
  seurat_obj <- .tiny_seurat(dat = dat)

  individual_vec <- as.character(seurat_obj@meta.data[, "Individual"])
  sparse_individual <- "indiv_5"
  drop_idx <- which(individual_vec == sparse_individual)[-(1:2)]

  seurat_obj[, SeuratObject::Cells(seurat_obj)[-drop_idx]]
}

.redo_helper_obj <- function(bool_diet, cap_multiplier){
  cache_key <- paste0("redo_helper_", bool_diet, "_", cap_multiplier)
  if(is.null(.fixture_cache[[cache_key]])){
    seurat_obj <- .redo_helper_seurat()
    .fixture_cache[[cache_key]] <- suppressWarnings(
      eSVD_helper(batch_var_prefix = NULL,
                  case_control_levels = c("0", "1"),
                  case_control_var = "CC",
                  categorical_vars = c("Sex"),
                  id_var = "Individual",
                  numerical_vars = "Age",
                  seurat_obj = seurat_obj,
                  bool_diet = bool_diet,
                  cap_multiplier = cap_multiplier,
                  k = 2)
    )
  }
  .fixture_cache[[cache_key]]
}

# ---- an object that still has its counts ------------------------------------

test_that("T-REDO-01: redoing at the cap the object was built with changes nothing", {
  esvd_obj <- .small_esvd_obj()
  expect_equal(esvd_obj$param$nuisance_cap_multiplier, 10)

  res <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 10)
  )

  .expect_same_results(res = res, reference = esvd_obj, label = "same cap")
  expect_identical(class(res), class(esvd_obj))
})

## [oracle] the stages after the fit, run again by hand at the new cap with
## the settings `.small_esvd_obj()` was built with.
##
## The fixture must be one on which the cap matters, or the comparison holds
## because nothing moved. That is asserted first, at each cap.
test_that("T-REDO-02: redoing at another cap equals running the later stages again", {
  esvd_obj <- .small_esvd_obj()
  dat_list <- .tiny_data()

  for(cap_multiplier in c(0.5, 1, Inf)){
    label <- paste0("cap_multiplier = ", cap_multiplier)

    reference <- suppressWarnings(
      estimate_nuisance(input_obj = esvd_obj,
                        bool_covariates_as_library = TRUE,
                        cap_multiplier = cap_multiplier)
    )
    reference <- compute_posterior(input_obj = reference,
                                   alpha_max = 2 * max(dat_list$dat),
                                   bool_covariates_as_library = TRUE,
                                   library_min = 0.1)
    reference <- compute_test_statistic(input_obj = reference)
    reference <- suppressWarnings(compute_pvalue(input_obj = reference))

    res <- .muffle_locfdr_fallback(
      recompute_pvalue(input_obj = esvd_obj, cap_multiplier = cap_multiplier)
    )

    expect_false(isTRUE(all.equal(res$teststat_vec, esvd_obj$teststat_vec)),
                 info = label)
    .expect_same_results(res = res, reference = reference, label = label)
  }
})

## Finding of the code review (2026-09-29). The redo repeats the posterior
## with the settings in `param`. `compute_posterior()` recorded them with
## `.combine_two_named_lists()`, which keeps an entry that is already there,
## so after a second call `param` described the FIRST call and the redo, at
## the same cap, returned statistics up to 0.99 away from the object's own.
##
## The settings of the second call must MOVE the statistic on this fixture,
## or a redo at the stale settings agrees with it and the test says nothing.
## That is asserted first. (`alpha_max = 5, library_min = 2` does not: every
## library size of the fixture is above 2.4.)
test_that("T-REDO-12: the redo uses the settings of the last compute_posterior call", {
  esvd_obj <- .small_esvd_obj()
  expect_equal(esvd_obj$param$posterior_library_min, 0.1)
  expect_equal(esvd_obj$param$posterior_pseudocount, 0)

  esvd_second <- compute_posterior(input_obj = esvd_obj,
                                   alpha_max = esvd_obj$param$posterior_alpha_max,
                                   bool_covariates_as_library = TRUE,
                                   library_min = 20,
                                   pseudocount = 1)
  esvd_second <- compute_test_statistic(input_obj = esvd_second)
  esvd_second <- suppressWarnings(compute_pvalue(input_obj = esvd_second))

  expect_true(max(abs(esvd_second$teststat_vec - esvd_obj$teststat_vec)) > 0.1)
  expect_equal(esvd_second$param$posterior_library_min, 20)
  expect_equal(esvd_second$param$posterior_pseudocount, 1)

  res <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = esvd_second, cap_multiplier = 10)
  )
  .expect_same_results(res = res, reference = esvd_second,
                       label = "after a second compute_posterior")
})

## The same for the per-gene path followed by the matrix path.
test_that("T-REDO-12b: compute_posterior after compute_test_per_gene records its own settings", {
  esvd_obj <- .small_esvd_obj()

  esvd_first <- .muffle_locfdr_fallback(
    compute_test_per_gene(input_obj = esvd_obj,
                          alpha_max = 1000,
                          library_min = 0.1)
  )
  expect_equal(esvd_first$param$posterior_library_min, 0.1)

  esvd_second <- compute_posterior(input_obj = esvd_first,
                                   alpha_max = 5,
                                   library_min = 20)
  expect_equal(esvd_second$param$posterior_alpha_max, 5)
  expect_equal(esvd_second$param$posterior_library_min, 20)
})

## The redo may be repeated: it reads the rates before the cap, which it
## never changes, so going to 1 and back to 10 returns the starting point.
test_that("T-REDO-03: redoing twice depends only on the last cap", {
  esvd_obj <- .small_esvd_obj()

  res <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1)
  )
  res <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = res, cap_multiplier = 10)
  )

  .expect_same_results(res = res, reference = esvd_obj, label = "1 then 10")
})

test_that("T-REDO-04: the redo leaves the fit and the uncapped rates untouched and records the new cap", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  res <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1)
  )

  for(element_name in c("x_mat", "y_mat", "z_mat", "nuisance_mle_vec",
                        "nuisance_library_median_vec", "gene_mean_count_vec",
                        "gene_sparsity_vec")){
    expect_identical(res[[latest_fit]][[element_name]],
                     esvd_obj[[latest_fit]][[element_name]],
                     info = element_name)
  }
  expect_identical(res$dat, esvd_obj$dat)
  expect_identical(res$covariates, esvd_obj$covariates)

  expect_equal(res$param$nuisance_cap_multiplier, 1)
  fit <- res[[latest_fit]]
  expect_equal(res$param$nuisance_num_capped,
               sum(fit$nuisance_mle_vec > fit$nuisance_library_median_vec))
  expect_equal(as.numeric(fit$nuisance_vec),
               as.numeric(pmin(fit$nuisance_mle_vec,
                               fit$nuisance_library_median_vec)),
               tolerance = 1e-12)

  # Every other setting is as it was.
  other_names <- setdiff(names(esvd_obj$param),
                         c("nuisance_cap_multiplier", "nuisance_num_capped",
                           "nuisance_num_boundary"))
  expect_identical(res$param[other_names], esvd_obj$param[other_names])
})

# ---- an object built with bool_diet = TRUE ----------------------------------

## [oracle] `eSVD()` from scratch at the new cap. Two runs of `eSVD()` on the
## same input give the same fit (the SVD start is deterministic), so the redo
## of a run at cap 10 must land on the run at cap 1.
test_that("T-REDO-05: a diet object is redone from the Seurat object and equals eSVD() from scratch", {
  skip_on_cran()
  seurat_obj <- .tiny_seurat()
  esvd_obj <- .redo_esvd_obj(bool_diet = TRUE, cap_multiplier = 10)
  reference <- .redo_esvd_obj(bool_diet = TRUE, cap_multiplier = 1)
  expect_null(esvd_obj$dat)
  expect_null(esvd_obj$covariates)

  res <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = esvd_obj,
                     cap_multiplier = 1,
                     seurat_obj = seurat_obj)
  )

  expect_false(isTRUE(all.equal(res$teststat_vec, esvd_obj$teststat_vec)))
  .expect_same_results(res = res, reference = reference, label = "diet")

  # A diet object comes back diet.
  expect_identical(names(res), names(reference))
  expect_identical(names(res[[res$latest_Fit]]),
                   names(reference[[reference$latest_Fit]]))
  expect_identical(res$param, reference$param)
})

test_that("T-REDO-05b: an object built with bool_diet = FALSE equals eSVD() from scratch, posterior matrices included", {
  skip_on_cran()
  esvd_obj <- .redo_esvd_obj(bool_diet = FALSE, cap_multiplier = 10)
  reference <- .redo_esvd_obj(bool_diet = FALSE, cap_multiplier = 1)
  expect_false(is.null(esvd_obj$dat))

  # No Seurat object: the counts are on the object.
  res <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1)
  )

  expect_false(is.null(res[[res$latest_Fit]]$posterior_mean_mat))
  .expect_same_results(res = res, reference = reference, label = "not diet")
  expect_identical(res$param, reference$param)
})

## [oracle] the covariates and counts of the same run with `bool_diet = FALSE`,
## which carries them. `Log_UMI` and the rescaled `Age` depend on which cells
## are present, so the comparison is made on the cohort the helper filtered:
## the Seurat object passed in still has the dropped individual's cells and
## the all-zero genes.
test_that("T-REDO-06: the counts and covariates rebuilt from the Seurat object are the fitted ones", {
  skip_on_cran()
  seurat_obj <- .redo_helper_seurat()
  reference <- .redo_helper_obj(bool_diet = FALSE, cap_multiplier = 10)
  esvd_obj <- .redo_helper_obj(bool_diet = TRUE, cap_multiplier = 10)
  analyzed_vec <- names(reference$gene_status)[reference$gene_status == "analyzed"]

  # The cohort was filtered: fewer cells and fewer genes than the Seurat object.
  expect_true(nrow(reference$covariates) < ncol(seurat_obj))
  expect_true(length(analyzed_vec) < nrow(seurat_obj))

  res <- .get_counts_and_covariates(input_obj = esvd_obj,
                                    seurat_obj = seurat_obj)

  expect_identical(dimnames(res$covariates), dimnames(reference$covariates))
  expect_equal(res$covariates, reference$covariates, tolerance = 1e-12)
  expect_identical(rownames(res$dat), rownames(reference$covariates))
  expect_identical(colnames(res$dat), analyzed_vec)
  expect_equal(as.matrix(res$dat),
               as.matrix(reference$dat[, analyzed_vec]),
               tolerance = 1e-12)
})

test_that("T-REDO-07: an eSVD_helper() object keeps its gene status and its padding through the redo", {
  skip_on_cran()
  seurat_obj <- .redo_helper_seurat()

  for(bool_diet in c(TRUE, FALSE)){
    label <- paste0("bool_diet = ", bool_diet)
    esvd_obj <- .redo_helper_obj(bool_diet = bool_diet, cap_multiplier = 10)
    reference <- .redo_helper_obj(bool_diet = bool_diet, cap_multiplier = 1)
    zero_idx <- which(esvd_obj$gene_status == "all_zero")
    expect_equal(unname(zero_idx), c(3, 12), info = label)

    res <- .muffle_locfdr_fallback(
      recompute_pvalue(input_obj = esvd_obj,
                       cap_multiplier = 1,
                       seurat_obj = seurat_obj)
    )

    expect_identical(res$gene_status, esvd_obj$gene_status, info = label)
    expect_true(all(is.na(res$teststat_vec[zero_idx])), info = label)
    expect_true(all(res$pvalue_list$fdr_vec[zero_idx] == 1), info = label)
    expect_identical(names(res$teststat_vec), names(esvd_obj$gene_status),
                     info = label)

    expect_false(isTRUE(all.equal(res$teststat_vec, esvd_obj$teststat_vec)),
                 info = label)
    .expect_same_results(res = res, reference = reference, label = label)
    expect_identical(names(res), names(reference), info = label)
    expect_identical(res$dat, reference$dat, info = label)
  }
})

## The redo repeats the test with the settings the object was built with, not
## with the defaults of the functions it calls. Each setting is tried on its
## own, away from its default, and is first shown to change the result: a
## setting that changes nothing on this fixture would let a redo that ignored
## it pass.
##
## [oracle] `eSVD()` from scratch with the same setting, at the new cap.
test_that("T-REDO-10: the redo uses the settings the object was built with", {
  skip_on_cran()
  seurat_obj <- .tiny_seurat()
  reference_default <- .redo_esvd_obj(bool_diet = TRUE, cap_multiplier = 1)

  setting_list <- list(
    alpha_max = list(alpha_max = 3),
    bool_stabilize_underdispersion = list(bool_stabilize_underdispersion = FALSE),
    library_min = list(library_min = 8),
    pseudocount = list(pseudocount = 1),
    bool_covariates_as_library = list(bool_covariates_as_library = FALSE),
    bool_adjust_covariates = list(bool_adjust_covariates = TRUE,
                                  bool_covariates_as_library = FALSE),
    bool_library_includes_interept = list(bool_library_includes_interept = FALSE)
  )
  # `bool_use_log` is not among them: it decides how the rate before the cap
  # is found, and the redo does not estimate it again.

  for(bool_diet in c(TRUE, FALSE)){
    for(i in seq_along(setting_list)){
      label <- paste0(names(setting_list)[i], ", bool_diet = ", bool_diet)
      arg_list <- c(list(dat = .tiny_counts(), bool_diet = bool_diet),
                    setting_list[[i]])

      esvd_obj <- suppressWarnings(
        do.call(.esvd_run, c(arg_list, list(cap_multiplier = 10)))
      )
      reference <- suppressWarnings(
        do.call(.esvd_run, c(arg_list, list(cap_multiplier = 1)))
      )
      expect_false(isTRUE(all.equal(reference$teststat_vec,
                                    reference_default$teststat_vec,
                                    tolerance = 1e-6)),
                   info = label)

      res <- .muffle_locfdr_fallback(
        recompute_pvalue(input_obj = esvd_obj,
                         cap_multiplier = 1,
                         seurat_obj = seurat_obj)
      )

      .expect_same_results(res = res, reference = reference, label = label)
      expect_identical(res$param, reference$param, info = label)
    }
  }
})

## `min_cells_per_individual = 0` is how a cohort with a two-cell individual
## is analyzed at all. A redo that fell back to the default of 3 would refuse
## the object it was given.
test_that("T-REDO-11: the redo keeps min_cells_per_individual", {
  skip_on_cran()
  seurat_obj <- .redo_helper_seurat()
  run_esvd <- function(cap_multiplier, bool_diet){
    suppressWarnings(
      eSVD_helper(batch_var_prefix = NULL,
                  case_control_levels = c("0", "1"),
                  case_control_var = "CC",
                  categorical_vars = c("Sex"),
                  id_var = "Individual",
                  numerical_vars = "Age",
                  seurat_obj = seurat_obj,
                  bool_diet = bool_diet,
                  cap_multiplier = cap_multiplier,
                  k = 2,
                  min_cells_per_id = 0)
    )
  }

  for(bool_diet in c(TRUE, FALSE)){
    label <- paste0("bool_diet = ", bool_diet)
    esvd_obj <- run_esvd(cap_multiplier = 10, bool_diet = bool_diet)
    reference <- run_esvd(cap_multiplier = 1, bool_diet = bool_diet)
    # The two-cell individual was kept.
    expect_true("indiv_5" %in% as.character(esvd_obj$individual), info = label)
    expect_equal(esvd_obj$param$test_min_cells_per_individual, 0, info = label)

    res <- .muffle_locfdr_fallback(
      recompute_pvalue(input_obj = esvd_obj,
                       cap_multiplier = 1,
                       seurat_obj = seurat_obj)
    )

    .expect_same_results(res = res, reference = reference, label = label)
  }
})

# ---- refusals ---------------------------------------------------------------

test_that("T-REDO-08: a diet object cannot be redone without the Seurat object it was built from", {
  skip_on_cran()
  seurat_obj <- .tiny_seurat()
  esvd_obj <- .redo_esvd_obj(bool_diet = TRUE, cap_multiplier = 10)

  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1),
               regexp = "seurat_obj")
  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1,
                                seurat_obj = "not a Seurat object"),
               regexp = "seurat_obj")

  # Cells of the fit that the Seurat object does not have.
  seurat_fewer_cells <- seurat_obj[, SeuratObject::Cells(seurat_obj)[-(1:5)]]
  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1,
                                seurat_obj = seurat_fewer_cells),
               regexp = "5 cell")

  # Genes of the fit that the Seurat object does not have.
  gene_vec <- rownames(seurat_obj)
  seurat_fewer_genes <- subset(seurat_obj, features = gene_vec[-(1:2)])
  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1,
                                seurat_obj = seurat_fewer_genes),
               regexp = "2 gene")

  # The right names and the wrong counts.
  dat_other <- .tiny_counts()
  dat_other[, 4] <- dat_other[, 4] + 1
  seurat_other_counts <- .tiny_seurat(dat = dat_other)
  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1,
                                seurat_obj = seurat_other_counts),
               regexp = "counts")

  # The right counts and another value of a covariate.
  seurat_other_age <- seurat_obj
  seurat_other_age@meta.data[1:20, "Age"] <- 90
  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1,
                                seurat_obj = seurat_other_age),
               regexp = "covariates")
})

## Finding of the code review (2026-09-29). The rebuilt covariates were
## compared with the fitted ones by their column sums, and a sum does not
## change when values move between cells. Exchanging a covariate between two
## individuals with the same number of cells was accepted, and the redo
## returned statistics 0.06 to 0.18 away from the right ones without a word.
##
## Individuals 1 and 2 are both controls with 20 cells each, of different
## `Sex` and `Age` (asserted, or the exchange would change nothing).
test_that("T-REDO-08b: a covariate exchanged between two individuals is refused", {
  skip_on_cran()
  seurat_obj <- .tiny_seurat()
  esvd_obj <- .redo_esvd_obj(bool_diet = TRUE, cap_multiplier = 10)

  metadata_df <- seurat_obj@meta.data
  idx1_vec <- which(metadata_df[, "Individual"] == "indiv_1")
  idx2_vec <- which(metadata_df[, "Individual"] == "indiv_2")
  expect_equal(length(idx1_vec), length(idx2_vec))

  for(variable in c("Age", "Sex")){
    value1 <- metadata_df[idx1_vec[1], variable]
    value2 <- metadata_df[idx2_vec[1], variable]
    expect_false(identical(as.character(value1), as.character(value2)),
                 info = variable)

    seurat_exchanged <- seurat_obj
    seurat_exchanged@meta.data[idx1_vec, variable] <- value2
    seurat_exchanged@meta.data[idx2_vec, variable] <- value1

    expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1,
                                  seurat_obj = seurat_exchanged),
                 regexp = "covariates", info = variable)
  }
})

## The counts of two cells exchanged: every gene's mean count is unchanged.
test_that("T-REDO-08c: counts exchanged between two cells are refused", {
  skip_on_cran()
  esvd_obj <- .redo_esvd_obj(bool_diet = TRUE, cap_multiplier = 10)

  dat <- .tiny_counts()
  dat_exchanged <- dat
  dat_exchanged[c(1, 100), ] <- dat[c(100, 1), ]
  expect_false(sum(dat[1, ]) == sum(dat[100, ]))
  expect_equal(colMeans(dat_exchanged), colMeans(dat))

  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = 1,
                                seurat_obj = .tiny_seurat(dat = dat_exchanged)),
               regexp = "counts|covariates")
})

test_that("T-REDO-09: an object that cannot be redone is refused with the reason", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  # Built before the rates were stored beside the cap.
  old_obj <- esvd_obj
  old_obj[[latest_fit]]$nuisance_mle_vec <- NULL
  expect_error(recompute_pvalue(input_obj = old_obj, cap_multiplier = 1),
               regexp = "nuisance_mle_vec")

  # Not tested yet, so there is nothing to redo.
  untested_obj <- esvd_obj
  untested_obj$teststat_vec <- NULL
  untested_obj$pvalue_list <- NULL
  expect_error(recompute_pvalue(input_obj = untested_obj, cap_multiplier = 1),
               regexp = "teststat_vec")

  # `initialize_esvd` accepts a count matrix without names, and everything
  # here is matched by name. `dat[, NULL]` selects no column and raises
  # nothing, so the missing names have to be refused up front.
  unnamed_obj <- esvd_obj
  rownames(unnamed_obj[[latest_fit]]$y_mat) <- NULL
  expect_error(recompute_pvalue(input_obj = unnamed_obj, cap_multiplier = 1),
               regexp = "gene names")
  unnamed_obj <- esvd_obj
  rownames(unnamed_obj[[latest_fit]]$x_mat) <- NULL
  expect_error(recompute_pvalue(input_obj = unnamed_obj, cap_multiplier = 1),
               regexp = "cell names")

  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = -1),
               regexp = "cap_multiplier")
  expect_error(recompute_pvalue(input_obj = list(), cap_multiplier = 1))
})

## The floor `min_val` recorded at the fit shares the units of the cap, so a
## redo at a cap at or below it would put every gene above its cap. It is
## refused before anything is recomputed.
test_that("T-REDO-13: a cap at or below the recorded min_val is refused", {
  skip_on_cran()
  esvd_obj <- .redo_esvd_obj(bool_diet = FALSE, cap_multiplier = 10)
  min_val <- esvd_obj$param$nuisance_min_val
  expect_equal(min_val, 1e-4)

  expect_error(recompute_pvalue(input_obj = esvd_obj, cap_multiplier = min_val),
               regexp = "`min_val`")
  expect_error(recompute_pvalue(input_obj = esvd_obj,
                                cap_multiplier = min_val / 2),
               regexp = "`min_val`")
  expect_error(recompute_pvalue(input_obj = esvd_obj,
                                cap_multiplier = 2 * min_val),
               regexp = NA)
})
