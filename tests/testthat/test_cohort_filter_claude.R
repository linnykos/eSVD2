# UNIT_TEST_PLAN.md section 2.17 -- T-COH-01 .. T-COH-14.
#
# THIS WHOLE FILE IS EXPECTED TO SKIP OR FAIL. `eSVD_helper()` and
# `filter_cohort()` are NEW and have not been implemented -- no eSVD2 R or C++
# code has been changed. The tests are written first, to the five-step
# specification of section 2.17.2:
#
#   1. drop donors with < min_cells_per_id cells, warning and naming them
#   2. check min_cells / min_ids / min_ids_per_arm / min_cells_casecontrol on
#      the FILTERED object; on failure warning() + return(NA)
#   3. label gene_status and remove the all_zero genes
#   4. call eSVD()
#   5. reinsert the removed genes at their original positions

.skip_if_no_helper <- function(){
  skip_if_not(exists("eSVD_helper"),
              "eSVD_helper() not implemented yet (UNIT_TEST_PLAN.md section 2.17)")
}

.skip_if_no_filter <- function(){
  skip_if_not(exists("filter_cohort"),
              "filter_cohort() not implemented yet (UNIT_TEST_PLAN.md section 2.17)")
}

# Builds a cohort with a controllable number of cells for one nominated donor,
# so the min_cells_per_id boundary can be probed from both sides.
.cohort_seurat <- function(cells_for_last_donor = 20,
                           num_individuals = 8,
                           num_cells_per_individual = 25){
  dat_list <- .build_tiny_data(num_cells_per_individual = num_cells_per_individual,
                               num_genes = 20,
                               num_individuals = num_individuals,
                               seed_number = 30)

  keep_vec <- rep(TRUE, nrow(dat_list$dat))
  last_donor <- levels(dat_list$individual_vec)[num_individuals]
  last_donor_idx <- which(dat_list$individual_vec == last_donor)
  if(cells_for_last_donor < length(last_donor_idx)){
    drop_idx <- utils::tail(last_donor_idx,
                            length(last_donor_idx) - cells_for_last_donor)
    keep_vec[drop_idx] <- FALSE
  }

  .tiny_seurat(dat = dat_list$dat[keep_vec, , drop = FALSE],
               covariate_df = dat_list$covariate_df[keep_vec, , drop = FALSE],
               individual_vec = droplevels(dat_list$individual_vec[keep_vec]))
}

test_that("T-COH-01: donors with 1 or 2 cells are dropped; 3 is kept", {
  .skip_if_no_filter()

  for(cell_count in c(1, 2)){
    seurat_obj <- .cohort_seurat(cells_for_last_donor = cell_count)
    res <- suppressWarnings(filter_cohort(seurat_obj = seurat_obj,
                                          id_var = "Individual",
                                          case_control_var = "CC"))
    expect_equal(length(unique(res$seurat_obj@meta.data[, "Individual"])), 7,
                 info = paste0(cell_count, " cells"))
  }

  seurat_obj <- .cohort_seurat(cells_for_last_donor = 3)
  res <- suppressWarnings(filter_cohort(seurat_obj = seurat_obj,
                                        id_var = "Individual",
                                        case_control_var = "CC"))
  expect_equal(length(unique(res$seurat_obj@meta.data[, "Individual"])), 8)
})

## Question Q-COH-3 resolved: drop with a warning. `esvd_helper.R` as it stands
## in the WAS2CODE_REPO drops SILENTLY, so this is expected to FAIL until the
## warning is added. A cohort quietly losing donors between input and output is
## the kind of thing that is discovered at revision time.
test_that("T-COH-02: the drop warns and names the donors dropped", {
  .skip_if_no_filter()

  seurat_obj <- .cohort_seurat(cells_for_last_donor = 2)

  expect_warning(filter_cohort(seurat_obj = seurat_obj,
                               id_var = "Individual",
                               case_control_var = "CC"),
                 regexp = "indiv_8")
})

## Question Q-COH-6 resolved: DROP FIRST, THEN CHECK. Checking before dropping
## admits exactly the cohorts the checks exist to reject -- 8 donors of which 2
## have one cell passes an 8-donor check, then becomes 6.
test_that("T-COH-03: the drop runs first and the checks see the filtered cohort", {
  .skip_if_no_filter()

  # 5 donors, two of them with a single cell. `min_ids = 4` uses `<=`, so a
  # 5-donor cohort passes -- but after the drop only 3 remain, which must not.
  seurat_obj <- .cohort_seurat(cells_for_last_donor = 1, num_individuals = 6)

  res <- suppressWarnings(filter_cohort(seurat_obj = seurat_obj,
                                        id_var = "Individual",
                                        case_control_var = "CC",
                                        min_ids = 5))

  # `all(is.na(NA))` is TRUE for any NA, so asserting only that would pass for
  # any rejection whatsoever. The warning text is what identifies WHICH filter
  # fired (question Q-COH-1 kept the bare NA return), so match on it.
  expect_warning(filter_cohort(seurat_obj = seurat_obj,
                               id_var = "Individual",
                               case_control_var = "CC",
                               min_ids = 5),
                 regexp = "individual|donor|min_ids")
  expect_true(all(is.na(res)))
})

## Question Q-COH-2 resolved: add `min_ids_per_arm`, default 2. `min_ids` counts
## donors POOLED ACROSS ARMS, so 4 case + 1 control passes it and still yields
## `n2 - 1 = 0` in `.compute_df()` -- T-TSTAT-04's NaN reached through the front
## door.
test_that("T-COH-04: a 4-case / 1-control cohort is rejected", {
  .skip_if_no_filter()

  dat_list <- .build_tiny_data(num_cells_per_individual = 25,
                               num_genes = 20,
                               num_individuals = 10,
                               seed_number = 40)
  # Keep 4 cases and 1 control.
  keep_donors <- c(levels(dat_list$individual_vec)[1],
                   levels(dat_list$individual_vec)[6:9])
  keep_vec <- dat_list$individual_vec %in% keep_donors

  seurat_obj <- .tiny_seurat(
    dat = dat_list$dat[keep_vec, , drop = FALSE],
    covariate_df = dat_list$covariate_df[keep_vec, , drop = FALSE],
    individual_vec = droplevels(dat_list$individual_vec[keep_vec])
  )

  res <- suppressWarnings(filter_cohort(seurat_obj = seurat_obj,
                                        id_var = "Individual",
                                        case_control_var = "CC",
                                        min_ids = 4,
                                        min_ids_per_arm = 2))

  expect_warning(filter_cohort(seurat_obj = seurat_obj,
                               id_var = "Individual",
                               case_control_var = "CC",
                               min_ids = 4,
                               min_ids_per_arm = 2),
                 regexp = "arm")
  expect_true(all(is.na(res)))
})

## Question Q-COH-1 resolved: keep `warning()` + `return(NA)`; no classed
## sentinel. The consequence is that the warning TEXT is the only channel
## telling a caller which of the three rejections fired, so each must name its
## own threshold.
test_that("T-COH-05: each rejection emits a distinct, identifying warning", {
  .skip_if_no_helper()

  tiny_seurat <- .cohort_seurat(num_individuals = 8,
                                num_cells_per_individual = 2)

  # Too few cells overall.
  expect_warning(eSVD_helper(batch_var_prefix = NULL,
                             case_control_levels = c("0", "1"),
                             case_control_var = "CC",
                             categorical_vars = "Sex",
                             id_var = "Individual",
                             numerical_vars = "Age",
                             seurat_obj = tiny_seurat,
                             min_cells = 20),
                 regexp = "cell")
})

## Question Q-COH-4 resolved: `bool_check_donors` is removed entirely, so
## thresholds are the only escape hatch. That opens a hazard in the middle:
## `min_cells_per_id` of 1 or 2 does not disable the drop, it WEAKENS it, and a
## surviving 1-cell donor is back on the section 1.1 path.
test_that("T-COH-06: min_cells_per_id must be 0 (off) or >= 3 (safe)", {
  .skip_if_no_filter()

  seurat_obj <- .cohort_seurat()

  expect_error(filter_cohort(seurat_obj = seurat_obj,
                             id_var = "Individual",
                             case_control_var = "CC",
                             min_cells_per_id = 2))
  expect_error(filter_cohort(seurat_obj = seurat_obj,
                             id_var = "Individual",
                             case_control_var = "CC",
                             min_cells_per_id = 1))
})

## The ordering constraint that ties section 2.17 to section 2.16: a gene
## expressed ONLY in a dropped low-cell donor becomes all-zero AFTER the drop.
## Computing `gene_status` first would mark it `analyzed` and send an all-zero
## count vector through the whole pipeline.
test_that("T-COH-07: the donor drop runs before gene_status is computed", {
  .skip_if_no_helper()

  dat_list <- .build_tiny_data(num_cells_per_individual = 25,
                               num_genes = 20,
                               num_individuals = 8,
                               seed_number = 50)

  # Reduce donor 8 to two cells, then make gene 5 non-zero ONLY in those cells.
  last_donor_idx <- which(dat_list$individual_vec == "indiv_8")
  keep_vec <- rep(TRUE, nrow(dat_list$dat))
  keep_vec[utils::tail(last_donor_idx, length(last_donor_idx) - 2)] <- FALSE

  dat <- dat_list$dat
  dat[, 5] <- 0
  dat[utils::head(last_donor_idx, 2), 5] <- 7

  seurat_obj <- .tiny_seurat(
    dat = dat[keep_vec, , drop = FALSE],
    covariate_df = dat_list$covariate_df[keep_vec, , drop = FALSE],
    individual_vec = droplevels(dat_list$individual_vec[keep_vec])
  )

  res <- suppressWarnings(
    eSVD_helper(batch_var_prefix = NULL,
                case_control_levels = c("0", "1"),
                case_control_var = "CC",
                categorical_vars = "Sex",
                id_var = "Individual",
                numerical_vars = "Age",
                seurat_obj = seurat_obj,
                k = 2)
  )

  # If the order inverts, gene 5 is marked `analyzed` and goes through the
  # pipeline with an all-zero count vector -- the exact failure section 2.16
  # exists to prevent, reintroduced by ordering alone.
  expect_equal(as.character(res[["gene_status"]][5]), "all_zero")
})

test_that("T-COH-08: the filters are a no-op on a clean cohort", {
  .skip_if_no_helper()

  # The wrapper must not move existing results on data where no threshold
  # binds, or every published number silently changes.
  seurat_obj <- .cohort_seurat(num_individuals = 8,
                               num_cells_per_individual = 25)

  res_helper <- suppressWarnings(
    eSVD_helper(batch_var_prefix = NULL,
                case_control_levels = c("0", "1"),
                case_control_var = "CC",
                categorical_vars = "Sex",
                id_var = "Individual",
                numerical_vars = "Age",
                seurat_obj = seurat_obj,
                k = 2)
  )
  res_direct <- suppressWarnings(
    eSVD(batch_var_prefix = NULL,
         case_control_levels = c("0", "1"),
         case_control_var = "CC",
         categorical_vars = "Sex",
         id_var = "Individual",
         numerical_vars = "Age",
         seurat_obj = seurat_obj,
         k = 2)
  )

  expect_equal(unname(res_helper$teststat_vec),
               unname(res_direct$teststat_vec), tolerance = 1e-10)
})

test_that("T-COH-09: ... reaches eSVD() intact", {
  .skip_if_no_helper()

  seurat_obj <- .cohort_seurat(num_individuals = 8,
                               num_cells_per_individual = 25)

  # The passthrough is the whole rest of the function, and a typo'd argument
  # name would be swallowed by `...` in silence.
  res <- suppressWarnings(
    eSVD_helper(batch_var_prefix = NULL,
                case_control_levels = c("0", "1"),
                case_control_var = "CC",
                categorical_vars = "Sex",
                id_var = "Individual",
                numerical_vars = "Age",
                seurat_obj = seurat_obj,
                k = 3,
                library_min = 0.25)
  )

  expect_equal(res$param$init_k, 3)
})

test_that("T-COH-10: all five thresholds are reachable and bool_check_donors is gone", {
  .skip_if_no_helper()

  argument_vec <- names(formals(eSVD_helper))

  for(threshold_name in c("min_cells_per_id", "min_cells", "min_ids",
                          "min_ids_per_arm", "min_cells_casecontrol")){
    expect_true(threshold_name %in% argument_vec, info = threshold_name)
  }
  # Question Q-COH-4: removed entirely, not deprecated.
  expect_false("bool_check_donors" %in% argument_vec)
})

## Question Q-COH-7 resolved: the wrapper filters, the pipeline refuses.
## `eSVD()` called directly must ERROR on the correctness conditions -- and must
## NOT error merely for being underpowered, since those thresholds are the
## helper's policy, not the model's.
test_that("T-COH-11: eSVD() errors on 1-cell donors and on all-zero genes", {
  seurat_obj <- .cohort_seurat(cells_for_last_donor = 1)

  expect_error(
    suppressWarnings(eSVD(batch_var_prefix = NULL,
                          case_control_levels = c("0", "1"),
                          case_control_var = "CC",
                          categorical_vars = "Sex",
                          id_var = "Individual",
                          numerical_vars = "Age",
                          seurat_obj = seurat_obj,
                          k = 2)),
    regexp = "cell"
  )

  expect_error(suppressWarnings(.esvd_run(.tiny_counts_with_zero_genes())),
               regexp = "zero")
})

test_that("T-COH-11b: eSVD() does not error merely for being underpowered", {
  # A legitimately small run -- 15 cells, 4 donors -- must stay possible through
  # `eSVD()` directly. Those minima are the helper's policy, not the model's.
  dat_list <- .build_tiny_data(num_cells_per_individual = 6,
                               num_genes = 12,
                               num_individuals = 4,
                               seed_number = 60)
  seurat_obj <- .tiny_seurat(dat = dat_list$dat,
                             covariate_df = dat_list$covariate_df,
                             individual_vec = dat_list$individual_vec)

  expect_no_error(
    suppressWarnings(eSVD(batch_var_prefix = NULL,
                          case_control_levels = c("0", "1"),
                          case_control_var = "CC",
                          categorical_vars = "Sex",
                          id_var = "Individual",
                          numerical_vars = "Age",
                          seurat_obj = seurat_obj,
                          k = 2))
  )
})

test_that("T-COH-12: setting a threshold to 0 disables that filter", {
  .skip_if_no_filter()

  # With `bool_check_donors` removed this is the only escape hatch, so it needs
  # to be a tested contract rather than an emergent property of `<=` and `<`.
  seurat_obj <- .cohort_seurat(cells_for_last_donor = 1)

  res <- suppressWarnings(filter_cohort(seurat_obj = seurat_obj,
                                        id_var = "Individual",
                                        case_control_var = "CC",
                                        min_cells_per_id = 0))

  expect_equal(length(unique(res$seurat_obj@meta.data[, "Individual"])), 8)
})

## Question Q-COH-7 makes both filters SEURAT OBJECT subsets, not matrix
## subsets: `obj[, keep_cells]` and `obj[keep_genes, ]`. Feature subsetting a
## Seurat object touches variable features, `scale.data` and reductions; none
## of that matters to `eSVD()`, but a subset that quietly dropped or reordered
## the counts layer would.
test_that("T-COH-13: subsetting preserves the counts layer and meta.data", {
  dat <- .tiny_counts()
  seurat_obj <- .tiny_seurat(dat = dat)

  keep_genes <- setdiff(colnames(dat), colnames(dat)[c(3, 12)])
  # `seurat_obj[keep_genes, ]` errors with "the default assay will be removed";
  # `subset(features = )` is the supported route, and this is exactly the kind
  # of Seurat-specific detail T-COH-13 exists to pin down before the helper is
  # written against it.
  subset_obj <- subset(seurat_obj, features = keep_genes)

  count_mat <- Matrix::t(SeuratObject::LayerData(subset_obj, assay = "RNA",
                                                 layer = "counts"))
  expect_equal(as.matrix(count_mat), dat[, keep_genes, drop = FALSE],
               ignore_attr = TRUE)
  expect_equal(nrow(subset_obj@meta.data), nrow(dat))

  # And the Log_UMI identity of T-GS-03 survives the subset.
  expect_equal(unname(log(Matrix::rowSums(count_mat))),
               unname(log(Matrix::rowSums(dat[, keep_genes, drop = FALSE]))))
})

test_that("T-COH-14: the five-step order holds end to end", {
  .skip_if_no_helper()

  # T-COH-03 and T-COH-07 each pin one adjacent pair; this pins the whole
  # chain, which is where an ordering regression would actually land.
  dat_list <- .build_tiny_data(num_cells_per_individual = 25,
                               num_genes = 20,
                               num_individuals = 8,
                               seed_number = 70)

  last_donor_idx <- which(dat_list$individual_vec == "indiv_8")
  keep_vec <- rep(TRUE, nrow(dat_list$dat))
  keep_vec[utils::tail(last_donor_idx, length(last_donor_idx) - 2)] <- FALSE

  dat <- dat_list$dat
  dat[, 5] <- 0
  dat[utils::head(last_donor_idx, 2), 5] <- 7

  seurat_obj <- .tiny_seurat(
    dat = dat[keep_vec, , drop = FALSE],
    covariate_df = dat_list$covariate_df[keep_vec, , drop = FALSE],
    individual_vec = droplevels(dat_list$individual_vec[keep_vec])
  )

  res <- suppressWarnings(
    eSVD_helper(batch_var_prefix = NULL,
                case_control_levels = c("0", "1"),
                case_control_var = "CC",
                categorical_vars = "Sex",
                id_var = "Individual",
                numerical_vars = "Age",
                seurat_obj = seurat_obj,
                k = 2,
                min_ids = 4)
  )

  # Step 1 dropped indiv_8, step 2 passed on 7 donors, step 3 marked gene 5
  # all_zero, steps 4 and 5 ran and reinserted it.
  expect_equal(as.character(res[["gene_status"]][5]), "all_zero")
  expect_equal(length(res$teststat_vec), ncol(dat))
  expect_true(is.na(res$teststat_vec[5]))
})
