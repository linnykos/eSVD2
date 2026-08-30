# UNIT_TEST_PLAN.md sections 2.2, 2.11, 3.1, 5 and 6.
#
# The principle of section 5: every exported function gets one test that
# malformed input produces a USEFUL error rather than a downstream NA.
# Question Q-VAL-1 / C-10 resolved: `expect_error(regexp = ...)` with a short
# stable fragment, never a snapshot. Many of these are expected to FAIL,
# because `stopifnot()` is used pervasively and its messages name nothing.

# ---- section 2.2, initialize_esvd -------------------------------------------

test_that("T-INIT-01: bool_intercept = TRUE gives a non-zero intercept column", {
  dat_list <- .tiny_data()
  covariates <- .tiny_covariates()

  res_intercept <- suppressWarnings(
    initialize_esvd(dat = dat_list$dat, covariates = covariates,
                    metadata_individual = dat_list$individual_vec,
                    bool_intercept = TRUE, case_control_variable = "CC_1",
                    k = 2, metadata_case_control = covariates[, "CC_1"],
                    verbose = 0)
  )
  res_no_intercept <- suppressWarnings(
    initialize_esvd(dat = dat_list$dat, covariates = covariates,
                    metadata_individual = dat_list$individual_vec,
                    bool_intercept = FALSE, case_control_variable = "CC_1",
                    k = 2, metadata_case_control = covariates[, "CC_1"],
                    verbose = 0)
  )

  # The `bool_intercept = TRUE` branch is completely untested today, and a
  # `c(0, ...)` versus `c(a0, ...)` slip in it would be invisible.
  z_intercept <- res_intercept[["fit_Init"]]$z_mat[, "Intercept"]
  z_no_intercept <- res_no_intercept[["fit_Init"]]$z_mat[, "Intercept"]

  expect_true(any(abs(z_intercept) > 1e-8))
  expect_true(all(abs(z_no_intercept) < 1e-8))
})

test_that("T-INIT-02/03: offset_variables pins the named coefficients to 1", {
  dat_list <- .tiny_data()
  covariates <- .tiny_covariates()

  res_single <- suppressWarnings(
    initialize_esvd(dat = dat_list$dat, covariates = covariates,
                    metadata_individual = dat_list$individual_vec,
                    bool_intercept = TRUE, case_control_variable = "CC_1",
                    k = 2, metadata_case_control = covariates[, "CC_1"],
                    offset_variables = "Log_UMI", verbose = 0)
  )
  expect_equal(unname(res_single[["fit_Init"]]$z_mat[, "Log_UMI"]),
               rep(1, ncol(dat_list$dat)))

  # The multi-offset path is reachable and untested. Note it forces BOTH
  # coefficients to 1, which is a modelling assumption worth stating out loud.
  res_multi <- suppressWarnings(
    initialize_esvd(dat = dat_list$dat, covariates = covariates,
                    metadata_individual = dat_list$individual_vec,
                    bool_intercept = TRUE, case_control_variable = "CC_1",
                    k = 2, metadata_case_control = covariates[, "CC_1"],
                    offset_variables = c("Log_UMI", "Age"), verbose = 0)
  )
  expect_equal(unname(res_multi[["fit_Init"]]$z_mat[, "Log_UMI"]),
               rep(1, ncol(dat_list$dat)))
  expect_equal(unname(res_multi[["fit_Init"]]$z_mat[, "Age"]),
               rep(1, ncol(dat_list$dat)))
})

## Question Q-INIT-1 resolved: the sparse path should zero NAs the way the
## dense path does. Expected to FAIL -- `dat[is.na(dat)] <- 0` sits under an
## `is.matrix(dat)` guard, so a `dgCMatrix` keeps its NAs.
test_that("T-INIT-05: dense and sparse both zero NA entries", {
  dat_list <- .tiny_data()
  covariates <- .tiny_covariates()

  dat_dense <- dat_list$dat
  dat_dense[1, 1] <- NA
  dat_sparse <- methods::as(methods::as(dat_dense * 1.0, "dMatrix"),
                            "CsparseMatrix")

  res_dense <- suppressWarnings(
    initialize_esvd(dat = dat_dense, covariates = covariates,
                    metadata_individual = dat_list$individual_vec,
                    bool_intercept = TRUE, case_control_variable = "CC_1",
                    k = 2, metadata_case_control = covariates[, "CC_1"],
                    verbose = 0)
  )
  # The sparse run currently ERRORS inside glmnet, because the NA survives:
  # "missing value where TRUE/FALSE needed" from the Poisson response check.
  # Caught here so the failure is reported as this assertion rather than
  # aborting the whole file.
  res_sparse <- try(suppressWarnings(
    initialize_esvd(dat = dat_sparse, covariates = covariates,
                    metadata_individual = dat_list$individual_vec,
                    bool_intercept = TRUE, case_control_variable = "CC_1",
                    k = 2, metadata_case_control = covariates[, "CC_1"],
                    verbose = 0)
  ), silent = TRUE)

  expect_false(inherits(res_sparse, "try-error"),
               info = if(inherits(res_sparse, "try-error"))
                 conditionMessage(attr(res_sparse, "condition")) else "ok")
  skip_if(inherits(res_sparse, "try-error"),
          "sparse path does not zero NAs (Q-INIT-1 fix not applied)")

  expect_equal(res_dense[["fit_Init"]]$z_mat, res_sparse[["fit_Init"]]$z_mat,
               tolerance = 1e-6)
})

# ---- section 2.11, compute_test_per_gene ------------------------------------

test_that("T-PG-01/02: the two pipelines agree element-wise on every output", {
  esvd_obj <- .small_esvd_obj()

  # Strip the matrix pipeline's results so the per-gene path starts clean.
  per_gene_input <- esvd_obj
  latest_fit <- per_gene_input[["latest_Fit"]]
  per_gene_input[[latest_fit]]$posterior_mean_mat <- NULL
  per_gene_input[[latest_fit]]$posterior_var_mat <- NULL
  per_gene_input$teststat_vec <- NULL
  per_gene_input$case_mean <- NULL
  per_gene_input$control_mean <- NULL
  per_gene_input$pvalue_list <- NULL

  res_per_gene <- suppressWarnings(
    compute_test_per_gene(input_obj = per_gene_input,
                          alpha_max = 2 * max(esvd_obj$dat),
                          bool_covariates_as_library = TRUE,
                          library_min = 0.1,
                          verbose = 0)
  )

  # The existing test is `abs(sum(a - b)) <= 1e-3`, which per-gene errors of
  # opposite sign cancel out of. This is element-wise.
  expect_equal(unname(res_per_gene$teststat_vec),
               unname(esvd_obj$teststat_vec), tolerance = 1e-8)
  expect_equal(unname(res_per_gene$case_mean), unname(esvd_obj$case_mean),
               tolerance = 1e-8)
  expect_equal(unname(res_per_gene$control_mean),
               unname(esvd_obj$control_mean), tolerance = 1e-8)

  for(element_name in c("df_vec", "gaussian_teststat", "log10pvalue",
                        "fdr_vec")){
    expect_equal(unname(res_per_gene$pvalue_list[[element_name]]),
                 unname(esvd_obj$pvalue_list[[element_name]]),
                 tolerance = 1e-8, info = element_name)
  }
  expect_equal(res_per_gene$pvalue_list$null_mean,
               esvd_obj$pvalue_list$null_mean, tolerance = 1e-8)
  expect_equal(res_per_gene$pvalue_list$null_sd,
               esvd_obj$pvalue_list$null_sd, tolerance = 1e-8)
})

## Question Q-PG-1 resolved: `library_min = 0.1` everywhere. Expected to FAIL --
## `compute_posterior.eSVD` defaults to 0.1 and `compute_test_per_gene` to 1e-2,
## so the two paths return different numbers when called with their own
## defaults. `eSVD()` passes the value explicitly, which is why this has gone
## unnoticed.
test_that("T-PG-03: the two default library_min values agree", {
  expect_equal(formals(compute_test_per_gene)$library_min,
               formals(eSVD2:::compute_posterior.eSVD)$library_min)
})

test_that("T-PG-05: compute_test_per_gene does not write the posterior matrices", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  per_gene_input <- esvd_obj
  per_gene_input[[latest_fit]]$posterior_mean_mat <- NULL
  per_gene_input[[latest_fit]]$posterior_var_mat <- NULL

  res <- suppressWarnings(
    compute_test_per_gene(input_obj = per_gene_input,
                          alpha_max = 2 * max(esvd_obj$dat),
                          bool_covariates_as_library = TRUE,
                          library_min = 0.1, verbose = 0)
  )

  # Its documented memory contract, and the whole reason the function exists.
  expect_null(res[[latest_fit]]$posterior_mean_mat)
  expect_null(res[[latest_fit]]$posterior_var_mat)
})

# ---- section 3.1, the C++ data loader ---------------------------------------

test_that("T-CPP-LOAD-01: dense double, dense integer and dgCMatrix agree", {
  point <- .feasible_point("poisson")
  family_obj <- esvd_family("poisson")
  num_factors <- ncol(point$x_mat)

  dat_double <- point$dat
  dat_integer <- matrix(as.integer(point$dat), nrow = nrow(point$dat),
                        ncol = ncol(point$dat))
  dat_sparse <- methods::as(methods::as(dat_double, "dMatrix"), "CsparseMatrix")

  objective_vec <- sapply(list(dat_double, dat_integer, dat_sparse),
                          function(dat){
    objfn_all_r(XC = point$x_mat, YZ = point$y_mat, k = num_factors,
                loader = data_loader(dat), family = family_obj,
                s = point$s_vec, gamma = point$gamma_vec,
                l2penx = 0, l2peny = 0, l2penz = 0)
  })

  # Three loader implementations of one interface; this is the cheapest way to
  # test all three at once.
  expect_equal(objective_vec[1], objective_vec[2], tolerance = 1e-10)
  expect_equal(objective_vec[1], objective_vec[3], tolerance = 1e-10)
})

test_that("T-CPP-LOAD-03: an all-zero sparse column iterates like its dense twin", {
  point <- .feasible_point("poisson")
  family_obj <- esvd_family("poisson")
  num_factors <- ncol(point$x_mat)

  dat_dense <- point$dat
  dat_dense[, 2] <- 0
  dat_sparse <- methods::as(methods::as(dat_dense, "dMatrix"), "CsparseMatrix")

  # An all-zero gene is common in real single-cell data and hits the
  # `m_nnz < 1` early return in `SparseVecIterator::value`, which nothing else
  # reaches.
  expect_equal(
    objfn_all_r(XC = point$x_mat, YZ = point$y_mat, k = num_factors,
                loader = data_loader(dat_dense), family = family_obj,
                s = point$s_vec, gamma = point$gamma_vec,
                l2penx = 0, l2peny = 0, l2penz = 0),
    objfn_all_r(XC = point$x_mat, YZ = point$y_mat, k = num_factors,
                loader = data_loader(dat_sparse), family = family_obj,
                s = point$s_vec, gamma = point$gamma_vec,
                l2penx = 0, l2peny = 0, l2penz = 0),
    tolerance = 1e-10
  )
})

## [verified] `data_loader()` on a `dgeMatrix` and on a character matrix returns
## WITHOUT error, leaving a null external pointer; the failure surfaces later as
## "external pointer is not valid". A `dgeMatrix` is entirely plausible user
## input. Expected to FAIL.
test_that("T-CPP-LOAD-05: data_loader errors at construction on an unsupported type", {
  set.seed(10)
  dense_mat <- matrix(stats::rpois(20, lambda = 3) * 1.0, nrow = 5, ncol = 4)
  dge_mat <- methods::as(dense_mat, "denseMatrix")

  expect_error(data_loader(dge_mat), regexp = "matrix")
  expect_error(data_loader(matrix("a", nrow = 5, ncol = 4)), regexp = "matrix")
})

# ---- section 5, input validation --------------------------------------------

test_that("T-VAL-01/03/04/05: initialize_esvd rejects malformed input", {
  dat_list <- .tiny_data()
  covariates <- .tiny_covariates()

  # T-VAL-01: k greater than the gene count.
  expect_error(initialize_esvd(dat = dat_list$dat, covariates = covariates,
                               metadata_individual = dat_list$individual_vec,
                               k = ncol(dat_list$dat) + 5),
               regexp = "k")

  # T-VAL-03: covariates lacking an Intercept column.
  expect_error(initialize_esvd(dat = dat_list$dat,
                               covariates = covariates[, -1, drop = FALSE],
                               metadata_individual = dat_list$individual_vec,
                               k = 2),
               regexp = "Intercept")

  # T-VAL-04: metadata_individual not a factor.
  expect_error(initialize_esvd(dat = dat_list$dat, covariates = covariates,
                               metadata_individual =
                                 as.character(dat_list$individual_vec),
                               k = 2),
               regexp = "metadata_individual")

  # T-VAL-05: lambda outside its documented bounds.
  expect_error(initialize_esvd(dat = dat_list$dat, covariates = covariates,
                               metadata_individual = dat_list$individual_vec,
                               k = 2, lambda = 1e6),
               regexp = "lambda")
})

test_that("T-VAL-09/10/11: compute_posterior rejects malformed input", {
  esvd_obj <- .small_esvd_obj()

  # T-VAL-11: the two experimental booleans are mutually exclusive.
  expect_error(compute_posterior(input_obj = esvd_obj,
                                 bool_adjust_covariates = TRUE,
                                 bool_covariates_as_library = TRUE))

  # T-VAL-10: a negative nuisance value.
  broken_obj <- esvd_obj
  latest_fit <- broken_obj[["latest_Fit"]]
  broken_obj[[latest_fit]]$nuisance_vec[1] <- -1
  expect_error(compute_posterior(input_obj = broken_obj,
                                 alpha_max = 100,
                                 bool_covariates_as_library = TRUE))
})

test_that("T-VAL-15/16: compute_pvalue names what is missing", {
  esvd_obj <- .small_esvd_obj()

  no_teststat_obj <- esvd_obj
  no_teststat_obj$teststat_vec <- NULL
  expect_error(compute_pvalue(input_obj = no_teststat_obj),
               regexp = "teststat")
})

test_that("T-VAL-21: esvd_family lists the seven valid names on an unknown one", {
  expect_error(esvd_family("gaussain"), regexp = "gaussian")
})

test_that("T-VAL-26/27: eSVD rejects malformed top-level arguments", {
  seurat_obj <- .tiny_seurat()

  # T-VAL-26: case_control_levels not length 2.
  expect_error(eSVD(batch_var_prefix = NULL,
                    case_control_levels = c("0"),
                    case_control_var = "CC",
                    categorical_vars = "Sex",
                    id_var = "Individual",
                    numerical_vars = "Age",
                    seurat_obj = seurat_obj))

  # T-VAL-27: duplicated categorical_vars.
  expect_error(eSVD(batch_var_prefix = NULL,
                    case_control_levels = c("0", "1"),
                    case_control_var = "CC",
                    categorical_vars = c("Sex", "Sex"),
                    id_var = "Individual",
                    numerical_vars = "Age",
                    seurat_obj = seurat_obj))
})

test_that("T-VAL-28: reparameterization names the available fits", {
  esvd_obj <- .small_esvd_obj()

  expect_error(reparameterization_esvd_covariates(input_obj = esvd_obj,
                                                  fit_name = "fit_Nonexistent"),
               regexp = "fit_Nonexistent|fit_First|fit_Init")
})

## NEW, from section 7: `eSVD()` calls `SeuratObject::LayerData()` while
## `SeuratObject` sits in `Suggests`, and there is NOT ONE `requireNamespace()`
## call anywhere in `R/`. CRAN requires conditional use of a suggested package.
test_that("T-VAL-34: R/ guards its use of the suggested SeuratObject", {
  r_dir <- system.file("R", package = "eSVD2")
  source_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(source_dir), "package source not reachable")

  source_text <- unlist(lapply(list.files(source_dir, pattern = "[.]R$",
                                          full.names = TRUE), readLines))
  uses_seurat <- any(grepl("SeuratObject::", source_text, fixed = TRUE))
  has_guard <- any(grepl("requireNamespace", source_text, fixed = TRUE))

  expect_true(!uses_seurat || has_guard,
              info = "SeuratObject is used from Suggests with no requireNamespace guard")
})

## NEW, from section 7: `NAMESPACE` has no `export(eSVD)`. The package's main
## user-facing function is exported ONLY by the blanket
## `exportPattern("^[[:alpha:]]+")`, which question Q-CPP-1 removes -- so that
## change would silently un-export `eSVD()`. One `@export` tag, but it has to
## land in the same commit.
test_that("T-ESVD-11: eSVD is explicitly exported", {
  namespace_path <- testthat::test_path("..", "..", "NAMESPACE")
  skip_if_not(file.exists(namespace_path), "NAMESPACE not reachable")

  namespace_text <- readLines(namespace_path)
  expect_true(any(grepl("^export(eSVD)$", namespace_text)) ||
                any(grepl("^export\\(eSVD\\)$", namespace_text)),
              info = "eSVD is exported only by the blanket exportPattern")
})

# ---- section 6, verbose branches --------------------------------------------

## Section 1.2, [verified]: `opt_esvd(verbose = 2)` throws. The per-iteration
## branch references something that is not in scope.
test_that("T-VERB-01: opt_esvd(verbose = 2) does not throw", {
  point <- .feasible_point("poisson", num_cells = 20, num_genes = 8)

  expect_no_error(suppressWarnings(utils::capture.output(
    opt_esvd(input_obj = point$dat, x_init = point$x_mat,
             y_init = point$y_mat, family = "poisson", max_iter = 2,
             nuisance_vec = point$gamma_vec, verbose = 2)
  )))
})

test_that("T-VERB-02: verbose = 1 prints and verbose = 0 is silent", {
  esvd_obj <- .small_esvd_obj()

  expect_silent(compute_test_statistic(input_obj = esvd_obj, verbose = 0))
  expect_output(compute_test_statistic(input_obj = esvd_obj, verbose = 1))
})
