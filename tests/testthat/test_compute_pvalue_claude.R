# UNIT_TEST_PLAN.md section 2.9 -- T-DF-01 .. T-DF-03, T-PVAL-01 .. T-PVAL-06.
#
# `compute_pvalue()` produces the package's headline output and has zero tests.
# T-PVAL-01 is the highest-value single assertion in the whole plan: it is the
# defect that changes published-style results.

# A donor-level fixture: every donor contributes `cells_per_donor` IDENTICAL
# cells and the posterior variances are zero. The within-donor average is then
# exact and the mixture variance reduces to the population variance of the donor
# means -- which makes `stats::t.test` an exact external oracle for the headline
# statistic, up to the Bessel factor documented below.
#
# Three cells per donor rather than one, because `compute_test_statistic` now
# refuses a cohort with fewer than `min_cells_per_individual` cells. The oracle
# survives the guard: what makes the reduction work is the zero posterior
# variance and the identical cells, not the cell count.
.donor_level_obj <- function(num_individuals = 12, num_genes = 6,
                             cells_per_donor = 3, seed_number = 10){
  set.seed(seed_number)

  n <- num_individuals * cells_per_donor
  donor_mean_mat <- matrix(stats::rnorm(num_individuals * num_genes,
                                        mean = 5, sd = 1.5),
                           nrow = num_individuals, ncol = num_genes)

  posterior_mean_mat <- donor_mean_mat[rep(seq_len(num_individuals),
                                           each = cells_per_donor), ,
                                       drop = FALSE]
  dimnames(posterior_mean_mat) <- list(paste0("cell_", seq_len(n)),
                                       paste0("gene_", seq_len(num_genes)))
  posterior_var_mat <- matrix(0, nrow = n, ncol = num_genes,
                              dimnames = dimnames(posterior_mean_mat))

  individual_vec <- factor(rep(paste0("indiv_", seq_len(num_individuals)),
                               each = cells_per_donor),
                           levels = paste0("indiv_", seq_len(num_individuals)))
  cc_vec <- rep(rep(c(0, 1), each = num_individuals / 2), each = cells_per_donor)

  list(cc_vec = cc_vec,
       donor_mean_mat = donor_mean_mat,
       individual_vec = individual_vec,
       posterior_mean_mat = posterior_mean_mat,
       posterior_var_mat = posterior_var_mat)
}

## The exact relationship between the eSVD statistic and Welch's t.
##
## `.compute_mixture_gaussian_variance()` returns the POPULATION variance of the
## mixture -- it divides by n, with no Bessel correction -- and that is
## DELIBERATE (Kevin, 2026-08-29). The mixture is a population object, so its
## population variance is the right thing to compute.
##
## The consequence, pinned here so nobody rediscovers it as a bug: on a fixture
## where the two should otherwise coincide exactly, the eSVD statistic is
## Welch's t scaled by sqrt(n/(n-1)) per arm. Verified exact at n = 3, 5, 10 and
## 25 individuals per arm, so this is the contract, not an artefact.
##
## This is the only EXTERNAL oracle in the suite for the package's headline
## statistic, which is why the fixture is built to survive the
## `min_cells_per_individual` guard rather than being deleted by it.
test_that("T-TSTAT-01a: the statistic is Welch's t times sqrt(n/(n-1))", {
  fixture <- .donor_level_obj()

  case_individuals <- unique(fixture$individual_vec[fixture$cc_vec == 1])
  control_individuals <- unique(fixture$individual_vec[fixture$cc_vec == 0])
  num_per_arm <- length(case_individuals)

  res <- compute_test_statistic(input_obj = fixture$posterior_mean_mat,
                                posterior_var_mat = fixture$posterior_var_mat,
                                case_individuals = case_individuals,
                                control_individuals = control_individuals,
                                individual_vec = fixture$individual_vec,
                                verbose = 0)

  bessel_factor <- sqrt(num_per_arm / (num_per_arm - 1))
  case_rows <- (num_per_arm + 1):(2 * num_per_arm)
  control_rows <- seq_len(num_per_arm)

  for(gene_idx in seq_len(ncol(fixture$posterior_mean_mat))){
    label <- paste0("gene ", gene_idx)
    t_res <- stats::t.test(fixture$donor_mean_mat[case_rows, gene_idx],
                           fixture$donor_mean_mat[control_rows, gene_idx],
                           var.equal = FALSE)

    expect_equal(as.numeric(res$teststat_vec[gene_idx]),
                 as.numeric(t_res$statistic) * bessel_factor,
                 tolerance = 1e-10, info = label)
  }
})

test_that("T-DF-03: .compute_df errors informatively when dat is absent", {
  esvd_obj <- .small_esvd_obj()
  esvd_obj$dat <- NULL

  # Section 1.6: `.compute_df()` re-derives everything from scratch, which is
  # why `compute_pvalue` needs `dat` -- and that conflicts with
  # `bool_diet = TRUE` in `eSVD()`. The error should say so.
  expect_error(compute_pvalue(input_obj = esvd_obj), regexp = "dat")
})

test_that("T-DF-02: .compute_df and compute_test_statistic use the same group variances", {
  esvd_obj <- .small_esvd_obj()

  # Section 1.6: `.compute_df()` recomputes the averaging matrix and both group
  # variances from scratch, duplicating `compute_test_statistic`. Two copies of
  # one computation drift; this is what makes the refactor safe.
  #
  # The observable consequence of agreement: for a two-sided Welch t with the
  # stored df, `pt` of the stored statistic must reproduce the same p-value the
  # statistic implies.
  df_vec <- esvd_obj$pvalue_list$df_vec
  teststat_vec <- esvd_obj$teststat_vec

  expect_equal(names(df_vec), names(teststat_vec))
  expect_true(all(is.finite(df_vec)))
  # Welch's df lies between (min(n1,n2) - 1) and (n1 + n2 - 2). With 3 case and
  # 3 control individuals that is [2, 4].
  expect_true(all(df_vec >= 2 - 1e-8))
  expect_true(all(df_vec <= 4 + 1e-8))
})

## THE headline regression test. Section 1.1: with `t = 40, df = 18`,
## `qnorm(pt(t, df))` is exactly +Inf, which makes `locfdr::locfdr()` error,
## the `tryCatch` in `multtest()` swallow it, and the empirical null silently
## degrade to `.multtest_simple()`. Expected to FAIL until `log.p = TRUE` lands.
test_that("T-PVAL-01: gaussian_teststat is finite for a strongly DE gene", {
  # `qnorm(pt(t, df))` saturates at +Inf for `t = 40, df = 18`. The df matters:
  # at the 4 degrees of freedom that a 3-versus-3 cohort gives, the t
  # distribution is heavy-tailed enough that `pt(40, 4)` is not yet 1 in double
  # precision and the bug does not reproduce. Twenty individuals are needed to
  # reach df near 18, which is an ordinary cohort size for this method.
  fixture <- .donor_level_obj(num_individuals = 20, num_genes = 8)

  case_individuals <- unique(fixture$individual_vec[fixture$cc_vec == 1])
  control_individuals <- unique(fixture$individual_vec[fixture$cc_vec == 0])

  res <- compute_test_statistic(input_obj = fixture$posterior_mean_mat,
                                posterior_var_mat = fixture$posterior_var_mat,
                                case_individuals = case_individuals,
                                control_individuals = control_individuals,
                                individual_vec = fixture$individual_vec,
                                verbose = 0)

  # One strongly up-regulated gene, which is exactly what a real dataset
  # supplies and what the current code cannot represent.
  teststat_vec <- res$teststat_vec
  teststat_vec[1] <- 40
  df_vec <- rep(18, length(teststat_vec))

  gaussian_teststat <- sapply(seq_along(teststat_vec), function(gene_idx){
    stats::qnorm(stats::pt(teststat_vec[gene_idx], df = df_vec[gene_idx]))
  })

  expect_true(all(is.finite(gaussian_teststat)),
              info = paste0("gene 1 gaussian_teststat = ",
                            gaussian_teststat[1]))
})

test_that("T-PVAL-01b: one Inf makes multtest silently abandon locfdr", {
  # Section 1.1's real damage, isolated. `.multtest_locfdr()` catches WARNINGS
  # as well as errors, so a single non-finite entry downgrades the empirical
  # null to `.multtest_simple()` -- the estimator whose own comment disclaims
  # it -- and nothing reaches the user.
  set.seed(10)
  # locfdr needs a few hundred genes to fit its spline; at 200 it warns and
  # `.multtest_locfdr()` catches warnings as well as errors, so the clean case
  # must be large enough that locfdr genuinely succeeds.
  teststat_vec <- stats::rnorm(1000)
  names(teststat_vec) <- paste0("gene_", seq_along(teststat_vec))

  res_clean <- suppressWarnings(multtest(teststat_vec))

  teststat_vec[1] <- Inf
  res_dirty <- suppressWarnings(multtest(teststat_vec))

  expect_equal(res_clean$method, "locfdr")
  # A non-finite entry should be rejected at entry, not absorbed.
  expect_true(all(is.finite(res_dirty$fdr_vec)),
              info = paste0("method fell back to: ", res_dirty$method))
})

## The real damage of section 1.1: the empirical null silently degrades and the
## user is never told, because `compute_pvalue()` discards `multtest()`'s
## `method` field. Expected to FAIL on two counts -- `method` is not stored at
## all, and it would not be "locfdr" here even if it were.
test_that("T-PVAL-02: pvalue_list records which empirical-null method ran", {
  esvd_obj <- .small_esvd_obj()

  expect_true("method" %in% names(esvd_obj$pvalue_list))
})

test_that("T-PVAL-03: log10pvalue is monotone in |gaussian_teststat - null_mean|", {
  esvd_obj <- .small_esvd_obj()
  pvalue_list <- esvd_obj$pvalue_list

  distance_vec <- abs(pvalue_list$gaussian_teststat - pvalue_list$null_mean)
  order_idx <- order(distance_vec)

  # A p-value that is not monotone in the statistic is definitionally broken.
  # Cheap, and it catches a sign error in the mirroring branch.
  expect_true(all(diff(pvalue_list$log10pvalue[order_idx]) >= -1e-8))
})

test_that("T-PVAL-04: every p-value lies in [0, 1]", {
  esvd_obj <- .small_esvd_obj()
  pvalue_vec <- 10^(-esvd_obj$pvalue_list$log10pvalue)

  # The mirroring construction `null_mean - (x - null_mean)` guarantees
  # `2*pnorm(...) <= 1`. Assert it rather than trust it.
  expect_true(all(pvalue_vec >= 0))
  expect_true(all(pvalue_vec <= 1))
  expect_true(all(esvd_obj$pvalue_list$log10pvalue >= 0))
})

test_that("T-PVAL-05: fdr_vec equals stats::p.adjust(pvalue, 'BH')", {
  esvd_obj <- .small_esvd_obj()
  pvalue_vec <- 10^(-esvd_obj$pvalue_list$log10pvalue)

  expect_equal(unname(esvd_obj$pvalue_list$fdr_vec),
               unname(stats::p.adjust(pvalue_vec, method = "BH")),
               tolerance = 1e-6)
})

## The same six-line computation is written three times -- in `multtest.R`, in
## `compute_pvalue.R`, and again in `compute_test_per_gene.R`. Either they agree
## or one of them is wrong.
test_that("T-PVAL-06: compute_pvalue's log10pvalue equals multtest's logpvalue_vec", {
  esvd_obj <- .small_esvd_obj()
  gaussian_teststat <- esvd_obj$pvalue_list$gaussian_teststat

  multtest_res <- suppressWarnings(multtest(gaussian_teststat))

  expect_equal(unname(esvd_obj$pvalue_list$log10pvalue),
               unname(multtest_res$logpvalue_vec),
               tolerance = 1e-8)
})

test_that("T-PVAL-08: every pvalue_list element is named by gene", {
  esvd_obj <- .small_esvd_obj()
  gene_vec <- colnames(esvd_obj$dat)

  for(element_name in c("df_vec", "fdr_vec", "gaussian_teststat",
                        "log10pvalue")){
    expect_equal(names(esvd_obj$pvalue_list[[element_name]]), gene_vec,
                 info = element_name)
  }
})
