# UNIT_TEST_PLAN.md section 2.8 -- T-TSTAT-02 .. T-TSTAT-07.
#
# T-TSTAT-01 and T-DF-01 live in `test_compute_pvalue_claude.R`, beside the
# `stats::t.test` oracle they share.

test_that("T-TSTAT-02: an individual in both arms errors", {
  fixture_n <- 18
  set.seed(10)
  posterior_mean_mat <- matrix(stats::rnorm(fixture_n * 3), nrow = fixture_n,
                               ncol = 3,
                               dimnames = list(paste0("cell_", 1:fixture_n),
                                               paste0("gene_", 1:3)))
  posterior_var_mat <- matrix(0.1, nrow = fixture_n, ncol = 3,
                              dimnames = dimnames(posterior_mean_mat))
  individual_vec <- factor(rep(paste0("indiv_", 1:6), each = 3))

  # `indiv_3` deliberately appears in both arms.
  expect_error(
    compute_test_statistic(input_obj = posterior_mean_mat,
                           posterior_var_mat = posterior_var_mat,
                           case_individuals = c("indiv_1", "indiv_2", "indiv_3"),
                           control_individuals = c("indiv_3", "indiv_4",
                                                   "indiv_5"),
                           individual_vec = individual_vec,
                           verbose = 0)
  )
})

## Resolved 2026-08-29: an arm with one individual is now a hard error rather
## than a NaN degrees-of-freedom. `.compute_df()`'s denominator contains
## `(v/n)^2/(n-1)`, so `n1 = 1` gives df = 0 and every `stats::pt()` downstream
## returns NaN -- which used to surface as "missing values and NaN's not
## allowed" from inside `multtest`, naming nothing.
test_that("T-TSTAT-04: an arm with one individual errors, naming the counts", {
  set.seed(10)
  n <- 12
  posterior_mean_mat <- matrix(stats::rnorm(n * 3), nrow = n, ncol = 3,
                               dimnames = list(paste0("cell_", 1:n),
                                               paste0("gene_", 1:3)))
  posterior_var_mat <- matrix(0.1, nrow = n, ncol = 3,
                              dimnames = dimnames(posterior_mean_mat))
  individual_vec <- factor(rep(paste0("indiv_", 1:4), each = 3))

  expect_error(
    compute_test_statistic(input_obj = posterior_mean_mat,
                           posterior_var_mat = posterior_var_mat,
                           case_individuals = c("indiv_1"),
                           control_individuals = c("indiv_2", "indiv_3",
                                                   "indiv_4"),
                           individual_vec = individual_vec,
                           verbose = 0),
    regexp = "at least 2 individuals"
  )
})

## The companion guard: individuals with too few cells are meant to be dropped
## by `eSVD_helper`, so reaching the test statistic with one is an error
## (Kevin, 2026-08-29). Note this is a POWER/stability rule, not a numerical
## one -- verified that such a cohort returns finite statistics; the refusal is
## a deliberate policy about what the method is willing to be asked.
test_that("T-TSTAT-04b: a donor with too few cells errors, and 0 disables it", {
  set.seed(10)
  cell_count_vec <- c(20, 20, 20, 20, 20, 1)
  n <- sum(cell_count_vec)
  individual_vec <- factor(rep(paste0("indiv_", 1:6), times = cell_count_vec),
                           levels = paste0("indiv_", 1:6))
  posterior_mean_mat <- matrix(stats::rnorm(n * 3, mean = 5), nrow = n, ncol = 3,
                               dimnames = list(paste0("cell_", 1:n),
                                               paste0("gene_", 1:3)))
  posterior_var_mat <- matrix(0.1, nrow = n, ncol = 3,
                              dimnames = dimnames(posterior_mean_mat))

  expect_error(
    compute_test_statistic(input_obj = posterior_mean_mat,
                           posterior_var_mat = posterior_var_mat,
                           case_individuals = paste0("indiv_", 1:3),
                           control_individuals = paste0("indiv_", 4:6),
                           individual_vec = individual_vec,
                           verbose = 0),
    regexp = "indiv_6"
  )

  # The escape hatch has to work, since it is the only one.
  res <- compute_test_statistic(input_obj = posterior_mean_mat,
                                posterior_var_mat = posterior_var_mat,
                                case_individuals = paste0("indiv_", 1:3),
                                control_individuals = paste0("indiv_", 4:6),
                                individual_vec = individual_vec,
                                min_cells_per_individual = 0,
                                verbose = 0)
  expect_true(all(is.finite(res$teststat_vec)))
})

test_that("T-TSTAT-05: every output carries colnames(posterior_mean_mat)", {
  esvd_obj <- .small_esvd_obj()
  gene_vec <- colnames(esvd_obj$dat)

  # `report_results` builds its `genes` column from `names(teststat_vec)`.
  expect_equal(names(esvd_obj$teststat_vec), gene_vec)
  expect_equal(names(esvd_obj$case_mean), gene_vec)
  expect_equal(names(esvd_obj$control_mean), gene_vec)
})

test_that("T-TSTAT-06: permuting rows and individual_vec together changes nothing", {
  set.seed(10)
  n <- 24
  posterior_mean_mat <- matrix(stats::rnorm(n * 4), nrow = n, ncol = 4,
                               dimnames = list(paste0("cell_", 1:n),
                                               paste0("gene_", 1:4)))
  posterior_var_mat <- matrix(stats::runif(n * 4, 0.05, 0.2), nrow = n, ncol = 4,
                              dimnames = dimnames(posterior_mean_mat))
  individual_vec <- factor(rep(paste0("indiv_", 1:6), each = 4))
  case_individuals <- paste0("indiv_", 1:3)
  control_individuals <- paste0("indiv_", 4:6)

  res_original <- compute_test_statistic(
    input_obj = posterior_mean_mat,
    posterior_var_mat = posterior_var_mat,
    case_individuals = case_individuals,
    control_individuals = control_individuals,
    individual_vec = individual_vec,
    verbose = 0
  )

  set.seed(20)
  permutation_idx <- sample(n)
  res_permuted <- compute_test_statistic(
    input_obj = posterior_mean_mat[permutation_idx, , drop = FALSE],
    posterior_var_mat = posterior_var_mat[permutation_idx, , drop = FALSE],
    case_individuals = case_individuals,
    control_individuals = control_individuals,
    individual_vec = individual_vec[permutation_idx],
    verbose = 0
  )

  # Exchangeability. Catches any accidental positional (rather than name-based)
  # indexing in the averaging matrix.
  expect_equal(res_original$teststat_vec, res_permuted$teststat_vec,
               tolerance = 1e-10)
})

## An `idx_list` entry of length 0 is reachable via an unused factor level. The
## resulting all-zero row of the averaging matrix silently averages to 0, which
## biases the group mean with no signal that anything went wrong.
## Expected to FAIL -- there is no guard.
test_that("T-TSTAT-07: .construct_averaging_matrix rejects an empty index set", {
  idx_list <- list(c(1, 2, 3), integer(0), c(4, 5))

  expect_error(.construct_averaging_matrix(idx_list = idx_list, n = 5))
})

test_that("T-TSTAT-08: .construct_averaging_matrix rows sum to 1", {
  idx_list <- list(c(1, 2, 3), c(4, 5))
  avg_mat <- .construct_averaging_matrix(idx_list = idx_list, n = 5)

  expect_equal(as.numeric(Matrix::rowSums(avg_mat)), rep(1, length(idx_list)))
  expect_equal(dim(avg_mat), c(length(idx_list), 5L))

  # And it does compute the within-individual mean.
  set.seed(10)
  mat <- matrix(stats::rnorm(5 * 2), nrow = 5, ncol = 2)
  expect_equal(as.matrix(avg_mat %*% mat),
               rbind(colMeans(mat[1:3, , drop = FALSE]),
                     colMeans(mat[4:5, , drop = FALSE])),
               tolerance = 1e-10)
})

test_that("T-TSTAT-09: .determine_individual_indices partitions the cells", {
  individual_vec <- factor(rep(paste0("indiv_", 1:6), each = 4))
  case_individuals <- paste0("indiv_", 1:3)
  control_individuals <- paste0("indiv_", 4:6)

  res <- .determine_individual_indices(case_individuals = case_individuals,
                                       control_individuals = control_individuals,
                                       individual_vec = individual_vec)

  all_idx <- unlist(c(res$case_indiv_idx, res$control_indiv_idx))
  expect_equal(sort(all_idx), seq_along(individual_vec))
  expect_equal(anyDuplicated(all_idx), 0L)
})

test_that("T-TSTAT-10: a length mismatch between individual_vec and the matrix errors", {
  set.seed(10)
  posterior_mean_mat <- matrix(stats::rnorm(24), nrow = 12, ncol = 2,
                               dimnames = list(paste0("cell_", 1:12),
                                               paste0("gene_", 1:2)))
  posterior_var_mat <- matrix(0.1, nrow = 12, ncol = 2,
                              dimnames = dimnames(posterior_mean_mat))

  expect_error(
    compute_test_statistic(input_obj = posterior_mean_mat,
                           posterior_var_mat = posterior_var_mat,
                           case_individuals = paste0("indiv_", 1:3),
                           control_individuals = paste0("indiv_", 4:6),
                           individual_vec = factor(rep(paste0("indiv_", 1:6),
                                                       each = 3)),
                           verbose = 0)
  )
})
