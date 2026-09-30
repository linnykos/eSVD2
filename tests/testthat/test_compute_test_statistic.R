context("Test compute_test_statistic")

## compute_test_statistic is correct

test_that("compute_test_statistic works", {
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

  res <- compute_test_statistic(input_obj = eSVD_obj$fit_First$posterior_mean_mat,
                                posterior_var_mat = eSVD_obj$fit_First$posterior_var_mat,
                                case_individuals = case_individuals,
                                control_individuals = control_individuals,
                                covariate_individual = "individual",
                                individual_vec = individual_vec,
                                metadata = metadata)

  mean_val <- mean(res$teststat_vec)
  sd_val <- sd(res$teststat_vec)
  quantile_vec <- sapply(res$teststat_vec, function(x){
    1-2*abs(stats::pnorm(x, mean = mean_val, sd = sd_val)-.5)
  })
  expect_true(mean(abs(quantile_vec[true_cc_status == 1])) > 2*mean(abs(quantile_vec[true_cc_status == 2])))
})

######################

## .determine_individual_indices is correct

test_that(".determine_individual_indices works", {
  n <- 100
  metadata <- data.frame(individual = factor(rep(1:4, each = n/4)))
  rownames(metadata) <- paste0("c", 1:n)
  res <- .determine_individual_indices(case_individuals = c("1", "2"),
                                       control_individuals = c("3", "4"),
                                       individual_vec = metadata[,"individual"])

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("case_indiv_idx", "control_indiv_idx"))))

  tmp <- c(res$case_indiv_idx, res$control_indiv_idx)
  expect_true(is.list(tmp))
  expect_true(length(tmp) == 4)
  expect_true(length(unique(unlist(tmp))) == n)
})

#######################

## .construct_averaging_matrix is correct

test_that(".construct_averaging_matrix works", {
  set.seed(10)
  vec <- numeric(0)
  for(i in 1:4){
    vec <- c(vec, rep(i, round(50*runif(1))))
  }
  metadata <- data.frame(individual = factor(vec))
  n <- length(vec)
  rownames(metadata) <- paste0("c", 1:n)
  tmp <- .determine_individual_indices(case_individuals = c("1", "2"),
                                       control_individuals = c("3", "4"),
                                       individual_vec = metadata[,"individual"])
  all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
  res <- .construct_averaging_matrix(idx_list = all_indiv_idx,
                                     n = n)

  expect_true(inherits(res, "dgCMatrix"))
  expect_true(all(dim(res) == c(4,n)))
  res2 <- as.matrix(res)
  for(i in 1:4){
    idx <- which(metadata[,"individual"] == as.character(i))
    expect_true(all(res2[i,-idx] == 0))
  }

  response_vec <- runif(n)
  avg1 <- as.numeric(res %*% response_vec)
  avg2 <- sapply(1:4, function(i){
    idx <- which(metadata[,"individual"] == as.character(i))
    mean(response_vec[idx])
  })
  expect_true(sum(abs(avg1 - avg2)) <= 1e-5)
})

###############################

## .compute_mixture_gaussian_variance is correct

test_that(".compute_mixture_gaussian_variance works", {
  set.seed(10)
  n <- 5
  p <- 10
  avg_posterior_mean_mat <- matrix(runif(n*p), nrow = n, ncol = p)
  avg_posterior_var_mat <- matrix(runif(n*p), nrow = n, ncol = p)

  res <- .compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat,
    avg_posterior_var_mat = avg_posterior_var_mat
  )
  res2 <- sapply(1:p, function(j){
    mean(avg_posterior_var_mat[,j]) + mean(avg_posterior_mean_mat[,j]^2) - mean(avg_posterior_mean_mat[,j])^2
  })

  expect_true(sum(abs(res - res2)) <= 1e-6)
})

################################################################################

# UNIT_TEST_PLAN.md section 2.8 -- T-TSTAT-02 .. T-TSTAT-07.
#
# T-TSTAT-01 and T-DF-01 live in `test_compute_pvalue.R`, beside the
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

  # The matrix is mean-zero Gaussian, so some arm means are negative. That is
  # fine for the statistic and undefined for a log fold change, which since
  # version 1.1.0 is returned as NA with a warning (T-LFC-08).
  expect_warning(
    res_original <- compute_test_statistic(
      input_obj = posterior_mean_mat,
      posterior_var_mat = posterior_var_mat,
      case_individuals = case_individuals,
      control_individuals = control_individuals,
      individual_vec = individual_vec,
      verbose = 0
    ),
    regexp = "not positive"
  )

  set.seed(20)
  permutation_idx <- sample(n)
  expect_warning(
    res_permuted <- compute_test_statistic(
      input_obj = posterior_mean_mat[permutation_idx, , drop = FALSE],
      posterior_var_mat = posterior_var_mat[permutation_idx, , drop = FALSE],
      case_individuals = case_individuals,
      control_individuals = control_individuals,
      individual_vec = individual_vec[permutation_idx],
      verbose = 0
    ),
    regexp = "not positive"
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
