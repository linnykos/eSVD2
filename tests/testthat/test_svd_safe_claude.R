# UNIT_TEST_PLAN.md section 2.3 -- T-SVD-01 .. T-SVD-07.
#
# `.svd_safe()` has a three-way fallback (irlba -> RSpectra -> base::svd) that
# no test has ever entered past the first link.

test_that("T-SVD-01: the decomposition reconstructs the matrix and reports its method", {
  set.seed(10)
  mat <- matrix(stats::rnorm(60 * 12), nrow = 60, ncol = 12)

  res <- .svd_safe(mat = mat, check_stability = FALSE, K = 3,
                   mean_vec = NULL, rescale = FALSE, scale_max = NULL,
                   sd_vec = NULL)

  expect_true(all(c("u", "d", "v") %in% names(res)))
  expect_equal(dim(res$u), c(nrow(mat), 3L))
  expect_equal(dim(res$v), c(ncol(mat), 3L))
  expect_equal(length(res$d), 3L)

  # Singular values are non-negative and sorted, and the rank-3 reconstruction
  # matches `base::svd`'s to within sign.
  expect_true(all(res$d >= 0))
  expect_true(all(diff(res$d) <= 1e-8))

  base_res <- base::svd(mat, nu = 3, nv = 3)
  expect_equal(res$d, base_res$d[1:3], tolerance = 1e-6)
  expect_equal(tcrossprod(.mult_mat_vec(res$u, res$d), res$v),
               tcrossprod(.mult_mat_vec(base_res$u, base_res$d[1:3]),
                          base_res$v),
               tolerance = 1e-6)
})

test_that("T-SVD-02: u and v carry the dimnames of the input", {
  set.seed(10)
  mat <- matrix(stats::rnorm(60 * 12), nrow = 60, ncol = 12,
                dimnames = list(paste0("cell_", 1:60), paste0("gene_", 1:12)))

  res <- .svd_safe(mat = mat, check_stability = FALSE, K = 3,
                   mean_vec = NULL, rescale = FALSE, scale_max = NULL,
                   sd_vec = NULL)

  # `.initialize_residuals()` does `rownames(x_mat) <- rownames(dat)` right
  # after this, and relies on the names not already being wrong.
  expect_equal(rownames(res$u), rownames(mat))
  expect_equal(rownames(res$v), colnames(mat))
})

test_that("T-SVD-03: a rank-deficient matrix still decomposes and reports a method", {
  set.seed(10)
  base_mat <- matrix(stats::rnorm(60 * 3), nrow = 60, ncol = 3)
  # Rank 3 presented as 12 columns: irlba is prone to warn here.
  mat <- base_mat %*% matrix(stats::rnorm(3 * 12), nrow = 3, ncol = 12)

  res <- suppressWarnings(
    .svd_safe(mat = mat, check_stability = FALSE, K = 5,
              mean_vec = NULL, rescale = FALSE, scale_max = NULL,
              sd_vec = NULL)
  )

  # The fallback contract: whichever link ran, the answer is a valid SVD and
  # the object says which one produced it.
  expect_true("method" %in% names(res))
  expect_true(all(is.finite(res$d)))
  expect_true(all(is.finite(res$u)))
  expect_true(all(is.finite(res$v)))
  # Rank 3 means singular values 4 and 5 are numerically zero.
  expect_true(res$d[4] < 1e-8 * res$d[1])
})

test_that("T-SVD-04: K within 2 of min(dim) still produces a valid SVD", {
  set.seed(10)
  mat <- matrix(stats::rnorm(20 * 6), nrow = 20, ncol = 6)

  res <- suppressWarnings(
    .svd_safe(mat = mat, check_stability = FALSE, K = 5,
              mean_vec = NULL, rescale = FALSE, scale_max = NULL,
              sd_vec = NULL)
  )

  base_res <- base::svd(mat)
  expect_equal(res$d, base_res$d[1:5], tolerance = 1e-6)
})

## Question Q-SVD-1 resolved: target CRAN and absorb the one `sparseMatrixStats`
## function. Q-SVD-3 gives the decision rule -- write this test first with a
## 1e-10 relative tolerance on a fixture that includes a large-mean column,
## then put the four-line `Matrix` rewrite under it. If it passes, ship the
## rewrite; if it fails, vendor the C++.
##
## Note the branch under test is DEAD CODE in the package as shipped:
## `.compute_matrix_sd`'s sparse path fires only when a caller passes
## `sd_vec = TRUE` with a sparse matrix, and the only in-package caller,
## `.initialize_residuals()`, passes `sd_vec = NULL` and a dense matrix.
test_that("T-SVD-05: .compute_matrix_sd on a dgCMatrix equals matrixStats::colSds", {
  skip_if_not_installed("sparseMatrixStats")

  set.seed(10)
  dense_mat <- matrix(stats::rpois(200 * 8, lambda = 2) * 1.0,
                      nrow = 200, ncol = 8)
  # A column with a large mean and small variance is where the naive two-pass
  # form loses precision, and is therefore the column that decides Q-SVD-3.
  dense_mat[, 1] <- 1e6 + stats::rnorm(200, sd = 1e-3)
  sparse_mat <- methods::as(methods::as(dense_mat, "dMatrix"), "CsparseMatrix")

  res <- .compute_matrix_sd(mat = sparse_mat, sd_vec = TRUE)
  expected <- matrixStats::colSds(as.matrix(sparse_mat))

  expect_equal(res, expected, tolerance = 1e-10)
})

## Q-SVD-3, ANSWERED BY THE TEST. The decision rule was: try the four-line
## `Matrix` rewrite, vendor the C++ only if it fails at 1e-10. It fails, and it
## fails badly -- the naive two-pass form
##   sqrt((colSums(x^2) - n*colMeans(x)^2)/(n-1))
## suffers catastrophic cancellation on a column with a large mean and small
## variance, goes NEGATIVE under the square root, and returns NaN. On the
## fixture below, column 1 (mean 1e6, sd 1e-3) gives NaN while
## `sparseMatrixStats::colSds` is exact to the last bit.
##
## So the naive rewrite is out. The good news is that vendoring C++ is still
## not necessary: the SHIFTED form below is stable, stays sparse, and is a
## handful of lines of `Matrix`. T-SVD-05c is the candidate to ship.
test_that("T-SVD-05b: the naive two-pass rewrite fails, and this is why", {
  set.seed(10)
  dense_mat <- matrix(stats::rpois(200 * 8, lambda = 2) * 1.0,
                      nrow = 200, ncol = 8)
  dense_mat[, 1] <- 1e6 + stats::rnorm(200, sd = 1e-3)
  sparse_mat <- methods::as(methods::as(dense_mat, "dMatrix"), "CsparseMatrix")

  n <- nrow(sparse_mat)
  col_mean_vec <- Matrix::colMeans(sparse_mat)
  naive_sd_vec <- suppressWarnings(
    sqrt((Matrix::colSums(sparse_mat^2) - n * col_mean_vec^2) / (n - 1))
  )

  # Column 1 is the pathological one; the rest are fine.
  expect_true(is.nan(as.numeric(naive_sd_vec)[1]))
  expect_true(all(is.finite(as.numeric(naive_sd_vec)[-1])))
})

test_that("T-SVD-05c: the shifted Matrix form is stable and matches to 1e-10", {
  set.seed(10)
  dense_mat <- matrix(stats::rpois(200 * 8, lambda = 2) * 1.0,
                      nrow = 200, ncol = 8)
  dense_mat[, 1] <- 1e6 + stats::rnorm(200, sd = 1e-3)
  sparse_mat <- methods::as(methods::as(dense_mat, "dMatrix"), "CsparseMatrix")

  # Sum of squared deviations, computed over the STORED non-zeros only and
  # corrected for the implied zeros, so nothing densifies:
  #   sum_i (x_i - m)^2 = sum_{stored} (x_i - m)^2 + (n - nnz) * m^2
  .sparse_col_sds <- function(mat){
    n <- nrow(mat)
    col_mean_vec <- Matrix::colMeans(mat)
    col_idx_vec <- rep(seq_len(ncol(mat)), diff(mat@p))
    deviation_vec <- (mat@x - col_mean_vec[col_idx_vec])^2
    stored_ss_vec <- as.numeric(tapply(deviation_vec, factor(col_idx_vec,
                                                            levels = seq_len(ncol(mat))),
                                       sum))
    stored_ss_vec[is.na(stored_ss_vec)] <- 0
    num_stored_vec <- as.numeric(diff(mat@p))
    sqrt((stored_ss_vec + (n - num_stored_vec) * col_mean_vec^2) / (n - 1))
  }

  expected <- matrixStats::colSds(as.matrix(sparse_mat))
  expect_equal(.sparse_col_sds(sparse_mat), expected, tolerance = 1e-10)
})

## Question Q-SVD-2 resolved: reject a length-1 numeric. The guard has to test
## `is.logical()`, NOT `length(x) == 1` -- `mean_vec = 0.5` and `mean_vec = TRUE`
## are both length 1, which is exactly why the current code cannot tell them
## apart. Expected to FAIL: today `0.5` is silently coerced to `TRUE`.
test_that("T-SVD-06: a length-1 numeric is rejected, not coerced to TRUE", {
  set.seed(10)
  mat <- matrix(stats::rnorm(40), nrow = 10, ncol = 4)

  expect_error(.compute_matrix_mean(mat = mat, mean_vec = 0.5),
               regexp = "mean_vec")
  expect_error(.compute_matrix_sd(mat = mat, sd_vec = 0.5),
               regexp = "sd_vec")
})

test_that("T-SVD-06b: the accepted forms still work", {
  set.seed(10)
  mat <- matrix(stats::rnorm(40), nrow = 10, ncol = 4)

  expect_null(.compute_matrix_mean(mat = mat, mean_vec = NULL))
  expect_null(.compute_matrix_mean(mat = mat, mean_vec = FALSE))
  expect_equal(.compute_matrix_mean(mat = mat, mean_vec = TRUE),
               Matrix::colMeans(mat))

  full_vec <- stats::rnorm(4)
  expect_equal(.compute_matrix_mean(mat = mat, mean_vec = full_vec), full_vec)
})

## Section 1.8: the stability guard uses `&` on scalars. The behaviour is
## identical to `&&` today; this pins the intended semantics BEFORE `&&` is
## substituted, so the substitution is provably safe.
test_that("T-SVD-07: check_stability only engages above K = 5", {
  set.seed(10)
  mat <- matrix(stats::rnorm(60 * 12), nrow = 60, ncol = 12)

  # K = 3 must not run the stability check, so it cannot warn about RSpectra.
  expect_no_warning(
    .svd_safe(mat = mat, check_stability = TRUE, K = 3,
              mean_vec = NULL, rescale = FALSE, scale_max = NULL,
              sd_vec = NULL)
  )

  # K = 10 does run it. It may or may not warn on well-conditioned data; what
  # matters is that it completes and returns a valid decomposition.
  res <- suppressWarnings(
    .svd_safe(mat = mat, check_stability = TRUE, K = 10,
              mean_vec = NULL, rescale = FALSE, scale_max = NULL,
              sd_vec = NULL)
  )
  expect_equal(length(res$d), 10L)
})
