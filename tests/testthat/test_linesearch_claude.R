# The Newton line search used to raise its failure as an Rcpp::warning() from
# inside a C++ frame holding Eigen objects (CRAN_READINESS.md 4.2). It now
# returns a flag that `opt_x`/`opt_yz` count and attach as an attribute, and
# `opt_esvd.default` turns a non-zero count into one R-level warning.

test_that("T-CN-01: opt_x and opt_yz attach a line-search failure count", {
  point <- .feasible_point("poisson", num_cells = 20, num_genes = 8)
  loader <- data_loader(point$dat)
  family <- esvd_family("poisson")
  k <- ncol(point$x_mat)

  xc_mat <- opt_x(XC_init = point$x_mat, YZ = point$y_mat, k = k,
                  loader = loader, family = family, s = point$s_vec,
                  gamma = point$gamma_vec, l2penx = 0.1, verbose = 0)
  num_failed_x <- attr(xc_mat, "num_linesearch_failed", exact = TRUE)
  expect_true(is.numeric(num_failed_x) && length(num_failed_x) == 1)
  expect_true(num_failed_x >= 0 && num_failed_x <= nrow(point$dat))

  yz_mat <- opt_yz(YZ_init = point$y_mat, XC = xc_mat, k = k,
                   fixed_cols = integer(0), loader = loader, family = family,
                   s = point$s_vec, gamma = point$gamma_vec, l2peny = 0.1,
                   l2penz = 0.1, verbose = 0)
  num_failed_yz <- attr(yz_mat, "num_linesearch_failed", exact = TRUE)
  expect_true(is.numeric(num_failed_yz) && length(num_failed_yz) == 1)
  expect_true(num_failed_yz >= 0 && num_failed_yz <= ncol(point$dat))

  # `inplace = FALSE` returns a fresh matrix that carries the count too.
  xc_copy <- opt_x(XC_init = point$x_mat, YZ = point$y_mat, k = k,
                   loader = loader, family = family, s = point$s_vec,
                   gamma = point$gamma_vec, l2penx = 0.1, verbose = 0,
                   inplace = FALSE)
  expect_false(is.null(attr(xc_copy, "num_linesearch_failed", exact = TRUE)))
})

test_that("T-CN-02: a well-posed opt_esvd run is silent and its fit carries no C++ attribute", {
  point <- .feasible_point("poisson", num_cells = 20, num_genes = 8)

  expect_no_warning(
    res <- opt_esvd(input_obj = point$dat, x_init = point$x_mat,
                    y_init = point$y_mat, family = "poisson", max_iter = 3,
                    nuisance_vec = point$gamma_vec, verbose = 0)
  )
  expect_null(attr(res$x_mat, "num_linesearch_failed", exact = TRUE))
  expect_null(attr(res$y_mat, "num_linesearch_failed", exact = TRUE))
})

test_that("T-CN-03: .attr_or_zero reads the count and defaults to 0", {
  mat <- matrix(0, 2, 2)
  expect_equal(.attr_or_zero(mat, "num_linesearch_failed"), 0)
  attr(mat, "num_linesearch_failed") <- 3L
  expect_equal(.attr_or_zero(mat, "num_linesearch_failed"), 3)
})
