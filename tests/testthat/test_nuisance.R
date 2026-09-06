context("Test nuisance")

test_that("estimate_nuisance works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  eSVD_obj$teststat_vec <- NULL
  eSVD_obj$fit_First$nuisance_vec <- NULL

  res <- estimate_nuisance(input_obj = eSVD_obj,
                           verbose = 0)
  # plot(res$fit_First$nuisance_vec)

  expect_true("nuisance_vec" %in% names(res$fit_First))
  expect_true(all(res$fit_First$nuisance_vec > 0))
})

test_that("estimate_nuisance.default returns finite positive values when mu contains extreme values", {
  set.seed(1)
  n <- 50; p <- 10

  # Construct a mean_mat where some entries are very large (exp overflow -> Inf),
  # which can cause gamma_rate to return NaN without the validity guard.
  mean_mat <- matrix(abs(rnorm(n * p, mean = 5, sd = 2)), nrow = n, ncol = p)
  mean_mat[, 1] <- Inf   # extreme column: triggers NaN in gamma_rate
  mean_mat[, 2] <- 1e300 # near-overflow column

  library_mat <- matrix(rep(1, n * p), nrow = n, ncol = p)
  dat <- matrix(rpois(n * p, lambda = 3), nrow = n, ncol = p)
  storage.mode(dat) <- "double"
  colnames(dat) <- paste0("gene", seq_len(p))

  # The `Inf` column cannot be estimated on either route, so it falls to
  # `min_val` and the function warns about it; that warning is the contract.
  expect_warning(
    res <- estimate_nuisance.default(
      input_obj   = dat,
      mean_mat    = mean_mat,
      library_mat = library_mat,
      verbose     = 0
    ),
    "nuisance estimation failed"
  )

  expect_true(length(res) == p)
  expect_true(all(is.finite(res)))
  expect_true(all(res > 0))
})
