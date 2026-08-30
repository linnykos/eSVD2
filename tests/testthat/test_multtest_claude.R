# UNIT_TEST_PLAN.md section 2.10 -- T-MT-01 .. T-MT-08.
#
# `multtest()` is the paper's Type-1-error control mechanism and has zero
# coverage. Its three estimators of one quantity are each other's oracle.

test_that("T-MT-01: .multtest_locfdr recovers the null on clean N(0,1) data", {
  set.seed(10)
  teststat_vec <- stats::rnorm(1000)
  names(teststat_vec) <- paste0("gene_", seq_along(teststat_vec))

  res <- suppressWarnings(multtest(teststat_vec))

  # Establish that the primary estimator works, before testing what happens
  # when it does not.
  expect_equal(res$method, "locfdr")
  expect_equal(res$null_mean, 0, tolerance = 0.15)
  expect_equal(res$null_sd, 1, tolerance = 0.15)
})

## NEW FINDING, not in UNIT_TEST_PLAN.md and not in CRAN_READINESS.md.
##
## The three estimators are each other's oracle, and on 1000 draws from exactly
## N(0, 1) they do NOT agree:
##
##   locfdr          mean = +0.0065   sd = 0.9866    <- accurate
##   truncated_mle   mean = +0.0671   sd = 1.1079    <- sd 12% too HIGH
##   simple          mean = +0.0104   sd = 0.7901    <- sd 21% too LOW
##
## `.multtest_simple()` is the LAST fallback, and it underestimates the null
## standard deviation by about a fifth. Dividing by a null sd that is 21% too
## small inflates every z-score by about 27%, so the p-values it produces are
## substantially ANTI-CONSERVATIVE.
##
## Combine that with section 1.1 and the consequence is concrete: one strongly
## DE gene produces an `Inf`, `locfdr` errors, the `tryCatch` swallows it,
## `.multtest_truncatedGauss` also fails on the `Inf`, and the run silently
## lands on `.multtest_simple` -- whose own source comment disclaims it -- for
## EVERY gene in the dataset. This is the mechanism by which section 1.1
## changes results, quantified.
##
## The first test pins the observed biases so a change to any estimator is
## caught. The second asserts what the plan's T-MT-02 asks -- that the three
## agree -- and is expected to FAIL.

test_that("T-MT-02a: the fallback estimators are biased, in opposite directions", {
  set.seed(10)
  teststat_vec <- stats::rnorm(1000)
  names(teststat_vec) <- paste0("gene_", seq_along(teststat_vec))

  res_locfdr <- .multtest_locfdr(teststat_vec)
  res_truncated <- .multtest_truncatedGauss(teststat_vec,
                                            observed_quantile = c(0.05, 0.95))
  res_simple <- .multtest_simple(teststat_vec,
                                 observed_quantile = c(0.05, 0.95))

  # locfdr is accurate on data drawn from the null it is estimating.
  expect_equal(res_locfdr$null_sd, 1, tolerance = 0.05)

  # The truncated MLE is conservative: too wide a null means p-values that are
  # too large. Survivable.
  expect_true(res_truncated$null_sd > res_locfdr$null_sd)

  # `.multtest_simple` is ANTI-conservative: too narrow a null means p-values
  # that are too small, across every gene at once.
  expect_true(res_simple$null_sd < 0.85,
              info = paste0("simple null_sd = ", res_simple$null_sd))
  expect_true(res_simple$null_sd < res_locfdr$null_sd)
})

test_that("T-MT-02: the three estimators agree on clean N(0,1) data", {
  set.seed(10)
  teststat_vec <- stats::rnorm(1000)
  names(teststat_vec) <- paste0("gene_", seq_along(teststat_vec))

  res_locfdr <- .multtest_locfdr(teststat_vec)
  res_truncated <- .multtest_truncatedGauss(teststat_vec,
                                            observed_quantile = c(0.05, 0.95))
  res_simple <- .multtest_simple(teststat_vec,
                                 observed_quantile = c(0.05, 0.95))

  for(res in list(res_truncated, res_simple)){
    expect_equal(res$null_sd, res_locfdr$null_sd, tolerance = 0.1,
                 info = res$method)
  }
})

## Section 1.1's silent-degradation mechanism, isolated. `.multtest_locfdr()`
## catches WARNINGS as well as errors, and `locfdr` warns whenever the gene
## count is modest -- so the fallback is entered far more often than "when
## locfdr fails" suggests. At 200 genes the run silently uses `truncated_mle`,
## which on data drawn from exactly N(0, 1) returns mean 0.70 and sd 2.57.
##
## Expected to FAIL: the point is that a silently-substituted estimator gives a
## materially wrong null, and nothing in `pvalue_list` records that it happened.
test_that("T-MT-03: a silently substituted estimator still gets the null right", {
  set.seed(10)
  teststat_small_vec <- stats::rnorm(200)
  names(teststat_small_vec) <- paste0("gene_", seq_along(teststat_small_vec))

  res_small <- suppressWarnings(multtest(teststat_small_vec))

  # It silently used a different estimator than the 1000-gene case did.
  expect_equal(res_small$method, "truncated_mle")

  # And that estimator is badly wrong on data drawn from exactly N(0, 1).
  expect_equal(res_small$null_mean, 0, tolerance = 0.2,
               info = paste0("null_mean = ", res_small$null_mean))
  expect_equal(res_small$null_sd, 1, tolerance = 0.2,
               info = paste0("null_sd = ", res_small$null_sd))
})

## Section 1.1 fix step 3: a non-finite entry should be rejected at entry rather
## than propagating. Today an `Inf` reaches `mean`/`sd` inside
## `.multtest_truncatedGauss` and yields `NA`. Expected to FAIL.
test_that("T-MT-04: a non-finite teststat errors at entry", {
  set.seed(10)
  teststat_vec <- stats::rnorm(1000)
  names(teststat_vec) <- paste0("gene_", seq_along(teststat_vec))
  teststat_vec[1] <- Inf

  expect_error(multtest(teststat_vec), regexp = "finite")
})

## Section 1.7: `.multtest_truncatedGauss` bounds only `theta`; Nelder-Mead is
## free to step `sigma0` negative, and `pnorm(0, 0, -1)` is `NaN`. The test
## asserts the GUARD by calling the objective directly, since the failure has
## not been reproduced in a real run.
test_that("T-MT-05: the truncated-Gaussian fit returns a usable positive sd", {
  # `pnorm` with a negative sd is NaN, and `dnorm` with a negative sd is NaN
  # too -- so an unguarded objective returns NaN rather than +Inf and the
  # optimizer has no reason to step back.
  expect_true(is.nan(stats::pnorm(0, mean = 0, sd = -1)))
  expect_true(is.nan(stats::dnorm(0, mean = 0, sd = -1)))

  # What this test CAN check without reaching into the closure: the fitted sd
  # is finite and positive on well-behaved input. It does not exercise the
  # sigma0 <= 0 path, because `.multtest_truncatedGauss` does not expose its
  # objective -- section 1.7's guard needs the objective factored out before it
  # can be tested directly. Recorded so the gap is visible rather than implied
  # by the test's name.
  set.seed(10)
  teststat_vec <- stats::rnorm(500)
  res <- .multtest_truncatedGauss(teststat_vec,
                                  observed_quantile = c(0.05, 0.95))
  expect_true(is.finite(res$null_sd))
  expect_true(res$null_sd > 0)
})

## `optim`'s `convergence` code is ignored entirely. Expected to FAIL.
test_that("T-MT-06: .multtest_truncatedGauss checks optim's convergence code", {
  set.seed(10)
  teststat_vec <- stats::rnorm(500)

  res <- .multtest_truncatedGauss(teststat_vec,
                                  observed_quantile = c(0.05, 0.95))

  expect_true("convergence" %in% names(res))
})

test_that("T-MT-07: a degenerate middle 90% gives null_sd 0, and the caller copes", {
  # `pnorm(x, sd = 0)` is a step function, so every p-value becomes 0 or 1.
  teststat_vec <- c(rep(0, 180), seq(-5, 5, length.out = 20))
  names(teststat_vec) <- paste0("gene_", seq_along(teststat_vec))

  res_simple <- .multtest_simple(teststat_vec,
                                 observed_quantile = c(0.05, 0.95))
  expect_equal(res_simple$null_sd, 0)

  # Whatever `multtest()` does with that, it must not return NaN p-values.
  res <- suppressWarnings(multtest(teststat_vec))
  expect_true(all(!is.na(res$fdr_vec)),
              info = paste0("method = ", res$method))
})

test_that("T-MT-08: observed_quantile is respected", {
  set.seed(10)
  teststat_vec <- stats::rnorm(1000)

  res_wide <- .multtest_simple(teststat_vec,
                               observed_quantile = c(0.05, 0.95))
  res_narrow <- .multtest_simple(teststat_vec,
                                 observed_quantile = c(0.25, 0.75))

  # A strictly smaller central subset must give a strictly smaller sd.
  expect_true(res_narrow$null_sd < res_wide$null_sd)
})

test_that("T-MT-09: multtest's outputs are all named and finite", {
  set.seed(10)
  teststat_vec <- stats::rnorm(1000)
  names(teststat_vec) <- paste0("gene_", seq_along(teststat_vec))

  res <- suppressWarnings(multtest(teststat_vec))

  for(element_name in c("fdr_vec", "logpvalue_vec", "pvalue_vec")){
    expect_equal(names(res[[element_name]]), names(teststat_vec),
                 info = element_name)
    expect_true(all(is.finite(res[[element_name]])), info = element_name)
  }
  expect_true(all(res$pvalue_vec >= 0 & res$pvalue_vec <= 1))
})
