# UNIT_TEST_PLAN.md section 2.15 -- the SHIPPED subset of T-MPFR-01 .. T-MPFR-08.
#
# T-MPFR-01..04 and the MPFR half of T-MPFR-07 do not live here: they require
# `Rmpfr`, which is the dependency they exist to remove (question Q-MPFR-1,
# resolved "outside of the package"). They run once as
# `additional_context/rmpfr_experiment_claude.R`, and their output is recorded
# in `additional_context/RMPFR_REPORT.md`.
#
# What remains needs no `Rmpfr` and pins the same conclusion permanently:
# double precision is adequate for every regime this package can reach, the
# `log.p = TRUE` fix of section 1.1 is what is actually needed, and the one
# place information is genuinely lost is `report_results()`.

# Mills-ratio asymptotic expansion of the standard normal log-tail:
#   log Phi(-z) = -z^2/2 - log(z) - log(2*pi)/2
#                 + log(1 - z^-2 + 3*z^-4 - 15*z^-6 + 105*z^-8 - ...)
# The series is asymptotic, so it is loose near z = -10 and tightens rapidly.
.asymptotic_log_tail <- function(z_vec){
  z_abs_vec <- abs(z_vec)
  -z_abs_vec^2 / 2 - log(z_abs_vec) - log(2 * pi) / 2 +
    log1p(-z_abs_vec^(-2) + 3 * z_abs_vec^(-4) - 15 * z_abs_vec^(-6) +
            105 * z_abs_vec^(-8))
}

test_that("T-MPFR-05: stats::pnorm(log.p = TRUE) matches the analytic tail expansion", {
  z_vec <- c(-20, -40, -100, -1000, -1e4)

  double_log_vec <- stats::pnorm(z_vec, mean = 0, sd = 1, log.p = TRUE)
  asymptotic_vec <- .asymptotic_log_tail(z_vec)

  relative_error_vec <- abs((double_log_vec - asymptotic_vec) / double_log_vec)

  # This is the dependency-free replacement for the 200-bit MPFR oracle of
  # T-MPFR-02. If a double ever stopped being accurate in the log tail, this
  # fails -- with nothing in Suggests and in microseconds.
  expect_true(all(relative_error_vec < 1e-12),
              info = paste0("max relative error = ", max(relative_error_vec)))

  # The log tail is strictly decreasing and finite across the whole range,
  # which is the property `compute_pvalue` actually relies on.
  expect_true(all(is.finite(double_log_vec)))
  expect_true(all(diff(double_log_vec) < 0))
})

test_that("T-MPFR-05b: the non-log branch underflows, and the log branch does not", {
  # This is the mechanism behind section 1.1, stated as a test so the reason
  # for the `log.p = TRUE` fix is recorded rather than remembered.
  expect_equal(2 * stats::pnorm(-40, mean = 0, sd = 1, log.p = FALSE), 0)
  expect_true(is.finite(stats::pnorm(-40, mean = 0, sd = 1, log.p = TRUE)))

  # The underflow edge sits just below |z| = 38.5 for the two-sided p-value.
  expect_true(2 * stats::pnorm(-38, mean = 0, sd = 1, log.p = FALSE) > 0)
  expect_equal(2 * stats::pnorm(-39, mean = 0, sd = 1, log.p = FALSE), 0)
})

test_that("T-MPFR-06: the test statistics an actual fit produces stay in the safe range", {
  esvd_obj <- .small_esvd_obj()
  pvalue_list <- esvd_obj$pvalue_list

  standardized_vec <- abs(pvalue_list$gaussian_teststat - pvalue_list$null_mean) /
    pvalue_list$null_sd

  # The precision argument above is conditional on the reachable range. This is
  # what makes it stay true: if a future change starts producing |z| = 1e6, this
  # fails and the argument gets re-examined rather than silently expiring.
  expect_true(all(is.finite(standardized_vec)))
  expect_true(max(standardized_vec) < 200,
              info = paste0("max |z| = ", max(standardized_vec)))
})

## Section 1.1: `compute_pvalue()` computes `qnorm(pt(teststat, df))` with no
## `log.p`, which saturates to +/-Inf for a strongly differentially expressed
## gene. This asserts the FIXED form is finite, and that the current form is
## not -- so the test documents the defect and will keep documenting the fix.
test_that("T-MPFR-07: log.p is the fix; the unlogged form saturates", {
  teststat_val <- 40
  df_val <- 18

  # What the code does today.
  expect_equal(stats::qnorm(stats::pt(teststat_val, df = df_val)), Inf)

  # What it should do.
  fixed_val <- stats::qnorm(stats::pt(-teststat_val, df = df_val, log.p = TRUE),
                            log.p = TRUE)
  expect_true(is.finite(fixed_val))
  expect_equal(fixed_val, -8.915293, tolerance = 1e-5)

  # And the resulting log10 p-value is finite and large, rather than Inf.
  log10pvalue <- -(stats::pnorm(fixed_val, mean = 0, sd = 1, log.p = TRUE) /
                     log(10) + log10(2))
  expect_true(is.finite(log10pvalue))
  expect_true(log10pvalue > 15)
})

## The one place information is genuinely lost -- and more precision is not the
## cure, because the column is a double either way. Expected to FAIL until
## `report_results()` exposes `log10pvalue`.
test_that("T-MPFR-08: report_results can distinguish two extremely significant genes", {
  esvd_obj <- .small_esvd_obj()

  # Two genes 400 orders of magnitude apart in significance.
  esvd_obj$pvalue_list$log10pvalue[1] <- 400
  esvd_obj$pvalue_list$log10pvalue[2] <- 800

  res <- report_results(esvd_obj)

  # `10^(-400)` and `10^(-800)` are both exactly 0 in double precision, so the
  # reported p-values tie and the genes cannot be ranked -- while
  # `pvalue_list$log10pvalue` has held the distinction the whole time.
  expect_false(res$pvalue[1] == res$pvalue[2])
})

test_that("T-MPFR-08b: log10pvalue itself does distinguish them", {
  # The companion to T-MPFR-08: the information exists, it is simply not
  # exposed. This half passes today and is what makes the one-line fix obvious.
  log10pvalue_vec <- c(400, 800)
  expect_false(log10pvalue_vec[1] == log10pvalue_vec[2])
  expect_equal(10^(-log10pvalue_vec), c(0, 0))
})
