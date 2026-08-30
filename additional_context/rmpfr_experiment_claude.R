# rmpfr_experiment_claude.R
#
# Decision experiment for UNIT_TEST_PLAN.md section 2.15 (T-MPFR-01 .. T-MPFR-04,
# T-MPFR-07). These tests answer, once, whether `eSVD2` needs the `Rmpfr`
# dependency. They deliberately do NOT ship in `tests/`: they require `Rmpfr`,
# which is the dependency they exist to remove (question Q-MPFR-1, resolved
# "outside of the package").
#
# Run:   Rscript additional_context/rmpfr_experiment_claude.R
# Output is pasted verbatim into additional_context/RMPFR_REPORT.md.
#
# The shipped counterparts are T-MPFR-05, T-MPFR-06 and T-MPFR-08 in
# tests/testthat/test_tail_precision.R, none of which need `Rmpfr`.

rm(list = ls())

library(Rmpfr)

# The grid straddles the double underflow edge for the NON-log normal tail,
# which sits at roughly z = -38.5: below that, 2*pnorm(z) is exactly 0 in double
# precision. eSVD2 reaches this regime whenever one gene is strongly
# differentially expressed (see UNIT_TEST_PLAN.md section 1.1).
z_vec <- c(-5, -10, -20, -37, -38.5, -40, -100, -1000, -1e4)

# 200 bits is ~60 decimal digits, far beyond anything a double can represent.
# It is the oracle, not a candidate implementation.
prec_bits <- 200

print("===================================================================")
print("T-MPFR-01: is Rmpfr::pnorm on plain doubles the same as stats::pnorm?")
print("===================================================================")

double_log_vec <- stats::pnorm(z_vec, mean = 0, sd = 1, log.p = TRUE)
rmpfr_log_vec <- Rmpfr::pnorm(z_vec, mean = 0, sd = 1, log.p = TRUE)

print(paste0("class(Rmpfr::pnorm(<double>)) = ",
             paste0(class(rmpfr_log_vec), collapse = ", ")))
print(paste0("identical(stats, Rmpfr) = ", identical(double_log_vec,
                                                     rmpfr_log_vec)))

print("===================================================================")
print("T-MPFR-02: how accurate is the double, against a 200-bit oracle?")
print("===================================================================")

# Rmpfr::pnorm dispatches to the true MPFR routine only when its ARGUMENT is an
# mpfr object. Converting the result of a double computation would measure
# nothing.
z_mpfr <- Rmpfr::mpfr(z_vec, precBits = prec_bits)
oracle_log_vec <- Rmpfr::pnorm(z_mpfr, mean = 0, sd = 1, log.p = TRUE)

relative_error_vec <- as.numeric(abs(
  (Rmpfr::mpfr(double_log_vec, precBits = prec_bits) - oracle_log_vec) /
    oracle_log_vec
))

accuracy_df <- data.frame(z = z_vec,
                          double_log_p = double_log_vec,
                          relative_error = relative_error_vec)
print(accuracy_df)
print(paste0("max relative error = ", max(relative_error_vec)))

print("===================================================================")
print("T-MPFR-03: does Rmpfr prevent the underflow in the NON-log branch?")
print("===================================================================")

# multtest() computes `pvalue_vec` with log.p = FALSE, and that is the value
# that reaches stats::p.adjust().
stats_tail_vec <- 2 * stats::pnorm(c(-38.5, -40), mean = 0, sd = 1,
                                   log.p = FALSE)
rmpfr_tail_vec <- 2 * Rmpfr::pnorm(c(-38.5, -40), mean = 0, sd = 1,
                                   log.p = FALSE)

print(paste0("2*stats::pnorm(-38.5) = ", stats_tail_vec[1],
             " ; 2*stats::pnorm(-40) = ", stats_tail_vec[2]))
print(paste0("2*Rmpfr::pnorm(-38.5) = ", rmpfr_tail_vec[1],
             " ; 2*Rmpfr::pnorm(-40) = ", rmpfr_tail_vec[2]))
print(paste0("both underflow to exactly 0: ",
             all(stats_tail_vec == 0) && all(rmpfr_tail_vec == 0)))

print("===================================================================")
print("T-MPFR-04: can an mpfr p-value survive stats::p.adjust?")
print("===================================================================")

# This is the decisive test. Even a correctly written MPFR pipeline ends at
# p.adjust(), which is double-only.
pvalue_mpfr <- 2 * Rmpfr::pnorm(Rmpfr::mpfr(c(-40, -30, -2),
                                            precBits = prec_bits),
                                mean = 0, sd = 1, log.p = FALSE)
print("the mpfr p-values, at full precision:")
print(pvalue_mpfr)

fdr_vec <- stats::p.adjust(pvalue_mpfr, method = "BH")
print(paste0("class after p.adjust = ", paste0(class(fdr_vec), collapse = ", ")))
print("stats::p.adjust(<mpfr>, 'BH') returns:")
print(fdr_vec)
print(paste0("the 1e-350 entry came back as exactly 0: ", fdr_vec[1] == 0))

print("===================================================================")
print("T-MPFR-07: log.p is the fix, not Rmpfr")
print("===================================================================")

# UNIT_TEST_PLAN.md section 1.1: compute_pvalue() computes
# qnorm(pt(teststat, df)) without log.p, which saturates.
teststat_val <- 40
df_val <- 18

print(paste0("qnorm(pt(40, 18))                        = ",
             stats::qnorm(stats::pt(teststat_val, df = df_val))))
print(paste0("qnorm(pt(-40, 18, log.p=T), log.p=T)     = ",
             stats::qnorm(stats::pt(-teststat_val, df = df_val, log.p = TRUE),
                          log.p = TRUE)))

# And the same computation carried out entirely in MPFR, for comparison.
gaussian_teststat <- stats::qnorm(
  stats::pt(-teststat_val, df = df_val, log.p = TRUE), log.p = TRUE
)
double_log10p <- -(stats::pnorm(gaussian_teststat, mean = 0, sd = 1,
                                log.p = TRUE) / log(10) + log10(2))
oracle_log10p <- as.numeric(
  -(Rmpfr::pnorm(Rmpfr::mpfr(gaussian_teststat, precBits = prec_bits),
                 mean = 0, sd = 1, log.p = TRUE) /
      log(Rmpfr::mpfr(10, precBits = prec_bits)) +
      log10(Rmpfr::mpfr(2, precBits = prec_bits)))
)
print(paste0("log10pvalue, double = ", double_log10p))
print(paste0("log10pvalue, mpfr   = ", oracle_log10p))
print(paste0("relative difference = ",
             abs(double_log10p - oracle_log10p) / abs(oracle_log10p)))

print("===================================================================")
print("T-MPFR-05 (shipped, no Rmpfr): the dependency-free oracle")
print("===================================================================")

# Mills-ratio asymptotic expansion of the standard normal log-tail:
#   log Phi(-z) = -z^2/2 - log(z) - log(2*pi)/2
#                 + log(1 - z^-2 + 3*z^-4 - 15*z^-6 + 105*z^-8 - ...)
# This is what lets the shipped test assert double adequacy with no dependency.
.asymptotic_log_tail <- function(z_vec){
  z_abs_vec <- abs(z_vec)
  -z_abs_vec^2 / 2 - log(z_abs_vec) - log(2 * pi) / 2 +
    log1p(-z_abs_vec^(-2) + 3 * z_abs_vec^(-4) - 15 * z_abs_vec^(-6) +
            105 * z_abs_vec^(-8))
}

z_asym_vec <- c(-10, -20, -40, -100, -1000, -1e4)
asym_df <- data.frame(
  z = z_asym_vec,
  relative_error = abs(
    (stats::pnorm(z_asym_vec, mean = 0, sd = 1, log.p = TRUE) -
       .asymptotic_log_tail(z_asym_vec)) /
      stats::pnorm(z_asym_vec, mean = 0, sd = 1, log.p = TRUE)
  )
)
print(asym_df)

print("===================================================================")
print("T-MPFR-08: where information IS lost, and Rmpfr is not the cure")
print("===================================================================")

# report_results() returns pvalue = 10^(-log10pvalue). Two genes with clearly
# distinct log10pvalue become indistinguishable once exponentiated.
log10pvalue_vec <- c(400, 500, 800)
reported_vec <- 10^(-log10pvalue_vec)
print(data.frame(log10pvalue = log10pvalue_vec, reported_pvalue = reported_vec))
print(paste0("all reported p-values are identical: ",
             length(unique(reported_vec)) == 1))
print("=> the fix is to expose log10pvalue in report_results(), not more bits.")

print("===================================================================")
print("Session")
print("===================================================================")
print(R.version.string)
print(paste0("Rmpfr ", as.character(utils::packageVersion("Rmpfr"))))
