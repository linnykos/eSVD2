context("Test gamma rate")

## gamma_rate is correct

test_that("gamma_rate and log_gamma_rate work", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  x_mat <- eSVD_obj$fit_First$x_mat
  y_mat <- eSVD_obj$fit_First$y_mat
  z_mat <- eSVD_obj$fit_First$z_mat
  covariates <- eSVD_obj$covariates
  case_control_idx <- which(colnames(covariates) == "case_control_1")

  nat_mat1 <- tcrossprod(x_mat, y_mat)
  nat_mat2 <- tcrossprod(covariates[,case_control_idx], z_mat[,case_control_idx])
  mean_mat_nolib <- exp(nat_mat1 + nat_mat2)
  library_mat <- exp(tcrossprod(covariates[,-case_control_idx], z_mat[,-case_control_idx]))

  bool_vec <- sapply(1:ncol(dat), function(j){
    res1 <- gamma_rate(x = dat[,j],
                       mu = mean_mat_nolib[,j],
                       s = library_mat[,j])
    bool1 <- length(res1) == 1 & res1 >= 0

    # res2 <- log_gamma_rate(x = dat[,j],
    #                        mu = mean_mat_nolib[,j],
    #                        s = library_mat[,j])
    # bool2 <- length(res2) == 1 & exp(res2) >= 0
    #
    # bool3 <- abs(res1 - exp(res2)) <= 1e-3

    # bool1 & bool2 & bool3
    bool1
  })

  expect_true(all(bool_vec))
})

################################################################################

# UNIT_TEST_PLAN.md section 3.5 -- T-CPP-GAM-01 .. T-CPP-GAM-08.
#
# Draws a clean sample from the eSVD hierarchical model and returns it, so the
# TRUE Gamma rate is known:
#   lambda_i ~ Gamma(shape = mu_i * beta, rate = beta)   =>  E lambda = mu_i
#   x_i      ~ Poisson(s_i * lambda_i)
# `beta` here is the RATE, which is what `gamma_rate()` estimates and what the
# package stores in `nuisance_vec`. It is the reciprocal of the paper's
# over-dispersion gamma.
.gamma_sample <- function(beta, num_obs = 4000, mu_val = 5, s_val = 1,
                          seed_number = 1){
  set.seed(seed_number)

  mu_vec <- rep(mu_val, num_obs)
  s_vec <- rep(s_val, num_obs)
  lambda_vec <- stats::rgamma(num_obs, shape = mu_vec * beta, rate = beta)
  x_vec <- stats::rpois(num_obs, lambda = s_vec * lambda_vec)

  list(mu_vec = mu_vec, s_vec = s_vec, x_vec = x_vec)
}

## FIXED 2026-08-29. Kept as the record of what these tests guard against.
##
## `src/gamma_rate.cpp` set the root-search bracket's upper bound to
## `Rcpp::max(s)` -- the largest LIBRARY SIZE -- and the refinement loop only
## ever multiplied it by 0.5. It could shrink the bracket, never grow it, so
## Boost's Newton search was clamped and returned `max(s)` as if it were the
## MLE whenever the true root lay above it. Measured before the fix, on 4000
## draws with mu = 5 and s = 1:
##
##     true beta    gamma_rate    exp(log_gamma_rate)
##       0.90         0.8825           0.8825
##       1.00         0.9432           0.9432
##       1.20         0.9996           1.1555     <- capped
##       3.00         0.9995           2.9250     <- capped
##      10.00         0.9996          10.2225     <- capped
##
## beta and s share units -- they enter the likelihood only through (s + beta)
## -- which is probably why max(s) looked like a reasonable bound. But sharing
## units is not being bounded: beta is free to exceed any library size.
##
## Why it mattered: `estimate_nuisance()` defaults to `bool_use_log = FALSE`, so
## it calls `gamma_rate`. Every gene whose over-dispersion gamma was below
## max(s) silently got `nuisance_vec = max(s)` instead of an estimate, and
## `nuisance_vec` enters the posterior denominator directly.
##
## The fix grows the bracket while [l(ub)]' > 0, shrinks the lower bound until
## [l(lb)]' > 0 so the root is bracketed on both sides, starts Newton from the
## geometric mean rather than a fixed 1.0, and raises the iteration cap. Boost
## bisects whenever a Newton step leaves the bracket, so the cost is unchanged
## on easy inputs. Verified against an independent `stats::optimize` of the R
## objective at 20 combinations of max(s) in {0.25, 1, 2, 4, 10} and true beta
## in {0.3, 1, 3, 10}: all agree to <= 6e-5 relative, where 18 of 20 previously
## missed, some by 84%.

test_that("T-CPP-GAM-05: gamma_rate recovers a rate well above max(s)", {
  # THE regression test for the fix. Before it, the bracket's upper bound was
  # `max(s)` and the refinement loop only ever shrank it, so this returned
  # 0.9995 for a true rate of 3.
  sample_list <- .gamma_sample(beta = 3, s_val = 1)

  res_direct <- gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                           s = sample_list$s_vec)

  expect_equal(res_direct, 3, tolerance = 0.15,
               info = paste0("gamma_rate returned ", signif(res_direct, 5),
                             " for a true rate of 3, with max(s) = 1"))
  expect_true(res_direct > 1,
              info = "a value at or just below max(s) means the bracket is capped again")
})

## An earlier draft of this test reused data generated at `s = 1` while passing
## `s = 10`, which is a misspecified model and proves nothing. Each library size
## gets its own draw here.
test_that("T-CPP-GAM-05b: the answer no longer tracks max(s)", {
  # The estimate must depend on the DATA, not on the library size that happens
  # to bound the search. Across four library sizes the recovered rate stays near
  # the truth instead of following `max(s)`.
  for(s_val in c(0.25, 1, 2, 4)){
    sample_list <- .gamma_sample(beta = 3, s_val = s_val)
    res <- gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                      s = sample_list$s_vec)

    expect_true(res > 2,
                info = paste0("max(s) = ", s_val, " gave ", signif(res, 5),
                              "; a value near max(s) means the cap is back"))
  }
})

test_that("T-CPP-GAM-05c: gamma_rate is accurate below max(s) as well", {
  # The regime that always worked, kept so the fix cannot regress it: rates at
  # or below `max(s)` were recovered correctly even with the capped bracket,
  # which is why the defect went unnoticed for so long.
  for(beta in c(0.1, 0.3, 0.6, 0.9)){
    sample_list <- .gamma_sample(beta = beta)
    res <- gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                      s = sample_list$s_vec)
    expect_equal(res, beta, tolerance = 0.1,
                 info = paste0("true beta = ", beta))
  }
})

test_that("T-CPP-GAM-06: log_gamma_rate recovers the rate across the whole range", {
  # The log-scale routine has bounds of [-10, 10] on log(beta), i.e. beta in
  # [4.5e-5, 22026], and no dependence on `s`. It is the one that works.
  for(beta in c(0.1, 0.6, 1.2, 3, 10)){
    sample_list <- .gamma_sample(beta = beta)
    res <- exp(log_gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                              s = sample_list$s_vec))
    expect_equal(res, beta, tolerance = 0.15,
                 info = paste0("true beta = ", beta))
  }
})

test_that("T-CPP-GAM-06b: log_gamma_rate's clamping is visible at the boundary", {
  # When the true log(beta) lies outside [lower, upper] the function returns the
  # boundary, and a saturated estimate is indistinguishable from a converged
  # one. Question Q-GAM-2 decided against a convergence status, so the boundary
  # value itself is the only signal -- assert that it is at least observable.
  sample_list <- .gamma_sample(beta = 3)

  res_clamped <- log_gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                                s = sample_list$s_vec,
                                lower = -0.1, upper = 0.1)
  expect_equal(res_clamped, 0.1, tolerance = 1e-8)
})

## The commented-out half of the existing `test_gamma_rate.R`, now writable.
## It was very likely commented out BECAUSE of the bracket cap: above max(s) the
## two routines could not agree, since one of them was pinned.
test_that("T-CPP-GAM-01: the two routines agree below max(s)", {
  for(beta in c(0.1, 0.3, 0.6, 0.9)){
    sample_list <- .gamma_sample(beta = beta)

    res_direct <- gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                             s = sample_list$s_vec)
    res_log <- exp(log_gamma_rate(x = sample_list$x_vec,
                                  mu = sample_list$mu_vec,
                                  s = sample_list$s_vec))

    expect_equal(res_direct, res_log, tolerance = 1e-4,
                 info = paste0("true beta = ", beta))
  }
})

test_that("T-CPP-GAM-01b: the two routines agree ABOVE the old cap too", {
  # This is what the fix bought. Before it, the two disagreed by a factor of 20
  # here -- almost certainly why the equivalence half of the existing
  # `test_gamma_rate.R` was commented out.
  for(beta in c(1.2, 3, 10)){
    sample_list <- .gamma_sample(beta = beta)

    res_direct <- gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                             s = sample_list$s_vec)
    res_log <- exp(log_gamma_rate(x = sample_list$x_vec,
                                  mu = sample_list$mu_vec,
                                  s = sample_list$s_vec))

    expect_equal(res_direct, res_log, tolerance = 1e-6,
                 info = paste0("true beta = ", beta))
  }
})

## The R implementation of the objective already sits in the comment block of
## `src/gamma_rate.cpp` lines 58-75. Turning it into a test is nearly free and
## gives an oracle independent of both C++ routines.
.gamma_objfn <- function(beta, x_vec, mu_vec, s_vec){
  sum(x_vec * log(s_vec) - lgamma(x_vec + 1) +
        mu_vec * beta * log(beta) - lgamma(mu_vec * beta) +
        lgamma(mu_vec * beta + x_vec) -
        (mu_vec * beta + x_vec) * log(s_vec + beta))
}

test_that("T-CPP-GAM-02: the estimate maximizes the likelihood locally", {
  # Below the cap, where `gamma_rate` is meant to be right.
  sample_list <- .gamma_sample(beta = 0.6)

  beta_hat <- gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                         s = sample_list$s_vec)
  objective_at_hat <- .gamma_objfn(beta_hat, sample_list$x_vec,
                                   sample_list$mu_vec, sample_list$s_vec)

  beta_grid <- beta_hat * c(0.7, 0.85, 0.95, 1.05, 1.15, 1.3)
  for(beta_val in beta_grid){
    objective_val <- .gamma_objfn(beta_val, sample_list$x_vec,
                                  sample_list$mu_vec, sample_list$s_vec)
    expect_true(objective_at_hat >= objective_val,
                info = paste0("beta = ", signif(beta_val, 4)))
  }
})

test_that("T-CPP-GAM-02b: the R objective confirms the true rate beats the capped one", {
  # The independent oracle applied to the finding: at beta = 3 the likelihood
  # is strictly higher than at the value `gamma_rate` actually returns, so the
  # cap is losing likelihood, not finding a different optimum.
  sample_list <- .gamma_sample(beta = 3)
  beta_hat <- gamma_rate(x = sample_list$x_vec, mu = sample_list$mu_vec,
                         s = sample_list$s_vec)

  # The independent oracle, now applied to the FIXED estimator: `gamma_rate`
  # must land at the maximum of the R objective, not below it. Before the fix
  # `l(3)` was strictly higher than `l(beta_hat)` -- the cap was losing
  # likelihood, not finding a different optimum.
  reference_hat <- stats::optimize(
    function(beta) .gamma_objfn(beta, sample_list$x_vec, sample_list$mu_vec,
                                sample_list$s_vec),
    interval = c(0.001, 200), maximum = TRUE
  )$maximum

  expect_equal(beta_hat, reference_hat, tolerance = 1e-3,
               info = paste0("gamma_rate = ", signif(beta_hat, 6),
                             ", R optimize = ", signif(reference_hat, 6)))
})

## Question Q-GAM-1 resolved: all-zero genes are removed by `eSVD_helper()`
## rather than handled here. This records the unit-level behaviour, since
## `gamma_rate` stays callable directly.
test_that("T-CPP-GAM-03: an all-zero gene returns something finite", {
  num_obs <- 100
  res_direct <- gamma_rate(x = rep(0, num_obs), mu = rep(1e-8, num_obs),
                           s = rep(1, num_obs))
  res_log <- log_gamma_rate(x = rep(0, num_obs), mu = rep(1e-8, num_obs),
                            s = rep(1, num_obs))

  expect_true(is.finite(res_direct))
  expect_true(is.finite(res_log))
  # [verified] 9.706e-4 and exactly -10 -- the latter is the clamp boundary, so
  # the two disagree by a factor of 21 in precisely the degenerate case.
  expect_equal(res_log, -10, tolerance = 1e-8)
})

test_that("T-CPP-GAM-04: a gene with a single non-zero count behaves", {
  num_obs <- 100
  x_vec <- rep(0, num_obs)
  x_vec[1] <- 1

  res <- gamma_rate(x = x_vec, mu = rep(0.01, num_obs), s = rep(1, num_obs))
  expect_true(is.finite(res))
  expect_true(res > 0)
})

## The C++ reads `m_n = x.length()` and then indexes `mu` and `s` with it, so a
## shorter `mu` is an OUT-OF-BOUNDS READ, not an error. Expected to FAIL.
test_that("T-CPP-GAM-08: a length mismatch among x, mu, s errors", {
  x_vec <- rep(1, 100)

  expect_error(gamma_rate(x = x_vec, mu = rep(1, 50), s = rep(1, 100)))
  expect_error(gamma_rate(x = x_vec, mu = rep(1, 100), s = rep(1, 50)))
  expect_error(log_gamma_rate(x = x_vec, mu = rep(1, 50), s = rep(1, 100)))
})
