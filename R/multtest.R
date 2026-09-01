#' Perform multiple-testing adjustment using Efron's empirical null
#'
#' Estimates the empirical null distribution \eqn{N(\mu_0, \sigma_0^2)} of the
#' Gaussianized test statistics and converts every statistic to a two-sided
#' p-value under it, followed by a Benjamini-Hochberg adjustment.
#'
#' Three estimators of the null are tried in order, each only when the
#' previous one fails: \code{locfdr::locfdr} (\code{method = "locfdr"}), the
#' truncated-Gaussian maximum-likelihood estimator of Efron (2007, Equation
#' 4.12) fitted to the statistics between the two \code{observed_quantile}
#' quantiles (\code{method = "truncated_mle"}), and a moment-matching
#' estimator on the same central subset (\code{method = "simple"}). A warning
#' is raised whenever the chain falls below \code{locfdr}, so a degraded null
#' is never silent; the estimator that ran is returned as \code{method}.
#'
#' @param teststat_vec       Named vector of finite Gaussian-like test
#'                           statistics, one per gene.
#' @param observed_quantile  Two numbers in \code{(0, 1)}, the lower and upper
#'                           quantiles of \code{teststat_vec} that delimit the
#'                           central subset used by the two fallback
#'                           estimators.
#'
#' @returns A list with \code{fdr_vec} (BH-adjusted p-values), \code{logpvalue_vec}
#' (\eqn{-\log_{10}} of the two-sided p-values), \code{method} (which estimator
#' produced the null), \code{null_mean}, \code{null_sd}, and \code{pvalue_vec}
#' (the two-sided p-values). All vectors carry \code{names(teststat_vec)}.
#' @export
multtest <- function(teststat_vec,
                     observed_quantile = c(0.05, 0.95)){
  if(!is.numeric(teststat_vec) || !all(is.finite(teststat_vec))){
    stop("`teststat_vec` must be a numeric vector of finite values; found ",
         sum(!is.finite(teststat_vec)), " non-finite entries. ",
         "A non-finite statistic upstream usually means `qnorm(pt())` ",
         "saturated; see `compute_pvalue`")
  }
  stopifnot(length(observed_quantile) == 2,
            all(observed_quantile > 0), all(observed_quantile < 1),
            observed_quantile[1] < observed_quantile[2])

  res <- .multtest_locfdr(teststat_vec)
  if(!.null_estimate_is_usable(res)){
    res <- .multtest_truncatedGauss(teststat_vec,
                                    observed_quantile = observed_quantile)
  }
  if(!.null_estimate_is_usable(res)){
    res <- .multtest_simple(teststat_vec,
                            observed_quantile = observed_quantile)
  }

  # A silently degraded null is the failure mode of CRAN_READINESS.md 1.1:
  # one bad gene used to move the whole dataset onto the crudest estimator
  # with nothing telling the user. The warning is the fix.
  if(res$method != "locfdr"){
    warning("locfdr could not estimate the empirical null; using the `",
            res$method, "` estimator instead (null_mean = ",
            signif(res$null_mean, 4), ", null_sd = ", signif(res$null_sd, 4),
            "). This usually means too few genes for locfdr's spline fit")
  }
  if(!is.finite(res$null_mean) || !is.finite(res$null_sd)){
    stop("no estimator could fit an empirical null to `teststat_vec` (",
         length(teststat_vec), " statistics; the central ",
         100 * diff(observed_quantile), "% window holds ",
         sum(teststat_vec >= stats::quantile(teststat_vec, observed_quantile[1]) &
               teststat_vec <= stats::quantile(teststat_vec, observed_quantile[2])),
         "). Too few genes; the empirical null needs at least a few dozen")
  }
  if(!.null_estimate_is_usable(res)){
    warning("the empirical null is degenerate (null_sd = ", res$null_sd,
            "): the central ", 100 * diff(observed_quantile),
            "% of `teststat_vec` has no spread, so every p-value is 0 or 1")
  }

  null_mean <- res$null_mean
  null_sd <- res$null_sd
  method <- res$method

  logpvalue_vec <- .two_sided_log_pvalue(teststat_vec = teststat_vec,
                                         null_mean = null_mean,
                                         null_sd = null_sd)
  pvalue_vec <- exp(logpvalue_vec)
  names(pvalue_vec) <- names(teststat_vec)
  fdr_vec <- stats::p.adjust(pvalue_vec, method = "BH")
  names(fdr_vec) <- names(teststat_vec)
  logpvalue_vec <- -logpvalue_vec / log(10)
  names(logpvalue_vec) <- names(teststat_vec)

  list(fdr_vec = fdr_vec,
       logpvalue_vec = logpvalue_vec,
       method = method,
       null_mean = null_mean,
       null_sd = null_sd,
       pvalue_vec = pvalue_vec)
}

####################

#' Two-sided log p-values under a Gaussian null
#'
#' Computed in log space so a statistic 40 null standard deviations out gives
#' a finite \code{log(p)} rather than an underflowed \code{0}. The result is
#' capped at \code{log(1)} so that a statistic sitting exactly on
#' \code{null_mean} reports \code{p = 1} even under a degenerate
#' \code{null_sd} of zero.
#'
#' @param teststat_vec  Numeric vector.
#' @param null_mean     Scalar.
#' @param null_sd       Non-negative scalar; \code{0} gives a step function.
#'
#' @returns Vector of natural-log two-sided p-values.
#' @noRd
.two_sided_log_pvalue <- function(teststat_vec, null_mean, null_sd){
  distance_vec <- abs(teststat_vec - null_mean)
  # 0/0 at a degenerate null would be NaN; a statistic on the null mean has
  # p = 1 regardless of the spread.
  z_vec <- ifelse(distance_vec == 0, 0, distance_vec / null_sd)
  logp_vec <- stats::pnorm(-z_vec, log.p = TRUE) + log(2)
  pmin(logp_vec, 0)
}

#' Is an estimated null usable?
#'
#' @param res List with \code{null_mean}, \code{null_sd} and optionally
#'            \code{convergence} (an \code{optim} code; non-zero means the
#'            fit did not converge).
#'
#' @returns Logical scalar.
#' @noRd
.null_estimate_is_usable <- function(res){
  converged <- is.null(res$convergence) || isTRUE(res$convergence == 0)
  converged &&
    length(res$null_mean) == 1 && is.finite(res$null_mean) &&
    length(res$null_sd) == 1 && is.finite(res$null_sd) && res$null_sd > 0
}

.multtest_locfdr <- function(teststat_vec){
  # Warnings are caught as well as errors: locfdr warns when its spline
  # misfits, typically at a few hundred genes, and a misfit null is not one
  # to use. The cost is that the fallback is entered more often than "when
  # locfdr errors", which is why `multtest()` warns whenever it happens.
  res <- tryCatch({
    locfdr_res <- locfdr::locfdr(teststat_vec, plot = 0)
    c(locfdr_res$fp0["mlest", "delta"],
      locfdr_res$fp0["mlest", "sigma"])
  }, warning = function(e){
    rep(NA, 2)
  }, error = function(e){
    rep(NA, 2)
  })

  list(method = "locfdr",
       null_mean = unname(res[1]),
       null_sd = unname(res[2]))
}

#' Mean and standard deviation of a standard normal truncated to \code{[a, b]}
#'
#' @param a Lower truncation point.
#' @param b Upper truncation point, \code{a < b}.
#'
#' @returns List with \code{mean} and \code{sd}.
#' @noRd
.truncated_normal_moments <- function(a, b){
  mass <- stats::pnorm(b) - stats::pnorm(a)
  trunc_mean <- (stats::dnorm(a) - stats::dnorm(b)) / mass
  trunc_var <- 1 + (a * stats::dnorm(a) - b * stats::dnorm(b)) / mass -
    trunc_mean^2

  list(mean = trunc_mean, sd = sqrt(trunc_var))
}

# observed_quantile is two numbers, a lower and and a upper
# see equation 4.8 onwards of "SIZE, POWER AND FALSE DISCOVERY RATES" by Efron
# TODO: If there is an obvious break in the values, hard-set the observed quantiles. We can do this by mixture-modeling
#
# Efron's Equation 4.12 has three parameters, (delta0, sigma0, p0), and the
# probability that a statistic lands in the central window is
# theta = p0 * H(delta0, sigma0), with H the null mass inside the window.
# An earlier version of this function optimized theta directly as a free
# fourth quantity, which drops the constraint p0 <= 1. Without it sigma0 is
# identified only by the shape of the truncated density, which for a window
# holding 90% of the mass is nearly flat towards sigma0 -> Inf: on 200 draws
# from exactly N(0, 1) that version returned null_sd = 2.57. Keeping p0 <= 1
# restores the estimator Efron describes, and with it the precision.
.multtest_truncatedGauss <- function(teststat_vec,
                                     observed_quantile){
  z_vec <- teststat_vec

  N <- length(z_vec)
  tmp <- stats::quantile(z_vec, probs = observed_quantile)
  lb <- tmp[1]; ub <- tmp[2]
  idx <- intersect(which(z_vec >= lb), which(z_vec <= ub))
  N0 <- length(idx)
  z0_vec <- z_vec[idx]

  fn <- function(param_vec){
    delta0 <- param_vec[1]
    sigma0 <- exp(param_vec[2])
    p0 <- param_vec[3]
    denom <- stats::pnorm(ub, mean = delta0, sd = sigma0) -
      stats::pnorm(lb, mean = delta0, sd = sigma0)
    theta <- p0 * denom
    # `optim` cannot take Inf, so return a finite penalty off the domain.
    if(!is.finite(denom) || denom <= 0 || theta <= 0 || theta >= 1) return(1e100)

    # full log-likelihood, Equation 4.12
    loglik <- N0*log(theta) + (N-N0)*log(1-theta) +
      sum(stats::dnorm(z0_vec, mean = delta0, sd = sigma0, log = TRUE)) -
      N0*log(denom)
    -loglik
  }

  # Initialize from the moment-matched estimate, which is already close.
  init <- .multtest_simple(teststat_vec = teststat_vec,
                           observed_quantile = observed_quantile)
  init_mass <- stats::pnorm(ub, mean = init$null_mean, sd = init$null_sd) -
    stats::pnorm(lb, mean = init$null_mean, sd = init$null_sd)
  init_p0 <- min(0.999, (N0 / N) / init_mass)

  res <- tryCatch({
    if(!is.finite(init$null_sd) || init$null_sd <= 0 || !is.finite(init_p0)){
      stop("degenerate initialization")
    }
    optim_res <- stats::optim(par = c(init$null_mean, log(init$null_sd), init_p0),
                              fn = fn,
                              method = "L-BFGS-B",
                              lower = c(-Inf, -Inf, 1e-4),
                              upper = c(Inf, Inf, 1),
                              control = list(maxit = 2000))
    c(optim_res$par[1], exp(optim_res$par[2]), optim_res$convergence)
  }, error = function(e){
    c(NA, NA, NA)
  })

  list(method = "truncated_mle",
       convergence = unname(res[3]),
       null_mean = unname(res[1]),
       null_sd = unname(res[2]))
}

# A deliberately simple estimator, used when everything else has failed.
# The sample between the two observed quantiles is a truncated sample, so its
# raw mean and standard deviation are biased -- by 21% for the sd at the
# default 5%/95% window. Matching its moments to those of a normal truncated
# at the same null quantiles removes the bias exactly when the window holds
# only null statistics, and leaves a mild upward bias otherwise.
.multtest_simple <- function(teststat_vec,
                             observed_quantile){
  tmp <- stats::quantile(teststat_vec, probs = observed_quantile)
  lb <- tmp[1]; ub <- tmp[2]
  idx <- intersect(which(teststat_vec >= lb), which(teststat_vec <= ub))
  vec <- teststat_vec[idx]

  moments <- .truncated_normal_moments(a = stats::qnorm(observed_quantile[1]),
                                       b = stats::qnorm(observed_quantile[2]))
  null_sd <- stats::sd(vec) / moments$sd
  null_mean <- mean(vec) - null_sd * moments$mean

  list(method = "simple",
       null_mean = unname(null_mean),
       null_sd = unname(null_sd))
}
