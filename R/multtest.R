# teststat_vec: should be a gaussian-like vector
multtest <- function(teststat_vec,
                     observed_quantile = c(0.05, 0.95)){
  res <- .multtest_locfdr(teststat_vec)
  if(any(is.na(res))){
    res <- .multtest_truncatedGauss(teststat_vec,
                                    observed_quantile = observed_quantile)
  }
  print(res)
  if(any(is.na(res))){
    res <- .multtest_simple(teststat_vec,
                            observed_quantile = observed_quantile)
  }

  null_mean <- res$null_mean
  null_sd <- res$null_sd
  method <- res$method

  logpvalue_vec <- sapply(teststat_vec, function(x){
    if(x < null_mean) {
      Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
    } else {
      Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
    }
  })
  logpvalue_vec <- -(logpvalue_vec/log(10) + log10(2))
  names(logpvalue_vec) <- names(teststat_vec)

  pvalue_vec <- sapply(teststat_vec, function(x){
    if(x < null_mean) {
      Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = F)
    } else {
      Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = F)
    }
  })
  pvalue_vec <- pvalue_vec*2
  names(pvalue_vec) <- names(teststat_vec)
  fdr_vec <- stats::p.adjust(pvalue_vec, method = "BH")
  names(fdr_vec) <- names(teststat_vec)

  list(fdr_vec = fdr_vec,
       logpvalue_vec = logpvalue_vec,
       method = method,
       null_mean = null_mean,
       null_sd = null_sd,
       pvalue_vec = pvalue_vec)
}

####################

.multtest_locfdr <- function(teststat_vec){
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
       null_mean = res[1],
       null_sd = res[2])
}

# observed_quantile is two numbers, a lower and and a upper
# see equation 4.8 onwards of "SIZE, POWER AND FALSE DISCOVERY RATES" by Efron
# TODO: If there is an obvious break in the values, hard-set the observed quantiles. We can do this by mixture-modeling
.multtest_truncatedGauss <- function(teststat_vec,
                                     observed_quantile){
  z_vec <- teststat_vec

  fn <- function(param_vec,
                 lb,
                 N,
                 N0,
                 z0_vec,
                 ub){
    delta0 <- param_vec[1]
    sigma0 <- param_vec[2]
    theta <-  param_vec[3]
    if(theta > 0.99 | theta < 0.01) return(Inf)
    denom <- stats::pnorm(ub, mean = delta0, sd = sigma0) - stats::pnorm(lb, mean = delta0, sd = sigma0)

    # compute each sample's likelihood
    sample_llvec <- sapply(z0_vec, function(z){
      stats::dnorm(z, mean = delta0, sd = sigma0, log = T)
    })

    # compute full log-likelihood, Equation 4.12
    loglik <- N0*log(theta) + (N-N0)*log(1-theta) + sum(sample_llvec) - N0*log(denom)
    -loglik
  }

  N <- length(z_vec)
  tmp <- stats::quantile(z_vec, probs = observed_quantile)
  lb <- tmp[1]; ub <- tmp[2]
  idx <- intersect(which(z_vec >= lb), which(z_vec <= ub))
  N0 <- length(idx)
  z0_vec <- z_vec[idx]

  init_delta0 <- mean(z0_vec)
  init_sigma0 <- stats::sd(z0_vec)
  init_theta <- 0.95*(stats::pnorm(ub, mean = init_delta0, sd = init_sigma0) - stats::pnorm(lb, mean = init_delta0, sd = init_sigma0))

  res <- tryCatch({
    optim_res <- stats::optim(par = c(init_delta0, init_sigma0, init_theta),
                        fn = fn,
                        method = "Nelder-Mead",
                        lb = lb,
                        N = N,
                        N0 = N0,
                        z0_vec = z0_vec,
                        ub = ub)
    c(optim_res$par[1], optim_res$par[2])
  }, error = function(e){
    rep(NA, 2)
  })

  list(method = "truncated_mle",
       null_mean = res[1],
       null_sd = res[2])
}

# this is a purposefully overly simple method, meant to use when everything else has failed. Realistically, it is not a good estimator
.multtest_simple <- function(teststat_vec,
                             observed_quantile){
  tmp <- stats::quantile(teststat_vec, probs = observed_quantile)
  lb <- tmp[1]; ub <- tmp[2]
  idx <- intersect(which(teststat_vec >= lb), which(teststat_vec <= ub))
  vec <- teststat_vec[idx]

  list(method = "simple",
       null_mean = mean(vec),
       null_sd = stats::sd(vec))
}
