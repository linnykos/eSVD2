# Candidate nuisance estimators for the overdispersion brainstorm
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Sourced by 01_cache_fits_claude.R and 02_run_candidates_claude.R. Nothing here
# changes the package: every candidate is a function that maps the cached fit
# to a vector of Gamma RATES (what the package stores in `nuisance_vec`), which
# is then written over the devel estimate before `compute_posterior()`.

# The likelihood ---------------------------------------------------------------

# Log-likelihood of one gene in rho = log(beta), dropping the terms free of
# beta. Marginally A ~ NB(shape = mu * beta, prob = s / (s + beta)); see the
# header of src/gamma_rate.cpp.
.nb_loglik <- function(rho,
                       x_vec,
                       mu_vec,
                       s_vec){
  beta_val <- exp(rho)
  shape_vec <- mu_vec * beta_val
  sum(lgamma(shape_vec + x_vec) - lgamma(shape_vec) -
        shape_vec * log1p(s_vec / beta_val) - x_vec * log(s_vec + beta_val))
}

# The limit of `.nb_loglik()` as beta -> Inf, with the same terms dropped, so
# the two are comparable.
.poisson_loglik <- function(x_vec,
                            mu_vec,
                            s_vec){
  sum(x_vec * log(mu_vec) - mu_vec * s_vec)
}

# Observed information of rho at `rho`, by a central difference.
.nb_information <- function(rho,
                            x_vec,
                            mu_vec,
                            s_vec,
                            step_val = 0.05){
  ll_vec <- sapply(rho + c(-1, 0, 1) * step_val, function(r){
    .nb_loglik(rho = r, x_vec = x_vec, mu_vec = mu_vec, s_vec = s_vec)
  })
  -(ll_vec[1] - 2 * ll_vec[2] + ll_vec[3]) / step_val^2
}

# The pieces every candidate shares -------------------------------------------

# The two matrices `estimate_nuisance()` builds, with the arguments the
# comparison pipeline uses (`bool_covariates_as_library = TRUE`, intercept in
# the library).
.fitted_matrices <- function(esvd_obj,
                             cc_var){
  fit <- esvd_obj[[esvd_obj[["latest_Fit"]]]]
  covariates <- esvd_obj$covariates
  library_idx <- which(colnames(covariates) != cc_var)
  nat_mat <- tcrossprod(fit$x_mat, fit$y_mat) +
    tcrossprod(covariates[, -library_idx, drop = FALSE],
               fit$z_mat[, -library_idx, drop = FALSE])

  list(library_mat = exp(tcrossprod(covariates[, library_idx, drop = FALSE],
                                    fit$z_mat[, library_idx, drop = FALSE])),
       mean_mat = exp(nat_mat))
}

# Per-gene summaries of the likelihood that the candidates are built from.
# `bool_boundary` marks the genes whose likelihood is maximized at the Poisson
# limit: the Poisson log-likelihood is at least that of the returned rate.
.gene_summaries <- function(dat,
                            mean_mat,
                            library_mat,
                            rate_mle_vec,
                            profile_level_vec = c(0.5, 0.9),
                            rho_lower = -10){
  p <- ncol(dat)
  res_list <- lapply(seq_len(p), function(j){
    x_vec <- as.numeric(dat[, j])
    mu_vec <- mean_mat[, j]
    s_vec <- library_mat[, j]

    rho_hat <- log(rate_mle_vec[j])
    ll_hat <- .nb_loglik(rho = rho_hat, x_vec = x_vec, mu_vec = mu_vec,
                         s_vec = s_vec)
    ll_poisson <- .poisson_loglik(x_vec = x_vec, mu_vec = mu_vec,
                                  s_vec = s_vec)
    bool_boundary <- ll_poisson >= ll_hat - 1e-6
    ll_max <- max(ll_hat, ll_poisson)

    # Lower end of the profile-likelihood interval of rho. The likelihood of
    # a boundary gene is flat to the right, but it still falls to the left,
    # so the lower end exists even where the maximizer does not.
    profile_lower_vec <- sapply(profile_level_vec, function(level){
      target_val <- ll_max - stats::qchisq(level, df = 1) / 2
      fn <- function(rho){
        .nb_loglik(rho = rho, x_vec = x_vec, mu_vec = mu_vec,
                   s_vec = s_vec) - target_val
      }
      rho_upper <- min(rho_hat, 25)
      if(fn(rho_lower) >= 0) return(rho_lower)
      if(fn(rho_upper) <= 0) return(rho_upper)
      stats::uniroot(fn, interval = c(rho_lower, rho_upper))$root
    })

    information_val <- .nb_information(rho = rho_hat, x_vec = x_vec,
                                       mu_vec = mu_vec, s_vec = s_vec)
    fitted_vec <- mu_vec * s_vec
    moment_val <- sum((x_vec - fitted_vec)^2 - fitted_vec) /
      sum(mu_vec * s_vec^2)

    c(bool_boundary = bool_boundary,
      information = information_val,
      ll_gain = ll_hat - ll_poisson,
      moment_excess = moment_val,
      pearson = mean((x_vec - fitted_vec)^2 / fitted_vec),
      profile_lower = exp(profile_lower_vec),
      s_max = max(s_vec),
      s_median = stats::median(s_vec))
  })
  res_df <- as.data.frame(do.call(rbind, res_list))
  colnames(res_df)[grep("^profile_lower", colnames(res_df))] <-
    paste0("profile_lower_", 100 * profile_level_vec)
  res_df$bool_boundary <- as.logical(res_df$bool_boundary)
  res_df
}

# Empirical-Bayes MAP rate under a Normal prior on the unit-free log rate,
# log(beta_j / median_i s_ji). The centre and spread are read off the INTERIOR
# genes only. In nearly Poisson data those are the genes whose residuals
# happened to look overdispersed, so the centre sits below the truth, which
# errs towards more overdispersion and a more cautious test. The spread uses
# the lower half of the interior genes, because near-divergence corrupts only
# the upper tail; the median sampling variance is subtracted, as DESeq2 does,
# and the result is floored at DESeq2's 0.25. A non-NULL `prior_var` replaces
# the estimate, for the sweep of 05_cap_sweep_claude.R.
#
# Two switches bring the prior closer to DESeq2's, for
# 07_deseq2_variants_claude.R. `bool_trend = TRUE` replaces the constant
# centre by a line in the log of the gene's mean count, fitted to the interior
# genes. A non-NULL `sampling_var` replaces the inverse observed information,
# e.g. by DESeq2's trigamma((m - p) / 2).
.map_rates <- function(dat,
                       mean_mat,
                       library_mat,
                       rate_mle_vec,
                       summary_df,
                       bool_trend = FALSE,
                       prior_var = NULL,
                       prior_var_floor = 0.25,
                       rho_range = c(-10, 25),
                       sampling_var = NULL){
  interior_vec <- !summary_df$bool_boundary
  eta_vec <- log(rate_mle_vec) - log(summary_df$s_median)
  log_count_vec <- log(colMeans(dat))

  slope_val <- 0
  if(bool_trend){
    trend_fit <- stats::lm(eta ~ log_count,
                           data = data.frame(eta = eta_vec[interior_vec],
                                             log_count = log_count_vec[interior_vec]))
    slope_val <- as.numeric(stats::coef(trend_fit)[2])
  }
  # The intercept is the median residual in both cases, so that the constant
  # centre is the special case of slope 0.
  residual_vec <- (eta_vec - slope_val * log_count_vec)[interior_vec]
  intercept_val <- stats::median(residual_vec)
  centre_vec <- intercept_val + slope_val * log_count_vec

  lower_sd_val <- (intercept_val -
                     stats::quantile(residual_vec, probs = 0.25)) /
    stats::qnorm(0.75)
  if(is.null(sampling_var)){
    interior_idx <- which(interior_vec & summary_df$information > 0)
    sampling_var <- stats::median(1 / summary_df$information[interior_idx])
  }
  if(is.null(prior_var)){
    prior_var <- max(lower_sd_val^2 - sampling_var, prior_var_floor)
  }

  rate_vec <- sapply(seq_len(ncol(dat)), function(j){
    x_vec <- as.numeric(dat[, j])
    prior_mean <- centre_vec[j] + log(summary_df$s_median[j])
    fn <- function(rho){
      .nb_loglik(rho = rho, x_vec = x_vec, mu_vec = mean_mat[, j],
                 s_vec = library_mat[, j]) -
        (rho - prior_mean)^2 / (2 * prior_var)
    }
    exp(stats::optimize(fn, interval = rho_range, maximum = TRUE)$maximum)
  })

  list(centre = stats::median(centre_vec),
       prior_var = prior_var,
       rate_vec = rate_vec,
       sampling_var = sampling_var,
       slope = slope_val)
}

# The candidates ---------------------------------------------------------------

# Every candidate as a named list of rate vectors. `rate_master_vec` and
# `rate_true_vec` may be NULL, in which case their candidates are left out.
.candidate_rates <- function(dat,
                             mean_mat,
                             library_mat,
                             rate_mle_vec,
                             summary_df,
                             rate_master_vec = NULL,
                             rate_true_vec = NULL){
  s_median_vec <- summary_df$s_median
  interior_vec <- !summary_df$bool_boundary
  # Medians over the interior genes: with half the genes at the boundary,
  # the median over all genes is itself a diverged value.
  median_val <- stats::median(rate_mle_vec[interior_vec])
  unit_free_median <- stats::median((rate_mle_vec / s_median_vec)[interior_vec])
  map_res <- .map_rates(dat = dat, mean_mat = mean_mat,
                        library_mat = library_mat,
                        rate_mle_vec = rate_mle_vec,
                        summary_df = summary_df)

  boundary_to_median_vec <- rate_mle_vec
  boundary_to_median_vec[!interior_vec] <-
    (unit_free_median * s_median_vec)[!interior_vec]

  candidate_list <- list(
    mle = rate_mle_vec,
    cap_max_s = pmin(rate_mle_vec, summary_df$s_max),
    cap_10s = pmin(rate_mle_vec, 10 * s_median_vec),
    cap_50s = pmin(rate_mle_vec, 50 * s_median_vec),
    cap_3median = pmin(rate_mle_vec, 3 * median_val),
    cap_10median = pmin(rate_mle_vec, 10 * median_val),
    winsor_95 = pmin(rate_mle_vec, stats::quantile(rate_mle_vec, probs = 0.95)),
    boundary_to_median = boundary_to_median_vec,
    common_raw = rep(median_val, length(rate_mle_vec)),
    common_unit = unit_free_median * s_median_vec,
    eb_map = map_res$rate_vec,
    profile_lower_50 = summary_df$profile_lower_50,
    profile_lower_90 = summary_df$profile_lower_90,
    moment_cap_50s = 1 / pmax(summary_df$moment_excess,
                              1 / (50 * s_median_vec))
  )
  if(!is.null(rate_master_vec)) candidate_list$master <- rate_master_vec
  # The generator's rate is in units where the library size averages 1; the
  # fit puts the gene intercept into the library, so the same overdispersion
  # is the rate times the fitted library size.
  if(!is.null(rate_true_vec)){
    candidate_list$oracle <- rate_true_vec * s_median_vec
  }

  attr(candidate_list, "map_centre") <- map_res$centre
  attr(candidate_list, "map_prior_var") <- map_res$prior_var
  candidate_list
}

# From a rate vector to p-values ----------------------------------------------

# The last three stages of `.run_pipeline()` in
# version_comparison/helpers_claude.R, with the same arguments. Warnings are
# muffled: the locfdr fallback warning fires on many of these small cohorts.
# `bool_median_stabilize` replaces the package's rescaling by the geometric
# mean with a rescaling by the median, done here before the posterior.
.rates_to_pvalues <- function(esvd_obj,
                              rate_vec,
                              bool_median_stabilize = FALSE){
  latest_fit <- esvd_obj[["latest_Fit"]]
  bool_stabilize <- TRUE
  if(bool_median_stabilize){
    if(stats::median(rate_vec) > 1){
      rate_vec <- rate_vec / stats::median(rate_vec)
    }
    bool_stabilize <- FALSE
  }
  esvd_obj[[latest_fit]]$nuisance_vec[] <- rate_vec

  withCallingHandlers({
    esvd_obj <- eSVD2::compute_posterior(
      input_obj = esvd_obj,
      alpha_max = 2 * max(esvd_obj$dat),
      bool_covariates_as_library = TRUE,
      bool_stabilize_underdispersion = bool_stabilize,
      library_min = 0.1
    )
    esvd_obj <- eSVD2::compute_test_statistic(input_obj = esvd_obj,
                                              verbose = 0)
    esvd_obj <- eSVD2::compute_pvalue(input_obj = esvd_obj)
  }, warning = function(w){invokeRestart("muffleWarning")})

  pvalue_list <- esvd_obj$pvalue_list
  list(gene_df = data.frame(fdr = as.numeric(pvalue_list$fdr_vec),
                            gaussian_teststat = as.numeric(pvalue_list$gaussian_teststat),
                            log10p = as.numeric(pvalue_list$log10pvalue),
                            logFC = as.numeric(esvd_obj$log2fc_vec),
                            logFC_se = as.numeric(esvd_obj$log2fc_se_vec),
                            teststat = as.numeric(esvd_obj$teststat_vec)),
       method = pvalue_list$method,
       null_mean = as.numeric(pvalue_list$null_mean),
       null_sd = as.numeric(pvalue_list$null_sd))
}
