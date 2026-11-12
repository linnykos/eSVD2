# `embed_prior`: the reporting rule kept by run A (`autoresearch_lfcse2/
# nov11`) of the second autoresearch loop: a log2 fold change (the log
# ratio of mean expression between the arms, cell-state differences
# included) and its standard error, from one fitted eSVD-DE cohort.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11.
#
# Mirrors commit 40e7e67 of `autoresearch_lfcse2/report_rule_claude.R` (the
# code is identical; only this header differs).
# Score: Gaussian negative log-likelihood of the true log2 fold change
# (ESVD2_TRUTH=marginal), nats per gene, lower is better, on the bank
# `lfcse-bank-2`.
#              baseline `null_rescaled`   `donor_wls`   this rule
#   dev            0.3066                   0.1508        0.1506
#   test           0.2969                   0.1319        0.1321
#
# Use: `res_df <- report_lfc_se(view)`; reads `dat`, `library_mat`,
# `covariates`, `individual_vec`, `case_control_vec`, `x_mat`, `y_mat` of
# the view and returns `gene`, `log2fc`, `se`. Needs `stats` and `eSVD2`.
# The embedding enters the SE only (the trend the residual variance is
# moderated toward); the logFC is that of `donor_wls`.
#
# The rule. Each donor is summarized, gene by gene, by l = log((summed
# counts + 1/2) / summed fitted library), and by its embedding offset
# e = log of the library-weighted mean over its cells of exp(x' y). The
# summaries l are regressed on the donor-level design D (intercept, arm,
# covariates constant within donors) with precision weights, which gives
# the arm coefficient, its variance factor a and the residual variance S.
# S is moderated across genes toward a trend in two covariates: the gene's
# level, and V, the weighted variance around D of the gene's embedding
# offset, which is the donor-to-donor variance that the embedding predicts
# for the gene. se^2 = a * S~ / (ln 2)^2. The logFC is shifted by the
# empirical null's mean times the median SE.

# Inverse of the trigamma function by Newton's method on its reciprocal
# (Smyth 2004, appendix); the start 0.5 + 1 / x is exact in both limits.
.trigamma_inverse <- function(x_val){
  y_val <- 0.5 + 1 / x_val
  for(iter in seq_len(50)){
    tri <- trigamma(y_val)
    step <- tri * (1 - tri / x_val) / psigamma(y_val, deriv = 2)
    y_val <- y_val + step
    # Relative step below single precision: converged.
    if(-step / y_val < 1e-8) break
  }
  y_val
}

# Moderate per-gene variances `var_vec`, each with `df` degrees of
# freedom, toward their trend in `covariate_vec` and, additively, in
# `second_vec` (Smyth 2004): the trend is a smoothing spline of the log
# variance on the first, plus a smoothing spline of what is left on the
# second, the excess of the residual
# variance over trigamma(df / 2) gives the prior degrees of freedom d0, and
# the moderated variance is (d0 * s0^2 + df * var) / (d0 + df).
.moderate_variance <- function(var_vec, df, covariate_vec, second_vec){
  idx <- which(is.finite(var_vec) & var_vec > 0 & is.finite(covariate_vec) &
                 is.finite(second_vec))
  log_vec <- log(var_vec[idx])
  spline_fit <- stats::smooth.spline(x = covariate_vec[idx], y = log_vec)
  trend_vec <- stats::predict(spline_fit, x = covariate_vec)$y
  second_fit <- stats::smooth.spline(x = second_vec[idx],
                                     y = log_vec - trend_vec[idx])
  trend_vec <- trend_vec + stats::predict(second_fit, x = second_vec)$y
  excess <- mean((log_vec - trend_vec[idx])^2) - trigamma(df / 2)
  if(excess <= 0){
    # No spread beyond sampling noise: every gene takes the trend.
    return(exp(trend_vec - digamma(df / 2) + log(df / 2)))
  }
  prior_df <- 2 * .trigamma_inverse(excess)
  prior_vec <- exp(trend_vec - digamma(df / 2) + log(df / 2) +
                     digamma(prior_df / 2) - log(prior_df / 2))
  (prior_df * prior_vec + df * var_vec) / (prior_df + df)
}

# The donor-level design: intercept, arm, and every column of `covariates`
# that is constant within donors (sex, age).
.donor_design <- function(covariates, individual_vec, arm_vec, donor_vec){
  individual_chr <- as.character(individual_vec)
  first_idx <- match(donor_vec, individual_chr)
  bool_donor_level <- apply(covariates, 2, function(col_vec){
    all(col_vec == col_vec[first_idx][match(individual_chr, donor_vec)])
  })
  design_mat <- cbind(1, arm_vec,
                      covariates[first_idx, bool_donor_level, drop = FALSE])
  qr_res <- qr(design_mat)
  # Columns 1 and 2 (intercept, arm) are pivoted first, so the constant
  # and arm columns of `covariates` itself drop out as redundant.
  design_mat[, qr_res$pivot[seq_len(qr_res$rank)], drop = FALSE]
}

.least_squares_fit <- function(design_mat, response_mat){
  design_mat %*% solve(crossprod(design_mat),
                       crossprod(design_mat, response_mat))
}

# Weighted least squares of each column of `log_mat` on `design_mat`, with
# precision weights from the two-part variance tau^2 + phi / (expected
# count) fitted over all donors and genes. Returns the arm coefficient
# (natural log), its variance factor, the residual variance on
# n - rank(design) degrees of freedom, the mean fitted value, and the
# weights.
.donor_regression <- function(log_mat, library_sum_mat, design_mat){
  fit_mat <- .least_squares_fit(design_mat, log_mat)
  # A donor's offset common to all genes: the median over genes of its
  # residuals, removed before the variance is estimated.
  log_mat <- log_mat - apply(log_mat - fit_mat, 1, stats::median)
  fit_mat <- .least_squares_fit(design_mat, log_mat)
  residual_df <- nrow(design_mat) - ncol(design_mat)
  inverse_count_vec <- as.numeric(exp(-(fit_mat + log(library_sum_mat))))
  squared_vec <- as.numeric((log_mat - fit_mat)^2) *
    nrow(design_mat) / residual_df
  slope <- stats::cov(squared_vec, inverse_count_vec) /
    stats::var(inverse_count_vec)
  # Neither part of a variance can be negative.
  variance_coef <- pmax(c(mean(squared_vec) -
                            slope * mean(inverse_count_vec), slope), 0)
  weight_mat <- matrix(1 / (variance_coef[1] +
                              variance_coef[2] * inverse_count_vec),
                       nrow = nrow(fit_mat))
  num_genes <- ncol(log_mat)
  coef_arm_vec <- numeric(num_genes)
  arm_var_vec <- numeric(num_genes)
  residual_var <- numeric(num_genes)
  for(j in seq_len(num_genes)){
    weighted_design <- design_mat * weight_mat[, j]
    gram_inv_mat <- solve(crossprod(design_mat, weighted_design))
    coef_vec <- gram_inv_mat %*% crossprod(weighted_design, log_mat[, j])
    coef_arm_vec[j] <- coef_vec[2]
    arm_var_vec[j] <- gram_inv_mat[2, 2]
    residual_var[j] <- sum(weight_mat[, j] *
                             (log_mat[, j] - design_mat %*% coef_vec)^2) /
      residual_df
  }
  list(arm_var_vec = arm_var_vec,
       coef_arm_vec = coef_arm_vec,
       level_vec = colMeans(fit_mat),
       residual_df = residual_df,
       residual_var = residual_var,
       weight_mat = weight_mat)
}

report_lfc_se <- function(view){
  individual_chr <- as.character(view$individual_vec)
  library_sum_mat <- rowsum(view$library_mat, group = individual_chr)
  donor_vec <- rownames(library_sum_mat)
  count_mat <- rowsum(view$dat, group = individual_chr)
  # Half a count keeps the logarithm of a donor with no count finite (the
  # continuity correction of a Poisson count).
  raw_mat <- log((count_mat + 0.5) / library_sum_mat)
  # The embedding offset of each donor and gene.
  embedding_mat <- exp(tcrossprod(view$x_mat, view$y_mat))
  offset_mat <- log(rowsum(view$library_mat * embedding_mat,
                           group = individual_chr) / library_sum_mat)
  arm_vec <- view$case_control_vec[match(donor_vec, individual_chr)]
  design_mat <- .donor_design(covariates = view$covariates,
                              individual_vec = view$individual_vec,
                              arm_vec = arm_vec, donor_vec = donor_vec)
  res <- .donor_regression(log_mat = raw_mat,
                           library_sum_mat = library_sum_mat,
                           design_mat = design_mat)
  lfc_vec <- res$coef_arm_vec / log(2)
  # The donor-to-donor variance the embedding predicts for each gene: the
  # weighted variance of its offset around the design, on the log scale.
  offset_res_mat <- offset_mat - .least_squares_fit(design_mat, offset_mat)
  offset_var <- colSums(res$weight_mat * offset_res_mat^2) / res$residual_df
  se_vec <- sqrt(res$arm_var_vec *
                   .moderate_variance(var_vec = res$residual_var,
                                      df = res$residual_df,
                                      covariate_vec = res$level_vec,
                                      second_vec = log(offset_var))) / log(2)

  z_vec <- lfc_vec / se_vec
  # `multtest()` warns when locfdr falls back on its cruder estimators; the
  # estimate is still used.
  null_res <- suppressWarnings(eSVD2::multtest(z_vec[is.finite(z_vec)]))

  data.frame(gene = colnames(view$dat),
             log2fc = lfc_vec - null_res$null_mean * stats::median(se_vec),
             se = se_vec,
             stringsAsFactors = FALSE)
}
