# The reporting rule on trial: a log2 fold change and its standard error
# from one fitted eSVD-DE cohort. Drafted by Claude for Kevin Z. Lin.
#
# This is the ONE file an experiment edits. Contract (notes/MEMO_claude.md,
# section 2): `report_lfc_se(view)` returns a data frame with one row per
# gene in the order of `colnames(view$dat)` and columns `gene`, `log2fc`,
# `se`. It may call base R, `stats`, `Matrix`, `locfdr` and eSVD2's exported
# functions, always with `::`; the top level holds function definitions and
# named constants only.
#
# The rule. For donor d and gene g let y_dg be the donor's summed counts
# and L_dg its summed fitted library sizes, and l_dg = log((y_dg + 1/2) /
# L_dg). With D the donor-level design (intercept, arm, and every covariate
# constant within donors), the log donor means are regressed on D by least
# squares, gene by gene:
#   lfc  = (arm coefficient) / ln 2,
#   se^2 = S~ * a / (ln 2)^2,
# where a = [(D'D)^-1]_arm and S~ is the residual variance on n - rank(D)
# degrees of freedom moderated toward its trend in mean expression
# (`.moderate_variance()`). The empirical null's mean d0 is fitted to
# z = lfc / se over the genes, and the rule reports
# (lfc - d0 * median(se), se).

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
# freedom, toward their trend in `covariate_vec` (Smyth 2004). On the log
# scale a scaled chi-square variance has mean log(s0^2) + digamma(df / 2) -
# log(df / 2) and variance trigamma(df / 2); the trend is a smoothing
# spline (smoothness by generalized cross-validation), and the excess of
# the residual variance over trigamma(df / 2) gives the prior degrees of
# freedom d0 through trigamma(d0 / 2). The moderated variance is
# (d0 * s0^2 + df * var) / (d0 + df).
.moderate_variance <- function(var_vec, df, covariate_vec){
  idx <- which(is.finite(var_vec) & var_vec > 0 & is.finite(covariate_vec))
  log_vec <- log(var_vec[idx])
  spline_fit <- stats::smooth.spline(x = covariate_vec[idx], y = log_vec)
  trend_vec <- stats::predict(spline_fit, x = covariate_vec)$y
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
# that is constant within donors (sex, age). The fitted library removes
# those covariates cell by cell; the donor-level regression removes what
# is left of them between donors and charges the arm coefficient for it.
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

.arm_of_donor <- function(case_control_vec, individual_vec, donor_vec){
  case_control_vec[match(donor_vec, as.character(individual_vec))]
}

report_lfc_se <- function(view){
  # Donor summary: the ratio of sums, sum of counts over sum of fitted
  # library sizes, which weights each cell by its depth.
  individual_chr <- as.character(view$individual_vec)
  count_mat <- rowsum(view$dat, group = individual_chr)
  library_sum_mat <- rowsum(view$library_mat, group = individual_chr)
  mean_mat <- count_mat / library_sum_mat
  # Half a count keeps the logarithm of a donor with no count finite (the
  # continuity correction of a Poisson count).
  log_mat <- log((count_mat + 0.5) / library_sum_mat)
  donor_vec <- rownames(mean_mat)
  arm_vec <- .arm_of_donor(case_control_vec = view$case_control_vec,
                           individual_vec = view$individual_vec,
                           donor_vec = donor_vec)

  # Least squares of the log donor means on the donor-level design.
  design_mat <- .donor_design(covariates = view$covariates,
                              individual_vec = view$individual_vec,
                              arm_vec = arm_vec, donor_vec = donor_vec)
  fit_mat <- design_mat %*% solve(crossprod(design_mat),
                                  crossprod(design_mat, log_mat))
  # A donor's offset common to all genes (depth the library missed): the
  # median over genes of its residuals, removed before the variance is
  # estimated. The arm's share of such offsets is what the centring removes
  # from the logFC, so it does not belong in the SE either.
  offset_vec <- apply(log_mat - fit_mat, 1, stats::median)
  log_mat <- log_mat - offset_vec
  fit_mat <- design_mat %*% solve(crossprod(design_mat),
                                  crossprod(design_mat, log_mat))
  residual_df <- nrow(design_mat) - ncol(design_mat)
  residual_sd <- sqrt(colSums((log_mat - fit_mat)^2) / residual_df)

  # Precision weights (voom; Law et al. 2014). The square root of a gene's
  # residual SD is smoothed against its mean log expected donor count; a
  # donor's observation is given the weight 1 / trend^4 at its own expected
  # count, so that a donor with few counts for a gene counts for less.
  expected_mat <- fit_mat + log(library_sum_mat)
  gene_size_vec <- colMeans(expected_mat)
  trend_fit <- stats::smooth.spline(x = gene_size_vec, y = sqrt(residual_sd))
  clamped_vec <- pmin(pmax(as.numeric(expected_mat), min(gene_size_vec)),
                      max(gene_size_vec))
  weight_mat <- matrix(1 / stats::predict(trend_fit, x = clamped_vec)$y^4,
                       nrow = nrow(expected_mat))

  # Weighted least squares, gene by gene.
  num_genes <- ncol(log_mat)
  lfc_vec <- numeric(num_genes)
  arm_var_vec <- numeric(num_genes)
  residual_var <- numeric(num_genes)
  for(j in seq_len(num_genes)){
    weighted_design <- design_mat * weight_mat[, j]
    gram_inv_mat <- solve(crossprod(design_mat, weighted_design))
    coef_vec <- gram_inv_mat %*% crossprod(weighted_design, log_mat[, j])
    lfc_vec[j] <- coef_vec[2] / log(2)
    arm_var_vec[j] <- gram_inv_mat[2, 2]
    residual_var[j] <- sum(weight_mat[, j] *
                             (log_mat[, j] - design_mat %*% coef_vec)^2) /
      residual_df
  }
  between_vec <- arm_var_vec *
    .moderate_variance(var_vec = residual_var, df = residual_df,
                       # Mean of the two arms' log mean ratios: unlike the
                       # log of the pooled mean, it does not move with the
                       # arm effect, so the SE cannot follow the logFC.
                       covariate_vec = (
                         log(colMeans(mean_mat[arm_vec == 1, , drop = FALSE])) +
                           log(colMeans(mean_mat[arm_vec == 0, , drop = FALSE]))
                       ) / 2)

  se_vec <- sqrt(between_vec) / log(2)

  z_vec <- lfc_vec / se_vec
  # `multtest()` warns when locfdr falls back on its cruder estimators; the
  # estimate is still used.
  null_res <- suppressWarnings(eSVD2::multtest(z_vec[is.finite(z_vec)]))

  data.frame(gene = colnames(view$dat),
             log2fc = lfc_vec - null_res$null_mean * stats::median(se_vec),
             se = se_vec,
             stringsAsFactors = FALSE)
}
