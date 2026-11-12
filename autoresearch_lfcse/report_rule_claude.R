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
# The rule. For donor d and gene g let m_dg be the ratio of sums, the
# donor's counts over its fitted library sizes. With arm means M1, M0 over
# n1, n0 donors and the squared coefficients of variation C1 = V1 / M1^2,
# C0 = V0 / M0^2 of the donor means (V the sample variance),
#   lfc = log2(M1 / M0),
#   se  = sqrt(C1~ / n1 + C0~ / n0) / ln 2,
# where C~ is C moderated toward its trend in mean expression
# (`.moderate_variance()`). The empirical null (mean d0, SD s0) is fitted
# to z = lfc / se over the genes, and the rule reports
# (lfc - d0 * median(se), s0 * se).

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

.arm_of_donor <- function(case_control_vec, individual_vec, donor_vec){
  case_control_vec[match(donor_vec, as.character(individual_vec))]
}

.lfc_from_donor_means <- function(mean_mat, arm_vec){
  case_mat <- mean_mat[arm_vec == 1, , drop = FALSE]
  control_mat <- mean_mat[arm_vec == 0, , drop = FALSE]
  n1 <- nrow(case_mat)
  n0 <- nrow(control_mat)
  case_mean <- colMeans(case_mat)
  control_mean <- colMeans(control_mat)
  case_var <- apply(case_mat, 2, stats::var)
  control_var <- apply(control_mat, 2, stats::var)
  covariate_vec <- log(colMeans(mean_mat))
  case_cv2 <- .moderate_variance(var_vec = case_var / case_mean^2,
                                 df = n1 - 1, covariate_vec = covariate_vec)
  control_cv2 <- .moderate_variance(var_vec = control_var / control_mean^2,
                                    df = n0 - 1,
                                    covariate_vec = covariate_vec)
  list(lfc_vec = log2(case_mean / control_mean),
       se_vec = sqrt(case_cv2 / n1 + control_cv2 / n0) / log(2))
}

report_lfc_se <- function(view){
  # Donor summary: the ratio of sums, sum of counts over sum of fitted
  # library sizes, which weights each cell by its depth.
  individual_chr <- as.character(view$individual_vec)
  mean_mat <- rowsum(view$dat, group = individual_chr) /
    rowsum(view$library_mat, group = individual_chr)
  arm_vec <- .arm_of_donor(case_control_vec = view$case_control_vec,
                           individual_vec = view$individual_vec,
                           donor_vec = rownames(mean_mat))
  raw <- .lfc_from_donor_means(mean_mat = mean_mat, arm_vec = arm_vec)

  z_vec <- raw$lfc_vec / raw$se_vec
  # `multtest()` warns when locfdr falls back on its cruder estimators; the
  # estimate is still used.
  null_res <- suppressWarnings(eSVD2::multtest(z_vec[is.finite(z_vec)]))

  data.frame(gene = colnames(view$dat),
             log2fc = raw$lfc_vec -
               null_res$null_mean * stats::median(raw$se_vec),
             se = null_res$null_sd * raw$se_vec,
             stringsAsFactors = FALSE)
}
