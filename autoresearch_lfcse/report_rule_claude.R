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
# (lfc - d0 * median(se), sqrt(F) * se), where F is the factor by which
# the donor-level covariates of the design widen the variance of the
# case-control coefficient (`.design_factor()`; 1 without covariates).

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

# Variance factor of the donor-level design. The fitted library removes
# every covariate that is constant within donors (sex, age) from the donor
# means, so the two-arm variance stands in for the variance of the
# case-control coefficient of a donor-level regression on those
# covariates. With D the donors x columns design (intercept, arm, and the
# donor-level covariates), that variance is sigma^2 [(D'D)^-1]_arm with
# sigma^2 on n - rank(D) degrees of freedom, against sigma^2 (1/n1 + 1/n0)
# on n - 2 for the two arms alone. The factor is the ratio of the two.
.design_factor <- function(covariates, individual_vec, arm_vec, donor_vec){
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
  keep_idx <- qr_res$pivot[seq_len(qr_res$rank)]
  design_mat <- design_mat[, keep_idx, drop = FALSE]
  num_donors <- nrow(design_mat)
  n1 <- sum(arm_vec == 1)
  n0 <- sum(arm_vec == 0)
  arm_var <- solve(crossprod(design_mat))[2, 2]
  list(design_mat = design_mat,
       factor = (arm_var / (1 / n1 + 1 / n0)) *
         ((num_donors - 2) / (num_donors - ncol(design_mat))))
}

.arm_of_donor <- function(case_control_vec, individual_vec, donor_vec){
  case_control_vec[match(donor_vec, as.character(individual_vec))]
}

.lfc_from_donor_means <- function(mean_mat, log_mat, arm_vec){
  case_mat <- mean_mat[arm_vec == 1, , drop = FALSE]
  control_mat <- mean_mat[arm_vec == 0, , drop = FALSE]
  n1 <- nrow(case_mat)
  n0 <- nrow(control_mat)
  case_mean <- colMeans(case_mat)
  control_mean <- colMeans(control_mat)
  # The squared CV of an arm's donor means is estimated by the variance of
  # their logarithms (equal to first order, and not driven by the largest
  # donor).
  case_var <- apply(log_mat[arm_vec == 1, , drop = FALSE], 2, stats::var)
  control_var <- apply(log_mat[arm_vec == 0, , drop = FALSE], 2, stats::var)
  covariate_vec <- log(colMeans(mean_mat))
  case_cv2 <- .moderate_variance(var_vec = case_var,
                                 df = n1 - 1, covariate_vec = covariate_vec)
  control_cv2 <- .moderate_variance(var_vec = control_var,
                                    df = n0 - 1,
                                    covariate_vec = covariate_vec)
  # The logFC is the difference of the arms' mean log donor means (the log
  # ratio of geometric means over donors), which a donor far above the
  # others moves less than it moves the arithmetic mean.
  list(lfc_vec = (colMeans(log_mat[arm_vec == 1, , drop = FALSE]) -
                    colMeans(log_mat[arm_vec == 0, , drop = FALSE])) / log(2),
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
  # Half a count keeps the logarithm of a donor with no count finite (the
  # continuity correction of a Poisson count).
  log_mat <- log((rowsum(view$dat, group = individual_chr) + 0.5) /
                   rowsum(view$library_mat, group = individual_chr))
  raw <- .lfc_from_donor_means(mean_mat = mean_mat, log_mat = log_mat,
                               arm_vec = arm_vec)

  # The within-donor term 1.2.0 stores, per donor mean: the donor's average
  # posterior variance over its number of cells, added to the between term.
  num_cells_vec <- as.numeric(table(individual_chr)[rownames(mean_mat)])
  within_mat <- rowsum(view$posterior_var_mat, group = individual_chr) /
    num_cells_vec^2
  bool_case <- arm_vec == 1
  within_vec <-
    colSums(within_mat[bool_case, , drop = FALSE]) /
    (sum(bool_case) * colMeans(mean_mat[bool_case, , drop = FALSE]))^2 +
    colSums(within_mat[!bool_case, , drop = FALSE]) /
    (sum(!bool_case) * colMeans(mean_mat[!bool_case, , drop = FALSE]))^2
  raw$se_vec <- sqrt(raw$se_vec^2 + within_vec / log(2)^2)

  design_res <- .design_factor(covariates = view$covariates,
                               individual_vec = view$individual_vec,
                               arm_vec = arm_vec,
                               donor_vec = rownames(mean_mat))

  # The logFC is the arm coefficient of the donor-level least-squares
  # regression of the log donor means on the design, the regression whose
  # variance the design factor describes. Without donor-level covariates
  # it is the difference of the arms' mean log donor means.
  coef_mat <- solve(crossprod(design_res$design_mat),
                    crossprod(design_res$design_mat, log_mat))
  raw$lfc_vec <- coef_mat[2, ] / log(2)

  z_vec <- raw$lfc_vec / raw$se_vec
  # `multtest()` warns when locfdr falls back on its cruder estimators; the
  # estimate is still used.
  null_res <- suppressWarnings(eSVD2::multtest(z_vec[is.finite(z_vec)]))

  data.frame(gene = colnames(view$dat),
             log2fc = raw$lfc_vec -
               null_res$null_mean * stats::median(raw$se_vec),
             se = sqrt(design_res$factor) * raw$se_vec,
             stringsAsFactors = FALSE)
}
