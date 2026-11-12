# The reporting rule on trial: a log2 fold change and its standard error
# from one fitted eSVD-DE cohort. Drafted by Claude for Kevin Z. Lin.
#
# This is the ONE file an experiment edits. Contract (notes/MEMO_claude.md,
# section 2): `report_lfc_se(view)` returns a data frame with one row per
# gene in the order of `colnames(view$dat)` and columns `gene`, `log2fc`,
# `se`. It may call base R, `stats`, `Matrix`, `locfdr` and eSVD2's exported
# functions, always with `::`; the top level holds function definitions and
# named constants only. The logFC or the SE must use the fitted embedding,
# explicitly (`x_mat`, `y_mat`) or through the posterior.
#
# The rule estimates the case-control effect at fixed cell state: the part
# of the arm contrast that runs along the embedding's loadings is treated
# as a difference in cell state between the arms' donors and removed.
#
# Each donor is summarized, gene by gene, by l = log((summed counts + 1/2)
# / summed fitted library), and by its embedding offset e = log of the
# library-weighted mean over its cells of exp(x' y). Let e_r be the
# residual of e on the donor-level design D (intercept, arm, covariates
# constant within donors) and b the one coefficient per cohort of the
# least-squares regression of the residuals of l on e_r over all donors
# and genes. Then l - b e_r is regressed on D with precision weights, which
# gives the arm coefficient theta_g (that of l, since e_r is orthogonal to
# the arm), its variance factor a and a residual variance S free of the
# donor-to-donor variation along the embedding.
#
# The fit's reparameterization makes the embedding orthogonal to the arm,
# so the arms' difference in cell state sits in theta as a term y_g' d.
# The vector d is estimated by least squares of theta_g on the loadings
# y_g over the genes (with an intercept), and
#   lfc_g = (theta_g - y_g' d) / ln 2,
#   se_g^2 = (a S~ + y_g' Cov(d) y_g) / (ln 2)^2,
# with S~ the moderated S. The embedding offset was fitted to the same
# counts, so S understates what is left in theta; the SE is rescaled by the
# SD of the empirical null fitted to lfc / se, and the logFC is shifted by
# the null's mean times the median SE.

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
# freedom, toward their trend in `covariate_vec` (Smyth 2004): the trend is
# a smoothing spline of the log variance, the excess of the residual
# variance over trigamma(df / 2) gives the prior degrees of freedom d0, and
# the moderated variance is (d0 * s0^2 + df * var) / (d0 + df).
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
  # The offset's residual on the design, and its one coefficient.
  offset_res_mat <- offset_mat - .least_squares_fit(design_mat, offset_mat)
  raw_res_mat <- raw_mat - .least_squares_fit(design_mat, raw_mat)
  offset_coef <- sum(offset_res_mat * raw_res_mat) / sum(offset_res_mat^2)

  res <- .donor_regression(log_mat = raw_mat - offset_coef * offset_res_mat,
                           library_sum_mat = library_sum_mat,
                           design_mat = design_mat)
  # The arms' difference in cell state: the arm coefficients regressed on
  # the loadings over the genes.
  loading_mat <- cbind(1, view$y_mat)
  gram_inv_mat <- solve(crossprod(loading_mat))
  state_coef_vec <- gram_inv_mat %*% crossprod(loading_mat, res$coef_arm_vec)
  # The intercept stays in the logFC; the empirical null centres it below.
  state_coef_vec[1] <- 0
  state_part_vec <- as.numeric(loading_mat %*% state_coef_vec)
  lfc_vec <- (res$coef_arm_vec - state_part_vec) / log(2)
  # The variance of the removed part: the residual variance of that
  # regression over the genes times the loadings' leverage.
  coef_res_vec <- res$coef_arm_vec -
    as.numeric(loading_mat %*% (gram_inv_mat %*%
                                  crossprod(loading_mat, res$coef_arm_vec)))
  state_var_vec <- sum(coef_res_vec^2) /
    (nrow(loading_mat) - ncol(loading_mat)) *
    rowSums((loading_mat[, -1, drop = FALSE] %*%
               gram_inv_mat[-1, -1, drop = FALSE]) *
              loading_mat[, -1, drop = FALSE])
  se_vec <- sqrt(res$arm_var_vec *
                   .moderate_variance(var_vec = res$residual_var,
                                      df = res$residual_df,
                                      covariate_vec = res$level_vec) +
                   state_var_vec) / log(2)

  z_vec <- lfc_vec / se_vec
  # `multtest()` warns when locfdr falls back on its cruder estimators; the
  # estimate is still used.
  null_res <- suppressWarnings(eSVD2::multtest(z_vec[is.finite(z_vec)]))

  data.frame(gene = colnames(view$dat),
             log2fc = lfc_vec - null_res$null_mean * stats::median(se_vec),
             se = null_res$null_sd * se_vec,
             stringsAsFactors = FALSE)
}
