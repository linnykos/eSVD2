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
# Baseline `raw_null`: the unshrunk donor means with the empirical null.
# For donor d and gene g let m_dg be the mean of A / library over the
# donor's cells. With arm means M1, M0 over n1, n0 donors and the sample
# variances V1, V0 of the donor means,
#   lfc = log2(M1 / M0),
#   se  = sqrt(V1 / (n1 M1^2) + V0 / (n0 M0^2)) / ln 2.
# The empirical null (mean d0, SD s0) is fitted to z = lfc / se over the
# genes, and the rule reports (lfc - d0 * se, s0 * se).

.donor_means <- function(mat, individual_vec){
  individual_chr <- as.character(individual_vec)
  num_cells_vec <- table(individual_chr)
  sum_mat <- rowsum(mat, group = individual_chr)
  sum_mat / as.numeric(num_cells_vec[rownames(sum_mat)])
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
  list(lfc_vec = log2(case_mean / control_mean),
       se_vec = sqrt(case_var / (n1 * case_mean^2) +
                       control_var / (n0 * control_mean^2)) / log(2))
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
