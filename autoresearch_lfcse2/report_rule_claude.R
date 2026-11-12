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
# Baseline `null_rescaled`: what eSVD2 1.2.0 stores, shifted and rescaled
# by its empirical null. 1.2.0 averages the posterior mean and variance of
# each cell within a donor, takes the mean M and the mixture variance V of
# the donor averages in each arm, and stores
#   lfc = log2(M_1 / M_0),
#   se  = sqrt(V_1 / (n_1 M_1^2) + V_0 / (n_0 M_0^2)) / ln 2.
# The posterior mean of a cell is (A + mu * gamma) / (library + gamma),
# with mu = exp(x' y + z_cc * cc), so the embedding enters both numbers.
# With (d0, s0) the mean and SD of the empirical null that 1.2.0 fitted to
# its Gaussianized statistics, the rule reports (lfc - d0 * se, s0 * se).

report_lfc_se <- function(view){
  null_mean <- attr(view$stored_df, "null_mean")
  null_sd <- attr(view$stored_df, "null_sd")
  data.frame(gene = colnames(view$dat),
             log2fc = view$stored_df$log2fc -
               null_mean * view$stored_df$log2fc_se,
             se = null_sd * view$stored_df$log2fc_se,
             stringsAsFactors = FALSE)
}
