# Reference candidate `null_rescaled`: the stored pair shifted and rescaled
# by the empirical null that 1.2.0 fitted to its Gaussianized statistics.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11. A drop-in body for
# `report_lfc_se()` (same contract as ../report_rule_claude.R), scored once
# at the build for the reference rows of notes/BANK.md.

report_lfc_se <- function(view){
  null_mean <- attr(view$stored_df, "null_mean")
  null_sd <- attr(view$stored_df, "null_sd")
  data.frame(gene = colnames(view$dat),
             log2fc = view$stored_df$log2fc -
               null_mean * view$stored_df$log2fc_se,
             se = null_sd * view$stored_df$log2fc_se,
             stringsAsFactors = FALSE)
}
