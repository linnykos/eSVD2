# Reference candidate `stored`: what eSVD2 1.2.0 reports.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11. A drop-in body for
# `report_lfc_se()` (same contract as ../report_rule_claude.R), scored once
# at the build for the reference rows of notes/BANK.md.

report_lfc_se <- function(view){
  data.frame(gene = colnames(view$dat),
             log2fc = view$stored_df$log2fc,
             se = view$stored_df$log2fc_se,
             stringsAsFactors = FALSE)
}
