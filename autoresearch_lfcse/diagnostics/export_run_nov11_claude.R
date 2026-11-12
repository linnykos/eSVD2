# Tables of the autoresearch run `autoresearch_lfcse/nov11`, for its writeup.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11.
#
# Copies the run's tables out of autoresearch_lfcse/output/ (not tracked)
# into autoresearch_lfcse/notes/csv/RUN_nov11/ (tracked), so that
# notes/RUN_nov11_claude.Rmd knits from a clone. Run from the repository
# root, on the machine that ran the loop:
#   Rscript autoresearch_lfcse/diagnostics/export_run_nov11_claude.R

rm(list = ls())

# Paths ------------------------------------------------------------------------

output_dir <- file.path("autoresearch_lfcse", "output")
csv_dir <- file.path("autoresearch_lfcse", "notes", "csv", "RUN_nov11")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# The two kept rows later withdrawn as invalid (their SE follows the
# estimated logFC), and the kept commit the branch was reset to.
withdrawn_vec <- c("da743dc", "b3fc556")
reset_commit <- "1de04e6"

# The log of experiments -------------------------------------------------------

results_df <- utils::read.delim(file.path(output_dir, "results.tsv"),
                                quote = "", comment.char = "",
                                stringsAsFactors = FALSE)
results_df$experiment <- seq_len(nrow(results_df)) - 1

# The change against the best commit at the time (negative is better). The
# experiment after the withdrawn pair was judged against the withdrawn best;
# from the reset onward the reference is the commit the branch was reset to.
results_df$delta <- NA_real_
best_score <- results_df$score[1]
for(i in seq_len(nrow(results_df))[-1]){
  results_df$delta[i] <- results_df$score[i] - best_score
  if(results_df$status[i] == "keep") best_score <- results_df$score[i]
  if(i > 1 && results_df$commit[i - 1] == withdrawn_vec[2]){
    best_score <- results_df$score[results_df$commit == reset_commit]
  }
}
results_df$outcome <- ifelse(results_df$status == "keep", "kept", "discarded")
results_df$outcome[results_df$commit %in% withdrawn_vec] <- "withdrawn"
results_df$outcome[1] <- "baseline"

# What the log's free-text descriptions say of each experiment, in columns:
# the memo's direction, the idea, the paired SE of `delta`, and the strata
# that moved. The paired SE of the withdrawn pair was not kept.
experiment_text <- "
experiment;commit;direction;idea;paired_se;strata
1;0afa880;D4;Drop the null-SD rescaling (before moderation);0.0103;10 donors 0.397 to 0.545
2;8c9ee26;D3;Shift by the constant $\\delta_0 \\cdot \\text{median}(s)$ instead of $\\delta_0 s_g$;0.0003;`sparse_up` 0.0054 to 0.0021
3;874fbb3;D2;Widen the SE by $t_{0.975,\\text{df}} / z_{0.975}$, Welch df;0.0002;absorbed by the null SD then in the rule
4;a9d2d81;D1;Donor summary as the ratio of sums;0.0004;20 and 40 donors, no covariates
5;5029a75;D4;Moderate each arm's squared CV toward a trend, prior df from the genes;0.0032;10 donors 0.397 to 0.305, every stratum
6;3df2a64;D4;Drop the null-SD rescaling (after moderation, before the design factor);0.0059;`sexage` at 10 donors, coverage 0.86
7;47fda43;D2;Replace the null-SD rescaling by the design's variance factor for donor-level covariates;0.0032;`sexage` -0.010 to -0.024, `dense` -0.010 to -0.055, 40 donors -0.373 to -0.402
8;99ace1a;D5;Add the stored within-donor term;0.0002;coverage 0.944 to 0.951
9;e0f5012;D2;Between term from the variance of log donor means;0.0002;donor SD 0.5
10;eea00eb;D1;Donors weighted by cell count;0.0054;unequal cells -0.073 to +0.040, RMSE 0.310
11;938139f;D5;The fit's case-control coefficient as the logFC;0.0063;unequal cells, RMSE 0.313
12;3b6c49c;D1;logFC from the mean log donor means;0.0017;RMSE 0.2934 to 0.2878, donor SD 0.5 0.369 to 0.341
13;d60190e;D1;logFC as the arm coefficient of the donor-level regression on donor-level covariates;0.0007;`sexage` -0.032 to -0.040, RMSE 0.2856
14;7c8fd08;D2;Residual variance pooled over the arms, moderated once (tie, 33 lines shorter);0.0002;none
15;b66d26b;D5;Drop the within-donor term again (20 lines shorter);0.0001;coverage 0.957 to 0.951
16;af1e2f2;D3;Centre at the median logFC;0.0002;`sparse_up`
17;d28cf26;D1;Precision weights from voom's mean-variance trend;0.0002;donor SD 0.2 -0.480 to -0.483
18;cca12c8;D5;$\\exp(x^\\top y)$ in the donor offset;0.0043;10 donors, donor SD 0.5, coverage 0.933
19;9e012ce;D3;No centring;0.0056;`sparse_up` -0.070 to +0.022, every effects level
20;a95c420;D3;Harmonic-mean SE in the shift;0.0000;none, no simpler
21;499cb2d;D3;Null mean by central matching;0.0015;every effects level
22;b7f94d1;D3;Biweight centre at unit scale;0.0001;`sparse_up`
23;bf5f87c;D4;Moderation toward one common value;0.0005;every stratum slightly
24;8f98e69;D5;Fitted Gamma rate as a second trend covariate;0.0006;every stratum
25;c346eca;D1;Half count scaled by the donor's library;0.0001;unequal cells
26;5014764;D4;Moderation trend against the expected count;0.0005;equal to no trend
27;1de04e6;D2;Remove each donor's offset common to all genes;0.0001;none individually
28;da743dc;D4;Moderation trend against the arithmetic-geometric gap of the donor ratios;NA;every stratum (withdrawn)
29;b3fc556;D4;Variance as a direct smooth of that gap;NA;none individually (withdrawn)
30;0982011;D4;Gap covariate computed on residuals;0.0059;10 donors 0.236 to 0.411 (against the withdrawn best)
31;bf6e2d4;D4;Arm-invariant covariate for the moderation trend;0.0001;none individually
32;4538139;D2;Donor-level covariates as random effects (ridge);0.0023;`sexage` -0.043 to -0.020, coverage 0.941
33;5b2cfc5;D1;Second pass of the precision weights;0.0002;none
34;80cdb6e;D1;Weights from $\\tau^2 + \\phi / \\mu$ in place of the voom spline;0.0002;unequal cells -0.0987 to -0.1007, RMSE 0.2855
35;e4fb126;D5;Gamma-Poisson sampling term in the weights;0.0002;`slope_weak` 0.989
36;d596ad6;D2;Drop the per-donor offset;0.0001;none individually
37;5df592c;D4;Moderation trend against the mean fitted log ratio;0.0004;`sparse_up` -0.0711 to -0.0766, 10 donors 0.2831 to 0.2800
38;c9a9547;D1;Reweighted second fit of $\\tau^2 + \\phi/\\mu$;0.0000;none, longer
"
experiment_df <- utils::read.table(text = experiment_text, sep = ";",
                                   header = TRUE, quote = "",
                                   comment.char = "", strip.white = TRUE,
                                   stringsAsFactors = FALSE)
# The text must describe the log's rows, in the log's order.
stopifnot(nrow(experiment_df) == nrow(results_df) - 1,
          all(experiment_df$commit == results_df$commit[-1]),
          all(experiment_df$experiment == results_df$experiment[-1]))
experiment_df$time <- results_df$time[-1]
experiment_df$score <- results_df$score[-1]
experiment_df$delta <- results_df$delta[-1]
experiment_df$status <- results_df$status[-1]
experiment_df$outcome <- results_df$outcome[-1]

utils::write.csv(results_df[, setdiff(colnames(results_df), "description")],
                 file.path(csv_dir, "results.csv"), row.names = FALSE)
utils::write.csv(experiment_df, file.path(csv_dir, "experiments.csv"),
                 row.names = FALSE)

# The validity check of the withdrawn commits ----------------------------------

# Computed during the run by a throwaway script on 12 dev cohorts and
# recorded here by hand, because the script was not kept: the ratio of the
# median SE of strong-effect genes to that of the other genes, and the
# Spearman correlation between the SE and the absolute reported logFC.
validity_df <- data.frame(
  rule = c("gap_trend", "before the gap covariate", "donor_wls"),
  commit = c("b3fc556", "1de04e6", "5df592c"),
  num_cohorts = 12,
  se_ratio = c(1.10, 1.00, 0.98),
  spearman = c(0.21, 0.03, 0.01),
  stringsAsFactors = FALSE
)
utils::write.csv(validity_df, file.path(csv_dir, "validity.csv"),
                 row.names = FALSE)

# Per-example scores -----------------------------------------------------------

score_file_list <- list(
  list(split = "dev", rule = "raw_null", file = "baseline_dev_scores.csv"),
  list(split = "test", rule = "raw_null", file = "baseline_test_scores.csv"),
  list(split = "dev", rule = "donor_wls", file = "best_scores.csv"),
  list(split = "test", rule = "donor_wls", file = "test_scores.csv")
)
score_list <- lapply(score_file_list, function(entry){
  score_df <- utils::read.csv(file.path(output_dir, entry$file),
                              stringsAsFactors = FALSE)
  stopifnot(all(score_df$split == entry$split),
            all(score_df$status == "ok"))
  score_df$rule <- entry$rule
  score_df[, setdiff(colnames(score_df), c("error", "bank", "seconds"))]
})
score_df <- do.call(rbind, score_list)
utils::write.csv(score_df, file.path(csv_dir, "scores.csv"),
                 row.names = FALSE)

# Summary blocks ---------------------------------------------------------------

block_file_list <- list(
  list(split = "dev", rule = "raw_null", file = "baseline_dev_block.txt"),
  list(split = "test", rule = "raw_null", file = "baseline_test_block.txt"),
  list(split = "dev", rule = "donor_wls", file = "best_block.txt"),
  list(split = "test", rule = "donor_wls", file = "run_test.log")
)
block_list <- lapply(block_file_list, function(entry){
  line_vec <- readLines(file.path(output_dir, entry$file), warn = FALSE)
  line_vec <- line_vec[grepl("^[A-Za-z0-9_.]+:", line_vec)]
  data.frame(split = entry$split, rule = entry$rule,
             key = sub(":.*$", "", line_vec),
             value = trimws(sub("^[^:]+:", "", line_vec)),
             stringsAsFactors = FALSE)
})
block_df <- do.call(rbind, block_list)
# The best commit's blocks must be the best commit's.
best_commit <- results_df$commit[max(which(results_df$outcome == "kept"))]
stopifnot(all(block_df$value[block_df$key == "commit" &
                               block_df$rule == "donor_wls"] == best_commit))
utils::write.csv(block_df, file.path(csv_dir, "blocks.csv"),
                 row.names = FALSE)

print(paste0("Wrote ", nrow(results_df), " result rows, ", nrow(score_df),
             " example scores and ", nrow(block_df), " block lines to ",
             csv_dir))
