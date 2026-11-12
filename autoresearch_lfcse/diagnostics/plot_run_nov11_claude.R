# Figures of the autoresearch run `autoresearch_lfcse/nov11`.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11.
#
# Reads the run's tables under autoresearch_lfcse/output/ and writes four
# PNGs to autoresearch_lfcse/notes/figures/RUN_nov11/. Run from the
# repository root:
#   Rscript autoresearch_lfcse/diagnostics/plot_run_nov11_claude.R

rm(list = ls())
library(ggplot2)

# Paths ------------------------------------------------------------------------

output_dir <- file.path("autoresearch_lfcse", "output")
figure_dir <- file.path("autoresearch_lfcse", "notes", "figures", "RUN_nov11")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# The two kept rows later withdrawn as invalid (their SE follows the
# estimated logFC); see RUN_nov11_claude.md.
withdrawn_vec <- c("da743dc", "b3fc556")
factor_vec <- c("donors", "cells", "rates", "effects", "covariates",
                "donor_sd")

# Tables -----------------------------------------------------------------------

results_df <- utils::read.delim(file.path(output_dir, "results.tsv"),
                                quote = "", comment.char = "",
                                stringsAsFactors = FALSE)
score_list <- list(
  dev = list(baseline = utils::read.csv(file.path(output_dir,
                                                  "baseline_dev_scores.csv")),
             best = utils::read.csv(file.path(output_dir,
                                              "best_scores.csv"))),
  test = list(baseline = utils::read.csv(file.path(output_dir,
                                                   "baseline_test_scores.csv")),
              best = utils::read.csv(file.path(output_dir,
                                               "test_scores.csv")))
)
for(split in names(score_list)){
  stopifnot(all(score_list[[split]]$baseline$id == score_list[[split]]$best$id))
}

# Figure 1: trajectory ---------------------------------------------------------

results_df$experiment <- seq_len(nrow(results_df)) - 1
results_df$outcome <- ifelse(results_df$status == "keep", "kept", "discarded")
results_df$outcome[results_df$commit %in% withdrawn_vec] <-
  "kept, then withdrawn"
results_df$outcome <- factor(results_df$outcome,
                             levels = c("kept", "discarded",
                                        "kept, then withdrawn"))
baseline_score <- results_df$score[1]
valid_idx <- which(results_df$outcome == "kept")
best_score <- results_df$score[valid_idx[length(valid_idx)]]

plot1 <- ggplot2::ggplot(results_df, ggplot2::aes(x = experiment, y = score))
plot1 <- plot1 + ggplot2::geom_hline(yintercept = baseline_score,
                                     linetype = "dashed", color = "grey50")
plot1 <- plot1 + ggplot2::geom_hline(yintercept = best_score,
                                     linetype = "dashed", color = "#0072B2")
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(shape = outcome,
                                                  color = outcome), size = 2.5)
plot1 <- plot1 + ggplot2::scale_shape_manual(values = c(16, 1, 4))
plot1 <- plot1 + ggplot2::scale_color_manual(values = c("#0072B2", "grey40",
                                                        "#D55E00"))
plot1 <- plot1 + ggplot2::labs(
  x = "Experiment (0 is the baseline)",
  y = "Dev score (NLL, nats per gene; lower is better)",
  title = "Dev score of every experiment",
  subtitle = "Dashed lines: baseline (grey) and final best (blue)"
)
plot1 <- plot1 + ggplot2::theme_bw()
ggplot2::ggsave(file.path(figure_dir, "trajectory.png"), plot1,
                width = 8, height = 4.5, dpi = 150)

# Figure 2: strata -------------------------------------------------------------

strata_list <- list()
for(split in names(score_list)){
  for(rule in c("baseline", "best")){
    score_df <- score_list[[split]][[rule]]
    for(factor_name in factor_vec){
      level_vec <- as.character(score_df[[factor_name]])
      for(level in unique(level_vec)){
        value_vec <- score_df$score[level_vec == level]
        strata_list[[length(strata_list) + 1]] <- data.frame(
          split = split, rule = rule,
          stratum = paste0(factor_name, ": ", level),
          mean = mean(value_vec),
          se = stats::sd(value_vec) / sqrt(length(value_vec)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
strata_df <- do.call(rbind, strata_list)
strata_df$stratum <- factor(strata_df$stratum,
                            levels = rev(unique(strata_df$stratum)))

plot2 <- ggplot2::ggplot(strata_df, ggplot2::aes(x = mean, y = stratum,
                                                 color = rule))
plot2 <- plot2 + ggplot2::geom_errorbar(
  ggplot2::aes(xmin = mean - se, xmax = mean + se), width = 0.3,
  orientation = "y", position = ggplot2::position_dodge(width = 0.5)
)
plot2 <- plot2 + ggplot2::geom_point(
  position = ggplot2::position_dodge(width = 0.5)
)
plot2 <- plot2 + ggplot2::facet_wrap(~ split)
plot2 <- plot2 + ggplot2::scale_color_manual(values = c("grey40", "#0072B2"))
plot2 <- plot2 + ggplot2::labs(
  x = "Mean score in the stratum (NLL, nats per gene; lower is better)",
  y = NULL, title = "Score by stratum, baseline against best",
  subtitle = "Bars: one standard error over the stratum's examples"
)
plot2 <- plot2 + ggplot2::theme_bw()
ggplot2::ggsave(file.path(figure_dir, "strata.png"), plot2,
                width = 9, height = 5.5, dpi = 150)

# Figure 3: paired differences on test -----------------------------------------

paired_df <- score_list$test$best[, c("id", factor_vec)]
paired_df$difference <- score_list$test$best$score -
  score_list$test$baseline$score
paired_df <- paired_df[order(paired_df$difference), ]
paired_df$rank <- seq_len(nrow(paired_df))
paired_df$stratum <- paste0(paired_df$donors, " donors, ",
                            paired_df$covariates)

plot3 <- ggplot2::ggplot(paired_df, ggplot2::aes(x = rank, y = difference,
                                                 color = stratum))
plot3 <- plot3 + ggplot2::geom_hline(yintercept = 0)
plot3 <- plot3 + ggplot2::geom_point(size = 1.8)
plot3 <- plot3 + ggplot2::labs(
  x = "Test example, sorted by the difference",
  y = "Best minus baseline (nats per gene; below 0 is better)",
  title = "Per-example change on the test split",
  color = "Donors, covariates"
)
plot3 <- plot3 + ggplot2::theme_bw()
ggplot2::ggsave(file.path(figure_dir, "paired.png"), plot3,
                width = 8, height = 4.5, dpi = 150)

# Figure 4: what the rule reports ----------------------------------------------

chosen_df <- score_list$test$best[, c("id", factor_vec)]
chosen_df$baseline_se <- score_list$test$baseline$median_se
chosen_df$best_se <- score_list$test$best$median_se
chosen_df$donors <- factor(chosen_df$donors)
chosen_df$panel <- paste0("donor SD ", chosen_df$donor_sd, ", covariates: ",
                          chosen_df$covariates)

plot4 <- ggplot2::ggplot(chosen_df, ggplot2::aes(x = baseline_se, y = best_se,
                                                 color = donors))
plot4 <- plot4 + ggplot2::geom_abline(slope = 1, intercept = 0,
                                      color = "grey50")
plot4 <- plot4 + ggplot2::geom_point(size = 1.8)
plot4 <- plot4 + ggplot2::facet_wrap(~ panel)
plot4 <- plot4 + ggplot2::labs(
  x = "Median reported SE, baseline (log2 units)",
  y = "Median reported SE, best (log2 units)",
  title = "The median standard error the rule reports, per test example",
  subtitle = "Grey line: equal"
)
plot4 <- plot4 + ggplot2::theme_bw()
ggplot2::ggsave(file.path(figure_dir, "chosen.png"), plot4,
                width = 8, height = 6, dpi = 150)

# Numbers the report quotes ----------------------------------------------------

for(split in names(score_list)){
  difference_vec <- score_list[[split]]$best$score -
    score_list[[split]]$baseline$score
  print(paste0(split, ": gain ", round(mean(difference_vec), 4),
               ", paired SE ",
               round(stats::sd(difference_vec) / sqrt(length(difference_vec)),
                     4),
               ", examples improved ", sum(difference_vec < 0), " of ",
               length(difference_vec)))
}
