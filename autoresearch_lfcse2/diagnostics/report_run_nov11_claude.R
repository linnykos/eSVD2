# Tables and figures of the two runs of the second lfcse loop, `nov11`
# (truth: marginal) and `nov11-direct` (truth: direct), for their reports.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11.
#
# Reads the runs' tables under autoresearch_lfcse2/output/run_marginal/ and
# run_direct/ (not tracked) and writes tracked copies to
# autoresearch_lfcse2/notes/csv/RUN_nov11/ and eight PNGs to
# autoresearch_lfcse2/notes/figures/RUN_nov11/. Run from the repository
# root, on the machine that ran the loop:
#   Rscript autoresearch_lfcse2/diagnostics/report_run_nov11_claude.R

rm(list = ls())
library(ggplot2)

# Paths ------------------------------------------------------------------------

output_dir <- file.path("autoresearch_lfcse2", "output")
csv_dir <- file.path("autoresearch_lfcse2", "notes", "csv", "RUN_nov11")
figure_dir <- file.path("autoresearch_lfcse2", "notes", "figures",
                        "RUN_nov11")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# The folder of each run, and the name of the rule it kept.
run_dir_vec <- c(marginal = "run_marginal", direct = "run_direct")
best_rule_vec <- c(marginal = "embed_prior", direct = "state_adjusted")
factor_vec <- c("donors", "cells", "expression", "effects", "covariates",
                "donor_sd", "state", "arm_sd", "depth")
rule_color_vec <- c(null_rescaled = "grey40", donor_wls = "#D55E00",
                    embed_prior = "#0072B2", state_adjusted = "#0072B2")

# Tables -----------------------------------------------------------------------

.read_block <- function(file){
  line_vec <- readLines(file, warn = FALSE)
  line_vec <- line_vec[grepl("^[A-Za-z0-9_.]+:", line_vec)]
  data.frame(key = sub(":.*$", "", line_vec),
             value = trimws(sub("^[^:]+:", "", line_vec)),
             stringsAsFactors = FALSE)
}

result_list <- list()
score_list <- list()
block_list <- list()
for(run in names(run_dir_vec)){
  run_dir <- file.path(output_dir, run_dir_vec[[run]])
  results_df <- utils::read.delim(file.path(run_dir, "results.tsv"),
                                  quote = "", comment.char = "",
                                  stringsAsFactors = FALSE)
  results_df$run <- run
  results_df$experiment <- seq_len(nrow(results_df)) - 1
  # The change against the best commit at the time (negative is better).
  results_df$delta <- NA_real_
  best_score <- results_df$score[1]
  for(i in seq_len(nrow(results_df))[-1]){
    results_df$delta[i] <- results_df$score[i] - best_score
    if(results_df$status[i] == "keep") best_score <- results_df$score[i]
  }
  result_list[[run]] <- results_df

  file_list <- list(
    list(split = "dev", rule = "null_rescaled",
         scores = "baseline_dev_scores.csv", block = "baseline_dev_block.txt"),
    list(split = "test", rule = "null_rescaled",
         scores = "baseline_test_scores.csv",
         block = "baseline_test_block.txt"),
    list(split = "dev", rule = best_rule_vec[[run]],
         scores = "best_scores.csv", block = "best_block.txt"),
    list(split = "test", rule = best_rule_vec[[run]],
         scores = "test_scores.csv", block = "run_test.log")
  )
  for(entry in file_list){
    score_df <- utils::read.csv(file.path(run_dir, entry$scores),
                                stringsAsFactors = FALSE)
    stopifnot(all(score_df$split == entry$split),
              all(score_df$status == "ok"))
    score_df$run <- run
    score_df$rule <- entry$rule
    score_df <- score_df[, setdiff(colnames(score_df),
                                   c("error", "bank", "seconds"))]
    key <- paste(run, entry$split, entry$rule)
    score_list[[key]] <- score_df
    # The stored `donor_wls` reference of the same examples, once per split:
    # its score is the rule's score less the stored difference.
    if(entry$rule != "null_rescaled"){
      wls_df <- score_df
      wls_df$rule <- "donor_wls"
      wls_df$score <- score_df$score - score_df$vs_wls
      keep_vec <- c("id", "split", "cell_idx", factor_vec, "status", "run",
                    "rule", "score")
      wls_df[, setdiff(colnames(wls_df), keep_vec)] <- NA
      score_list[[paste(run, entry$split, "donor_wls")]] <- wls_df
    }
    block_df <- .read_block(file.path(run_dir, entry$block))
    block_df$run <- run
    block_df$split <- entry$split
    block_df$rule <- entry$rule
    block_list[[key]] <- block_df
  }
}
results_df <- do.call(rbind, result_list)
# The marginal run's baseline was scored on test before the scorer printed
# the score against both truths, so its table lacks those two columns.
column_vec <- unique(unlist(lapply(score_list, colnames)))
score_list <- lapply(score_list, function(sub_df){
  sub_df[setdiff(column_vec, colnames(sub_df))] <- NA
  sub_df[, column_vec, drop = FALSE]
})
score_df <- do.call(rbind, score_list)
block_df <- do.call(rbind, block_list)
rownames(results_df) <- NULL
rownames(score_df) <- NULL
rownames(block_df) <- NULL
utils::write.csv(results_df, file.path(csv_dir, "results.csv"),
                 row.names = FALSE)
utils::write.csv(score_df, file.path(csv_dir, "scores.csv"),
                 row.names = FALSE)
utils::write.csv(block_df, file.path(csv_dir, "blocks.csv"),
                 row.names = FALSE)
utils::write.csv(utils::read.csv(file.path(output_dir, "diagnostics",
                                           "oracle_dev_scores.csv")),
                 file.path(csv_dir, "oracle_dev.csv"), row.names = FALSE)
utils::write.csv(utils::read.csv(file.path(output_dir, "diagnostics",
                                           "bank_checks.csv")),
                 file.path(csv_dir, "bank_checks.csv"), row.names = FALSE)

# Figures ----------------------------------------------------------------------

for(run in names(run_dir_vec)){
  best_rule <- best_rule_vec[[run]]
  truth_label <- paste0("truth: ", run)

  ## Figure: trajectory
  run_df <- results_df[results_df$run == run, ]
  run_df$outcome <- factor(ifelse(run_df$status == "keep", "kept",
                                  "discarded"),
                           levels = c("kept", "discarded"))
  kept_idx <- which(run_df$status == "keep")
  plot1 <- ggplot2::ggplot(run_df, ggplot2::aes(x = experiment, y = score))
  plot1 <- plot1 + ggplot2::geom_hline(yintercept = run_df$score[1],
                                       linetype = "dashed", color = "grey50")
  plot1 <- plot1 + ggplot2::geom_hline(
    yintercept = run_df$score[kept_idx[length(kept_idx)]],
    linetype = "dashed", color = "#0072B2"
  )
  plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(shape = outcome,
                                                    color = outcome),
                                       size = 2.5)
  plot1 <- plot1 + ggplot2::scale_shape_manual(values = c(16, 1))
  plot1 <- plot1 + ggplot2::scale_color_manual(values = c("#0072B2",
                                                          "grey40"))
  plot1 <- plot1 + ggplot2::scale_x_continuous(
    breaks = seq(0, max(run_df$experiment))
  )
  plot1 <- plot1 + ggplot2::labs(
    x = "Experiment (0 is the baseline)",
    y = "Dev score (NLL, nats per gene; lower is better)",
    title = paste0("Dev score of every experiment (", truth_label, ")"),
    subtitle = "Dashed lines: baseline (grey) and final best (blue)"
  )
  plot1 <- plot1 + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figure_dir, paste0(run, "_trajectory.png")),
                  plot1, width = 8, height = 4.5, dpi = 150)

  ## Figure: strata
  strata_list <- list()
  for(split in c("dev", "test")){
    for(rule in c("null_rescaled", "donor_wls", best_rule)){
      sub_df <- score_df[score_df$run == run & score_df$split == split &
                           score_df$rule == rule, ]
      for(factor_name in factor_vec){
        level_vec <- as.character(sub_df[[factor_name]])
        for(level in unique(level_vec)){
          value_vec <- sub_df$score[level_vec == level]
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
  strata_df$rule <- factor(strata_df$rule,
                           levels = c("null_rescaled", "donor_wls",
                                      best_rule))
  plot2 <- ggplot2::ggplot(strata_df, ggplot2::aes(x = mean, y = stratum,
                                                   color = rule))
  plot2 <- plot2 + ggplot2::geom_errorbar(
    ggplot2::aes(xmin = mean - se, xmax = mean + se), width = 0.3,
    orientation = "y", position = ggplot2::position_dodge(width = 0.6)
  )
  plot2 <- plot2 + ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.6)
  )
  plot2 <- plot2 + ggplot2::facet_wrap(~ split)
  plot2 <- plot2 + ggplot2::scale_color_manual(
    values = unname(rule_color_vec[levels(strata_df$rule)])
  )
  plot2 <- plot2 + ggplot2::labs(
    x = "Mean score in the stratum (NLL, nats per gene; lower is better)",
    y = NULL,
    title = paste0("Score by stratum (", truth_label, ")"),
    subtitle = "Bars: one standard error over the stratum's examples"
  )
  plot2 <- plot2 + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figure_dir, paste0(run, "_strata.png")), plot2,
                  width = 9.5, height = 7, dpi = 150)

  ## Figure: paired differences on test
  test_best_df <- score_df[score_df$run == run & score_df$split == "test" &
                             score_df$rule == best_rule, ]
  test_base_df <- score_df[score_df$run == run & score_df$split == "test" &
                             score_df$rule == "null_rescaled", ]
  stopifnot(all(test_best_df$id == test_base_df$id))
  paired_df <- test_best_df[, c("id", factor_vec)]
  paired_df$difference <- test_best_df$score - test_base_df$score
  paired_df <- paired_df[order(paired_df$difference), ]
  paired_df$rank <- seq_len(nrow(paired_df))
  plot3 <- ggplot2::ggplot(paired_df, ggplot2::aes(x = rank, y = difference,
                                                   color = state))
  plot3 <- plot3 + ggplot2::geom_hline(yintercept = 0)
  plot3 <- plot3 + ggplot2::geom_point(size = 1.6)
  plot3 <- plot3 + ggplot2::scale_color_manual(values = c("#0072B2",
                                                          "#D55E00",
                                                          "grey40"))
  plot3 <- plot3 + ggplot2::labs(
    x = "Test example, sorted by the difference",
    y = "Best minus baseline (nats per gene; below 0 is better)",
    title = paste0("Per-example change on the test split (", truth_label,
                   ")"),
    color = "state"
  )
  plot3 <- plot3 + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figure_dir, paste0(run, "_paired.png")), plot3,
                  width = 8, height = 4.5, dpi = 150)

  ## Figure: what the rule reports
  chosen_df <- test_best_df[, c("id", factor_vec)]
  chosen_df$baseline_se <- test_base_df$median_se
  chosen_df$best_se <- test_best_df$median_se
  chosen_df$donors <- factor(chosen_df$donors)
  chosen_df$panel <- paste0("state: ", chosen_df$state, ", donor SD ",
                            chosen_df$donor_sd)
  plot4 <- ggplot2::ggplot(chosen_df, ggplot2::aes(x = baseline_se,
                                                   y = best_se,
                                                   color = donors))
  plot4 <- plot4 + ggplot2::geom_abline(slope = 1, intercept = 0,
                                        color = "grey50")
  plot4 <- plot4 + ggplot2::geom_point(size = 1.5)
  plot4 <- plot4 + ggplot2::facet_wrap(~ panel)
  plot4 <- plot4 + ggplot2::scale_color_manual(values = c("#D55E00",
                                                          "#009E73",
                                                          "#0072B2"))
  plot4 <- plot4 + ggplot2::labs(
    x = "Median reported SE, baseline (log2 units)",
    y = "Median reported SE, best (log2 units)",
    title = paste0("The median standard error per test example (",
                   truth_label, ")"),
    subtitle = "Grey line: equal"
  )
  plot4 <- plot4 + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figure_dir, paste0(run, "_chosen.png")), plot4,
                  width = 8.5, height = 6, dpi = 150)
}

# Numbers the reports quote ----------------------------------------------------

for(run in names(run_dir_vec)){
  for(split in c("dev", "test")){
    best_df <- score_df[score_df$run == run & score_df$split == split &
                          score_df$rule == best_rule_vec[[run]], ]
    for(reference in c("null_rescaled", "donor_wls")){
      reference_df <- score_df[score_df$run == run & score_df$split == split &
                                 score_df$rule == reference, ]
      difference_vec <- best_df$score -
        reference_df$score[match(best_df$id, reference_df$id)]
      print(paste0(run, ", ", split, ", against ", reference, ": ",
                   round(mean(difference_vec), 4), ", paired SE ",
                   round(stats::sd(difference_vec) /
                           sqrt(length(difference_vec)), 4),
                   ", examples improved ", sum(difference_vec < 0), " of ",
                   length(difference_vec)))
    }
  }
}
print("Done! :)")
