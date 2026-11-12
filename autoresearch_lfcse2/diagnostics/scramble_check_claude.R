# How much the kept rules depend on the embedding: each rule on a cohort's
# view, and on the same view with the rows of the embedding scores permuted
# across cells (the posterior recomputed from the permuted fit), and with
# the rows of the loadings permuted across genes as well.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11.
#
# The memo (section 1) requires the logFC or the SE to use the embedding;
# this is the check it names. The script reads the truth to score both
# versions, so none of its code may enter the trial file. Run from the
# repository root:
#   Rscript autoresearch_lfcse2/diagnostics/scramble_check_claude.R
# Writes autoresearch_lfcse2/notes/csv/RUN_nov11/scramble.csv, one row per
# rule and example (every fourth dev example).

rm(list = ls())

# Paths ------------------------------------------------------------------------

repo_dir <- getwd()
loop_dir <- file.path(repo_dir, "autoresearch_lfcse2")
source(file.path(loop_dir, "harness", "generator_claude.R"))
.load_esvd2(repo_dir)
bank_dir <- .bank_dir()
csv_dir <- file.path(loop_dir, "notes", "csv", "RUN_nov11")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
num_cores <- as.integer(Sys.getenv("ESVD2_CORES", unset = "8"))
rule_file_vec <- c(state_adjusted = "RULE_nov11_direct_claude.R",
                   embed_prior = "RULE_nov11_marginal_claude.R")

# Each rule file alone in its own environment: the files must need nothing
# else of the project.
rule_list <- lapply(rule_file_vec, function(rule_file){
  rule_env <- new.env()
  sys.source(file.path(loop_dir, "notes", rule_file), envir = rule_env)
  rule_env$report_lfc_se
})

.nll <- function(value, truth_vec){
  mean(0.5 * log(2 * pi * value$se^2) +
         (value$log2fc - truth_vec)^2 / (2 * value$se^2))
}

# The check ---------------------------------------------------------------------

manifest_df <- utils::read.csv(file.path(bank_dir, "manifest_dev.csv"),
                               stringsAsFactors = FALSE)
manifest_df <- manifest_df[seq(1, nrow(manifest_df), by = 4), ]

row_list <- parallel::mclapply(seq_len(nrow(manifest_df)), function(i){
  slim <- readRDS(file.path(bank_dir, manifest_df$file_public[i]))
  truth <- readRDS(file.path(bank_dir, manifest_df$file_truth[i]))
  direct_vec <- as.numeric(truth$population$cc_effect_vec) / log(2)
  marginal_vec <- as.numeric(truth$truth_log2fc_vec)
  view <- .build_view(slim)
  # The same fit with each cell given another cell's embedding scores.
  set.seed(manifest_df$seed[i] + 60)
  scrambled_slim <- slim
  scrambled_slim$x_mat <- slim$x_mat[sample(nrow(slim$x_mat)), , drop = FALSE]
  rownames(scrambled_slim$x_mat) <- rownames(slim$x_mat)
  scrambled_view <- .build_view(scrambled_slim)
  # And with each gene given another gene's loadings as well.
  both_slim <- scrambled_slim
  both_slim$y_mat <- slim$y_mat[sample(nrow(slim$y_mat)), , drop = FALSE]
  rownames(both_slim$y_mat) <- rownames(slim$y_mat)
  both_view <- .build_view(both_slim)

  do.call(rbind, lapply(names(rule_list), function(rule){
    value <- rule_list[[rule]](view)
    scrambled_value <- rule_list[[rule]](scrambled_view)
    both_value <- rule_list[[rule]](both_view)
    data.frame(
      id = manifest_df$id[i], state = manifest_df$state[i],
      donors = manifest_df$donors[i], rule = rule,
      change_in_se = stats::median(abs(scrambled_value$log2fc -
                                         value$log2fc) / value$se),
      se_ratio = stats::median(scrambled_value$se / value$se),
      score_direct = .nll(value, direct_vec),
      score_direct_scrambled = .nll(scrambled_value, direct_vec),
      score_marginal = .nll(value, marginal_vec),
      score_marginal_scrambled = .nll(scrambled_value, marginal_vec),
      change_in_se_both = stats::median(abs(both_value$log2fc -
                                              value$log2fc) / value$se),
      score_direct_both = .nll(both_value, direct_vec),
      score_marginal_both = .nll(both_value, marginal_vec),
      stringsAsFactors = FALSE
    )
  }))
}, mc.cores = num_cores, mc.preschedule = FALSE)
scramble_df <- do.call(rbind, row_list)
utils::write.csv(scramble_df, file.path(csv_dir, "scramble.csv"),
                 row.names = FALSE)

value_vec <- c("change_in_se", "se_ratio", "score_direct",
               "score_direct_scrambled", "score_marginal",
               "score_marginal_scrambled", "change_in_se_both",
               "score_direct_both", "score_marginal_both")
print(stats::aggregate(scramble_df[, value_vec],
                       by = list(rule = scramble_df$rule,
                                 state = scramble_df$state),
                       FUN = mean), digits = 3)
print("Done! :)")
