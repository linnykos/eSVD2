# Gabitto split-half study, step 4: how each rule replicates between the halves
# Drafted by Claude for Kevin Z. Lin, 2026-11-15

library(ggplot2)

rm(list = ls())

# Reads the four `gabitto_<label>_rules.RData` files of step 3 and, for each
# split and each rule, compares the two halves' log2 fold changes and SEs:
# the replicate score, the correlation and the sign agreement over all genes
# and over the genes that at least one half nominates, and the
# discovery-to-replication statistics. Fits nothing.
#
# Run from the repository root:
#   Rscript additional_context/2026-11-15_gabitto/04_evaluate_claude.R
#
# Writes <out_dir>/gabitto_evaluation.RData and <out_dir>/figures/*.png, and
# the tables the writeup reads to csv/ beside this script.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context", "2026-11-15_gabitto")
stopifnot(dir.exists(study_dir))
source(file.path(study_dir, "paths_claude.R"))
source(file.path(study_dir, "helpers_claude.R"))
csv_dir <- file.path(study_dir, "csv")
if(Sys.getenv("GABITTO_OUT_DIR") != "") csv_dir <- file.path(out_dir, "csv")
figure_dir <- file.path(out_dir, "figures")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

split_vec <- c("split1", "split2")
fdr_level <- 0.05

# Loading ----------------------------------------------------------------------

rules_list <- lapply(fit_label_vec, function(fit_label){
  rules_env <- new.env()
  load(file.path(out_dir, paste0("gabitto_", fit_label, "_rules.RData")),
       envir = rules_env)
  as.list(rules_env)
})
names(rules_list) <- fit_label_vec
method_vec <- colnames(rules_list[[1]]$lfc_mat)
stopifnot(all(sapply(rules_list, function(rules){
  identical(colnames(rules$lfc_mat), method_vec)
})))

# Each rule's own nomination: BH over the half's genes of the Wald test of
# the pair the rule reports.
for(fit_label in fit_label_vec){
  rules <- rules_list[[fit_label]]
  wald_p_mat <- .wald_pvalue(rules$lfc_mat, rules$se_mat)
  wald_p_mat[!is.finite(wald_p_mat)] <- NA_real_
  rules_list[[fit_label]]$wald_p_mat <- wald_p_mat
  rules_list[[fit_label]]$wald_q_mat <- apply(wald_p_mat, 2, stats::p.adjust,
                                              method = "BH")
}

# One row per fit --------------------------------------------------------------

fit_df <- do.call(rbind, lapply(fit_label_vec, function(fit_label){
  rules <- rules_list[[fit_label]]
  data.frame(fit = fit_label,
             num_donors = nrow(rules$donor_summary_df),
             num_case = sum(rules$donor_summary_df$case_control == 1),
             num_control = sum(rules$donor_summary_df$case_control == 0),
             num_cells = rules$num_cells,
             num_genes = nrow(rules$gene_df),
             fit_minutes = round(rules$fit_seconds / 60),
             null_mean = rules$null_mean,
             null_sd = rules$null_sd,
             num_esvd_deg = sum(rules$gene_df$esvd_fdr < fdr_level),
             num_deseq2_deg = sum(rules$gene_df$deseq2_padj < fdr_level,
                                  na.rm = TRUE),
             share_capped = mean(rules$gene_df$status %in%
                                   c("capped", "boundary")),
             num_rule_errors = sum(rules$error_vec != ""),
             stringsAsFactors = FALSE)
}))
print(fit_df)

# One row per split, gene and method -------------------------------------------

pair_df <- do.call(rbind, lapply(split_vec, function(split){
  rules1 <- rules_list[[paste0(split, "_half1")]]
  rules2 <- rules_list[[paste0(split, "_half2")]]
  gene_vec <- intersect(rules1$gene_df$gene, rules2$gene_df$gene)
  idx1 <- match(gene_vec, rules1$gene_df$gene)
  idx2 <- match(gene_vec, rules2$gene_df$gene)
  do.call(rbind, lapply(method_vec, function(method){
    data.frame(
      split = split,
      method = method,
      gene = gene_vec,
      lfc1 = rules1$lfc_mat[idx1, method],
      se1 = rules1$se_mat[idx1, method],
      lfc2 = rules2$lfc_mat[idx2, method],
      se2 = rules2$se_mat[idx2, method],
      wald_p1 = rules1$wald_p_mat[idx1, method],
      wald_p2 = rules2$wald_p_mat[idx2, method],
      own1 = rules1$wald_q_mat[idx1, method] < fdr_level,
      own2 = rules2$wald_q_mat[idx2, method] < fdr_level,
      esvd1 = rules1$gene_df$esvd_fdr[idx1] < fdr_level,
      esvd2 = rules2$gene_df$esvd_fdr[idx2] < fdr_level,
      deseq1 = rules1$gene_df$deseq2_padj[idx1] < fdr_level,
      deseq2 = rules2$gene_df$deseq2_padj[idx2] < fdr_level,
      zero_fraction = (rules1$gene_df$zero_fraction[idx1] +
                         rules2$gene_df$zero_fraction[idx2]) / 2,
      mean_count = (rules1$gene_df$mean_count[idx1] +
                      rules2$gene_df$mean_count[idx2]) / 2,
      capped = rules1$gene_df$status[idx1] != "estimated" |
        rules2$gene_df$status[idx2] != "estimated",
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(pair_df) <- NULL
for(name in c("own1", "own2", "esvd1", "esvd2", "deseq1", "deseq2")){
  pair_df[[name]][is.na(pair_df[[name]])] <- FALSE
}

# Agreement over gene sets -----------------------------------------------------

# The gene sets. `nominated_esvd` and `nominated_deseq2` are the same genes
# for every rule (the package's test, or DESeq2's, calls the gene at 5% FDR
# in at least one half); `nominated_own` uses the rule's own Wald test and
# so differs between rules.
.gene_sets <- function(sub_df){
  list(all = rep(TRUE, nrow(sub_df)),
       nominated_esvd = sub_df$esvd1 | sub_df$esvd2,
       nominated_deseq2 = sub_df$deseq1 | sub_df$deseq2,
       nominated_own = sub_df$own1 | sub_df$own2,
       zeros_below_80 = sub_df$zero_fraction < 0.8,
       zeros_80_to_95 = sub_df$zero_fraction >= 0.8 &
         sub_df$zero_fraction < 0.95,
       zeros_above_95 = sub_df$zero_fraction >= 0.95,
       rate_estimated = !sub_df$capped,
       rate_capped = sub_df$capped)
}

agreement_df <- do.call(rbind, lapply(split_vec, function(split){
  do.call(rbind, lapply(method_vec, function(method){
    sub_df <- pair_df[pair_df$split == split & pair_df$method == method, ]
    set_list <- .gene_sets(sub_df)
    do.call(rbind, lapply(names(set_list), function(gene_set){
      cbind(data.frame(split = split, method = method, gene_set = gene_set,
                       stringsAsFactors = FALSE),
            .agreement_row(sub_df, which(set_list[[gene_set]])))
    }))
  }))
}))

discovery_df <- do.call(rbind, lapply(split_vec, function(split){
  do.call(rbind, lapply(method_vec, function(method){
    sub_df <- pair_df[pair_df$split == split & pair_df$method == method, ]
    nomination_list <- list(esvd = list(sub_df$esvd1, sub_df$esvd2),
                            deseq2 = list(sub_df$deseq1, sub_df$deseq2),
                            own = list(sub_df$own1, sub_df$own2))
    do.call(rbind, lapply(names(nomination_list), function(nomination){
      cbind(data.frame(split = split, method = method, nomination = nomination,
                       stringsAsFactors = FALSE),
            .discovery_row(sub_df, nomination_list[[nomination]][[1]],
                           nomination_list[[nomination]][[2]]))
    }))
  }))
}))

# Which of eSVD-DE's calls replicate ------------------------------------------

# Each gene the package's test calls in one half, with the `stored` logFC
# of that half and of the other. The calls are grouped by whether the gene's
# nuisance rate is at the cap, and by whether the Wald test of
# `donor_wls@esvd_nofloor` in the discovering half has nominal p < 0.05.
calls_df <- do.call(rbind, lapply(split_vec, function(split){
  stored_df <- pair_df[pair_df$split == split & pair_df$method == "stored", ]
  wls_df <- pair_df[pair_df$split == split &
                      pair_df$method == "donor_wls@esvd_nofloor", ]
  stopifnot(identical(stored_df$gene, wls_df$gene))
  long_df <- rbind(
    data.frame(disc = stored_df$lfc1, rep = stored_df$lfc2,
               wls_disc = wls_df$lfc1, wls_rep = wls_df$lfc2,
               wls_p = wls_df$wald_p1, capped = stored_df$capped,
               called = stored_df$esvd1),
    data.frame(disc = stored_df$lfc2, rep = stored_df$lfc1,
               wls_disc = wls_df$lfc2, wls_rep = wls_df$lfc1,
               wls_p = wls_df$wald_p2, capped = stored_df$capped,
               called = stored_df$esvd2)
  )
  long_df <- long_df[long_df$called, ]
  group_list <- list(
    "all calls" = rep(TRUE, nrow(long_df)),
    "rate estimated" = !long_df$capped,
    "rate at the cap" = long_df$capped,
    "donor_wls p < 0.05" = long_df$wls_p < 0.05,
    "donor_wls p >= 0.05" = long_df$wls_p >= 0.05,
    "donor_wls p < 0.01" = long_df$wls_p < 0.01
  )
  do.call(rbind, lapply(names(group_list), function(group){
    sub_df <- long_df[which(group_list[[group]]), ]
    data.frame(split = split, group = group, num_calls = nrow(sub_df),
               sign_rep = mean(sign(sub_df$disc) == sign(sub_df$rep)),
               sign_rep_wls = mean(sign(sub_df$wls_disc) ==
                                     sign(sub_df$wls_rep)),
               slope_rep = sum(sub_df$disc * sub_df$rep) / sum(sub_df$disc^2),
               median_abs_lfc = stats::median(abs(sub_df$disc)),
               stringsAsFactors = FALSE)
  }))
}))
print(calls_df)

# The score against one common target ------------------------------------------

# Each rule's pair in one half scored against the other half's pseudobulk
# DESeq2 estimate, which is independent of the rule's half and the same
# target for every rule: the NLL of (rule's logFC - DESeq2's logFC of the
# other half) under N(0, rule's SE^2 + DESeq2's SE^2), averaged over the
# two directions. A rule cannot gain here by a bias its two halves share.
cross_df <- do.call(rbind, lapply(split_vec, function(split){
  deseq2_df <- pair_df[pair_df$split == split & pair_df$method == "deseq2", ]
  do.call(rbind, lapply(method_vec, function(method){
    sub_df <- pair_df[pair_df$split == split & pair_df$method == method, ]
    stopifnot(identical(sub_df$gene, deseq2_df$gene))
    bool_ok <- is.finite(sub_df$lfc1) & is.finite(sub_df$lfc2) &
      is.finite(sub_df$se1) & is.finite(sub_df$se2) &
      is.finite(deseq2_df$lfc1) & is.finite(deseq2_df$lfc2) &
      is.finite(deseq2_df$se1) & is.finite(deseq2_df$se2)
    nll_vec <- c(.nll_vec(sub_df$lfc1[bool_ok] - deseq2_df$lfc2[bool_ok],
                          sub_df$se1[bool_ok]^2 + deseq2_df$se2[bool_ok]^2),
                 .nll_vec(sub_df$lfc2[bool_ok] - deseq2_df$lfc1[bool_ok],
                          sub_df$se2[bool_ok]^2 + deseq2_df$se1[bool_ok]^2))
    cor_val <- (stats::cor(sub_df$lfc1[bool_ok], deseq2_df$lfc2[bool_ok]) +
                  stats::cor(sub_df$lfc2[bool_ok], deseq2_df$lfc1[bool_ok])) / 2
    data.frame(split = split, method = method, num_genes = sum(bool_ok),
               score_cross = mean(nll_vec), pearson_cross = cor_val,
               stringsAsFactors = FALSE)
  }))
}))

# The summary table of the writeup ---------------------------------------------

.take <- function(df, key_name, key, column, prefix){
  sub_df <- df[df[[key_name]] == key, c("split", "method", column)]
  colnames(sub_df)[-(1:2)] <- paste0(prefix, column)
  sub_df
}
summary_df <- Reduce(function(x, y){merge(x, y, by = c("split", "method"))}, list(
  .take(agreement_df, "gene_set", "all",
        c("num_genes", "score", "z_rms", "coverage", "pearson", "sign_agree"),
        "all_"),
  .take(agreement_df, "gene_set", "nominated_esvd",
        c("num_genes", "score", "pearson", "sign_agree"), "esvd_"),
  .take(agreement_df, "gene_set", "nominated_deseq2",
        c("num_genes", "score", "pearson", "sign_agree"), "deseq2_"),
  .take(agreement_df, "gene_set", "nominated_own",
        c("num_genes", "score", "pearson", "sign_agree"), "own_"),
  .take(discovery_df, "nomination", "esvd",
        c("sign_rep", "slope_rep", "rep_p05"), "esvd_"),
  .take(discovery_df, "nomination", "deseq2",
        c("sign_rep", "slope_rep", "rep_p05"), "deseq2_"),
  .take(discovery_df, "nomination", "own",
        c("num_half1", "num_half2", "num_both", "sign_rep", "slope_rep",
          "rep_p05"), "own_"),
  cross_df[, c("split", "method", "score_cross", "pearson_cross")]
))
summary_df <- summary_df[order(summary_df$split,
                               match(summary_df$method, method_vec)), ]
rownames(summary_df) <- NULL
print(summary_df[, c("split", "method", "all_score", "all_z_rms",
                     "all_pearson", "esvd_pearson", "esvd_sign_rep",
                     "own_num_genes", "own_pearson", "own_sign_rep",
                     "score_cross")], digits = 3)

# Writing ----------------------------------------------------------------------

.round_df <- function(df, digits = 4){
  for(name in colnames(df)){
    if(is.double(df[[name]])) df[[name]] <- signif(df[[name]], digits)
  }
  df
}
utils::write.csv(.round_df(fit_df), file = file.path(csv_dir, "fits.csv"),
                 row.names = FALSE)
utils::write.csv(.round_df(agreement_df),
                 file = file.path(csv_dir, "agreement.csv"), row.names = FALSE)
utils::write.csv(.round_df(discovery_df),
                 file = file.path(csv_dir, "discovery.csv"), row.names = FALSE)
utils::write.csv(.round_df(cross_df),
                 file = file.path(csv_dir, "cross_score.csv"),
                 row.names = FALSE)
utils::write.csv(.round_df(calls_df),
                 file = file.path(csv_dir, "esvd_calls.csv"), row.names = FALSE)
utils::write.csv(.round_df(summary_df),
                 file = file.path(csv_dir, "summary.csv"), row.names = FALSE)
error_df <- do.call(rbind, lapply(fit_label_vec, function(fit_label){
  data.frame(fit = fit_label, method = method_vec,
             error = rules_list[[fit_label]]$error_vec,
             seconds = round(rules_list[[fit_label]]$seconds_vec, 1),
             stringsAsFactors = FALSE)
}))
utils::write.csv(error_df, file = file.path(csv_dir, "rule_runs.csv"),
                 row.names = FALSE)
# The per-gene pairs, so that the writeup knits from a clone.
gz_con <- gzfile(file.path(csv_dir, "gene_pairs.csv.gz"), open = "w")
utils::write.csv(.round_df(pair_df, digits = 5), file = gz_con,
                 row.names = FALSE)
close(gz_con)

date_of_run <- Sys.time()
session_info <- utils::sessionInfo()
note <- paste("Step 4 of the Gabitto split-half study. pair_df: one row per",
              "split, method and gene with the two halves' log2FC and SE and",
              "the nominations (own = BH 5% on the rule's Wald test, esvd =",
              "the package's FDR, deseq = pseudobulk DESeq2's padj).",
              "agreement_df, discovery_df, cross_df, summary_df: the tables",
              "of csv/. fit_df: one row per fit.")
save(pair_df, agreement_df, discovery_df, cross_df, calls_df, summary_df,
     fit_df,
     error_df, method_vec, fdr_level, date_of_run, session_info, note,
     file = file.path(out_dir, "gabitto_evaluation.RData"))

# Figures ----------------------------------------------------------------------

# The panel text of one method under one nomination.
.panel_labels <- function(method, gene_set, nomination){
  do.call(rbind, lapply(split_vec, function(split){
    agreement_row <- agreement_df[agreement_df$split == split &
                                    agreement_df$method == method &
                                    agreement_df$gene_set == gene_set, ]
    discovery_row <- discovery_df[discovery_df$split == split &
                                    discovery_df$method == method &
                                    discovery_df$nomination == nomination, ]
    data.frame(split = split,
               label = paste0(agreement_row$num_genes, " nominated genes\n",
                              "r = ", round(agreement_row$pearson, 2),
                              ", same sign in the other half ",
                              round(100 * discovery_row$sign_rep), "%\n",
                              "replication slope ",
                              round(discovery_row$slope_rep, 2)),
               stringsAsFactors = FALSE)
  }))
}

for(method in method_vec){
  sub_df <- pair_df[pair_df$method == method, ]
  file_stub <- gsub("@", "-", method)

  sub_df$nomination <- .nomination_factor(sub_df$own1, sub_df$own2)
  plot1 <- .plot_replication(
    plot_df = sub_df,
    label_df = .panel_labels(method, "nominated_own", "own"),
    title = paste0(method, ": all genes, coloured by the rule's own Wald ",
                   "test (BH 5%)"),
    bool_only_nominated = FALSE
  )
  ggplot2::ggsave(file.path(figure_dir, paste0("replication_", file_stub,
                                               "_own.png")),
                  plot = plot1, width = 8, height = 4.8, dpi = 200)

  sub_df$nomination <- .nomination_factor(sub_df$esvd1, sub_df$esvd2)
  plot1 <- .plot_replication(
    plot_df = sub_df,
    label_df = .panel_labels(method, "nominated_esvd", "esvd"),
    title = paste0(method, ": genes nominated by eSVD-DE's test (FDR 5%) in ",
                   "at least one half"),
    bool_only_nominated = TRUE
  )
  ggplot2::ggsave(file.path(figure_dir, paste0("replication_", file_stub,
                                               "_esvd.png")),
                  plot = plot1, width = 8, height = 4.8, dpi = 200)
}

summary_plot_list <- list(
  all_score = list("Replicate score over all genes (nats per gene, lower is better)", NULL),
  score_cross = list("Score against the other half's DESeq2 (nats per gene, lower is better)", NULL),
  all_z_rms = list("RMS of (lfc1 - lfc2) / sqrt(se1^2 + se2^2); 1 is calibrated", 1),
  all_pearson = list("Correlation between the halves, all genes", NULL),
  esvd_pearson = list("Correlation between the halves, genes nominated by eSVD-DE's test", NULL),
  own_pearson = list("Correlation between the halves, genes nominated by the rule's own test", NULL),
  esvd_sign_rep = list("Same sign in the other half, genes nominated by eSVD-DE's test", NULL),
  own_sign_rep = list("Same sign in the other half, genes nominated by the rule's own test", NULL),
  own_slope_rep = list("Replication slope, genes nominated by the rule's own test", 1)
)
for(column in names(summary_plot_list)){
  plot1 <- .plot_summary(summary_df = summary_df, column = column,
                         xlab = summary_plot_list[[column]][[1]],
                         method_vec = method_vec,
                         reference = summary_plot_list[[column]][[2]])
  ggplot2::ggsave(file.path(figure_dir, paste0("summary_", column, ".png")),
                  plot = plot1, width = 7, height = 5.2, dpi = 200)
}

print("Done! :)")
