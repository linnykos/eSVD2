# Beyond-score study, step 3: score what step 2 stored
# Drafted by Claude for Kevin Z. Lin, 2026-11-12

rm(list = ls())

# Reads the per-gene tables of step 2 for the three banks and writes the
# tables the writeup reads, all under csv/ beside this script:
#   cohort_scores.csv.gz one row per bank x cohort x method x truth kind
#   strata.csv           mean score and the head-to-head with DESeq2, per
#                        bank x method x truth kind x factor level
#   gene_bins.csv        the same per gene bin (sparsity, nuisance status,
#                        effect class)
#   library.csv          each library variant against the fitted library,
#                        paired within cohorts, per factor level
#   significance.csv     significant genes against the +/- 2 SE interval
#   replicates.csv.gz    agreement of the logFC between the two replicate
#                        cohorts of an example (third bank)
#   replicate_genes.csv.gz the genes of a few examples, for the scatter plots
#   bank3_cohorts.csv    realized values of the third bank
#   bank_genes.csv       genes per sparsity bin and the share at the
#                        library floor, per bank
#   failures.csv         non-finite returns, per bank x method
#
# Run from the repository root:
#   Rscript additional_context/2026-11-12_beyond-score/03_summarize_claude.R
# Cores: `ESVD2_CORES`, default 5.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context",
                       "2026-11-12_beyond-score")
stopifnot(dir.exists(study_dir))
csv_dir <- file.path(study_dir, "csv")
dir.create(csv_dir, showWarnings = FALSE)
num_cores <- as.integer(Sys.getenv("ESVD2_CORES", unset = "5"))
root_dir <- Sys.getenv("ESVD2_BANK_DIR", unset = "")
if(root_dir == "") root_dir <- tools::R_user_dir("eSVD2", which = "cache")
eval_root <- file.path(root_dir, "lfcse-beyond-score")
bank_vec <- c("bank1", "bank2", "bank3")
bank_vec <- bank_vec[dir.exists(file.path(eval_root, bank_vec))]
truth_kind_vec <- c("marginal", "direct")
# A gene is called at this false discovery rate, by any test.
fdr_level <- 0.05
# Upper edges of the gene bins by the share of zero counts.
sparsity_break_vec <- c(0, 0.5, 0.8, 0.95, 1)
sparsity_label_vec <- c("under 50%", "50 to 80%", "80 to 95%", "over 95%")
# The examples whose genes are exported for the replicate scatter plots:
# 20 donors, 50 cells each, cell-level confounders, donor SD 0.5 where the
# design has it.
scatter_method_vec <- c("stored", "donor_wls", "donor_wls@esvd_nofloor",
                        "donor_wls@total_cov", "state_adjusted", "deseq2")

# Scoring one cohort -----------------------------------------------------------

.nll_vec <- function(lfc_vec, se_vec, truth_vec){
  0.5 * log(2 * pi * se_vec^2) + (lfc_vec - truth_vec)^2 / (2 * se_vec^2)
}

.gene_class <- function(truth_vec){
  class_vec <- rep("strong", length(truth_vec))
  class_vec[abs(truth_vec) < 0.5] <- "weak"
  class_vec[truth_vec == 0] <- "null"
  class_vec
}

.sparsity_bin <- function(zero_fraction_vec){
  as.character(cut(zero_fraction_vec, breaks = sparsity_break_vec,
                   labels = sparsity_label_vec, include.lowest = TRUE,
                   right = FALSE))
}

.mean_or_na <- function(vec){
  if(length(vec) > 0) mean(vec) else NA_real_
}

# The cohort-level row of one method under one truth.
.score_method <- function(lfc_vec, se_vec, truth_vec, direct_vec, gene_df,
                          deseq2_nll_vec, own_padj_vec){
  bool_ok <- is.finite(lfc_vec) & is.finite(se_vec) & se_vec > 0
  nll_vec <- .nll_vec(lfc_vec, se_vec, truth_vec)
  class_vec <- .gene_class(truth_vec)
  bool_null <- class_vec == "null" & bool_ok
  bool_weak <- class_vec == "weak" & bool_ok
  z_vec <- lfc_vec / se_vec
  bool_common <- bool_ok & is.finite(deseq2_nll_vec)
  # The +/- 2 SE interval leaves out zero.
  bool_excludes <- bool_ok & abs(lfc_vec) > 2 * se_vec
  bool_own <- !is.na(own_padj_vec) & own_padj_vec <= fdr_level
  bool_esvd <- !is.na(gene_df$esvd_fdr) & gene_df$esvd_fdr <= fdr_level
  bool_truly_null <- truth_vec == 0

  # The second gate's sums (autoresearch_lfcse2/harness/evaluate_claude.R):
  # the log SE of each strong-effect gene less the log SE a smoothing
  # spline in the log mean count gives the other genes at its level.
  bool_strong <- abs(direct_vec) >= 0.5 & bool_ok
  bool_other <- abs(direct_vec) < 0.5 & bool_ok
  strong_sum <- 0
  strong_num <- 0
  if(sum(bool_strong) >= 5 && sum(bool_other) >= 50){
    level_vec <- log(gene_df$mean_count)
    spline_fit <- stats::smooth.spline(x = level_vec[bool_other],
                                       y = log(se_vec[bool_other]))
    expected_vec <- stats::predict(spline_fit, x = level_vec[bool_strong])$y
    strong_sum <- sum(log(se_vec[bool_strong]) - expected_vec)
    strong_num <- sum(bool_strong)
  }

  data.frame(
    num_finite = sum(bool_ok),
    score = mean(nll_vec[bool_ok]),
    mse = mean((lfc_vec[bool_ok] - truth_vec[bool_ok])^2),
    coverage_95 = mean(abs(lfc_vec[bool_ok] - truth_vec[bool_ok]) <=
                         stats::qnorm(0.975) * se_vec[bool_ok]),
    z_sd_null = if(sum(bool_null) > 1) stats::sd(z_vec[bool_null]) else
      NA_real_,
    bias_null = if(any(bool_null)) stats::median(lfc_vec[bool_null]) else
      NA_real_,
    weak_cross = sum(lfc_vec[bool_weak] * truth_vec[bool_weak]),
    weak_square = sum(truth_vec[bool_weak]^2),
    strong_sum = strong_sum,
    strong_num = strong_num,
    median_se = stats::median(se_vec[bool_ok]),
    num_common = sum(bool_common),
    score_common = .mean_or_na(nll_vec[bool_common]),
    deseq2_common = .mean_or_na(deseq2_nll_vec[bool_common]),
    num_null = sum(bool_truly_null),
    num_nonnull = sum(!bool_truly_null),
    num_excludes = sum(bool_excludes),
    num_own = sum(bool_own),
    num_own_null = sum(bool_own & bool_truly_null),
    num_own_excludes = sum(bool_own & bool_excludes),
    num_esvd = sum(bool_esvd),
    num_esvd_null = sum(bool_esvd & bool_truly_null),
    num_esvd_excludes = sum(bool_esvd & bool_excludes),
    num_esvd_covers_truth = sum(bool_esvd & bool_ok &
                                  abs(lfc_vec - truth_vec) <= 2 * se_vec),
    stringsAsFactors = FALSE
  )
}

# The gene-bin rows of one method under one truth.
.score_bins <- function(lfc_vec, se_vec, truth_vec, gene_df, deseq2_nll_vec){
  bool_ok <- is.finite(lfc_vec) & is.finite(se_vec) & se_vec > 0
  nll_vec <- .nll_vec(lfc_vec, se_vec, truth_vec)
  bin_list <- list(
    sparsity = .sparsity_bin(gene_df$zero_fraction),
    status = ifelse(gene_df$status == "estimated", "estimated", "capped"),
    class = .gene_class(truth_vec)
  )
  do.call(rbind, lapply(names(bin_list), function(bin_type){
    bin_vec <- bin_list[[bin_type]]
    do.call(rbind, lapply(unique(bin_vec), function(bin){
      idx <- which(bin_vec == bin & bool_ok)
      idx_common <- idx[is.finite(deseq2_nll_vec[idx])]
      data.frame(bin_type = bin_type,
                 bin = bin,
                 num_genes = length(idx),
                 score = .mean_or_na(nll_vec[idx]),
                 mse = .mean_or_na((lfc_vec[idx] - truth_vec[idx])^2),
                 coverage_95 = .mean_or_na(abs(lfc_vec[idx] - truth_vec[idx]) <=
                                             stats::qnorm(0.975) * se_vec[idx]),
                 num_common = length(idx_common),
                 score_common = .mean_or_na(nll_vec[idx_common]),
                 deseq2_common = .mean_or_na(deseq2_nll_vec[idx_common]),
                 # Sums for the slope of the logFC on the truth over the
                 # bin's genes with a nonzero truth.
                 slope_cross = sum(lfc_vec[idx] * truth_vec[idx]),
                 slope_square = sum(truth_vec[idx]^2),
                 median_se = if(length(idx) > 0) stats::median(se_vec[idx]) else
                   NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }))
}

.read_cohort <- function(file){
  out <- readRDS(file)
  gene_df <- out$gene_df
  method_vec <- colnames(out$lfc_mat)
  meta_df <- as.data.frame(lapply(out$meta, as.character),
                           stringsAsFactors = FALSE)
  meta_df <- meta_df[, setdiff(colnames(meta_df),
                               c("id", "split", "cell_idx", "seed")),
                     drop = FALSE]
  key_df <- cbind(data.frame(bank = out$bank, cohort = out$cohort,
                             id = out$id, replicate = out$replicate,
                             stringsAsFactors = FALSE), meta_df)
  cohort_list <- list()
  bin_list <- list()
  for(truth_kind in truth_kind_vec){
    truth_vec <- gene_df[[paste0("truth_", truth_kind)]]
    deseq2_nll_vec <- .nll_vec(out$lfc_mat[, "deseq2"], out$se_mat[, "deseq2"],
                               truth_vec)
    deseq2_nll_vec[!(out$se_mat[, "deseq2"] > 0)] <- NA
    for(method in method_vec){
      lfc_vec <- out$lfc_mat[, method]
      se_vec <- out$se_mat[, method]
      # A rule's own test is the Wald test of its pair with BH; DESeq2's
      # own is its adjusted p-value.
      own_padj_vec <- if(method == "deseq2"){
        gene_df$deseq2_padj
      } else if(method == "deseq2_nocell"){
        gene_df$deseq2_nocell_padj
      } else {
        stats::p.adjust(2 * stats::pnorm(-abs(lfc_vec / se_vec)),
                        method = "BH")
      }
      label_df <- data.frame(method = method, truth_kind = truth_kind,
                             stringsAsFactors = FALSE)
      cohort_list[[length(cohort_list) + 1]] <- cbind(
        key_df, label_df,
        .score_method(lfc_vec = lfc_vec, se_vec = se_vec,
                      truth_vec = truth_vec,
                      direct_vec = gene_df$truth_direct, gene_df = gene_df,
                      deseq2_nll_vec = deseq2_nll_vec,
                      own_padj_vec = own_padj_vec)
      )
      bin_list[[length(bin_list) + 1]] <- cbind(
        key_df[, c("bank", "cohort")], label_df,
        .score_bins(lfc_vec = lfc_vec, se_vec = se_vec, truth_vec = truth_vec,
                    gene_df = gene_df, deseq2_nll_vec = deseq2_nll_vec),
        row.names = NULL
      )
    }
  }
  gene_summary_df <- data.frame(
    bank = out$bank,
    cohort = out$cohort,
    bin = sparsity_label_vec,
    num_genes = as.numeric(table(factor(.sparsity_bin(gene_df$zero_fraction),
                                        levels = sparsity_label_vec))),
    floor_sum = as.numeric(tapply(gene_df$floor_fraction,
                                  factor(.sparsity_bin(gene_df$zero_fraction),
                                         levels = sparsity_label_vec), sum)),
    capped_sum = as.numeric(tapply(gene_df$status != "estimated",
                                   factor(.sparsity_bin(gene_df$zero_fraction),
                                          levels = sparsity_label_vec), sum)),
    stringsAsFactors = FALSE
  )
  gene_summary_df[is.na(gene_summary_df)] <- 0
  failure_df <- data.frame(bank = out$bank, cohort = out$cohort,
                           method = method_vec,
                           num_nonfinite = colSums(!is.finite(out$lfc_mat) |
                                                     !is.finite(out$se_mat)),
                           error = as.character(out$error_vec[method_vec]),
                           stringsAsFactors = FALSE)
  list(cohort_df = do.call(rbind, cohort_list),
       bin_df = do.call(rbind, bin_list),
       gene_summary_df = gene_summary_df,
       failure_df = failure_df)
}

# Replicates -------------------------------------------------------------------

# Agreement of one method's logFC between the two replicate cohorts of an
# example, over a set of genes.
.replicate_row <- function(lfc_a_vec, lfc_b_vec, idx, min_genes = 8){
  bool_ok <- is.finite(lfc_a_vec[idx]) & is.finite(lfc_b_vec[idx])
  a_vec <- lfc_a_vec[idx][bool_ok]
  b_vec <- lfc_b_vec[idx][bool_ok]
  if(length(a_vec) < min_genes){
    return(data.frame(num_genes = length(a_vec), pearson = NA_real_,
                      spearman = NA_real_, mean_lfc = NA_real_,
                      sd_lfc = NA_real_, sign_agree = NA_real_))
  }
  data.frame(num_genes = length(a_vec),
             pearson = stats::cor(a_vec, b_vec),
             spearman = stats::cor(a_vec, b_vec, method = "spearman"),
             mean_lfc = mean(c(a_vec, b_vec)),
             sd_lfc = stats::sd(c(a_vec, b_vec)),
             sign_agree = mean(sign(a_vec) == sign(b_vec)))
}

.read_replicates <- function(id){
  out_a <- readRDS(file.path(eval_root, "bank3", paste0(id, "_a.rds")))
  out_b <- readRDS(file.path(eval_root, "bank3", paste0(id, "_b.rds")))
  direct_vec <- out_a$gene_df$truth_direct
  marginal_vec <- out_a$gene_df$truth_marginal
  zero_vec <- (out_a$gene_df$zero_fraction + out_b$gene_df$zero_fraction) / 2
  sparsity_vec <- .sparsity_bin(zero_vec)
  meta_df <- as.data.frame(lapply(out_a$meta, as.character),
                           stringsAsFactors = FALSE)
  meta_df <- meta_df[, setdiff(colnames(meta_df), c("cell_idx", "seed")),
                     drop = FALSE]
  # A true DEG has a case-control effect of its own; at `dense` every gene
  # does and the strong ones stand in for the DEGs.
  set_vec <- if(out_a$meta$deg == "dense"){
    ifelse(abs(direct_vec) >= 0.5, "deg", "other")
  } else {
    ifelse(direct_vec != 0, "deg", "null")
  }
  method_vec <- colnames(out_a$lfc_mat)
  row_list <- list()
  for(method in method_vec){
    lfc_a_vec <- out_a$lfc_mat[, method]
    lfc_b_vec <- out_b$lfc_mat[, method]
    for(gene_set in unique(set_vec)){
      for(sparsity in c("all", sparsity_label_vec)){
        idx <- which(set_vec == gene_set &
                       (sparsity == "all" | sparsity_vec == sparsity))
        row_list[[length(row_list) + 1]] <- cbind(
          meta_df,
          data.frame(method = method, gene_set = gene_set,
                     sparsity = sparsity, stringsAsFactors = FALSE),
          .replicate_row(lfc_a_vec, lfc_b_vec, idx)
        )
      }
    }
  }
  gene_df <- NULL
  if(out_a$meta$donors == 20 && out_a$meta$cells == "equal" &&
     out_a$meta$deg != "dense"){
    gene_df <- do.call(rbind, lapply(scatter_method_vec, function(method){
      data.frame(id = id, deg = out_a$meta$deg, state = out_a$meta$state,
                 confound = out_a$meta$confound,
                 donor_sd = out_a$meta$donor_sd,
                 method = method, gene = out_a$gene_df$gene,
                 gene_set = set_vec, sparsity = sparsity_vec,
                 truth_direct = direct_vec, truth_marginal = marginal_vec,
                 lfc_a = out_a$lfc_mat[, method],
                 lfc_b = out_b$lfc_mat[, method],
                 stringsAsFactors = FALSE)
    }))
  }
  list(replicate_df = do.call(rbind, row_list), gene_df = gene_df)
}

# Reading every cohort ---------------------------------------------------------

.bind_fill <- function(df_list){
  column_vec <- unique(unlist(lapply(df_list, colnames)))
  do.call(rbind, lapply(df_list, function(df){
    df[setdiff(column_vec, colnames(df))] <- NA
    df[, column_vec, drop = FALSE]
  }))
}

file_vec <- unlist(lapply(bank_vec, function(bank){
  list.files(file.path(eval_root, bank), pattern = "\\.rds$",
             full.names = TRUE)
}))
print(paste0("Scoring ", length(file_vec), " cohorts of ",
             paste(bank_vec, collapse = ", ")))
res_list <- parallel::mclapply(file_vec, .read_cohort, mc.cores = num_cores)
stopifnot(all(sapply(res_list, is.list)))
cohort_df <- .bind_fill(lapply(res_list, function(res){res$cohort_df}))
bin_df <- do.call(rbind, lapply(res_list, function(res){res$bin_df}))
gene_summary_df <- do.call(rbind, lapply(res_list, function(res){
  res$gene_summary_df
}))
failure_df <- do.call(rbind, lapply(res_list, function(res){res$failure_df}))
rm(res_list)
cohort_df$vs_deseq2 <- cohort_df$score_common - cohort_df$deseq2_common
factor_list <- list(bank1 = c("donors", "cells", "rates", "effects",
                              "covariates", "donor_sd"),
                    bank2 = c("donors", "cells", "expression", "effects",
                              "covariates", "donor_sd", "state", "arm_sd",
                              "depth"),
                    bank3 = c("donors", "confound", "deg", "state", "cells",
                              "donor_sd", "arm_sd", "depth"))

# Compressed, and with the columns the writeup's paired comparisons need:
# the full table is about 25 MB as text.
keep_vec <- c("bank", "cohort", "id", "replicate",
              unique(unlist(factor_list)), "method", "truth_kind",
              "num_finite", "score", "mse", "coverage_95", "vs_deseq2")
small_df <- cohort_df[, intersect(keep_vec, colnames(cohort_df))]
for(name in c("score", "mse", "coverage_95", "vs_deseq2")){
  small_df[[name]] <- signif(small_df[[name]], 6)
}
gz_con <- gzfile(file.path(csv_dir, "cohort_scores.csv.gz"), open = "w")
utils::write.csv(small_df, file = gz_con, row.names = FALSE)
close(gz_con)

# Strata -----------------------------------------------------------------------

.se_of_mean <- function(vec){
  vec <- vec[is.finite(vec)]
  if(length(vec) < 2) return(NA_real_)
  stats::sd(vec) / sqrt(length(vec))
}

# One summary row of a set of cohort rows of one method.
.summary_row <- function(sub_df){
  data.frame(
    num_cohorts = nrow(sub_df),
    score = mean(sub_df$score),
    score_se = .se_of_mean(sub_df$score),
    rmse = sqrt(mean(sub_df$mse)),
    coverage_95 = mean(sub_df$coverage_95),
    z_sd_null = mean(sub_df$z_sd_null, na.rm = TRUE),
    bias_null = mean(sub_df$bias_null, na.rm = TRUE),
    slope_weak = sum(sub_df$weak_cross) / sum(sub_df$weak_square),
    se_strong = exp(sum(sub_df$strong_sum) / sum(sub_df$strong_num)),
    median_se = mean(sub_df$median_se),
    vs_deseq2 = mean(sub_df$vs_deseq2, na.rm = TRUE),
    vs_deseq2_se = .se_of_mean(sub_df$vs_deseq2),
    win_deseq2 = mean(sub_df$vs_deseq2 < 0, na.rm = TRUE),
    num_nonfinite_cohorts = sum(sub_df$num_finite < max(sub_df$num_finite)),
    stringsAsFactors = FALSE
  )
}

group_df <- unique(cohort_df[, c("bank", "method", "truth_kind")])
strata_list <- list()
for(i in seq_len(nrow(group_df))){
  bool_group <- cohort_df$bank == group_df$bank[i] &
    cohort_df$method == group_df$method[i] &
    cohort_df$truth_kind == group_df$truth_kind[i]
  sub_df <- cohort_df[bool_group, ]
  strata_list[[length(strata_list) + 1]] <- cbind(
    group_df[i, ], factor = "all", level = "all", .summary_row(sub_df),
    row.names = NULL
  )
  for(factor_name in factor_list[[group_df$bank[i]]]){
    for(level in sort(unique(sub_df[[factor_name]]))){
      strata_list[[length(strata_list) + 1]] <- cbind(
        group_df[i, ], factor = factor_name, level = level,
        .summary_row(sub_df[sub_df[[factor_name]] == level, ]),
        row.names = NULL
      )
    }
  }
}
strata_df <- do.call(rbind, strata_list)
utils::write.csv(strata_df, file = file.path(csv_dir, "strata.csv"),
                 row.names = FALSE)

# Gene bins --------------------------------------------------------------------

# Per bank x method x truth kind x bin: scores averaged over genes (each
# cohort weighted by its genes in the bin); the head-to-head with DESeq2 is
# the mean over cohorts of the per-cohort difference, with its paired SE.
bin_df$vs_deseq2 <- bin_df$score_common - bin_df$deseq2_common
bin_key <- paste(bin_df$bank, bin_df$method, bin_df$truth_kind,
                 bin_df$bin_type, bin_df$bin, sep = "\r")
gene_bin_df <- do.call(rbind, lapply(split(bin_df, bin_key), function(sub_df){
  sub_df <- sub_df[sub_df$num_genes > 0, ]
  weight_vec <- sub_df$num_genes
  data.frame(bank = sub_df$bank[1], method = sub_df$method[1],
             truth_kind = sub_df$truth_kind[1], bin_type = sub_df$bin_type[1],
             bin = sub_df$bin[1],
             num_genes = sum(weight_vec),
             num_cohorts = nrow(sub_df),
             score = sum(weight_vec * sub_df$score) / sum(weight_vec),
             rmse = sqrt(sum(weight_vec * sub_df$mse) / sum(weight_vec)),
             coverage_95 = sum(weight_vec * sub_df$coverage_95) /
               sum(weight_vec),
             slope = sum(sub_df$slope_cross) / sum(sub_df$slope_square),
             median_se = stats::median(sub_df$median_se),
             vs_deseq2 = mean(sub_df$vs_deseq2, na.rm = TRUE),
             vs_deseq2_se = .se_of_mean(sub_df$vs_deseq2),
             win_deseq2 = mean(sub_df$vs_deseq2 < 0, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
rownames(gene_bin_df) <- NULL
utils::write.csv(gene_bin_df, file = file.path(csv_dir, "gene_bins.csv"),
                 row.names = FALSE)

# The sparsity bins crossed with the confounding level of the third bank,
# for the rules the recommendation turns on.
if("bank3" %in% bank_vec){
  confound_df <- unique(cohort_df[cohort_df$bank == "bank3",
                                  c("cohort", "confound")])
  cross_df <- bin_df[bin_df$bank == "bank3" & bin_df$bin_type == "sparsity" &
                       bin_df$num_genes > 0, ]
  cross_df$confound <- confound_df$confound[match(cross_df$cohort,
                                                  confound_df$cohort)]
  cross_key <- paste(cross_df$method, cross_df$truth_kind, cross_df$bin,
                     cross_df$confound, sep = "\r")
  cross_out_df <- do.call(rbind, lapply(split(cross_df, cross_key),
                                        function(sub_df){
    weight_vec <- sub_df$num_genes
    data.frame(method = sub_df$method[1], truth_kind = sub_df$truth_kind[1],
               bin = sub_df$bin[1], confound = sub_df$confound[1],
               num_genes = sum(weight_vec),
               score = sum(weight_vec * sub_df$score) / sum(weight_vec),
               rmse = sqrt(sum(weight_vec * sub_df$mse) / sum(weight_vec)),
               coverage_95 = sum(weight_vec * sub_df$coverage_95) /
                 sum(weight_vec),
               slope = sum(sub_df$slope_cross) / sum(sub_df$slope_square),
               vs_deseq2 = mean(sub_df$vs_deseq2, na.rm = TRUE),
               vs_deseq2_se = .se_of_mean(sub_df$vs_deseq2),
               stringsAsFactors = FALSE)
  }))
  rownames(cross_out_df) <- NULL
  utils::write.csv(cross_out_df,
                   file = file.path(csv_dir, "sparsity_by_confound.csv"),
                   row.names = FALSE)
}

# Library variants -------------------------------------------------------------

# Each variant `rule@library` against `rule` (the fitted library), paired
# within cohorts: a positive difference means the fitted library is better.
variant_vec <- grep("@", unique(cohort_df$method), value = TRUE)
library_list <- list()
for(bank in bank_vec){
  for(truth_kind in truth_kind_vec){
    bool_base <- cohort_df$bank == bank & cohort_df$truth_kind == truth_kind
    for(variant in variant_vec){
      rule <- sub("@.*", "", variant)
      variant_df <- cohort_df[bool_base & cohort_df$method == variant, ]
      if(nrow(variant_df) == 0) next
      base_df <- cohort_df[bool_base & cohort_df$method == rule, ]
      base_df <- base_df[match(variant_df$cohort, base_df$cohort), ]
      diff_vec <- variant_df$score - base_df$score
      mse_ratio_vec <- variant_df$mse / base_df$mse
      for(factor_name in c("all", factor_list[[bank]])){
        level_vec <- if(factor_name == "all") rep("all", nrow(variant_df)) else
          variant_df[[factor_name]]
        for(level in sort(unique(level_vec))){
          idx <- which(level_vec == level)
          library_list[[length(library_list) + 1]] <- data.frame(
            bank = bank, truth_kind = truth_kind, rule = rule,
            library = sub(".*@", "", variant), factor = factor_name,
            level = level, num_cohorts = length(idx),
            score_variant = mean(variant_df$score[idx]),
            score_fitted = mean(base_df$score[idx]),
            diff = mean(diff_vec[idx]),
            diff_se = .se_of_mean(diff_vec[idx]),
            fitted_better = mean(diff_vec[idx] > 0),
            rmse_variant = sqrt(mean(variant_df$mse[idx])),
            rmse_fitted = sqrt(mean(base_df$mse[idx])),
            coverage_variant = mean(variant_df$coverage_95[idx]),
            coverage_fitted = mean(base_df$coverage_95[idx]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}
library_df <- do.call(rbind, library_list)
utils::write.csv(library_df, file = file.path(csv_dir, "library.csv"),
                 row.names = FALSE)

# Significance against the interval --------------------------------------------

# Pooled over the cohorts of a bank (and, for the third bank, of each
# level of `deg`): of the genes a test calls, the share whose +/- 2 SE
# interval leaves out zero.
.significance_row <- function(sub_df){
  data.frame(
    num_cohorts = nrow(sub_df),
    num_own = sum(sub_df$num_own),
    own_excludes = sum(sub_df$num_own_excludes) / sum(sub_df$num_own),
    own_fdp = sum(sub_df$num_own_null) / sum(sub_df$num_own),
    own_power = sum(sub_df$num_own - sub_df$num_own_null) /
      sum(sub_df$num_nonnull),
    num_esvd = sum(sub_df$num_esvd),
    esvd_excludes = sum(sub_df$num_esvd_excludes) / sum(sub_df$num_esvd),
    esvd_fdp = sum(sub_df$num_esvd_null) / sum(sub_df$num_esvd),
    esvd_power = sum(sub_df$num_esvd - sub_df$num_esvd_null) /
      sum(sub_df$num_nonnull),
    esvd_covers_truth = sum(sub_df$num_esvd_covers_truth) /
      sum(sub_df$num_esvd),
    excludes_per_cohort = mean(sub_df$num_excludes),
    stringsAsFactors = FALSE
  )
}
significance_list <- list()
for(i in seq_len(nrow(group_df))){
  bool_group <- cohort_df$bank == group_df$bank[i] &
    cohort_df$method == group_df$method[i] &
    cohort_df$truth_kind == group_df$truth_kind[i]
  sub_df <- cohort_df[bool_group, ]
  effect_name <- if(group_df$bank[i] == "bank3") "deg" else "effects"
  significance_list[[length(significance_list) + 1]] <- cbind(
    group_df[i, ], effects = "all", .significance_row(sub_df),
    row.names = NULL
  )
  for(level in sort(unique(sub_df[[effect_name]]))){
    significance_list[[length(significance_list) + 1]] <- cbind(
      group_df[i, ], effects = level,
      .significance_row(sub_df[sub_df[[effect_name]] == level, ]),
      row.names = NULL
    )
  }
}
significance_df <- do.call(rbind, significance_list)
utils::write.csv(significance_df,
                 file = file.path(csv_dir, "significance.csv"),
                 row.names = FALSE)

# Replicates -------------------------------------------------------------------

if("bank3" %in% bank_vec){
  id_vec <- unique(cohort_df$id[cohort_df$bank == "bank3"])
  bool_pair <- file.exists(file.path(eval_root, "bank3",
                                     paste0(id_vec, "_a.rds"))) &
    file.exists(file.path(eval_root, "bank3", paste0(id_vec, "_b.rds")))
  replicate_list <- parallel::mclapply(id_vec[bool_pair], .read_replicates,
                                       mc.cores = num_cores)
  stopifnot(all(sapply(replicate_list, is.list)))
  replicate_df <- do.call(rbind, lapply(replicate_list, function(res){
    res$replicate_df
  }))
  # Compressed and rounded: the table has 150,000 rows.
  replicate_df <- replicate_df[, c("id", "donors", "confound", "deg", "state",
                                   "method", "gene_set", "sparsity",
                                   "num_genes", "pearson", "spearman",
                                   "mean_lfc", "sd_lfc", "sign_agree")]
  for(name in c("pearson", "spearman", "mean_lfc", "sd_lfc", "sign_agree")){
    replicate_df[[name]] <- round(replicate_df[[name]], 4)
  }
  gz_con <- gzfile(file.path(csv_dir, "replicates.csv.gz"), open = "w")
  utils::write.csv(replicate_df, file = gz_con, row.names = FALSE)
  close(gz_con)
  replicate_gene_df <- do.call(rbind, lapply(replicate_list, function(res){
    res$gene_df
  }))
  # One example per level of `deg` x `state` at cell-level confounding (or
  # the first available), to keep the table small.
  example_df <- unique(replicate_gene_df[, c("id", "deg", "state", "confound",
                                             "donor_sd")])
  example_df <- example_df[order(example_df$confound != "cell",
                                 -as.numeric(example_df$donor_sd),
                                 example_df$id), ]
  example_df <- example_df[!duplicated(example_df[, c("deg", "state")]), ]
  replicate_gene_df <- replicate_gene_df[replicate_gene_df$id %in%
                                           example_df$id, ]
  for(name in c("truth_direct", "truth_marginal", "lfc_a", "lfc_b")){
    replicate_gene_df[[name]] <- round(replicate_gene_df[[name]], 4)
  }
  gz_con <- gzfile(file.path(csv_dir, "replicate_genes.csv.gz"), open = "w")
  utils::write.csv(replicate_gene_df, file = gz_con, row.names = FALSE)
  close(gz_con)

  manifest_df <- utils::read.csv(file.path(root_dir, "lfcse-bank-3",
                                           "manifest.csv"),
                                 stringsAsFactors = FALSE)
  utils::write.csv(manifest_df[, setdiff(colnames(manifest_df),
                                         c("file_truth", "error"))],
                   file = file.path(csv_dir, "bank3_cohorts.csv"),
                   row.names = FALSE)
}

# Bank descriptives and failures -----------------------------------------------

bank_gene_df <- stats::aggregate(cbind(num_genes, floor_sum, capped_sum) ~
                                   bank + bin, data = gene_summary_df,
                                 FUN = sum)
bank_gene_df$floor_share <- bank_gene_df$floor_sum / bank_gene_df$num_genes
bank_gene_df$capped_share <- bank_gene_df$capped_sum / bank_gene_df$num_genes
utils::write.csv(bank_gene_df, file = file.path(csv_dir, "bank_genes.csv"),
                 row.names = FALSE)

failure_out_df <- stats::aggregate(
  cbind(num_genes = num_nonfinite, num_cohorts = num_nonfinite > 0,
        num_errors = error != "") ~ bank + method,
  data = failure_df, FUN = sum
)
utils::write.csv(failure_out_df, file = file.path(csv_dir, "failures.csv"),
                 row.names = FALSE)

print(strata_df[strata_df$factor == "all" & strata_df$truth_kind == "marginal",
                c("bank", "method", "score", "rmse", "coverage_95",
                  "vs_deseq2", "win_deseq2")])
print("Done! :)")
