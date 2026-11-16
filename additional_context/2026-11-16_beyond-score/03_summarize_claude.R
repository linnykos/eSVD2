# Second beyond-score study, step 3: score what step 2 stored
# Drafted by Claude for Kevin Z. Lin, 2026-11-16

rm(list = ls())

# Reads the per-gene tables of step 2 for the fourth bank and writes the
# tables the writeup reads, all under csv/ beside this script:
#   bank4_cohorts.csv        realized values of the bank, one row per example
#   bank4_cohort_scores.csv.gz one row per fitted cohort x method
#   bank4_strata.csv         one row per method x factor level: the score
#                            against the truth and its diagnostics
#   bank4_gene_bins.csv      the same per gene bin (sparsity, nuisance status,
#                            fitted depth coefficient within 1 of 1 or not)
#   bank4_calls.csv          one row per method x test x factor level x
#                            nuisance status: calls, false discoveries, sign
#                            errors, power
#   bank4_replicates.csv     one row per example x rank x method: what a
#                            cohort without a truth can measure between two
#                            replicates (the statistics of the Gabitto and
#                            Green studies) beside the score against the
#                            truth
#   bank4_replicate_strata.csv the same averaged per method x factor level
#   bank4_failures.csv       non-finite returns and errors per method
#
# Run from the repository root:
#   Rscript additional_context/2026-11-16_beyond-score/03_summarize_claude.R
# Cores: `ESVD2_CORES`, default 5.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context",
                       "2026-11-16_beyond-score")
stopifnot(dir.exists(study_dir))
csv_dir <- file.path(study_dir, "csv")
dir.create(csv_dir, showWarnings = FALSE)
num_cores <- as.integer(Sys.getenv("ESVD2_CORES", unset = "5"))
root_dir <- Sys.getenv("ESVD2_BANK_DIR", unset = "")
if(root_dir == "") root_dir <- tools::R_user_dir("eSVD2", which = "cache")
eval_dir <- file.path(root_dir, "lfcse-beyond-score", "bank4")
bank_dir <- file.path(root_dir, "lfcse-bank-4")
factor_vec <- c("rank", "donors", "balance", "deg", "state", "donor_sd", "tail",
                "capped")
# A gene is called at this false discovery rate, by any test.
fdr_level <- 0.05
sparsity_break_vec <- c(0, 0.5, 0.8, 0.95, 1)
sparsity_label_vec <- c("under 50%", "50 to 80%", "80 to 95%", "over 95%")
# A gene counts as a true DEG of the bank when its true log2 fold change is
# at least this large; the weak effects below it are neither DEGs nor nulls.
deg_threshold <- 0.25

# Scoring one cohort -----------------------------------------------------------

.nll_vec <- function(lfc_vec, se_vec, truth_vec){
  0.5 * log(2 * pi * se_vec^2) + (lfc_vec - truth_vec)^2 / (2 * se_vec^2)
}

.sparsity_bin <- function(zero_fraction_vec){
  as.character(cut(zero_fraction_vec, breaks = sparsity_break_vec,
                   labels = sparsity_label_vec, include.lowest = TRUE,
                   right = FALSE))
}

.mean_or_na <- function(vec){
  if(length(vec) > 0) mean(vec) else NA_real_
}

# The Wald test of a reported pair against the standard normal, with BH
# over the cohort's genes: each rule's own test, as in the Gabitto and
# Green studies.
.own_padj <- function(lfc_vec, se_vec){
  p_vec <- 2 * stats::pnorm(-abs(lfc_vec / se_vec))
  p_vec[!is.finite(p_vec)] <- NA_real_
  stats::p.adjust(p_vec, method = "BH")
}

# Counts of one test's calls over the genes `bool_set`. A sign error is a
# called gene with a nonzero truth whose logFC has the other sign.
.call_counts <- function(bool_call, bool_set, lfc_vec, truth_vec){
  bool_call <- bool_call & bool_set
  bool_nonnull <- truth_vec != 0
  data.frame(num_genes = sum(bool_set),
             num_nonnull = sum(bool_set & bool_nonnull),
             num_deg = sum(bool_set & abs(truth_vec) >= deg_threshold),
             num_calls = sum(bool_call),
             num_null_calls = sum(bool_call & !bool_nonnull),
             num_sign_errors = sum(bool_call & bool_nonnull &
                                     sign(lfc_vec) != sign(truth_vec),
                                   na.rm = TRUE),
             num_deg_calls = sum(bool_call & abs(truth_vec) >= deg_threshold),
             stringsAsFactors = FALSE)
}

.read_cohort <- function(file){
  out <- readRDS(file)
  gene_df <- out$gene_df
  truth_vec <- gene_df$truth_marginal
  method_vec <- colnames(out$lfc_mat)
  meta_df <- as.data.frame(lapply(out$meta, as.character),
                           stringsAsFactors = FALSE)
  key_df <- cbind(data.frame(cohort = out$cohort, replicate = out$replicate,
                             stringsAsFactors = FALSE),
                  meta_df[, c("id", factor_vec, "draw")])
  bool_capped <- gene_df$status != "estimated"
  sparsity_vec <- .sparsity_bin(gene_df$zero_fraction)
  bool_esvd <- !is.na(gene_df$esvd_fdr) & gene_df$esvd_fdr <= fdr_level
  deseq2_nll_vec <- .nll_vec(out$lfc_mat[, "deseq2"], out$se_mat[, "deseq2"],
                             truth_vec)
  deseq2_nll_vec[!(out$se_mat[, "deseq2"] > 0)] <- NA
  status_list <- list(all = rep(TRUE, nrow(gene_df)), estimated = !bool_capped,
                      capped = bool_capped)
  # A fitted depth coefficient more than 1 away from 1 makes the unfloored
  # library range over orders of magnitude between the cells of a gene.
  bool_depth_off <- abs(gene_df$depth_coef - 1) > 1
  bin_list <- c(
    lapply(stats::setNames(sparsity_label_vec, sparsity_label_vec),
           function(bin){sparsity_vec == bin}),
    list(estimated = !bool_capped, capped = bool_capped,
         depth_within = !bool_depth_off, depth_beyond = bool_depth_off)
  )
  bin_type_vec <- c(rep("sparsity", length(sparsity_label_vec)), "status",
                    "status", "depth", "depth")

  cohort_list <- list()
  bin_row_list <- list()
  call_list <- list()
  for(method in method_vec){
    lfc_vec <- out$lfc_mat[, method]
    se_vec <- out$se_mat[, method]
    bool_ok <- is.finite(lfc_vec) & is.finite(se_vec) & se_vec > 0
    nll_vec <- .nll_vec(lfc_vec, se_vec, truth_vec)
    z_vec <- lfc_vec / se_vec
    bool_null <- truth_vec == 0 & bool_ok
    bool_weak <- truth_vec != 0 & abs(truth_vec) < 0.5 & bool_ok
    bool_deg <- abs(truth_vec) >= deg_threshold & bool_ok
    bool_common <- bool_ok & is.finite(deseq2_nll_vec)

    # The second gate's sums (autoresearch_lfcse2/harness/evaluate_claude.R):
    # the log SE of each strong-effect gene less the log SE a smoothing
    # spline in the log mean count gives the other genes at its level.
    bool_strong <- abs(truth_vec) >= 0.5 & bool_ok
    bool_other <- abs(truth_vec) < 0.5 & bool_ok
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

    cohort_list[[method]] <- cbind(key_df, data.frame(
      method = method,
      num_finite = sum(bool_ok),
      score = mean(nll_vec[bool_ok]),
      mse = mean((lfc_vec[bool_ok] - truth_vec[bool_ok])^2),
      coverage_95 = mean(abs(lfc_vec[bool_ok] - truth_vec[bool_ok]) <=
                           stats::qnorm(0.975) * se_vec[bool_ok]),
      z_sd_null = if(sum(bool_null) > 1) stats::sd(z_vec[bool_null]) else
        NA_real_,
      weak_cross = sum(lfc_vec[bool_weak] * truth_vec[bool_weak]),
      weak_square = sum(truth_vec[bool_weak]^2),
      # The gate on retained signal, with a truth: the slope of the logFC on
      # the truth over the true DEGs.
      deg_cross = sum(lfc_vec[bool_deg] * truth_vec[bool_deg]),
      deg_square = sum(truth_vec[bool_deg]^2),
      strong_sum = strong_sum,
      strong_num = strong_num,
      median_se = stats::median(se_vec[bool_ok]),
      sd_lfc = stats::sd(lfc_vec[bool_ok]),
      score_common = .mean_or_na(nll_vec[bool_common]),
      deseq2_common = .mean_or_na(deseq2_nll_vec[bool_common]),
      stringsAsFactors = FALSE
    ))

    bin_row_list[[method]] <- do.call(rbind, lapply(seq_along(bin_list),
                                                    function(j){
      idx <- which(bin_list[[j]] & bool_ok)
      data.frame(cohort = out$cohort, method = method,
                 bin_type = bin_type_vec[j], bin = names(bin_list)[j],
                 num_genes = length(idx),
                 score = .mean_or_na(nll_vec[idx]),
                 mse = .mean_or_na((lfc_vec[idx] - truth_vec[idx])^2),
                 coverage_95 = .mean_or_na(abs(lfc_vec[idx] - truth_vec[idx]) <=
                                             stats::qnorm(0.975) * se_vec[idx]),
                 slope_cross = sum(lfc_vec[idx] * truth_vec[idx]),
                 slope_square = sum(truth_vec[idx]^2),
                 stringsAsFactors = FALSE)
    }))

    own_padj_vec <- if(method == "deseq2"){
      gene_df$deseq2_padj
    } else {
      .own_padj(lfc_vec, se_vec)
    }
    bool_own <- !is.na(own_padj_vec) & own_padj_vec <= fdr_level
    test_list <- list(own = bool_own, esvd = bool_esvd)
    call_list[[method]] <- do.call(rbind, lapply(names(test_list),
                                                 function(test){
      do.call(rbind, lapply(names(status_list), function(status){
        cbind(data.frame(cohort = out$cohort, method = method, test = test,
                         status = status, stringsAsFactors = FALSE),
              .call_counts(bool_call = test_list[[test]],
                           bool_set = status_list[[status]],
                           lfc_vec = lfc_vec, truth_vec = truth_vec))
      }))
    }))
  }
  failure_df <- data.frame(cohort = out$cohort, method = method_vec,
                           num_nonfinite = colSums(!is.finite(out$lfc_mat) |
                                                     !is.finite(out$se_mat)),
                           error = as.character(out$error_vec[method_vec]),
                           stringsAsFactors = FALSE)
  list(cohort_df = do.call(rbind, cohort_list),
       bin_df = do.call(rbind, bin_row_list),
       call_df = do.call(rbind, call_list),
       failure_df = failure_df)
}

# Replicates -------------------------------------------------------------------

# What the Gabitto and Green studies computed between two halves, between
# the two replicate cohorts of an example at one rank, beside the truth.
.read_replicates <- function(id, rank){
  out_a <- readRDS(file.path(eval_dir, paste0(id, "_a_k", rank, ".rds")))
  out_b <- readRDS(file.path(eval_dir, paste0(id, "_b_k", rank, ".rds")))
  stopifnot(identical(out_a$gene_df$gene, out_b$gene_df$gene))
  truth_vec <- out_a$gene_df$truth_marginal
  meta_df <- as.data.frame(lapply(out_a$meta, as.character),
                           stringsAsFactors = FALSE)
  key_df <- meta_df[, c("id", factor_vec, "draw")]
  bool_deg <- abs(truth_vec) >= deg_threshold
  bool_null <- truth_vec == 0
  esvd_a <- !is.na(out_a$gene_df$esvd_fdr) & out_a$gene_df$esvd_fdr <= fdr_level
  esvd_b <- !is.na(out_b$gene_df$esvd_fdr) & out_b$gene_df$esvd_fdr <= fdr_level
  .cor_or_na <- function(x_vec, y_vec, idx){
    idx <- idx[is.finite(x_vec[idx]) & is.finite(y_vec[idx])]
    if(length(idx) < 8) return(NA_real_)
    stats::cor(x_vec[idx], y_vec[idx])
  }
  # Discovery in one replicate, the sign in the other, pooled over the two
  # directions; and whether the discovered sign is the true one.
  .sign_rows <- function(lfc_a, lfc_b, call_a, call_b){
    disc_vec <- c(lfc_a[call_a], lfc_b[call_b])
    rep_vec <- c(lfc_b[call_a], lfc_a[call_b])
    true_vec <- c(truth_vec[call_a], truth_vec[call_b])
    bool_ok <- is.finite(disc_vec) & is.finite(rep_vec)
    disc_vec <- disc_vec[bool_ok]
    rep_vec <- rep_vec[bool_ok]
    true_vec <- true_vec[bool_ok]
    c(num = length(disc_vec),
      sign_rep = if(length(disc_vec) >= 5) mean(sign(disc_vec) ==
                                                  sign(rep_vec)) else NA_real_,
      slope_rep = if(length(disc_vec) >= 5) sum(disc_vec * rep_vec) /
        sum(disc_vec^2) else NA_real_,
      # Of the calls, the share that are exact nulls or have the wrong sign.
      wrong = if(length(disc_vec) >= 5) mean(true_vec == 0 |
                                               sign(disc_vec) !=
                                               sign(true_vec)) else NA_real_)
  }
  deseq2_a <- list(lfc = out_a$lfc_mat[, "deseq2"], se = out_a$se_mat[, "deseq2"])
  deseq2_b <- list(lfc = out_b$lfc_mat[, "deseq2"], se = out_b$se_mat[, "deseq2"])
  row_list <- lapply(colnames(out_a$lfc_mat), function(method){
    lfc_a <- out_a$lfc_mat[, method]
    lfc_b <- out_b$lfc_mat[, method]
    se_a <- out_a$se_mat[, method]
    se_b <- out_b$se_mat[, method]
    bool_ok <- is.finite(lfc_a) & is.finite(lfc_b) & is.finite(se_a) &
      is.finite(se_b) & se_a > 0 & se_b > 0
    idx <- which(bool_ok)
    diff_vec <- lfc_a[idx] - lfc_b[idx]
    var_vec <- se_a[idx]^2 + se_b[idx]^2
    bool_cross <- bool_ok & is.finite(deseq2_a$lfc) & is.finite(deseq2_b$lfc) &
      is.finite(deseq2_a$se) & is.finite(deseq2_b$se)
    cross_vec <- c(
      .nll_vec(lfc_a[bool_cross], sqrt(se_a[bool_cross]^2 +
                                         deseq2_b$se[bool_cross]^2),
               deseq2_b$lfc[bool_cross]),
      .nll_vec(lfc_b[bool_cross], sqrt(se_b[bool_cross]^2 +
                                         deseq2_a$se[bool_cross]^2),
               deseq2_a$lfc[bool_cross])
    )
    own_a <- if(method == "deseq2") out_a$gene_df$deseq2_padj else
      .own_padj(lfc_a, se_a)
    own_b <- if(method == "deseq2") out_b$gene_df$deseq2_padj else
      .own_padj(lfc_b, se_b)
    own_a <- !is.na(own_a) & own_a <= fdr_level
    own_b <- !is.na(own_b) & own_b <= fdr_level
    esvd_row <- .sign_rows(lfc_a, lfc_b, esvd_a, esvd_b)
    own_row <- .sign_rows(lfc_a, lfc_b, own_a, own_b)
    cbind(key_df, data.frame(
      method = method,
      # Against the truth, mean of the two replicates.
      score_truth = mean(c(.nll_vec(lfc_a[idx], se_a[idx], truth_vec[idx]),
                           .nll_vec(lfc_b[idx], se_b[idx], truth_vec[idx]))),
      rmse_truth = sqrt(mean(c((lfc_a[idx] - truth_vec[idx])^2,
                               (lfc_b[idx] - truth_vec[idx])^2))),
      # Without the truth.
      score_replicate = mean(.nll_vec(diff_vec, sqrt(var_vec), 0)),
      z_rms = sqrt(mean(diff_vec^2 / var_vec)),
      score_cross = mean(cross_vec),
      pearson_all = .cor_or_na(lfc_a, lfc_b, idx),
      pearson_deg = .cor_or_na(lfc_a, lfc_b, which(bool_deg)),
      pearson_null = .cor_or_na(lfc_a, lfc_b, which(bool_null)),
      pearson_esvd = .cor_or_na(lfc_a, lfc_b, which(esvd_a | esvd_b)),
      num_esvd = sum(esvd_a | esvd_b),
      esvd_sign_rep = esvd_row[["sign_rep"]],
      esvd_slope_rep = esvd_row[["slope_rep"]],
      esvd_wrong = esvd_row[["wrong"]],
      num_own = sum(own_a | own_b),
      num_own_both = sum(own_a & own_b),
      own_sign_rep = own_row[["sign_rep"]],
      own_wrong = own_row[["wrong"]],
      stringsAsFactors = FALSE
    ))
  })
  do.call(rbind, row_list)
}

# Reading every cohort ---------------------------------------------------------

# Only the cohorts of the manifest: the results folder may hold cohorts of
# examples that a later design no longer has.
manifest_df <- utils::read.csv(file.path(bank_dir, "manifest.csv"),
                               stringsAsFactors = FALSE)
file_vec <- list.files(eval_dir, pattern = "\\.rds$", full.names = TRUE)
file_vec <- file_vec[sub("_[ab]_k[0-9]+\\.rds$", "", basename(file_vec)) %in%
                       manifest_df$id[manifest_df$status == "ok"]]
print(paste0("Scoring ", length(file_vec), " fitted cohorts of bank 4"))
res_list <- parallel::mclapply(file_vec, .read_cohort, mc.cores = num_cores)
stopifnot(all(sapply(res_list, is.list)))
cohort_df <- do.call(rbind, lapply(res_list, function(res){res$cohort_df}))
bin_df <- do.call(rbind, lapply(res_list, function(res){res$bin_df}))
call_df <- do.call(rbind, lapply(res_list, function(res){res$call_df}))
failure_df <- do.call(rbind, lapply(res_list, function(res){res$failure_df}))
rm(res_list)
rownames(cohort_df) <- NULL
cohort_df$vs_deseq2 <- cohort_df$score_common - cohort_df$deseq2_common
method_vec <- unique(cohort_df$method)

# The shapes the recommendation is weighted by. A cohort is Gabitto-like
# when it has 40 or 80 donors, a third of them controls, and is fitted at
# the default rank of 30; Green-like when it has 200 donors of few cells at
# rank 30.
# The second block (outlying donors, the higher capped share) is kept
# apart as `green_block`, and it is left out of every other shape.
bool_first <- cohort_df$tail == "gaussian" & cohort_df$capped == "low"
cohort_df$shape <- "other"
cohort_df$shape[bool_first & cohort_df$rank == "30" &
                  cohort_df$balance == "third" &
                  cohort_df$donors %in% c("40", "80")] <- "gabitto_like"
cohort_df$shape[bool_first & cohort_df$rank == "30" &
                  cohort_df$balance == "third" &
                  cohort_df$donors == "200"] <- "green_like"
cohort_df$shape[!bool_first] <- "green_block"
strata_factor_vec <- c(factor_vec, "shape")

small_df <- cohort_df[, c("cohort", "id", "replicate", strata_factor_vec,
                          "method", "num_finite", "score", "mse",
                          "coverage_95", "vs_deseq2")]
for(name in c("score", "mse", "coverage_95", "vs_deseq2")){
  small_df[[name]] <- signif(small_df[[name]], 6)
}
gz_con <- gzfile(file.path(csv_dir, "bank4_cohort_scores.csv.gz"), open = "w")
utils::write.csv(small_df, file = gz_con, row.names = FALSE)
close(gz_con)

# Strata -----------------------------------------------------------------------

.se_of_mean <- function(vec){
  vec <- vec[is.finite(vec)]
  if(length(vec) < 2) return(NA_real_)
  stats::sd(vec) / sqrt(length(vec))
}

.summary_row <- function(sub_df){
  data.frame(
    num_cohorts = nrow(sub_df),
    score = mean(sub_df$score),
    score_se = .se_of_mean(sub_df$score),
    rmse = sqrt(mean(sub_df$mse)),
    coverage_95 = mean(sub_df$coverage_95),
    z_sd_null = mean(sub_df$z_sd_null, na.rm = TRUE),
    slope_weak = sum(sub_df$weak_cross) / sum(sub_df$weak_square),
    slope_deg = sum(sub_df$deg_cross) / sum(sub_df$deg_square),
    se_strong = exp(sum(sub_df$strong_sum) / sum(sub_df$strong_num)),
    median_se = mean(sub_df$median_se),
    sd_lfc = mean(sub_df$sd_lfc),
    vs_deseq2 = mean(sub_df$vs_deseq2, na.rm = TRUE),
    vs_deseq2_se = .se_of_mean(sub_df$vs_deseq2),
    win_deseq2 = mean(sub_df$vs_deseq2 < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# The strata of every table. The first block (Gaussian donors, the lower
# capped share) gives `all`, every level of its factors, the shapes, and
# rank crossed with the donors. The second block is reported under the
# factor `block`, by what it adds (`outlier`, `capped_high`), the rank and
# the donors, beside the first block's examples of the same donors, a
# third of them controls (`matched`).
.strata_sets <- function(sub_df){
  bool_first <- sub_df$tail == "gaussian" & sub_df$capped == "low"
  set_list <- list(list(factor = "all", level = "all", bool = bool_first))
  for(factor_name in setdiff(strata_factor_vec, c("tail", "capped"))){
    for(level in sort(unique(sub_df[[factor_name]][bool_first]))){
      set_list[[length(set_list) + 1]] <- list(
        factor = factor_name, level = level,
        bool = bool_first & sub_df[[factor_name]] == level
      )
    }
  }
  for(rank in sort(unique(sub_df$rank))){
    for(donors in sort(unique(sub_df$donors))){
      bool_cell <- sub_df$rank == rank & sub_df$donors == donors
      set_list[[length(set_list) + 1]] <- list(
        factor = "rank_donors", level = paste0("k", rank, "_n", donors),
        bool = bool_first & bool_cell
      )
      if(!donors %in% c("40", "200")) next
      block_list <- list(
        outlier = sub_df$tail == "outlier",
        capped_high = sub_df$capped == "high",
        matched = bool_first & sub_df$balance == "third" & sub_df$draw == "1",
        # The first block's examples that the `capped_high` ones repeat.
        matched_state = bool_first & sub_df$balance == "third" &
          sub_df$draw == "1" & sub_df$state == "donor"
      )
      for(block in names(block_list)){
        set_list[[length(set_list) + 1]] <- list(
          factor = "block",
          level = paste0(block, "_k", rank, "_n", donors),
          bool = block_list[[block]] & bool_cell
        )
      }
    }
  }
  set_list
}

strata_df <- do.call(rbind, lapply(method_vec, function(method){
  sub_df <- cohort_df[cohort_df$method == method, ]
  do.call(rbind, lapply(.strata_sets(sub_df), function(set){
    cbind(data.frame(method = method, factor = set$factor, level = set$level,
                     stringsAsFactors = FALSE),
          .summary_row(sub_df[set$bool, ]))
  }))
}))
utils::write.csv(strata_df, file = file.path(csv_dir, "bank4_strata.csv"),
                 row.names = FALSE)

# Gene bins --------------------------------------------------------------------

# Per method x rank x bin, each cohort weighted by its genes in the bin.
bin_df$rank <- cohort_df$rank[match(bin_df$cohort, cohort_df$cohort)]
bin_key <- paste(bin_df$method, bin_df$rank, bin_df$bin_type, bin_df$bin,
                 sep = "\r")
gene_bin_df <- do.call(rbind, lapply(split(bin_df, bin_key), function(sub_df){
  sub_df <- sub_df[sub_df$num_genes > 0, ]
  # A bin with no gene in any cohort (no depth coefficient far from 1 at
  # rank 5, say) has no row.
  if(nrow(sub_df) == 0) return(NULL)
  weight_vec <- sub_df$num_genes
  data.frame(method = sub_df$method[1], rank = sub_df$rank[1],
             bin_type = sub_df$bin_type[1], bin = sub_df$bin[1],
             num_genes = sum(weight_vec),
             score = sum(weight_vec * sub_df$score) / sum(weight_vec),
             rmse = sqrt(sum(weight_vec * sub_df$mse) / sum(weight_vec)),
             coverage_95 = sum(weight_vec * sub_df$coverage_95) /
               sum(weight_vec),
             slope = sum(sub_df$slope_cross) / sum(sub_df$slope_square),
             stringsAsFactors = FALSE)
}))
rownames(gene_bin_df) <- NULL
utils::write.csv(gene_bin_df, file = file.path(csv_dir, "bank4_gene_bins.csv"),
                 row.names = FALSE)

# Calls ------------------------------------------------------------------------

# Pooled over the cohorts of a stratum. The false discovery proportion
# counts exact nulls, so it exists at `deg = half_null` only; the sign
# error share is over the calls with a nonzero truth.
key_df <- cohort_df[match(call_df$cohort, cohort_df$cohort),
                    c(strata_factor_vec, "draw")]
call_df <- cbind(call_df, key_df)
.call_row <- function(sub_df){
  num_cohorts <- length(unique(sub_df$cohort))
  null_df <- sub_df[sub_df$deg == "half_null", ]
  data.frame(
    num_cohorts = num_cohorts,
    calls_per_cohort = sum(sub_df$num_calls) / num_cohorts,
    share_called = sum(sub_df$num_calls) / sum(sub_df$num_genes),
    fdp = if(nrow(null_df) > 0) sum(null_df$num_null_calls) /
      sum(null_df$num_calls) else NA_real_,
    sign_error = sum(sub_df$num_sign_errors) /
      sum(sub_df$num_calls - sub_df$num_null_calls),
    # Exact nulls and wrong signs together, over all calls.
    wrong = sum(sub_df$num_null_calls + sub_df$num_sign_errors) /
      sum(sub_df$num_calls),
    power_deg = sum(sub_df$num_deg_calls) / sum(sub_df$num_deg),
    stringsAsFactors = FALSE
  )
}
call_group_df <- unique(call_df[, c("method", "test", "status")])
call_out_df <- do.call(rbind, lapply(seq_len(nrow(call_group_df)), function(i){
  sub_df <- call_df[call_df$method == call_group_df$method[i] &
                      call_df$test == call_group_df$test[i] &
                      call_df$status == call_group_df$status[i], ]
  do.call(rbind, lapply(.strata_sets(sub_df), function(set){
    cbind(call_group_df[i, ], factor = set$factor, level = set$level,
          .call_row(sub_df[set$bool, ]), row.names = NULL)
  }))
}))
utils::write.csv(call_out_df, file = file.path(csv_dir, "bank4_calls.csv"),
                 row.names = FALSE)

# Replicates -------------------------------------------------------------------

pair_df <- unique(cohort_df[, c("id", "rank")])
replicate_list <- parallel::mclapply(seq_len(nrow(pair_df)), function(i){
  .read_replicates(id = pair_df$id[i], rank = pair_df$rank[i])
}, mc.cores = num_cores)
stopifnot(all(sapply(replicate_list, is.data.frame)))
replicate_df <- do.call(rbind, replicate_list)
rownames(replicate_df) <- NULL
replicate_df$shape <- cohort_df$shape[match(paste(replicate_df$id,
                                                  replicate_df$rank),
                                            paste(cohort_df$id,
                                                  cohort_df$rank))]
value_vec <- setdiff(colnames(replicate_df),
                     c("id", strata_factor_vec, "draw", "method"))
out_df <- replicate_df
for(name in value_vec) out_df[[name]] <- signif(out_df[[name]], 5)
utils::write.csv(out_df, file = file.path(csv_dir, "bank4_replicates.csv"),
                 row.names = FALSE)

replicate_strata_df <- do.call(rbind, lapply(method_vec, function(method){
  sub_df <- replicate_df[replicate_df$method == method, ]
  do.call(rbind, lapply(.strata_sets(sub_df), function(set){
    cbind(data.frame(method = method, factor = set$factor, level = set$level,
                     num_examples = sum(set$bool), stringsAsFactors = FALSE),
          as.data.frame(as.list(colMeans(sub_df[set$bool, value_vec,
                                                drop = FALSE],
                                         na.rm = TRUE))))
  }))
}))
utils::write.csv(replicate_strata_df,
                 file = file.path(csv_dir, "bank4_replicate_strata.csv"),
                 row.names = FALSE)

# Bank descriptives and failures -----------------------------------------------

utils::write.csv(manifest_df[, setdiff(colnames(manifest_df),
                                       c("file_truth", "error"))],
                 file = file.path(csv_dir, "bank4_cohorts.csv"),
                 row.names = FALSE)
failure_out_df <- stats::aggregate(
  cbind(num_genes = num_nonfinite, num_cohorts = num_nonfinite > 0,
        num_errors = error != "") ~ method,
  data = failure_df, FUN = sum
)
utils::write.csv(failure_out_df,
                 file = file.path(csv_dir, "bank4_failures.csv"),
                 row.names = FALSE)

print(strata_df[strata_df$factor %in% c("all", "rank", "shape", "block"),
                c("method", "factor", "level", "score", "rmse", "coverage_95",
                  "slope_deg", "vs_deseq2")], digits = 3)
print("Done! :)")
