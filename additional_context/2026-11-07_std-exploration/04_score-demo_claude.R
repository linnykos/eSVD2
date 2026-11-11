rm(list = ls())

# Score demo for the autoresearch memo: no fitting, reads cache/ only
# Drafted by Claude for Kevin Z. Lin, 2026-11-11
#
# The autoresearch memo (additional_context/AUTORESEARCH_LFC-SE_MEMO_claude.md)
# proposes scoring a reported (logFC, SE) by the Gaussian negative
# log-likelihood of the true log2 fold change, averaged over the genes of a
# cohort. This script computes that score for the nine candidates of the
# exploration on every cached fit of scripts 02 and 03, to see (i) whether it
# orders the candidates as the exploration's separate diagnostics did,
# (ii) how large one paired standard error is, (iii) what winsorizing the
# squared z-score changes, and (iv) how much a thresholding rule gains, which
# is the exploit the memo's attenuation gate is there to stop.
# Writes to output/:
#   04_score-by-cohort.csv   one row per cached cohort x candidate
#   04_score-summary.csv     mean score per scenario x candidate
#   04_score-overall.csv     mean over all cohorts, paired against raw_null
# About 1 minute on Kevin's laptop.

exploration_dir <- getwd()
source(file.path(exploration_dir, "helpers_claude.R"))
.load_esvd2_devel(exploration_dir)
cache_dir <- file.path(exploration_dir, "cache")
output_dir <- file.path(exploration_dir, "output")

# Settings ---------------------------------------------------------------------

candidate_vec <- c("stored", "null_rescaled", "between", "between_plus",
                   "jackknife", "deshrunk", "raw", "raw_null", "deseq2",
                   "threshold")
# The squared z-score (logFC - truth)^2 / SE^2 of one gene is capped here.
# Inf is the plain score; 5 SE is the cap the memo recommends.
z_cap_vec <- c(Inf, 5)
# A weak effect, for the attenuation slope: nonzero and below this in |log2|.
weak_effect_max <- 0.5
# The thresholding rule: a gene within this many SE of zero is reported as 0.
threshold_z <- 2

# The score --------------------------------------------------------------------

# Mean over genes of 0.5 * log(2 * pi * se^2) + min(z^2, cap^2) / 2,
# z = (lfc - truth) / se, in nats per gene. Lower is better. With no cap this
# is the negative log density of the truth under N(lfc, se^2); in expectation
# over cohorts it is smallest, gene by gene, when se^2 is the mean squared
# error of lfc.
.gaussian_nll <- function(lfc_vec, se_vec, truth_vec, z_cap = Inf){
  bool_ok <- is.finite(lfc_vec) & is.finite(se_vec) & se_vec > 0
  z_sq_vec <- ((lfc_vec - truth_vec) / se_vec)^2
  nll_vec <- 0.5 * log(2 * pi * se_vec^2) + pmin(z_sq_vec, z_cap^2) / 2
  mean(nll_vec[bool_ok])
}

# The exploit a proper score invites: report 0 with a small SE for every gene
# that is not significant, and leave the significant genes alone. The small SE
# is a tenth of the gene's own, so no constant is tuned to the bank.
.add_threshold_candidate <- function(candidate_df){
  z_vec <- candidate_df$raw_null_lfc / candidate_df$raw_null_se
  bool_small <- abs(z_vec) < threshold_z
  candidate_df$threshold_lfc <- ifelse(bool_small, 0, candidate_df$raw_null_lfc)
  candidate_df$threshold_se <- ifelse(bool_small,
                                      candidate_df$raw_null_se / 10,
                                      candidate_df$raw_null_se)
  candidate_df
}

.score_cohort <- function(candidate_df, cohort_label, scenario, script){
  candidate_df <- .add_derived_candidates(candidate_df)
  candidate_df <- .add_threshold_candidate(candidate_df)
  truth_vec <- candidate_df$truth
  bool_weak <- truth_vec != 0 & abs(truth_vec) < weak_effect_max
  row_list <- lapply(candidate_vec, function(candidate){
    lfc_vec <- candidate_df[[paste0(candidate, "_lfc")]]
    se_vec <- candidate_df[[paste0(candidate, "_se")]]
    slope_weak <- NA_real_
    if(sum(bool_weak) >= 10){
      # Through the origin: the fraction of a weak effect that is reported.
      slope_weak <- sum(lfc_vec[bool_weak] * truth_vec[bool_weak]) /
        sum(truth_vec[bool_weak]^2)
    }
    data.frame(script = script,
               scenario = scenario,
               cohort = cohort_label,
               candidate = candidate,
               score = .gaussian_nll(lfc_vec, se_vec, truth_vec,
                                     z_cap = z_cap_vec[1]),
               score_capped = .gaussian_nll(lfc_vec, se_vec, truth_vec,
                                            z_cap = z_cap_vec[2]),
               rmse = sqrt(mean((lfc_vec - truth_vec)^2, na.rm = TRUE)),
               slope_weak = slope_weak,
               num_weak = sum(bool_weak))
  })
  do.call(rbind, row_list)
}

# Every cached cohort ----------------------------------------------------------

file_vec <- sort(list.files(cache_dir, pattern = "\\.rds$"))
score_list <- vector("list", length(file_vec))
for(file_idx in seq_along(file_vec)){
  out <- readRDS(file.path(cache_dir, file_vec[file_idx]))
  if(!is.null(out$error_message)) next
  script <- substr(file_vec[file_idx], 1, 2)
  label <- sub("\\.rds$", "", file_vec[file_idx])
  if(script == "02"){
    score_list[[file_idx]] <- .score_cohort(out$gene_df, label, out$scenario,
                                            script)
  } else {
    score_list[[file_idx]] <- rbind(
      .score_cohort(out$gene_df_a, paste0(label, "_a"), out$scenario, script),
      .score_cohort(out$gene_df_b, paste0(label, "_b"), out$scenario, script)
    )
  }
}
score_df <- do.call(rbind, score_list)
utils::write.csv(score_df, file.path(output_dir, "04_score-by-cohort.csv"),
                 row.names = FALSE)

# Summaries --------------------------------------------------------------------

score_df$cell <- paste(score_df$script, score_df$scenario)
summary_df <- stats::aggregate(
  cbind(score, score_capped, rmse) ~ cell + candidate,
  data = score_df, FUN = mean)
utils::write.csv(summary_df, file.path(output_dir, "04_score-summary.csv"),
                 row.names = FALSE)

# Paired against raw_null over cohorts, the comparison the keep rule makes.
reference_df <- score_df[score_df$candidate == "raw_null", ]
overall_list <- lapply(candidate_vec, function(candidate){
  df <- score_df[score_df$candidate == candidate, ]
  stopifnot(identical(df$cohort, reference_df$cohort))
  diff_vec <- df$score - reference_df$score
  diff_capped_vec <- df$score_capped - reference_df$score_capped
  num_cohorts <- nrow(df)
  bool_weak <- !is.na(df$slope_weak)
  data.frame(candidate = candidate,
             num_cohorts = num_cohorts,
             score = mean(df$score),
             delta = mean(diff_vec),
             paired_se = stats::sd(diff_vec) / sqrt(num_cohorts),
             score_capped = mean(df$score_capped),
             delta_capped = mean(diff_capped_vec),
             paired_se_capped = stats::sd(diff_capped_vec) / sqrt(num_cohorts),
             rmse = mean(df$rmse),
             slope_weak = stats::weighted.mean(df$slope_weak[bool_weak],
                                               df$num_weak[bool_weak]))
  })
overall_df <- do.call(rbind, overall_list)
utils::write.csv(overall_df, file.path(output_dir, "04_score-overall.csv"),
                 row.names = FALSE)

options(width = 200)
print(overall_df, digits = 3, row.names = FALSE)
wide_df <- stats::reshape(summary_df[, c("cell", "candidate", "score_capped")],
                          idvar = "candidate", timevar = "cell",
                          direction = "wide")
print(wide_df, digits = 3, row.names = FALSE)

print("Done! :)")
