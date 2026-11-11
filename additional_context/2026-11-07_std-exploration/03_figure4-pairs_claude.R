rm(list = ls())

# Pairs of cohorts: does each candidate triple recover Figure 4's terms?
# Drafted by Claude for Kevin Z. Lin, 2026-11-07
#
# Two populations A and B share everything but their case-control effects,
# which have a known size ratio `k_true` = S_A / S_B and a known correlation
# `r_true` on the DE genes. For each replicate pair a cohort is drawn from
# each, both are refitted, and the Writeup50 decomposition
# (slope, rho_A, k, r) is computed under every candidate triple and under
# DESeq2, on the Rmd's 95%-SE-trimmed universe, then compared with the truth.
# The sign agreement of A's DEGs in B (Figure 4a / the sign-agreement figure)
# is computed alongside.
# Writes to output/:
#   03_decomposition.csv   one row per scenario x replicate x candidate
#   03_summary.csv         mean and SD over replicates, with the truth
# Per-pair fits are cached in cache/03_<scenario>_rep<r>.rds.
# About 10 minutes on Kevin's laptop.

library(DESeq2)

exploration_dir <- getwd()
source(file.path(exploration_dir, "helpers_claude.R"))
.load_esvd2_devel(exploration_dir)
cache_dir <- file.path(exploration_dir, "cache")
output_dir <- file.path(exploration_dir, "output")
progress_file <- file.path(output_dir, "03_progress.txt")

# Settings ---------------------------------------------------------------------

num_replicates <- 20
se_trim_quantile <- 0.95
fdr_threshold <- 0.05
k_true <- 2
r_true <- 0.8

# A: 20 donors. B: 20 donors (equal), 10 donors (Prater-like), or 20 donors
# with rates spanning 0.5 to 200 (many capped genes, as in the real cohorts).
# A fourth scenario weakens the effects (see below).
# `weak_effects` halves the effect SD so that the realized reliability rho_A
# is near the real cohorts' 0.4 to 0.65 rather than the 0.9 of the others.
scenario_list <- list(
  equal = list(num_individuals_b = 20),
  small_b = list(num_individuals_b = 10),
  wide_nuisance = list(num_individuals_b = 20,
                       nuisance_range = c(0.5, 200)),
  weak_effects = list(num_individuals_b = 20, effect_sd_a = 0.3)
)

# One replicate pair -----------------------------------------------------------

.fit_cohort <- function(population, seed_number, num_individuals){
  sim <- .simulate_cohort(population = population,
                          num_individuals = num_individuals,
                          seed_number = seed_number)
  res <- .run_esvd_pipeline(dat = sim$dat,
                            covariate_df = sim$covariate_df,
                            individual_vec = sim$individual_vec,
                            k = 3)
  if(!is.null(res$error_message)) return(res)
  candidate_df <- .compute_candidates(res$esvd_obj)
  deseq2_df <- .run_deseq2_pseudobulk(dat = sim$dat,
                                      cc_vec = sim$cc_vec,
                                      individual_vec = sim$individual_vec)
  gene_df <- merge(candidate_df, deseq2_df, by = "gene", all.x = TRUE,
                   sort = FALSE)
  gene_df <- gene_df[match(sim$gene_vec, gene_df$gene), ]
  gene_df$truth <- as.numeric(sim$truth_log2fc_vec[gene_df$gene])
  list(error_message = NULL, gene_df = gene_df)
}

.run_pair <- function(scenario_name, scenario_idx, rep_idx){
  cache_file <- file.path(cache_dir, paste0("03_", scenario_name, "_rep",
                                            rep_idx, ".rds"))
  if(file.exists(cache_file)) return(readRDS(cache_file))
  settings <- scenario_list[[scenario_name]]
  pair_args <- list(k_true = k_true, r_true = r_true,
                    seed_number = 200 + scenario_idx)
  if(!is.null(settings$nuisance_range)){
    pair_args$nuisance_range <- settings$nuisance_range
  }
  if(!is.null(settings$effect_sd_a)) pair_args$effect_sd_a <- settings$effect_sd_a
  populations <- do.call(.simulate_effect_pair, pair_args)
  seed_base <- 100000 * scenario_idx + 10 * rep_idx
  fit_a <- .fit_cohort(populations$population_a, seed_base + 1, 20)
  fit_b <- .fit_cohort(populations$population_b, seed_base + 2,
                       settings$num_individuals_b)
  out <- list(scenario = scenario_name, rep = rep_idx,
              error_message = c(fit_a$error_message, fit_b$error_message),
              gene_df_a = fit_a$gene_df, gene_df_b = fit_b$gene_df)
  saveRDS(out, cache_file)
  out
}

# Running ----------------------------------------------------------------------

.append_progress(progress_file, "03 started")
pair_list <- list()
for(scenario_idx in seq_along(scenario_list)){
  scenario_name <- names(scenario_list)[scenario_idx]
  for(rep_idx in seq_len(num_replicates)){
    out <- .run_pair(scenario_name, scenario_idx, rep_idx)
    pair_list[[paste0(scenario_name, "_", rep_idx)]] <- out
    if(length(out$error_message) > 0){
      .append_progress(progress_file, paste0(scenario_name, " rep ", rep_idx,
                                             " FAILED: ",
                                             paste(out$error_message,
                                                   collapse = " | ")))
    } else if(rep_idx %% 5 == 0){
      .append_progress(progress_file, paste0(scenario_name, " rep ", rep_idx,
                                             " done"))
    }
  }
}
ok_list <- Filter(function(out){length(out$error_message) == 0}, pair_list)
.append_progress(progress_file, paste0(length(ok_list), " of ",
                                       length(pair_list), " pairs fitted"))

# Decomposing ------------------------------------------------------------------

candidate_vec <- .candidate_vec()
decomposition_row_list <- list()
for(out in ok_list){
  # Candidates derived from the stored columns, for caches written before
  # they were added.
  a_df <- .add_derived_candidates(out$gene_df_a)
  b_df <- .add_derived_candidates(out$gene_df_b)
  stopifnot(identical(a_df$gene, b_df$gene))
  for(candidate in candidate_vec){
    lfc_col <- paste0(candidate, "_lfc")
    se_col <- paste0(candidate, "_se")
    pvalue_col <- paste0(candidate, "_pvalue")
    bool_universe <- .moments_universe(a_df[[lfc_col]], b_df[[lfc_col]],
                                       a_df[[se_col]], b_df[[se_col]],
                                       quantile_level = se_trim_quantile)
    decomposition <- .decompose(lfc_a_vec = a_df[[lfc_col]][bool_universe],
                                lfc_b_vec = b_df[[lfc_col]][bool_universe],
                                se_a_vec = a_df[[se_col]][bool_universe],
                                se_b_vec = b_df[[se_col]][bool_universe])
    truth <- .decompose_truth(truth_a_vec = a_df$truth[bool_universe],
                              truth_b_vec = b_df$truth[bool_universe],
                              lfc_a_vec = a_df[[lfc_col]][bool_universe])

    # A's DEGs at BH 0.05 on the candidate's own p-value, within the
    # untrimmed converged universe, as the sign-agreement figure does.
    bool_converged <- is.finite(a_df[[lfc_col]]) & is.finite(b_df[[lfc_col]]) &
      is.finite(a_df[[pvalue_col]])
    padj_vec <- rep(NA_real_, nrow(a_df))
    padj_vec[bool_converged] <- stats::p.adjust(
      a_df[[pvalue_col]][bool_converged], method = "BH")
    idx_discovery <- which(!is.na(padj_vec) & padj_vec <= fdr_threshold)
    bool_true_de <- a_df$truth[idx_discovery] != 0

    decomposition_row_list[[length(decomposition_row_list) + 1]] <- cbind(
      data.frame(scenario = out$scenario, rep = out$rep,
                 candidate = candidate, stringsAsFactors = FALSE),
      decomposition,
      data.frame(
        truth_signal_rms_a = truth$signal_rms_a,
        truth_signal_rms_b = truth$signal_rms_b,
        truth_rho_a = truth$rho_a,
        truth_k = truth$k,
        truth_r = truth$r,
        truth_slope = truth$slope_truth,
        num_discovery = length(idx_discovery),
        frac_discovery_true = mean(bool_true_de),
        frac_same_sign = if(length(idx_discovery) > 0){
          mean(sign(a_df[[lfc_col]][idx_discovery]) ==
                 sign(b_df[[lfc_col]][idx_discovery]))
        } else NA_real_,
        # The same sign check against the TRUTH of B, which is what the
        # replicate's sign is a noisy version of.
        frac_same_sign_truth_b = if(length(idx_discovery) > 0){
          mean(sign(a_df[[lfc_col]][idx_discovery]) ==
                 sign(b_df$truth[idx_discovery]))
        } else NA_real_))
  }
}
decomposition_df <- do.call(rbind, decomposition_row_list)
rownames(decomposition_df) <- NULL
utils::write.csv(decomposition_df,
                 file = file.path(output_dir, "03_decomposition.csv"),
                 row.names = FALSE)

# Summarising ------------------------------------------------------------------

term_vec <- c("slope", "rho_a", "k", "r", "num_discovery",
              "frac_discovery_true", "frac_same_sign",
              "frac_same_sign_truth_b", "noise_rms_a", "signal_rms_a")
summary_row_list <- list()
for(scenario_name in names(scenario_list)){
  for(candidate in candidate_vec){
    df <- decomposition_df[decomposition_df$scenario == scenario_name &
                             decomposition_df$candidate == candidate, ]
    if(nrow(df) == 0) next
    row <- data.frame(scenario = scenario_name, candidate = candidate,
                      num_replicates = nrow(df),
                      frac_valid = mean(df$bool_valid),
                      stringsAsFactors = FALSE)
    for(term in term_vec){
      row[[paste0(term, "_mean")]] <- mean(df[[term]], na.rm = TRUE)
      row[[paste0(term, "_sd")]] <- stats::sd(df[[term]], na.rm = TRUE)
    }
    row$truth_rho_a <- mean(df$truth_rho_a)
    row$truth_k <- mean(df$truth_k)
    row$truth_r <- mean(df$truth_r)
    row$truth_slope <- mean(df$truth_slope)
    row$truth_signal_rms_a <- mean(df$truth_signal_rms_a)
    summary_row_list[[length(summary_row_list) + 1]] <- row
  }
}
summary_df <- do.call(rbind, summary_row_list)
rownames(summary_df) <- NULL
utils::write.csv(summary_df,
                 file = file.path(output_dir, "03_summary.csv"),
                 row.names = FALSE)

print(summary_df[, c("scenario", "candidate", "frac_valid", "slope_mean",
                     "truth_slope", "rho_a_mean", "truth_rho_a", "k_mean",
                     "truth_k", "r_mean", "truth_r", "frac_same_sign_mean")],
      digits = 3)
.append_progress(progress_file, "03 finished")
print("Done! :)")
