rm(list = ls())

# Replicate cohorts: the sampling SE of eSVD-DE's log fold change
# Drafted by Claude for Kevin Z. Lin, 2026-11-07
#
# For each scenario, `num_replicates` cohorts are drawn from ONE population
# truth and each is refitted with eSVD2 1.2.0. The SD of a gene's logFC over
# the replicates is the sampling SE that Figure 4 of Writeup50 assumes a
# method's SE to be. Every candidate SE (README.md) is compared with it, by
# nuisance status, together with DESeq2 on pseudobulk as the reference.
# Writes to output/:
#   02_gene-level.csv   one row per scenario x candidate x gene
#   02_summary.csv      one row per scenario x candidate x stratum
#   02_scenario.csv     one row per scenario x replicate (null parameters,
#                       logFC slope on the truth, timing)
# Per-replicate fits are cached in cache/02_<scenario>_rep<r>.rds and reused.
# About 25 minutes on Kevin's laptop (4 s per eSVD fit, 3 s per DESeq2 fit).

library(DESeq2)

exploration_dir <- getwd()
source(file.path(exploration_dir, "helpers_claude.R"))
.load_esvd2_devel(exploration_dir)
cache_dir <- file.path(exploration_dir, "cache")
output_dir <- file.path(exploration_dir, "output")
progress_file <- file.path(output_dir, "02_progress.txt")

# Settings ---------------------------------------------------------------------

num_replicates <- 30
se_trim_quantile <- 0.95
fdr_threshold <- 0.05

# Each scenario changes one thing from the baseline: 20 donors (10 vs 10) with
# 50 cells each, 300 genes of which 30 carry +/- 0.8 (natural log), rates
# between 2 and 8, true rank 3 fitted with k = 3.
scenario_list <- list(
  baseline = list(),
  wide_nuisance = list(nuisance_range = c(0.5, 200)),
  small_cohort = list(num_individuals = 10),
  many_cells = list(num_cells_per_individual = 200),
  asymmetric = list(bool_asymmetric = TRUE),
  overfit_k = list(k_fit = 10),
  large_cohort = list(num_individuals = 40)
)
population_arg_vec <- c("bool_asymmetric", "nuisance_range")
cohort_arg_vec <- c("num_individuals", "num_cells_per_individual")

# One replicate ----------------------------------------------------------------

.run_replicate <- function(scenario_name, scenario_idx, rep_idx){
  cache_file <- file.path(cache_dir, paste0("02_", scenario_name, "_rep",
                                            rep_idx, ".rds"))
  if(file.exists(cache_file)) return(readRDS(cache_file))

  settings <- scenario_list[[scenario_name]]
  population_args <- settings[intersect(names(settings), population_arg_vec)]
  cohort_args <- settings[intersect(names(settings), cohort_arg_vec)]
  k_fit <- if(is.null(settings$k_fit)) 3 else settings$k_fit

  # The population is fixed within a scenario (its seed does not depend on
  # the replicate); only the cohort changes.
  population <- do.call(.simulate_population,
                        c(list(seed_number = 100 + scenario_idx),
                          population_args))
  sim <- do.call(.simulate_cohort,
                 c(list(population = population,
                        seed_number = 1000 * scenario_idx + rep_idx),
                   cohort_args))

  start_time <- Sys.time()
  res <- .run_esvd_pipeline(dat = sim$dat,
                            covariate_df = sim$covariate_df,
                            individual_vec = sim$individual_vec,
                            k = k_fit)
  esvd_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if(!is.null(res$error_message)){
    out <- list(error_message = res$error_message, stage = res$stage,
                scenario = scenario_name, rep = rep_idx)
    saveRDS(out, cache_file)
    return(out)
  }
  candidate_df <- .compute_candidates(res$esvd_obj)
  deseq2_df <- .run_deseq2_pseudobulk(dat = sim$dat,
                                      cc_vec = sim$cc_vec,
                                      individual_vec = sim$individual_vec)
  gene_df <- merge(candidate_df, deseq2_df, by = "gene", all.x = TRUE,
                   sort = FALSE)
  gene_df <- gene_df[match(sim$gene_vec, gene_df$gene), ]
  gene_df$truth <- as.numeric(sim$truth_log2fc_vec[gene_df$gene])
  gene_df$is_de <- as.logical(sim$is_de_vec[gene_df$gene])
  gene_df$scenario <- scenario_name
  gene_df$rep <- rep_idx

  out <- list(error_message = NULL,
              esvd_secs = esvd_secs,
              gene_df = gene_df,
              rep = rep_idx,
              scenario = scenario_name,
              warning_vec = unique(res$warning_vec))
  saveRDS(out, cache_file)
  out
}

# Running ----------------------------------------------------------------------

.append_progress(progress_file, "02 started")
replicate_list <- list()
for(scenario_idx in seq_along(scenario_list)){
  scenario_name <- names(scenario_list)[scenario_idx]
  for(rep_idx in seq_len(num_replicates)){
    out <- .run_replicate(scenario_name = scenario_name,
                          scenario_idx = scenario_idx,
                          rep_idx = rep_idx)
    replicate_list[[paste0(scenario_name, "_", rep_idx)]] <- out
    if(!is.null(out$error_message)){
      .append_progress(progress_file, paste0(scenario_name, " rep ", rep_idx,
                                             " FAILED at ", out$stage, ": ",
                                             out$error_message))
    } else if(rep_idx %% 10 == 0){
      .append_progress(progress_file, paste0(scenario_name, " rep ", rep_idx,
                                             " done (", round(out$esvd_secs, 1),
                                             " s per eSVD fit)"))
    }
  }
}

ok_list <- Filter(function(out){is.null(out$error_message)}, replicate_list)
.append_progress(progress_file, paste0(length(ok_list), " of ",
                                       length(replicate_list),
                                       " replicates fitted"))
# Candidates derived from the stored columns, for caches written before
# they were added.
all_df <- do.call(rbind, lapply(ok_list, function(out){
  .add_derived_candidates(out$gene_df)
}))
rownames(all_df) <- NULL

# Aggregating ------------------------------------------------------------------

candidate_vec <- .candidate_vec()

# Gene level: the oracle SD across replicates, the mean candidate SE, and
# the status a gene usually gets.
gene_row_list <- list()
for(scenario_name in names(scenario_list)){
  scenario_df <- all_df[all_df$scenario == scenario_name, ]
  if(nrow(scenario_df) == 0) next
  gene_split <- split(scenario_df, scenario_df$gene)
  for(candidate in candidate_vec){
    lfc_col <- paste0(candidate, "_lfc")
    se_col <- paste0(candidate, "_se")
    gene_row_list[[length(gene_row_list) + 1]] <- do.call(rbind, lapply(
      gene_split, function(df){
        lfc_vec <- df[[lfc_col]]
        se_vec <- df[[se_col]]
        bool_ok <- is.finite(lfc_vec) & is.finite(se_vec)
        data.frame(
          scenario = scenario_name,
          candidate = candidate,
          gene = df$gene[1],
          truth = df$truth[1],
          is_de = df$is_de[1],
          num_replicates = sum(bool_ok),
          frac_capped = mean(df$nuisance_status %in% c("capped", "boundary")),
          mean_shrink_weight = mean(df$shrink_weight),
          mean_count = mean(df$mean_count),
          mean_lfc = mean(lfc_vec[bool_ok]),
          oracle_sd = stats::sd(lfc_vec[bool_ok]),
          rmse = sqrt(mean((lfc_vec[bool_ok] - df$truth[1])^2)),
          mean_se = mean(se_vec[bool_ok]),
          rms_se = sqrt(mean(se_vec[bool_ok]^2)),
          stringsAsFactors = FALSE)
      }))
  }
}
gene_level_df <- do.call(rbind, gene_row_list)
rownames(gene_level_df) <- NULL
utils::write.csv(gene_level_df,
                 file = file.path(output_dir, "02_gene-level.csv"),
                 row.names = FALSE)

# Summary: per scenario x candidate x stratum. The strata use the
# per-replicate status for the pooled quantities (coverage, calibration,
# discoveries) and the usual status (`frac_capped >= 0.5`) for the per-gene
# ones (the SE-to-oracle ratio and rank correlation).
.stratum_bool <- function(status_vec, stratum){
  if(stratum == "all") return(rep(TRUE, length(status_vec)))
  if(stratum == "estimated") return(status_vec == "estimated")
  status_vec %in% c("capped", "boundary")
}
stratum_vec <- c("all", "estimated", "capped_or_boundary")

summary_row_list <- list()
for(scenario_name in names(scenario_list)){
  scenario_df <- all_df[all_df$scenario == scenario_name, ]
  if(nrow(scenario_df) == 0) next
  rep_vec <- sort(unique(scenario_df$rep))
  for(candidate in candidate_vec){
    lfc_col <- paste0(candidate, "_lfc")
    se_col <- paste0(candidate, "_se")
    pvalue_col <- paste0(candidate, "_pvalue")
    gene_df <- gene_level_df[gene_level_df$scenario == scenario_name &
                               gene_level_df$candidate == candidate, ]
    for(stratum in stratum_vec){
      bool_pooled <- .stratum_bool(scenario_df$nuisance_status, stratum) &
        is.finite(scenario_df[[lfc_col]]) & is.finite(scenario_df[[se_col]])
      pooled_df <- scenario_df[bool_pooled, ]
      usual_status <- ifelse(gene_df$frac_capped >= 0.5,
                             "capped", "estimated")
      bool_gene <- .stratum_bool(usual_status, stratum) &
        is.finite(gene_df$oracle_sd) & gene_df$oracle_sd > 0
      g_df <- gene_df[bool_gene, ]
      if(nrow(pooled_df) < 10 || nrow(g_df) < 5) next

      error_vec <- pooled_df[[lfc_col]] - pooled_df$truth
      se_vec <- pooled_df[[se_col]]
      z_null_vec <- (pooled_df[[lfc_col]] / se_vec)[!pooled_df$is_de]

      # The moments identity, replicate by replicate on the Rmd's trimmed
      # universe: rho from the SE against rho from the truth.
      moment_mat <- t(sapply(rep_vec, function(rep_idx){
        df <- pooled_df[pooled_df$rep == rep_idx, ]
        if(nrow(df) < 10) return(c(NA, NA, NA))
        bool_trim <- df[[se_col]] <= stats::quantile(df[[se_col]],
                                                     se_trim_quantile)
        df <- df[bool_trim, ]
        total_ms <- mean(df[[lfc_col]]^2)
        c(rho_moments = 1 - mean(df[[se_col]]^2) / total_ms,
          rho_truth = mean(df$truth^2) / total_ms,
          noise_ratio = mean(df[[se_col]]^2) /
            mean((df[[lfc_col]] - df$truth)^2))
      }))

      # Discoveries at BH 0.05 on the candidate's own p-value, per replicate.
      # BH is applied within the stratum here only for the `all` stratum to
      # be the method's own list; the strata then read off that list.
      discovery_mat <- t(sapply(rep_vec, function(rep_idx){
        df_all <- scenario_df[scenario_df$rep == rep_idx, ]
        bool_fin <- is.finite(df_all[[pvalue_col]])
        padj_vec <- rep(NA_real_, nrow(df_all))
        padj_vec[bool_fin] <- stats::p.adjust(df_all[[pvalue_col]][bool_fin],
                                              method = "BH")
        bool_deg <- !is.na(padj_vec) & padj_vec <= fdr_threshold
        bool_in <- .stratum_bool(df_all$nuisance_status, stratum)
        c(num_deg = sum(bool_deg & bool_in),
          num_false = sum(bool_deg & bool_in & !df_all$is_de),
          num_true = sum(bool_deg & bool_in & df_all$is_de),
          num_de = sum(bool_in & df_all$is_de))
      }))

      summary_row_list[[length(summary_row_list) + 1]] <- data.frame(
        scenario = scenario_name,
        candidate = candidate,
        stratum = stratum,
        num_genes_usual = nrow(g_df),
        frac_stratum = mean(.stratum_bool(scenario_df$nuisance_status,
                                          stratum)),
        median_se_over_oracle = stats::median(g_df$mean_se / g_df$oracle_sd),
        spearman_se_oracle = stats::cor(g_df$mean_se, g_df$oracle_sd,
                                        method = "spearman"),
        median_oracle_sd = stats::median(g_df$oracle_sd),
        median_se = stats::median(g_df$mean_se),
        calibration = sqrt(mean(error_vec^2) / mean(se_vec^2)),
        coverage_2se = mean(abs(error_vec) <= 2 * se_vec),
        sd_z_null = stats::sd(z_null_vec),
        rho_moments = mean(moment_mat[, 1], na.rm = TRUE),
        rho_truth = mean(moment_mat[, 2], na.rm = TRUE),
        noise_ratio = mean(moment_mat[, 3], na.rm = TRUE),
        mean_num_deg = mean(discovery_mat[, "num_deg"]),
        mean_num_false = mean(discovery_mat[, "num_false"]),
        mean_power = mean(discovery_mat[, "num_true"] /
                            pmax(discovery_mat[, "num_de"], 1)),
        stringsAsFactors = FALSE)
    }
  }
}
summary_df <- do.call(rbind, summary_row_list)
rownames(summary_df) <- NULL
utils::write.csv(summary_df,
                 file = file.path(output_dir, "02_summary.csv"),
                 row.names = FALSE)

# Scenario level: per replicate, the empirical null, the regression slope of
# each logFC on the truth (shrinkage toward zero would show as a slope below
# one), the median logFC of the null genes (the compositional shift), and
# the timing.
scenario_row_list <- list()
for(out in ok_list){
  df <- .add_derived_candidates(out$gene_df)
  slope <- function(lfc_vec){
    bool_ok <- is.finite(lfc_vec)
    unname(stats::coef(stats::lm(lfc_vec[bool_ok] ~ df$truth[bool_ok]))[2])
  }
  null_median <- function(lfc_vec){
    stats::median(lfc_vec[!df$is_de], na.rm = TRUE)
  }
  scenario_row_list[[length(scenario_row_list) + 1]] <- data.frame(
    scenario = out$scenario,
    rep = out$rep,
    esvd_secs = out$esvd_secs,
    num_warnings = length(out$warning_vec),
    null_mean = df$null_mean[1],
    null_sd = df$null_sd[1],
    raw_null_mean = df$raw_null_mean[1],
    raw_null_sd = df$raw_null_sd[1],
    frac_capped = mean(df$nuisance_status %in% c("capped", "boundary")),
    median_shrink_weight = stats::median(df$shrink_weight),
    slope_stored = slope(df$stored_lfc),
    slope_raw = slope(df$raw_lfc),
    slope_deseq2 = slope(df$deseq2_lfc),
    null_median_stored = null_median(df$stored_lfc),
    null_median_raw = null_median(df$raw_lfc),
    null_median_null_rescaled = null_median(df$null_rescaled_lfc),
    null_median_raw_null = null_median(df$raw_null_lfc),
    null_median_deseq2 = null_median(df$deseq2_lfc),
    rmse_stored = sqrt(mean((df$stored_lfc - df$truth)^2, na.rm = TRUE)),
    rmse_raw = sqrt(mean((df$raw_lfc - df$truth)^2, na.rm = TRUE)),
    rmse_null_rescaled = sqrt(mean((df$null_rescaled_lfc - df$truth)^2,
                                   na.rm = TRUE)),
    rmse_raw_null = sqrt(mean((df$raw_null_lfc - df$truth)^2, na.rm = TRUE)),
    rmse_deseq2 = sqrt(mean((df$deseq2_lfc - df$truth)^2, na.rm = TRUE)),
    stringsAsFactors = FALSE)
}
scenario_df <- do.call(rbind, scenario_row_list)
rownames(scenario_df) <- NULL
utils::write.csv(scenario_df,
                 file = file.path(output_dir, "02_scenario.csv"),
                 row.names = FALSE)

print(summary_df[summary_df$stratum == "all",
                 c("scenario", "candidate", "median_se_over_oracle",
                   "spearman_se_oracle", "coverage_2se", "rho_moments",
                   "rho_truth", "mean_num_false", "mean_power")],
      digits = 3)
.append_progress(progress_file, "02 finished")
print("Done! :)")
