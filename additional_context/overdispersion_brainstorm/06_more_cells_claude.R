# Does the shrinkage still calibrate the test with five times the cells?
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# The sampling variance of a gene's log rate falls with the number of cells,
# and the empirical-Bayes prior then moves the estimate less. A cap does not
# depend on the number of cells. This script repeats three regimes with 150
# cells per individual in place of 30 (3000 cells in place of 600) and scores
# a few candidates. Nothing is cached: each data set is simulated, fitted and
# scored in memory.
#
# Run from the package root:
#   Rscript additional_context/overdispersion_brainstorm/06_more_cells_claude.R
#
# Writes output/summary_more_cells.csv.

library(parallel)

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
work_dir <- file.path("additional_context", "overdispersion_brainstorm")
out_dir <- file.path(work_dir, "output")

library(eSVD2, lib.loc = file.path(cmp_dir, "lib", "devel_1.1.0"))
source(file.path(cmp_dir, "helpers_claude.R"))
source(file.path(work_dir, "helpers_candidates_claude.R"))

# The regimes ------------------------------------------------------------------

num_cells_per_individual <- 150
num_genes <- 300
num_reps <- 10
weak_effect_vec <- c(rep(c(0.3, -0.3), length.out = 30),
                     rep(0, num_genes - 30))
regime_spec_list <- list(
  null = list(num_de = 0),
  weak_de = list(cc_effect_vec = weak_effect_vec),
  wide_rate = list(cc_effect_vec = weak_effect_vec,
                   nuisance_range = c(0.5, 200))
)
candidate_vec <- c("mle", "cap_max_s", "cap_10s", "eb_map", "eb_map_cap_10s",
                   "profile_lower_90", "oracle")

job_df <- expand.grid(regime = names(regime_spec_list),
                      rep = seq_len(num_reps),
                      stringsAsFactors = FALSE)

# Scoring ----------------------------------------------------------------------

res_list <- parallel::mclapply(seq_len(nrow(job_df)), function(i){
  regime <- job_df$regime[i]
  rep_idx <- job_df$rep[i]
  # Offset by 20 so that no seed repeats one of the other scripts.
  seed_number <- 1000 * rep_idx + 20 + match(regime, names(regime_spec_list))
  sim <- do.call(.simulate_cohort,
                 c(regime_spec_list[[regime]],
                   list(num_cells_per_individual = num_cells_per_individual,
                        num_genes = num_genes,
                        seed_number = seed_number)))
  covariates <- eSVD2::format_covariates(dat = sim$dat,
                                         covariate_df = sim$covariate_df,
                                         rescale_numeric_variables = "Age")
  fit_res <- .run_pipeline(dat = sim$dat,
                           covariates = covariates,
                           cc_var = "CC_1",
                           individual_vec = sim$individual_vec)
  stopifnot(is.null(fit_res$error_message))
  esvd_obj <- fit_res$esvd_obj

  mat_list <- .fitted_matrices(esvd_obj = esvd_obj, cc_var = "CC_1")
  rate_mle_vec <- as.numeric(esvd_obj[[esvd_obj[["latest_Fit"]]]]$nuisance_vec)
  summary_df <- .gene_summaries(dat = sim$dat,
                                mean_mat = mat_list$mean_mat,
                                library_mat = mat_list$library_mat,
                                rate_mle_vec = rate_mle_vec)
  candidate_list <- .candidate_rates(dat = sim$dat,
                                     mean_mat = mat_list$mean_mat,
                                     library_mat = mat_list$library_mat,
                                     rate_mle_vec = rate_mle_vec,
                                     summary_df = summary_df,
                                     rate_true_vec = sim$nuisance_true_vec)
  # The shrunken rate with the cap as a backstop.
  candidate_list$eb_map_cap_10s <- pmin(candidate_list$eb_map,
                                        10 * summary_df$s_median)
  interior_idx <- which(!summary_df$bool_boundary & summary_df$information > 0)

  do.call(rbind, lapply(candidate_vec, function(candidate){
    rate_vec <- candidate_list[[candidate]]
    res <- .rates_to_pvalues(esvd_obj = esvd_obj, rate_vec = rate_vec)
    data.frame(candidate = candidate,
               regime = regime,
               rep = rep_idx,
               fp = sum(res$gene_df$fdr < 0.05 & !sim$is_de_vec),
               tp = sum(res$gene_df$fdr < 0.05 & sim$is_de_vec),
               type1 = mean(10^(-res$gene_df$log10p[!sim$is_de_vec]) < 0.05),
               frac_boundary = mean(summary_df$bool_boundary),
               frac_within_2fold = mean(abs(log(rate_vec / summary_df$s_median) -
                                              log(sim$nuisance_true_vec)) <
                                          log(2)),
               max_unit_free_rate = max(rate_vec / summary_df$s_median),
               prior_var = attr(candidate_list, "map_prior_var"),
               sampling_var = stats::median(1 / summary_df$information[interior_idx]))
  }))
}, mc.cores = 6)
stopifnot(!any(sapply(res_list, function(x){inherits(x, "try-error")})))
more_cells_df <- do.call(rbind, res_list)

# Tables -----------------------------------------------------------------------

summary_df <- stats::aggregate(cbind(fp, tp, type1, frac_boundary,
                                     frac_within_2fold, max_unit_free_rate, prior_var,
                                     sampling_var) ~ candidate + regime,
                               data = more_cells_df, FUN = mean)
se_df <- stats::aggregate(fp ~ candidate + regime, data = more_cells_df,
                          FUN = function(v){stats::sd(v) / sqrt(length(v))})
summary_df$fp_se <- se_df$fp
summary_df <- summary_df[order(summary_df$regime,
                               match(summary_df$candidate, candidate_vec)), ]
num_col_vec <- sapply(summary_df, is.numeric)
summary_df[num_col_vec] <- lapply(summary_df[num_col_vec], signif, digits = 3)
utils::write.csv(summary_df, file.path(out_dir, "summary_more_cells.csv"),
                 row.names = FALSE)

options(width = 200)
rownames(summary_df) <- NULL
print(summary_df)
