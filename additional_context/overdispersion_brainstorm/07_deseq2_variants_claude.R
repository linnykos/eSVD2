# Is each candidate's rate a better ESTIMATE than master's, and do DESeq2's
# two design choices help the shrinkage?
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Part 1 scores the rates themselves against the truth, in the fit's units
# (the rate divided by the gene's median fitted library size), for the
# candidates of 02_run_candidates_claude.R. Part 2 adds three variants of the
# empirical-Bayes rate: a centre that is a fitted trend in the gene's mean
# count, DESeq2's trigamma((m - p) / 2) as the sampling variance, and both.
#
# Run from the package root, after 02_run_candidates_claude.R:
#   Rscript additional_context/overdispersion_brainstorm/07_deseq2_variants_claude.R
#
# Writes output/summary_accuracy.csv and output/summary_deseq2_variants.csv.

library(parallel)

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
work_dir <- file.path("additional_context", "overdispersion_brainstorm")
cache_dir <- file.path(work_dir, "cache")
out_dir <- file.path(work_dir, "output")

library(eSVD2, lib.loc = file.path(cmp_dir, "lib", "devel"))
source(file.path(work_dir, "helpers_candidates_claude.R"))

gene_df <- utils::read.csv(file.path(out_dir, "genes.csv"))
summary_gene_df <- utils::read.csv(file.path(out_dir, "gene_summaries.csv"))
regime_vec <- c("null", "baseline", "weak_de", "low_count", "wide_rate",
                "trend", "near_poisson", "generate_null")

# Part 2 first: the variants, so that part 1 can score them too ---------------

stem_vec <- sub("\\.rds$", "", list.files(cache_dir, pattern = "\\.rds$"))
stem_vec <- stem_vec[sub("__rep.*$", "", stem_vec) %in% regime_vec]

res_list <- parallel::mclapply(stem_vec, function(stem){
  regime <- sub("__rep.*$", "", stem)
  rep_idx <- as.integer(sub("^.*__rep", "", stem))
  cache <- readRDS(file.path(cache_dir, paste0(stem, ".rds")))
  esvd_obj <- cache$esvd_obj
  dat <- as.matrix(esvd_obj$dat)
  x <- summary_gene_df[summary_gene_df$regime == regime &
                         summary_gene_df$rep == rep_idx, ]
  stopifnot(all(x$gene == cache$gene_vec))
  mat_list <- .fitted_matrices(esvd_obj = esvd_obj, cc_var = cache$cc_var)

  # m cells, and p coefficients per gene: the latent dimensions and the
  # covariates.
  fit <- esvd_obj[[esvd_obj[["latest_Fit"]]]]
  num_coefficients <- ncol(fit$y_mat) + ncol(fit$z_mat)
  trigamma_val <- trigamma((nrow(dat) - num_coefficients) / 2)

  variant_df <- data.frame(
    candidate = c("eb_map", "eb_trend", "eb_trigamma", "eb_trend_trigamma"),
    bool_trend = c(FALSE, TRUE, FALSE, TRUE),
    bool_trigamma = c(FALSE, FALSE, TRUE, TRUE)
  )
  do.call(rbind, lapply(seq_len(nrow(variant_df)), function(i){
    sampling_var <- NULL
    if(variant_df$bool_trigamma[i]) sampling_var <- trigamma_val
    map_res <- .map_rates(dat = dat,
                          mean_mat = mat_list$mean_mat,
                          library_mat = mat_list$library_mat,
                          rate_mle_vec = x$rate_mle,
                          summary_df = x,
                          bool_trend = variant_df$bool_trend[i],
                          sampling_var = sampling_var)
    res <- .rates_to_pvalues(esvd_obj = esvd_obj, rate_vec = map_res$rate_vec)
    data.frame(candidate = variant_df$candidate[i],
               regime = regime,
               rep = rep_idx,
               gene = cache$gene_vec,
               is_de = cache$is_de_vec,
               nuisance = map_res$rate_vec,
               fdr = res$gene_df$fdr,
               prior_var = map_res$prior_var,
               sampling_var = map_res$sampling_var,
               slope = map_res$slope)
  }))
}, mc.cores = 8)
stopifnot(!any(sapply(res_list, function(x){inherits(x, "try-error")})))
variant_gene_df <- do.call(rbind, res_list)

# The prototype must reproduce itself through the new code path.
check_df <- merge(variant_gene_df[variant_gene_df$candidate == "eb_map",
                                  c("regime", "rep", "gene", "nuisance")],
                  gene_df[gene_df$candidate == "eb_map",
                          c("regime", "rep", "gene", "nuisance")],
                  by = c("regime", "rep", "gene"))
stopifnot(nrow(check_df) == sum(variant_gene_df$candidate == "eb_map"),
          max(abs(check_df$nuisance.x / check_df$nuisance.y - 1)) < 1e-3)

num_reps <- length(unique(variant_gene_df$rep))
variant_key_df <- unique(variant_gene_df[, c("candidate", "regime")])
variant_summary_df <- do.call(rbind, lapply(seq_len(nrow(variant_key_df)),
                                            function(i){
  x <- variant_gene_df[variant_gene_df$candidate == variant_key_df$candidate[i] &
                         variant_gene_df$regime == variant_key_df$regime[i], ]
  fp_vec <- tapply(x$fdr < 0.05 & !x$is_de, x$rep, sum)
  data.frame(variant_key_df[i, ],
             fp = mean(fp_vec),
             fp_se = stats::sd(fp_vec) / sqrt(length(fp_vec)),
             tp = sum(x$fdr < 0.05 & x$is_de) / num_reps,
             prior_var = mean(x$prior_var),
             sampling_var = mean(x$sampling_var),
             slope = mean(x$slope))
}))
utils::write.csv(variant_summary_df,
                 file.path(out_dir, "summary_deseq2_variants.csv"),
                 row.names = FALSE)

# Part 1: the rate as an estimate ----------------------------------------------

candidate_vec <- c("master", "cap_max_s", "mle", "cap_10s", "eb_map",
                   "eb_trend", "eb_trigamma", "eb_trend_trigamma",
                   "profile_lower_90", "common_unit")
all_gene_df <- rbind(gene_df[gene_df$candidate %in% candidate_vec,
                             c("candidate", "regime", "rep", "gene",
                               "nuisance")],
                     variant_gene_df[variant_gene_df$candidate != "eb_map",
                                     c("candidate", "regime", "rep", "gene",
                                       "nuisance")])
all_gene_df <- merge(all_gene_df,
                     summary_gene_df[, c("regime", "rep", "gene", "s_median",
                                         "nuisance_true")],
                     by = c("regime", "rep", "gene"))
all_gene_df <- all_gene_df[!is.na(all_gene_df$nuisance_true), ]
all_gene_df$log_error <- log(all_gene_df$nuisance / all_gene_df$s_median) -
  log(all_gene_df$nuisance_true)

accuracy_key_df <- unique(all_gene_df[, c("candidate", "regime")])
accuracy_df <- do.call(rbind, lapply(seq_len(nrow(accuracy_key_df)),
                                     function(i){
  x <- all_gene_df[all_gene_df$candidate == accuracy_key_df$candidate[i] &
                     all_gene_df$regime == accuracy_key_df$regime[i], ]
  data.frame(accuracy_key_df[i, ],
             spearman = mean(sapply(split(x, x$rep), function(y){
               stats::cor(y$nuisance / y$s_median, y$nuisance_true,
                          method = "spearman")
             })),
             median_ratio = exp(stats::median(x$log_error)),
             median_abs_log_error = stats::median(abs(x$log_error)),
             frac_within_2fold = mean(abs(x$log_error) < log(2)))
}))
utils::write.csv(accuracy_df, file.path(out_dir, "summary_accuracy.csv"),
                 row.names = FALSE)

# Tables -----------------------------------------------------------------------

options(width = 200)
.print_wide <- function(df, column, row_vec, digits){
  regime_present_vec <- intersect(regime_vec, df$regime)
  res <- sapply(regime_present_vec, function(regime){
    sapply(row_vec, function(candidate){
      val <- df[df$candidate == candidate & df$regime == regime, column]
      if(length(val) == 0) return(NA)
      val
    })
  })
  print(column)
  print(round(res, digits))
}
for(column in c("spearman", "median_ratio", "median_abs_log_error",
                "frac_within_2fold")){
  .print_wide(accuracy_df, column, candidate_vec, 2)
}
variant_vec <- c("eb_map", "eb_trend", "eb_trigamma", "eb_trend_trigamma")
for(column in c("fp", "fp_se", "tp", "prior_var", "sampling_var", "slope")){
  .print_wide(variant_summary_df, column, variant_vec, 3)
}
