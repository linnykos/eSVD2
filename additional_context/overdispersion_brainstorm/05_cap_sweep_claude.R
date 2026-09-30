# How tight does the cap on the rate have to be?
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Sweeps the multiplier c of the cap beta_j <= c * median_i(s_ji) on the devel
# maximum-likelihood rates, from c = 1 to no cap, and the prior variance of
# the empirical-Bayes MAP rate, from strong to weak shrinkage.
#
# Run from the package root:
#   Rscript additional_context/overdispersion_brainstorm/05_cap_sweep_claude.R
#
# Writes output/summary_cap_sweep.csv.

library(parallel)

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
work_dir <- file.path("additional_context", "overdispersion_brainstorm")
cache_dir <- file.path(work_dir, "cache")
out_dir <- file.path(work_dir, "output")

library(eSVD2, lib.loc = file.path(cmp_dir, "lib", "devel_1.1.0"))
source(file.path(work_dir, "helpers_candidates_claude.R"))

multiplier_vec <- c(1, 2, 5, 10, 20, 50, 200, 1000, Inf)
prior_var_vec <- c(0.05, 0.25, 1, 4)
summary_gene_df <- utils::read.csv(file.path(out_dir, "gene_summaries.csv"))

# Scoring ----------------------------------------------------------------------

stem_vec <- sub("\\.rds$", "", list.files(cache_dir, pattern = "\\.rds$"))

res_list <- parallel::mclapply(stem_vec, function(stem){
  regime <- sub("__rep.*$", "", stem)
  rep_idx <- as.integer(sub("^.*__rep", "", stem))
  cache <- readRDS(file.path(cache_dir, paste0(stem, ".rds")))
  esvd_obj <- cache$esvd_obj
  dat <- as.matrix(esvd_obj$dat)
  x <- summary_gene_df[summary_gene_df$regime == regime &
                         summary_gene_df$rep == rep_idx, ]
  stopifnot(all(x$gene == cache$gene_vec))

  rate_list <- lapply(multiplier_vec, function(multiplier){
    pmin(x$rate_mle, multiplier * x$s_median)
  })
  names(rate_list) <- paste0("cap_", multiplier_vec)

  mat_list <- .fitted_matrices(esvd_obj = esvd_obj, cc_var = cache$cc_var)
  for(prior_var in prior_var_vec){
    map_res <- .map_rates(dat = dat,
                          mean_mat = mat_list$mean_mat,
                          library_mat = mat_list$library_mat,
                          rate_mle_vec = x$rate_mle,
                          summary_df = x,
                          prior_var = prior_var)
    rate_list[[paste0("map_", prior_var)]] <- map_res$rate_vec
  }

  do.call(rbind, lapply(names(rate_list), function(candidate){
    res <- .rates_to_pvalues(esvd_obj = esvd_obj,
                             rate_vec = rate_list[[candidate]])
    data.frame(candidate = candidate,
               regime = regime,
               rep = rep_idx,
               fp = sum(res$gene_df$fdr < 0.05 & !cache$is_de_vec),
               tp = sum(res$gene_df$fdr < 0.05 & cache$is_de_vec),
               type1 = mean(10^(-res$gene_df$log10p[!cache$is_de_vec]) < 0.05),
               frac_capped = mean(rate_list[[candidate]] < x$rate_mle))
  }))
}, mc.cores = 8)
stopifnot(!any(sapply(res_list, function(x){inherits(x, "try-error")})))
sweep_df <- do.call(rbind, res_list)

# Tables -----------------------------------------------------------------------

summary_df <- stats::aggregate(cbind(fp, tp, type1, frac_capped) ~
                                 candidate + regime,
                               data = sweep_df, FUN = mean)
se_df <- stats::aggregate(cbind(fp, tp) ~ candidate + regime, data = sweep_df,
                          FUN = function(v){stats::sd(v) / sqrt(length(v))})
stopifnot(all(se_df$candidate == summary_df$candidate),
          all(se_df$regime == summary_df$regime))
summary_df$fp_se <- se_df$fp
summary_df$tp_se <- se_df$tp
utils::write.csv(summary_df, file.path(out_dir, "summary_cap_sweep.csv"),
                 row.names = FALSE)

options(width = 200)
candidate_vec <- unique(sweep_df$candidate)
for(column in c("fp", "fp_se", "tp", "tp_se", "type1", "frac_capped")){
  res <- stats::reshape(summary_df[, c("candidate", "regime", column)],
                        idvar = "candidate", timevar = "regime",
                        direction = "wide")
  colnames(res) <- sub(paste0("^", column, "\\."), "", colnames(res))
  res <- res[match(candidate_vec, res$candidate), ]
  res[, -1] <- round(res[, -1], 2)
  rownames(res) <- NULL
  print(column)
  print(res)
}
