# Scores every candidate nuisance estimator on the cached fits
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# For each cached fit and each candidate: overwrite the nuisance rate, rerun
# compute_posterior -> compute_test_statistic -> compute_pvalue, and record the
# per-gene results. Needs 01_cache_fits_claude.R, 01b_run_master_claude.R and
# the master output of additional_context/version_comparison.
#
# Run from the package root:
#   Rscript additional_context/overdispersion_brainstorm/02_run_candidates_claude.R
#
# Writes output/genes.csv (one row per gene per data set per candidate),
# output/runs.csv (one row per data set per candidate) and
# output/gene_summaries.csv (the likelihood summaries the candidates use).

library(parallel)

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
work_dir <- file.path("additional_context", "overdispersion_brainstorm")
cache_dir <- file.path(work_dir, "cache")
out_dir <- file.path(work_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(eSVD2, lib.loc = file.path(cmp_dir, "lib", "devel"))
source(file.path(work_dir, "helpers_candidates_claude.R"))

# master's own output: the six regimes of the version comparison, and the
# four added here (01b_run_master_claude.R).
column_vec <- c("regime", "rep", "gene", "nuisance")
master_gene_df <- rbind(
  utils::read.csv(file.path(cmp_dir, "output", "master",
                            "genes.csv"))[, column_vec],
  utils::read.csv(file.path(out_dir, "master_extra_genes.csv"))[, column_vec]
)

# Scoring ----------------------------------------------------------------------

stem_vec <- sub("\\.rds$", "", list.files(cache_dir, pattern = "\\.rds$"))

res_list <- parallel::mclapply(stem_vec, function(stem){
  regime <- sub("__rep.*$", "", stem)
  rep_idx <- as.integer(sub("^.*__rep", "", stem))
  cache <- readRDS(file.path(cache_dir, paste0(stem, ".rds")))
  esvd_obj <- cache$esvd_obj
  dat <- as.matrix(esvd_obj$dat)

  mat_list <- .fitted_matrices(esvd_obj = esvd_obj, cc_var = cache$cc_var)
  rate_mle_vec <- as.numeric(esvd_obj[[esvd_obj[["latest_Fit"]]]]$nuisance_vec)
  summary_df <- .gene_summaries(dat = dat,
                                mean_mat = mat_list$mean_mat,
                                library_mat = mat_list$library_mat,
                                rate_mle_vec = rate_mle_vec)

  idx_vec <- which(master_gene_df$regime == regime &
                     master_gene_df$rep == rep_idx)
  rate_master_vec <- NULL
  if(length(idx_vec) > 0){
    stopifnot(all(master_gene_df$gene[idx_vec] == cache$gene_vec))
    rate_master_vec <- master_gene_df$nuisance[idx_vec]
  }
  candidate_list <- .candidate_rates(dat = dat,
                                     mean_mat = mat_list$mean_mat,
                                     library_mat = mat_list$library_mat,
                                     rate_mle_vec = rate_mle_vec,
                                     summary_df = summary_df,
                                     rate_master_vec = rate_master_vec,
                                     rate_true_vec = cache$nuisance_true_vec)

  # Two candidates change the rescaling and not the rate.
  stabilize_vec <- stats::setNames(rep(FALSE, length(candidate_list)),
                                   names(candidate_list))
  candidate_list$mle_median_stabilize <- candidate_list$mle
  candidate_list$cap_50s_median_stabilize <- candidate_list$cap_50s
  stabilize_vec <- c(stabilize_vec,
                     mle_median_stabilize = TRUE,
                     cap_50s_median_stabilize = TRUE)

  gene_list <- list()
  run_list <- list()
  for(candidate in names(candidate_list)){
    res <- .rates_to_pvalues(esvd_obj = esvd_obj,
                             rate_vec = candidate_list[[candidate]],
                             bool_median_stabilize = stabilize_vec[[candidate]])
    gene_list[[candidate]] <- cbind(
      data.frame(candidate = candidate,
                 regime = regime,
                 rep = rep_idx,
                 gene = cache$gene_vec,
                 is_de = cache$is_de_vec,
                 nuisance = candidate_list[[candidate]]),
      res$gene_df
    )
    run_list[[candidate]] <- data.frame(candidate = candidate,
                                        regime = regime,
                                        rep = rep_idx,
                                        method = res$method,
                                        null_mean = res$null_mean,
                                        null_sd = res$null_sd)
  }

  nuisance_true_vec <- cache$nuisance_true_vec
  if(is.null(nuisance_true_vec)) nuisance_true_vec <- rep(NA, ncol(dat))
  summary_df <- cbind(data.frame(regime = regime,
                                 rep = rep_idx,
                                 gene = cache$gene_vec,
                                 is_de = cache$is_de_vec,
                                 mean_count = colMeans(dat),
                                 nuisance_true = as.numeric(nuisance_true_vec),
                                 rate_mle = rate_mle_vec,
                                 truth_log2fc = as.numeric(cache$truth_log2fc_vec),
                                 map_centre = attr(candidate_list, "map_centre"),
                                 map_prior_var = attr(candidate_list,
                                                      "map_prior_var")),
                      summary_df)

  list(gene_df = do.call(rbind, gene_list),
       run_df = do.call(rbind, run_list),
       summary_df = summary_df)
}, mc.cores = 8)

bool_failed_vec <- sapply(res_list, function(x){inherits(x, "try-error")})
if(any(bool_failed_vec)){
  print(res_list[[which(bool_failed_vec)[1]]])
  stop(sum(bool_failed_vec), " data set(s) failed, the first being ",
       stem_vec[which(bool_failed_vec)[1]])
}

.write_rounded <- function(df, file){
  num_col_vec <- sapply(df, is.numeric)
  df[num_col_vec] <- lapply(df[num_col_vec], signif, digits = 6)
  utils::write.csv(df, file, row.names = FALSE)
}
.write_rounded(do.call(rbind, lapply(res_list, function(x){x$gene_df})),
               file.path(out_dir, "genes.csv"))
.write_rounded(do.call(rbind, lapply(res_list, function(x){x$run_df})),
               file.path(out_dir, "runs.csv"))
.write_rounded(do.call(rbind, lapply(res_list, function(x){x$summary_df})),
               file.path(out_dir, "gene_summaries.csv"))
print(paste0("Scored ", length(stem_vec), " data sets"))
