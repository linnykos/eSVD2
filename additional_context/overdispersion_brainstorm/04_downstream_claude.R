# Two fixes that act on the test statistic and leave the rate alone
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# 02_run_candidates_claude.R shows that even the TRUE rates inflate the false
# discoveries when the rates differ a lot between genes, so part of the
# problem is downstream of the estimate. This script tries two statistics on a
# few of the rate vectors that script wrote:
#
#   between:      Welch on the individuals' mean posterior, with the variance
#                 BETWEEN individuals only (Bessel-corrected), in place of the
#                 mixture variance, which adds the per-cell posterior variance
#   scale_trend:  the package's statistic, divided by a robust scale computed
#                 within tenths of the rate, so that the single empirical null
#                 sees statistics of one scale
#
# Run from the package root:
#   Rscript additional_context/overdispersion_brainstorm/04_downstream_claude.R
#
# Writes output/summary_downstream.csv.

library(parallel)

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
work_dir <- file.path("additional_context", "overdispersion_brainstorm")
cache_dir <- file.path(work_dir, "cache")
out_dir <- file.path(work_dir, "output")

library(eSVD2, lib.loc = file.path(cmp_dir, "lib", "devel_1.1.0"))

candidate_gene_df <- utils::read.csv(file.path(out_dir, "genes.csv"))
rate_candidate_vec <- c("mle", "cap_max_s", "cap_50s", "common_unit", "oracle")
candidate_gene_df <- candidate_gene_df[candidate_gene_df$candidate %in%
                                         rate_candidate_vec, ]

# The statistics ---------------------------------------------------------------

.muffle <- function(expr){
  withCallingHandlers(expr,
                      warning = function(w){invokeRestart("muffleWarning")})
}

.between_statistic <- function(esvd_obj){
  latest_fit <- esvd_obj[["latest_Fit"]]
  posterior_mean_mat <- esvd_obj[[latest_fit]]$posterior_mean_mat
  individual_vec <- esvd_obj$individual
  indiv_mean_mat <- rowsum(posterior_mean_mat, group = individual_vec) /
    as.numeric(table(individual_vec)[levels(individual_vec)])
  cc_by_indiv_vec <- tapply(esvd_obj$case_control, individual_vec, unique)
  cc_by_indiv_vec <- cc_by_indiv_vec[rownames(indiv_mean_mat)]

  case_mat <- indiv_mean_mat[cc_by_indiv_vec == 1, , drop = FALSE]
  control_mat <- indiv_mean_mat[cc_by_indiv_vec == 0, , drop = FALSE]
  n1 <- nrow(case_mat)
  n0 <- nrow(control_mat)
  case_term_vec <- apply(case_mat, 2, stats::var) / n1
  control_term_vec <- apply(control_mat, 2, stats::var) / n0
  teststat_vec <- (colMeans(case_mat) - colMeans(control_mat)) /
    sqrt(case_term_vec + control_term_vec)
  df_vec <- (case_term_vec + control_term_vec)^2 /
    (case_term_vec^2 / (n1 - 1) + control_term_vec^2 / (n0 - 1))

  eSVD2:::.t_to_gaussian(teststat_vec = teststat_vec, df_vec = df_vec)
}

.scale_trend_statistic <- function(gaussian_vec,
                                   rate_vec,
                                   num_groups = 10){
  group_vec <- ceiling(num_groups * rank(rate_vec, ties.method = "first") /
                         length(rate_vec))
  scale_vec <- tapply(gaussian_vec, group_vec, function(v){
    stats::IQR(v) / (2 * stats::qnorm(0.75))
  })
  centre_vec <- tapply(gaussian_vec, group_vec, stats::median)
  (gaussian_vec - centre_vec[group_vec]) / scale_vec[group_vec] *
    stats::median(scale_vec)
}

# Scoring ----------------------------------------------------------------------

stem_vec <- sub("\\.rds$", "", list.files(cache_dir, pattern = "\\.rds$"))

res_list <- parallel::mclapply(stem_vec, function(stem){
  regime <- sub("__rep.*$", "", stem)
  rep_idx <- as.integer(sub("^.*__rep", "", stem))
  cache <- readRDS(file.path(cache_dir, paste0(stem, ".rds")))
  esvd_obj <- cache$esvd_obj
  latest_fit <- esvd_obj[["latest_Fit"]]

  res <- list()
  for(candidate in rate_candidate_vec){
    x <- candidate_gene_df[candidate_gene_df$candidate == candidate &
                             candidate_gene_df$regime == regime &
                             candidate_gene_df$rep == rep_idx, ]
    if(nrow(x) == 0) next()
    stopifnot(all(x$gene == cache$gene_vec))

    esvd_obj[[latest_fit]]$nuisance_vec[] <- x$nuisance
    esvd_obj <- .muffle(eSVD2::compute_posterior(
      input_obj = esvd_obj,
      alpha_max = 2 * max(esvd_obj$dat),
      bool_covariates_as_library = TRUE,
      library_min = 0.1
    ))

    statistic_list <- list(
      between = .between_statistic(esvd_obj),
      scale_trend = .scale_trend_statistic(gaussian_vec = x$gaussian_teststat,
                                           rate_vec = x$nuisance)
    )
    for(statistic in names(statistic_list)){
      z_vec <- as.numeric(statistic_list[[statistic]])
      multtest_res <- .muffle(eSVD2::multtest(z_vec))
      theory_p_vec <- 2 * stats::pnorm(-abs(z_vec))
      res[[paste(candidate, statistic)]] <- data.frame(
        candidate = candidate,
        statistic = statistic,
        regime = regime,
        rep = rep_idx,
        is_de = cache$is_de_vec,
        z = z_vec,
        p_empirical = as.numeric(multtest_res$pvalue_vec),
        fdr_empirical = as.numeric(multtest_res$fdr_vec),
        p_theory = theory_p_vec,
        fdr_theory = stats::p.adjust(theory_p_vec, method = "BH")
      )
    }
  }
  do.call(rbind, res)
}, mc.cores = 8)
stopifnot(!any(sapply(res_list, function(x){inherits(x, "try-error")})))
downstream_df <- do.call(rbind, res_list)

# Tables -----------------------------------------------------------------------

num_reps <- length(unique(downstream_df$rep))
key_df <- unique(downstream_df[, c("candidate", "statistic", "regime")])
summary_df <- do.call(rbind, lapply(seq_len(nrow(key_df)), function(i){
  x <- downstream_df[downstream_df$candidate == key_df$candidate[i] &
                       downstream_df$statistic == key_df$statistic[i] &
                       downstream_df$regime == key_df$regime[i], ]
  data.frame(key_df[i, ],
             null_sd = stats::sd(x$z[!x$is_de]),
             type1_empirical = mean(x$p_empirical[!x$is_de] < 0.05),
             fp_empirical = sum(x$fdr_empirical < 0.05 & !x$is_de) / num_reps,
             tp_empirical = sum(x$fdr_empirical < 0.05 & x$is_de) / num_reps,
             type1_theory = mean(x$p_theory[!x$is_de] < 0.05),
             fp_theory = sum(x$fdr_theory < 0.05 & !x$is_de) / num_reps,
             tp_theory = sum(x$fdr_theory < 0.05 & x$is_de) / num_reps)
}))
utils::write.csv(summary_df, file.path(out_dir, "summary_downstream.csv"),
                 row.names = FALSE)

options(width = 200)
for(regime in c("null", "baseline", "near_poisson", "weak_de", "low_count",
                "wide_rate", "generate_null")){
  print(regime)
  x <- summary_df[summary_df$regime == regime, -3]
  x[, -(1:2)] <- round(x[, -(1:2)], 3)
  rownames(x) <- NULL
  print(x[order(x$statistic), ])
}
