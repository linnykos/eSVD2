# Compact tables for overdispersion_cap_claude.Rmd
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# The report reads only small, tracked CSVs. This script makes them from the
# large per-gene tables of 02_run_candidates_claude.R (which are not tracked)
# and from the cached fits.
#
# Run from the package root, after 02_run_candidates_claude.R:
#   Rscript additional_context/overdispersion_brainstorm/08_report_tables_claude.R
#
# Writes output/report_*.csv.

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
work_dir <- file.path("additional_context", "overdispersion_brainstorm")
cache_dir <- file.path(work_dir, "cache")
out_dir <- file.path(work_dir, "output")

source(file.path(work_dir, "helpers_candidates_claude.R"))

gene_df <- utils::read.csv(file.path(out_dir, "genes.csv"))
summary_gene_df <- utils::read.csv(file.path(out_dir, "gene_summaries.csv"))
gene_df$p <- 10^(-gene_df$log10p)

method_vec <- c("master", "mle", "cap_10s")
regime_vec <- c("null", "baseline", "weak_de", "low_count", "wide_rate",
                "trend", "near_poisson", "generate_null")

.write_rounded <- function(df, file_name, digits = 5){
  num_col_vec <- sapply(df, is.numeric)
  df[num_col_vec] <- lapply(df[num_col_vec], signif, digits = digits)
  utils::write.csv(df, file.path(out_dir, file_name), row.names = FALSE)
}

# Likelihood curves of four example genes --------------------------------------

stem <- "baseline__rep01"
cache <- readRDS(file.path(cache_dir, paste0(stem, ".rds")))
esvd_obj <- cache$esvd_obj
dat <- as.matrix(esvd_obj$dat)
mat_list <- .fitted_matrices(esvd_obj = esvd_obj, cc_var = cache$cc_var)
x <- summary_gene_df[summary_gene_df$regime == "baseline" &
                       summary_gene_df$rep == 1, ]
stopifnot(all(x$gene == colnames(dat)))

unit_free_vec <- x$rate_mle / x$s_median
interior_idx <- which(!x$bool_boundary)
.closest_interior <- function(prob){
  target_val <- stats::quantile(unit_free_vec[interior_idx], probs = prob)
  interior_idx[which.min(abs(unit_free_vec[interior_idx] - target_val))]
}
example_idx_vec <- c(
  "strongly overdispersed" = .closest_interior(0.02),
  "typical" = .closest_interior(0.5),
  "weakly overdispersed" = .closest_interior(0.98),
  "at the boundary" = which(x$bool_boundary)[1]
)

rate_grid_vec <- exp(seq(log(0.2), log(3e7), length.out = 200))
profile_df <- do.call(rbind, lapply(names(example_idx_vec), function(label){
  j <- example_idx_vec[[label]]
  x_vec <- as.numeric(dat[, j])
  mu_vec <- mat_list$mean_mat[, j]
  s_vec <- mat_list$library_mat[, j]
  ll_vec <- sapply(log(rate_grid_vec), function(rho){
    .nb_loglik(rho = rho, x_vec = x_vec, mu_vec = mu_vec, s_vec = s_vec)
  })
  ll_poisson <- .poisson_loglik(x_vec = x_vec, mu_vec = mu_vec, s_vec = s_vec)
  data.frame(example = label,
             gene = x$gene[j],
             rate = rate_grid_vec,
             loglik_drop = ll_vec - max(c(ll_vec, ll_poisson)),
             poisson_drop = ll_poisson - max(c(ll_vec, ll_poisson)),
             rate_mle = x$rate_mle[j],
             legacy_cap = x$s_max[j],
             proposed_cap = 10 * x$s_median[j],
             rate_true_fit_units = x$nuisance_true[j] * x$s_median[j],
             pearson = x$pearson[j])
}))
.write_rounded(profile_df, "report_profile_curves.csv")

# Rates against the truth, three replicates ------------------------------------

rate_df <- gene_df[gene_df$candidate %in% method_vec & gene_df$rep <= 3 &
                     gene_df$regime %in% c("baseline", "low_count",
                                           "wide_rate", "trend",
                                           "near_poisson"),
                   c("candidate", "regime", "rep", "gene", "nuisance")]
rate_df <- merge(rate_df,
                 summary_gene_df[, c("regime", "rep", "gene", "bool_boundary",
                                     "mean_count", "nuisance_true",
                                     "s_median")],
                 by = c("regime", "rep", "gene"))
rate_df$rate_unit_free <- rate_df$nuisance / rate_df$s_median
.write_rounded(rate_df[, c("candidate", "regime", "rep", "gene",
                           "bool_boundary", "mean_count", "nuisance",
                           "nuisance_true", "rate_unit_free", "s_median")],
               "report_rates.csv", digits = 4)

# The Pearson statistic and the boundary ---------------------------------------

pearson_df <- summary_gene_df[summary_gene_df$regime %in%
                                c("baseline", "low_count", "trend",
                                  "near_poisson"),
                              c("regime", "rep", "gene", "bool_boundary",
                                "pearson", "ll_gain", "mean_count")]

# The first-order term of the log-likelihood around the Poisson limit:
#   l(beta) = l_Poisson + score / (2 * beta) + O(beta^-2),
#   score = sum_i [(A_i - m_i)^2 - A_i] / mu_i,   m_i = mu_i * s_i.
# A negative score means the likelihood approaches its limit from below, so
# no finite rate does better than Poisson nearby.
pearson_df$score_at_poisson <- NA
for(stem in unique(paste0(pearson_df$regime, "__rep",
                          sprintf("%02d", pearson_df$rep)))){
  stem_cache <- readRDS(file.path(cache_dir, paste0(stem, ".rds")))
  stem_dat <- as.matrix(stem_cache$esvd_obj$dat)
  stem_mat_list <- .fitted_matrices(esvd_obj = stem_cache$esvd_obj,
                                    cc_var = stem_cache$cc_var)
  fitted_mat <- stem_mat_list$mean_mat * stem_mat_list$library_mat
  score_vec <- colSums(((stem_dat - fitted_mat)^2 - stem_dat) /
                         stem_mat_list$mean_mat)
  idx_vec <- which(pearson_df$regime == sub("__rep.*$", "", stem) &
                     pearson_df$rep == as.integer(sub("^.*__rep", "", stem)))
  stopifnot(all(pearson_df$gene[idx_vec] == colnames(stem_dat)))
  pearson_df$score_at_poisson[idx_vec] <- score_vec
}
stopifnot(!anyNA(pearson_df$score_at_poisson))
.write_rounded(pearson_df, "report_pearson.csv", digits = 4)

# Discoveries per replicate ----------------------------------------------------

discovery_df <- do.call(rbind, lapply(method_vec, function(method){
  do.call(rbind, lapply(regime_vec, function(regime){
    y <- gene_df[gene_df$candidate == method & gene_df$regime == regime, ]
    data.frame(candidate = method,
               regime = regime,
               rep = sort(unique(y$rep)),
               fp = as.numeric(tapply(y$fdr < 0.05 & !y$is_de, y$rep, sum)),
               tp = as.numeric(tapply(y$fdr < 0.05 & y$is_de, y$rep, sum)),
               num_de = as.numeric(tapply(y$is_de, y$rep, sum)))
  }))
}))
.write_rounded(discovery_df, "report_discovery_by_rep.csv")

# Quantiles of the null p-values -----------------------------------------------

prob_vec <- stats::ppoints(150)
qq_df <- do.call(rbind, lapply(method_vec, function(method){
  do.call(rbind, lapply(regime_vec, function(regime){
    p_vec <- gene_df$p[gene_df$candidate == method &
                         gene_df$regime == regime & !gene_df$is_de]
    # Evenly spaced on the -log10 scale of the expected p-value, so that the
    # tail, where FDR control acts, is not thinned away.
    expected_vec <- 10^(-seq(0, log10(length(p_vec)), length.out = 150))
    data.frame(candidate = method,
               regime = regime,
               expected = expected_vec,
               observed = as.numeric(stats::quantile(p_vec,
                                                     probs = expected_vec,
                                                     type = 1)))
  }))
}))
.write_rounded(qq_df, "report_null_qq.csv")

# The geometric mean that every rate is divided by -----------------------------

rescale_df <- do.call(rbind, lapply(method_vec, function(method){
  do.call(rbind, lapply(regime_vec, function(regime){
    y <- gene_df[gene_df$candidate == method & gene_df$regime == regime, ]
    geometric_mean_vec <- tapply(y$nuisance, y$rep, function(v){
      10^mean(log10(v))
    })
    rescaled_vec <- y$nuisance / geometric_mean_vec[as.character(y$rep)]
    data.frame(candidate = method,
               regime = regime,
               geometric_mean = mean(geometric_mean_vec),
               rescaled_q01 = as.numeric(stats::quantile(rescaled_vec, 0.01)),
               rescaled_q50 = stats::median(rescaled_vec),
               rescaled_q99 = as.numeric(stats::quantile(rescaled_vec, 0.99)),
               rescaled_max = max(rescaled_vec))
  }))
}))
.write_rounded(rescale_df, "report_rescaling.csv")

# The trend regime, by third of the expression ---------------------------------

trend_df <- gene_df[gene_df$regime == "trend" &
                      gene_df$candidate %in% c(method_vec, "oracle"), ]
trend_df <- merge(trend_df,
                  summary_gene_df[summary_gene_df$regime == "trend",
                                  c("rep", "gene", "mean_count",
                                    "nuisance_true", "rate_mle", "s_median")],
                  by = c("rep", "gene"))
third_vec <- rep(NA, nrow(trend_df))
for(rep_idx in unique(trend_df$rep)){
  idx_vec <- which(trend_df$rep == rep_idx)
  cut_vec <- stats::quantile(trend_df$mean_count[idx_vec],
                             probs = c(1 / 3, 2 / 3))
  third_vec[idx_vec] <- 1 + (trend_df$mean_count[idx_vec] > cut_vec[1]) +
    (trend_df$mean_count[idx_vec] > cut_vec[2])
}
trend_df$third <- c("low", "middle", "high")[third_vec]
num_reps <- length(unique(trend_df$rep))
key_df <- unique(trend_df[, c("candidate", "third")])
trend_summary_df <- do.call(rbind, lapply(seq_len(nrow(key_df)), function(i){
  y <- trend_df[trend_df$candidate == key_df$candidate[i] &
                  trend_df$third == key_df$third[i], ]
  data.frame(candidate = key_df$candidate[i],
             third = key_df$third[i],
             mean_count = stats::median(y$mean_count),
             rate_true_median = stats::median(y$nuisance_true),
             fp = sum(y$fdr < 0.05 & !y$is_de) / num_reps,
             tp = sum(y$fdr < 0.05 & y$is_de) / num_reps,
             num_de = sum(y$is_de) / num_reps,
             frac_changed = mean(abs(y$nuisance / y$rate_mle - 1) > 1e-6),
             frac_within_2fold = mean(abs(log(y$nuisance / y$s_median) -
                                            log(y$nuisance_true)) < log(2)))
}))
.write_rounded(trend_summary_df, "report_trend_by_third.csv")

# Where the true rates' false discoveries come from ---------------------------

# The oracle's null genes in the two regimes whose true rates vary most, cut
# by the size of the true rate (in units of the gene's library size).
oracle_df <- gene_df[gene_df$candidate == "oracle" &
                       gene_df$regime %in% c("wide_rate", "trend") &
                       !gene_df$is_de, ]
oracle_df <- merge(oracle_df,
                   summary_gene_df[, c("regime", "rep", "gene",
                                       "nuisance_true")],
                   by = c("regime", "rep", "gene"))
oracle_df$rate_bin <- cut(oracle_df$nuisance_true,
                          breaks = c(0, 10, 30, 100, Inf),
                          labels = c("below 10", "10 to 30", "30 to 100",
                                     "above 100"))
key_df <- unique(oracle_df[, c("regime", "rate_bin")])
oracle_summary_df <- do.call(rbind, lapply(seq_len(nrow(key_df)), function(i){
  y <- oracle_df[oracle_df$regime == key_df$regime[i] &
                   oracle_df$rate_bin == key_df$rate_bin[i], ]
  data.frame(regime = key_df$regime[i],
             rate_bin = key_df$rate_bin[i],
             num_null_genes = nrow(y),
             num_false_discoveries = sum(y$fdr < 0.05),
             frac_false_discoveries = mean(y$fdr < 0.05),
             sd_gaussian_teststat = stats::sd(y$gaussian_teststat))
}))
oracle_summary_df <- oracle_summary_df[order(oracle_summary_df$regime,
                                             oracle_summary_df$rate_bin), ]
.write_rounded(oracle_summary_df, "report_oracle_by_rate.csv")

print(paste0("Wrote ", length(list.files(out_dir, pattern = "^report_")),
             " report tables to ", out_dir))
