# Tables for OVERDISPERSION_BRAINSTORM.md
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Reads the output of 02_run_candidates_claude.R and writes the summary tables
# to output/summary_*.csv.
#
# Run from the package root:
#   Rscript additional_context/overdispersion_brainstorm/03_summarize_claude.R

rm(list = ls())

work_dir <- file.path("additional_context", "overdispersion_brainstorm")
out_dir <- file.path(work_dir, "output")

gene_df <- utils::read.csv(file.path(out_dir, "genes.csv"))
summary_df <- utils::read.csv(file.path(out_dir, "gene_summaries.csv"))
gene_df$p <- 10^(-gene_df$log10p)

regime_vec <- c("baseline", "null", "strong_de", "near_poisson",
                "minimal_design", "generate_null", "weak_de", "low_count",
                "wide_rate", "trend")
candidate_vec <- unique(gene_df$candidate)
num_reps <- length(unique(gene_df$rep))

# Discoveries ------------------------------------------------------------------

# Area under the ROC curve of |Gaussian statistic| for DE against null, which
# measures the ranking alone, whatever the empirical null did to the scale.
.compute_auc <- function(score_vec, bool_positive_vec){
  rank_vec <- rank(score_vec)
  num_positive <- sum(bool_positive_vec)
  num_negative <- sum(!bool_positive_vec)
  (sum(rank_vec[bool_positive_vec]) - num_positive * (num_positive + 1) / 2) /
    (num_positive * num_negative)
}

discovery_list <- list()
for(candidate in candidate_vec){
  for(regime in regime_vec){
    x <- gene_df[gene_df$candidate == candidate & gene_df$regime == regime, ]
    if(nrow(x) == 0) next()
    bool_sig_vec <- x$fdr < 0.05
    auc_val <- NA
    if(any(x$is_de)){
      auc_val <- mean(sapply(split(x, x$rep), function(y){
        .compute_auc(score_vec = abs(y$gaussian_teststat),
                     bool_positive_vec = y$is_de)
      }))
    }
    discovery_list[[paste(candidate, regime)]] <- data.frame(
      candidate = candidate,
      regime = regime,
      type1 = mean(x$p[!x$is_de] < 0.05),
      fp = sum(bool_sig_vec & !x$is_de) / num_reps,
      tp = sum(bool_sig_vec & x$is_de) / num_reps,
      fdp = sum(bool_sig_vec & !x$is_de) / max(1, sum(bool_sig_vec)),
      auc = auc_val,
      max_abs_teststat = max(abs(x$teststat))
    )
  }
}
discovery_df <- do.call(rbind, discovery_list)
utils::write.csv(discovery_df, file.path(out_dir, "summary_discovery.csv"),
                 row.names = FALSE)

.wide <- function(column, digits){
  res <- stats::reshape(discovery_df[, c("candidate", "regime", column)],
                        idvar = "candidate", timevar = "regime",
                        direction = "wide")
  colnames(res) <- sub(paste0("^", column, "\\."), "", colnames(res))
  res[, -1] <- round(res[, -1], digits)
  rownames(res) <- NULL
  res
}
print("False discoveries per data set at FDR < 0.05")
options(width = 200)
print(.wide("fp", 1))
print("True discoveries per data set at FDR < 0.05")
print(.wide("tp", 1))
print("Type-I error at p < 0.05, null genes")
print(.wide("type1", 3))
print("Pooled false discovery proportion")
print(.wide("fdp", 3))
print("AUC of |z|, DE against null")
print(.wide("auc", 3))
print("Largest |Welch statistic|")
print(.wide("max_abs_teststat", 1))

# The boundary -----------------------------------------------------------------

boundary_df <- do.call(rbind, lapply(regime_vec, function(regime){
  x <- summary_df[summary_df$regime == regime, ]
  interior_vec <- !x$bool_boundary
  unit_free_vec <- x$rate_mle / x$s_median
  data.frame(regime = regime,
             frac_boundary = mean(x$bool_boundary),
             frac_above_1e4 = mean(x$rate_mle > 1e4),
             pearson_boundary = mean(x$pearson[!interior_vec]),
             pearson_interior = mean(x$pearson[interior_vec]),
             pearson_all = mean(x$pearson),
             ratio_raw = stats::median((x$rate_mle / x$nuisance_true)[interior_vec]),
             ratio_unit_free = stats::median((unit_free_vec /
                                                x$nuisance_true)[interior_vec]),
             spearman_raw = stats::cor(x$rate_mle, x$nuisance_true,
                                       method = "spearman"),
             spearman_unit_free = stats::cor(unit_free_vec, x$nuisance_true,
                                             method = "spearman"),
             frac_at_legacy_cap = mean(x$rate_mle > x$s_max),
             map_prior_var = mean(x$map_prior_var))
}))
utils::write.csv(boundary_df, file.path(out_dir, "summary_boundary.csv"),
                 row.names = FALSE)
print("The boundary, and the units of the rate")
print(format(boundary_df, digits = 3))

# How the null scale of the statistic depends on the rate ---------------------

# Within each data set, null genes are ranked by the rate that entered the
# posterior and cut into five groups; the SD of the Gaussian statistic per
# group says whether genes with larger rates have wider null statistics.
scale_list <- list()
for(candidate in c("mle", "master", "cap_10s", "cap_50s", "eb_map",
                   "common_unit", "oracle")){
  for(regime in c("null", "baseline", "near_poisson", "wide_rate", "trend")){
    x <- gene_df[gene_df$candidate == candidate & gene_df$regime == regime &
                   !gene_df$is_de, ]
    if(nrow(x) == 0) next()
    group_vec <- unlist(lapply(split(x$nuisance, x$rep), function(v){
      ceiling(5 * rank(v, ties.method = "first") / length(v))
    }))
    x <- x[order(x$rep), ]
    sd_vec <- tapply(x$gaussian_teststat, group_vec, stats::sd)
    scale_list[[paste(candidate, regime)]] <- data.frame(
      candidate = candidate,
      regime = regime,
      t(round(as.numeric(sd_vec), 2))
    )
  }
}
scale_df <- do.call(rbind, scale_list)
colnames(scale_df)[3:7] <- paste0("rate_group_", 1:5)
rownames(scale_df) <- NULL
utils::write.csv(scale_df, file.path(out_dir, "summary_null_scale.csv"),
                 row.names = FALSE)
print("SD of the null Gaussian statistic, by fifth of the rate (1 = smallest)")
print(scale_df)
