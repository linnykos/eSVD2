rm(list = ls())

# What is a realistic donor SD? A check on the fitted ASD cohort
# Drafted by Claude for Kevin Z. Lin, 2026-11-11
#
# Question 6 of the alignment memo. The bank's generator gives every donor
# a random effect on the log scale with SD `donor_sd`; the exploration used
# 0.1 throughout. This script estimates that SD gene by gene on the
# Velmeshev L2/3 cohort from the fitted eSVD object, and the bank's two
# levels are set at the lower and upper quartiles across genes.
#
# For a gene, let m_d be donor d's mean of A / library over its cells and
# v_d the sampling variance of that mean (the within-donor variance of
# A / library over the donor's cell count). Within an arm, the variance of
# m_d across donors is the between-donor variance plus the mean of v_d, so
#   CV^2 = (var_d(m_d) - mean_d(v_d)) / mean_d(m_d)^2
# pooled over the two arms, and a log-normal donor effect with that CV has
# log-scale SD sqrt(log(1 + CV^2)). Working on the raw scale keeps donors
# whose mean is zero. A negative CV^2 is reported as SD 0.
#
# Run from the repository root:
#   Rscript autoresearch_lfcse/diagnostics/donor_sd_asd_claude.R <asd_esvd.RData>
# Writes autoresearch_lfcse/output/diagnostics/donor_sd_asd.csv. About five
# minutes, nearly all of it loading the 1.4 GB object.

# Paths ------------------------------------------------------------------------

arg_vec <- commandArgs(trailingOnly = TRUE)
if(length(arg_vec) != 1) stop("usage: donor_sd_asd_claude.R <asd_esvd.RData>")
asd_file <- arg_vec[1]
repo_dir <- getwd()
output_dir <- file.path(repo_dir, "autoresearch_lfcse", "output",
                        "diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Loading ----------------------------------------------------------------------

load(asd_file)
dat <- as.matrix(eSVD_obj$dat)
covariates <- eSVD_obj$covariates
fit <- eSVD_obj[[eSVD_obj[["latest_Fit"]]]]
cc_var <- eSVD_obj$param$init_case_control_variable
individual_vec <- as.character(eSVD_obj$individual)
cc_vec <- covariates[, cc_var]
print(paste0("Cells: ", nrow(dat), ", genes: ", ncol(dat), ", donors: ",
             length(unique(individual_vec)), ", case-control column: ",
             cc_var))

# The library as 1.2.0's default forms it: every covariate but the
# case-control one, the intercept included.
library_idx <- which(colnames(covariates) != cc_var)
library_mat <- exp(tcrossprod(covariates[, library_idx, drop = FALSE],
                              fit$z_mat[, library_idx, drop = FALSE]))
library_mat <- pmax(library_mat, 0.1)
raw_mat <- dat / library_mat
rm(library_mat)

# Per-donor moments ------------------------------------------------------------

donor_vec <- unique(individual_vec)
num_cells_vec <- as.numeric(table(individual_vec)[donor_vec])
sum_mat <- rowsum(raw_mat, group = individual_vec)[donor_vec, , drop = FALSE]
sum_sq_mat <- rowsum(raw_mat^2, group = individual_vec)[donor_vec, ,
                                                        drop = FALSE]
mean_mat <- sum_mat / num_cells_vec
within_var_mat <- (sum_sq_mat - num_cells_vec * mean_mat^2) /
  (num_cells_vec - 1)
sampling_var_mat <- within_var_mat / num_cells_vec
arm_vec <- cc_vec[match(donor_vec, individual_vec)]
print(paste0("Cells per donor: median ", stats::median(num_cells_vec),
             ", range ", min(num_cells_vec), " to ", max(num_cells_vec)))

# The pooled estimate ----------------------------------------------------------

between_ss_vec <- rep(0, ncol(dat))
sampling_sum_vec <- rep(0, ncol(dat))
df_val <- 0
for(arm in sort(unique(arm_vec))){
  idx_vec <- which(arm_vec == arm)
  arm_mean_mat <- mean_mat[idx_vec, , drop = FALSE]
  centre_vec <- colMeans(arm_mean_mat)
  between_ss_vec <- between_ss_vec +
    colSums(sweep(arm_mean_mat, MARGIN = 2, STATS = centre_vec)^2)
  # Each arm's sampling term is weighted by its degrees of freedom, as its
  # between-donor sum of squares is.
  sampling_sum_vec <- sampling_sum_vec + (length(idx_vec) - 1) *
    colMeans(sampling_var_mat[idx_vec, , drop = FALSE])
  df_val <- df_val + length(idx_vec) - 1
}
total_var_vec <- between_ss_vec / df_val
sampling_var_vec <- sampling_sum_vec / df_val
grand_mean_vec <- colMeans(mean_mat)
cv_sq_vec <- (total_var_vec - sampling_var_vec) / grand_mean_vec^2
donor_sd_vec <- sqrt(log1p(pmax(cv_sq_vec, 0)))

gene_df <- data.frame(gene = colnames(dat),
                      mean_count = colMeans(dat),
                      cv_sq = cv_sq_vec,
                      donor_sd = donor_sd_vec,
                      sampling_share = sampling_var_vec / total_var_vec,
                      stringsAsFactors = FALSE)
utils::write.csv(gene_df, file = file.path(output_dir, "donor_sd_asd.csv"),
                 row.names = FALSE)

# Reporting --------------------------------------------------------------------

prob_vec <- c(0.1, 0.25, 0.5, 0.75, 0.9)
print("Donor SD (log scale) across genes, quantiles:")
print(round(stats::quantile(donor_sd_vec, prob_vec, na.rm = TRUE), 3))
print(paste0("Share of genes with a non-positive estimate: ",
             round(mean(cv_sq_vec <= 0, na.rm = TRUE), 3)))
tercile_vec <- cut(gene_df$mean_count,
                   breaks = stats::quantile(gene_df$mean_count,
                                            c(0, 1 / 3, 2 / 3, 1)),
                   labels = c("low", "middle", "high"),
                   include.lowest = TRUE)
print("By tercile of mean count:")
print(round(do.call(rbind, tapply(donor_sd_vec, tercile_vec,
                                  stats::quantile, probs = prob_vec,
                                  na.rm = TRUE)), 3))
print("Sampling share of the across-donor variance, quantiles:")
print(round(stats::quantile(gene_df$sampling_share, prob_vec, na.rm = TRUE),
            3))
print("Done! :)")
