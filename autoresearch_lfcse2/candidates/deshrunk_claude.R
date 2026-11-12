# Reference candidate `deshrunk`: the stored logFC with the between-donor
# SE of the posterior donor means divided by the gene's mean shrinkage weight
# `library / (library + rate)`.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11. A drop-in body for
# `report_lfc_se()` (same contract as ../report_rule_claude.R), scored once
# at the build for the reference rows of notes/BANK.md.

.donor_means <- function(mat, individual_vec){
  individual_chr <- as.character(individual_vec)
  num_cells_vec <- table(individual_chr)
  sum_mat <- rowsum(mat, group = individual_chr)
  sum_mat / as.numeric(num_cells_vec[rownames(sum_mat)])
}

.arm_of_donor <- function(case_control_vec, individual_vec, donor_vec){
  case_control_vec[match(donor_vec, as.character(individual_vec))]
}

.lfc_from_donor_means <- function(mean_mat, arm_vec){
  case_mat <- mean_mat[arm_vec == 1, , drop = FALSE]
  control_mat <- mean_mat[arm_vec == 0, , drop = FALSE]
  n1 <- nrow(case_mat)
  n0 <- nrow(control_mat)
  case_mean <- colMeans(case_mat)
  control_mean <- colMeans(control_mat)
  case_var <- apply(case_mat, 2, stats::var)
  control_var <- apply(control_mat, 2, stats::var)
  list(lfc_vec = log2(case_mean / control_mean),
       se_vec = sqrt(case_var / (n1 * case_mean^2) +
                       control_var / (n0 * control_mean^2)) / log(2))
}

report_lfc_se <- function(view){
  mean_mat <- .donor_means(mat = view$posterior_mean_mat,
                           individual_vec = view$individual_vec)
  arm_vec <- .arm_of_donor(case_control_vec = view$case_control_vec,
                           individual_vec = view$individual_vec,
                           donor_vec = rownames(mean_mat))
  between <- .lfc_from_donor_means(mean_mat = mean_mat, arm_vec = arm_vec)
  denominator_mat <- sweep(view$library_mat, MARGIN = 2,
                           STATS = view$nuisance_vec, FUN = "+")
  weight_mat <- .donor_means(mat = view$library_mat / denominator_mat,
                             individual_vec = view$individual_vec)
  data.frame(gene = colnames(view$dat),
             log2fc = view$stored_df$log2fc,
             se = between$se_vec / colMeans(weight_mat),
             stringsAsFactors = FALSE)
}
