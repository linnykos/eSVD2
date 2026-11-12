# Checks of the bank `lfcse-bank-2` as built, for notes/BANK.md.
# Drafted by Claude for Kevin Z. Lin, 2026-11-11.
#
# Reads the manifests and, for the oracle rows, the truth of every dev
# example. This script reads the truth: none of its code may enter the
# trial file. Run from the repository root:
#   Rscript autoresearch_lfcse2/diagnostics/bank_checks_claude.R
# Writes autoresearch_lfcse2/output/diagnostics/bank_checks.csv (realized
# values by factor level) and oracle_dev_scores.csv (one row per example).

rm(list = ls())

# Paths ------------------------------------------------------------------------

repo_dir <- getwd()
loop_dir <- file.path(repo_dir, "autoresearch_lfcse2")
source(file.path(loop_dir, "harness", "generator_claude.R"))
.load_esvd2(repo_dir)
bank_dir <- .bank_dir()
output_dir <- file.path(loop_dir, "output", "diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
num_cores <- as.integer(Sys.getenv("ESVD2_CORES", unset = "8"))

# Realized values by factor level ----------------------------------------------

manifest_df <- do.call(rbind, lapply(c("dev", "test"), function(split){
  utils::read.csv(file.path(bank_dir, paste0("manifest_", split, ".csv")),
                  stringsAsFactors = FALSE)
}))
manifest_df$frac_cap_or_boundary <- manifest_df$frac_capped +
  manifest_df$frac_boundary
value_vec <- c("frac_cap_or_boundary", "frac_capped", "mean_count",
               "num_cells", "num_exact_null", "num_warnings", "fit_seconds")
check_list <- list()
for(factor_name in .factor_names()){
  for(level in sort(unique(manifest_df[[factor_name]]))){
    bool_vec <- manifest_df[[factor_name]] == level
    check_list[[length(check_list) + 1]] <- data.frame(
      factor = factor_name, level = as.character(level),
      num_examples = sum(bool_vec),
      t(colMeans(manifest_df[bool_vec, value_vec])),
      stringsAsFactors = FALSE
    )
  }
}
check_df <- do.call(rbind, check_list)
utils::write.csv(check_df, file.path(output_dir, "bank_checks.csv"),
                 row.names = FALSE)
print(check_df, digits = 3)
print(table(manifest_df$split, manifest_df$status))
print(paste0("Fit minutes summed over examples: ",
             round(sum(manifest_df$fit_seconds) / 60)))

# Oracle rows ------------------------------------------------------------------

# `donor_wls` (the first loop's rule) with a donor x gene offset subtracted
# from the log donor means before the regression. With the offset equal to
# the donor's true state effect (and block shift) it is the rule an ideal
# use of the embedding could reach; with the donor's whole true effect it
# is the floor that sampling noise alone leaves.
rule_env <- new.env()
sys.source(file.path(loop_dir, "candidates", "donor_wls_claude.R"),
           envir = rule_env)

.wls_with_offset <- function(view, offset_mat){
  individual_chr <- as.character(view$individual_vec)
  count_mat <- rowsum(view$dat, group = individual_chr)
  library_sum_mat <- rowsum(view$library_mat, group = individual_chr)
  donor_vec <- rownames(count_mat)
  log_mat <- log((count_mat + 0.5) / library_sum_mat) -
    offset_mat[donor_vec, , drop = FALSE]
  arm_vec <- rule_env$.arm_of_donor(case_control_vec = view$case_control_vec,
                                    individual_vec = view$individual_vec,
                                    donor_vec = donor_vec)
  design_mat <- rule_env$.donor_design(covariates = view$covariates,
                                       individual_vec = view$individual_vec,
                                       arm_vec = arm_vec,
                                       donor_vec = donor_vec)
  .least_squares <- function(response_mat){
    design_mat %*% solve(crossprod(design_mat),
                         crossprod(design_mat, response_mat))
  }
  fit_mat <- .least_squares(log_mat)
  log_mat <- log_mat - apply(log_mat - fit_mat, 1, stats::median)
  fit_mat <- .least_squares(log_mat)
  residual_df <- nrow(design_mat) - ncol(design_mat)
  inverse_count_vec <- as.numeric(exp(-(fit_mat + log(library_sum_mat))))
  squared_vec <- as.numeric((log_mat - fit_mat)^2) * nrow(design_mat) /
    residual_df
  slope <- stats::cov(squared_vec, inverse_count_vec) /
    stats::var(inverse_count_vec)
  variance_coef <- pmax(c(mean(squared_vec) -
                            slope * mean(inverse_count_vec), slope), 0)
  # With the whole donor effect removed both parts can be 0.
  if(all(variance_coef == 0)) variance_coef[1] <- mean(squared_vec)
  weight_mat <- matrix(1 / (variance_coef[1] +
                              variance_coef[2] * inverse_count_vec),
                       nrow = nrow(fit_mat))
  num_genes <- ncol(log_mat)
  lfc_vec <- numeric(num_genes)
  arm_var_vec <- numeric(num_genes)
  residual_var <- numeric(num_genes)
  for(j in seq_len(num_genes)){
    weighted_design <- design_mat * weight_mat[, j]
    gram_inv_mat <- solve(crossprod(design_mat, weighted_design))
    coef_vec <- gram_inv_mat %*% crossprod(weighted_design, log_mat[, j])
    lfc_vec[j] <- coef_vec[2] / log(2)
    arm_var_vec[j] <- gram_inv_mat[2, 2]
    residual_var[j] <- sum(weight_mat[, j] *
                             (log_mat[, j] - design_mat %*% coef_vec)^2) /
      residual_df
  }
  se_vec <- sqrt(arm_var_vec *
                   rule_env$.moderate_variance(var_vec = residual_var,
                                               df = residual_df,
                                               covariate_vec =
                                                 colMeans(fit_mat))) / log(2)
  z_vec <- lfc_vec / se_vec
  null_res <- suppressWarnings(eSVD2::multtest(z_vec[is.finite(z_vec)]))
  list(lfc_vec = lfc_vec - null_res$null_mean * stats::median(se_vec),
       se_vec = se_vec)
}

.nll <- function(res, truth_vec){
  mean(0.5 * log(2 * pi * res$se_vec^2) +
         (res$lfc_vec - truth_vec)^2 / (2 * res$se_vec^2))
}

dev_df <- manifest_df[manifest_df$split == "dev", ]
row_list <- parallel::mclapply(seq_len(nrow(dev_df)), function(i){
  view <- .build_view(readRDS(file.path(bank_dir, dev_df$file_public[i])))
  truth <- readRDS(file.path(bank_dir, dev_df$file_truth[i]))
  truth_vec <- as.numeric(truth$truth_log2fc_vec)
  state_mat <- tcrossprod(truth$donor_state_mat, truth$population$y_mat) +
    outer(truth$donor_df$block_shift,
          as.numeric(truth$population$bool_block_vec))
  dimnames(state_mat) <- dimnames(truth$donor_effect_mat)
  zero_mat <- 0 * truth$donor_effect_mat
  res_list <- list(
    wls = tryCatch(.wls_with_offset(view, zero_mat), error = function(e){NULL}),
    oracle_state = tryCatch(.wls_with_offset(view, state_mat),
                            error = function(e){NULL}),
    oracle_all = tryCatch(.wls_with_offset(view, truth$donor_effect_mat),
                          error = function(e){NULL})
  )
  score_vec <- sapply(res_list, function(res){
    if(is.null(res)) NA_real_ else .nll(res, truth_vec)
  })
  cbind(dev_df[i, c("id", .factor_names())], t(score_vec))
}, mc.cores = num_cores, mc.preschedule = FALSE)
oracle_df <- do.call(rbind, row_list)
utils::write.csv(oracle_df, file.path(output_dir, "oracle_dev_scores.csv"),
                 row.names = FALSE)
print(colMeans(oracle_df[, c("wls", "oracle_state", "oracle_all")],
               na.rm = TRUE))
print(stats::aggregate(cbind(wls, oracle_state, oracle_all) ~ state,
                       data = oracle_df, FUN = mean))
print("Done! :)")
