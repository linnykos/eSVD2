# Green split-half study, step 5: which step of `donor_wls` loses the
# agreement between the halves
# Drafted by Claude for Kevin Z. Lin, 2026-11-16

library(Seurat)

rm(list = ls())

# On the Green halves the logFC of `donor_wls@total` does not correlate
# between the halves for well-expressed genes, while `raw@total` (the same
# library) and DESeq2 do. `donor_wls` reads the counts and the library only
# through their sums per donor, so it can be rerun from the donor
# pseudobulk in seconds. This script reruns it with one step changed at a
# time (the donor offset, the precision weights, the donor-level
# covariates, each covariate alone, and the log donor means winsorized gene
# by gene) and records the correlation between the halves for each version. The summed-count library is used throughout,
# so no fit is needed.
#
# Run from the repository root, with the cohort to read as the argument:
#   Rscript additional_context/2026-11-16_green/05_diagnose-donor-wls_claude.R green
#   Rscript additional_context/2026-11-16_green/05_diagnose-donor-wls_claude.R gabitto
#
# Writes csv/donor_wls_diagnosis_<cohort>.csv and, for the share of donors
# that are outliers on the log scale, csv/outlying_donors_<cohort>.csv,
# beside this script.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context", "2026-11-16_green")
stopifnot(dir.exists(study_dir))
cohort <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(cohort %in% c("green", "gabitto"))
cohort_dir <- if(cohort == "green") study_dir else
  file.path(repo_dir, "additional_context", "2026-11-15_gabitto")
source(file.path(cohort_dir, "paths_claude.R"))
csv_dir <- file.path(study_dir, "csv")
rule_env <- new.env(parent = globalenv())
sys.source(file.path(repo_dir, "autoresearch_lfcse", "notes",
                     "RULE_nov11_claude.R"), envir = rule_env)

# The donor pseudobulk ---------------------------------------------------------

load(data_file)
load(file.path(out_dir, paste0(cohort, "_splits.RData")))
count_mat <- SeuratObject::LayerData(seurat_obj, assay = "RNA",
                                     layer = "counts")
metadata_df <- seurat_obj@meta.data
rm(seurat_obj)
donor_chr <- as.character(metadata_df[[id_var]])
# Genes x donors, then donors x genes.
indicator_mat <- Matrix::sparse.model.matrix(~ 0 + factor(donor_chr))
colnames(indicator_mat) <- levels(factor(donor_chr))
pseudobulk_mat <- as.matrix(Matrix::t(count_mat %*% indicator_mat))
cell_total_vec <- Matrix::colSums(count_mat)
# The `total` library of step 3: each cell's summed count over the mean
# summed count of the half's cells, summed within the donor. The mean is a
# constant per half and only moves the intercept.
donor_total_vec <- tapply(cell_total_vec, donor_chr, sum)[rownames(pseudobulk_mat)]
num_cells_vec <- table(donor_chr)[rownames(pseudobulk_mat)]
zero_fraction_vec <- 1 - Matrix::rowSums(count_mat > 0) / ncol(count_mat)
rm(count_mat, indicator_mat)

donor_meta_df <- unique(metadata_df[, c(id_var, case_control_var,
                                        categorical_vars, numerical_vars)])
donor_meta_df <- donor_meta_df[match(rownames(pseudobulk_mat),
                                     as.character(donor_meta_df[[id_var]])), ]
arm_vec <- as.numeric(as.character(donor_meta_df[[case_control_var]]) ==
                        case_control_levels[2])
covariate_mat <- cbind(
  sapply(categorical_vars, function(name){
    value_vec <- as.character(donor_meta_df[[name]])
    as.numeric(value_vec == sort(unique(value_vec))[2])
  }),
  sapply(numerical_vars, function(name){
    as.numeric(scale(donor_meta_df[[name]]))
  })
)
colnames(covariate_mat) <- c(categorical_vars, numerical_vars)

# Outlying donors --------------------------------------------------------------

# How often a donor's log mean of a gene lies far from the other donors',
# over the whole cohort and the genes with under 80% zeros: each gene's
# log donor means (summed counts over the donor's summed total) are
# standardized by the gene's median and MAD. Under a Gaussian 0.13% of the
# values are more than 3 MADs below the median.
dense_idx <- which(zero_fraction_vec < 0.8)
log_donor_mat <- log((pseudobulk_mat[, dense_idx] + 0.5) /
                       as.numeric(donor_total_vec / mean(donor_total_vec)))
median_vec <- apply(log_donor_mat, 2, stats::median)
deviation_mat <- t(t(log_donor_mat) - median_vec)
z_mat <- t(t(deviation_mat) / apply(log_donor_mat, 2, stats::mad))
bool_low_mat <- z_mat < -3
donor_share_vec <- rowMeans(bool_low_mat)
outlier_df <- data.frame(
  cohort = cohort,
  num_donors = nrow(log_donor_mat),
  num_genes = ncol(log_donor_mat),
  share_below_3 = mean(bool_low_mat),
  share_below_5 = mean(z_mat < -5),
  share_above_3 = mean(z_mat > 3),
  mean_deficit = mean(deviation_mat[bool_low_mat]),
  median_deficit = stats::median(deviation_mat[bool_low_mat]),
  # Donors with at least a tenth of their genes below the threshold, and
  # the largest share of any donor.
  share_donors_outlying = mean(donor_share_vec >= 0.1),
  max_donor_share = max(donor_share_vec),
  cor_with_log_depth = stats::cor(donor_share_vec,
                                  log(as.numeric(donor_total_vec) /
                                        as.numeric(num_cells_vec))),
  stringsAsFactors = FALSE
)
print(outlier_df)
utils::write.csv(outlier_df, file = file.path(csv_dir,
                                              paste0("outlying_donors_",
                                                     cohort, ".csv")),
                 row.names = FALSE)
rm(log_donor_mat, deviation_mat, z_mat, bool_low_mat)

# `donor_wls` with switches ----------------------------------------------------

# The rule's own code (autoresearch_lfcse/notes/RULE_nov11_claude.R) from
# the donor sums on, with three switches. `bool_winsorize` pulls each
# gene's log donor means to within 3 MADs of the gene's median before
# anything is fitted, which bounds what one outlying donor can contribute
# to a mean of logarithms. With every switch at its default
# and all covariates the result equals the rule's, which is checked below.
.donor_wls_version <- function(count_mat, library_sum_mat, arm_vec,
                               covariate_mat,
                               bool_offset = TRUE,
                               bool_weights = TRUE,
                               bool_winsorize = FALSE){
  log_mat <- log((count_mat + 0.5) / library_sum_mat)
  if(bool_winsorize){
    log_mat <- apply(log_mat, 2, function(value_vec){
      centre <- stats::median(value_vec)
      spread <- stats::mad(value_vec)
      pmin(pmax(value_vec, centre - 3 * spread), centre + 3 * spread)
    })
  }
  design_mat <- cbind(1, arm_vec, covariate_mat)
  .fit <- function(response_mat){
    design_mat %*% solve(crossprod(design_mat),
                         crossprod(design_mat, response_mat))
  }
  fit_mat <- .fit(log_mat)
  if(bool_offset){
    offset_vec <- apply(log_mat - fit_mat, 1, stats::median)
    log_mat <- log_mat - offset_vec
    fit_mat <- .fit(log_mat)
  }
  residual_df <- nrow(design_mat) - ncol(design_mat)
  weight_mat <- matrix(1, nrow = nrow(fit_mat), ncol = ncol(fit_mat))
  if(bool_weights){
    inverse_count_vec <- as.numeric(exp(-(fit_mat + log(library_sum_mat))))
    squared_vec <- as.numeric((log_mat - fit_mat)^2) *
      nrow(design_mat) / residual_df
    slope <- stats::cov(squared_vec, inverse_count_vec) /
      stats::var(inverse_count_vec)
    variance_coef <- pmax(c(mean(squared_vec) -
                              slope * mean(inverse_count_vec), slope), 0)
    weight_mat <- matrix(1 / (variance_coef[1] +
                                variance_coef[2] * inverse_count_vec),
                         nrow = nrow(fit_mat))
  }
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
                                               covariate_vec = colMeans(fit_mat))) /
    log(2)
  list(lfc_vec = lfc_vec, se_vec = se_vec)
}

# `raw` from the donor sums: the log ratio of the arms' arithmetic means of
# the donors' ratios of sums. (The rule itself averages the cells' ratios
# within a donor; with the summed-count library the two are close, and the
# stored `raw@total` is reported beside this one.)
.raw_version <- function(count_mat, library_sum_mat, arm_vec){
  mean_mat <- count_mat / library_sum_mat
  log2(colMeans(mean_mat[arm_vec == 1, , drop = FALSE]) /
         colMeans(mean_mat[arm_vec == 0, , drop = FALSE]))
}

version_list <- list(
  "as the rule" = list(),
  "no donor offset" = list(bool_offset = FALSE),
  "no weights" = list(bool_weights = FALSE),
  "no covariates" = list(covariates = character(0)),
  "no covariates, no offset" = list(covariates = character(0),
                                    bool_offset = FALSE),
  "no covariates, no offset, no weights" = list(covariates = character(0),
                                                bool_offset = FALSE,
                                                bool_weights = FALSE),
  "as the rule, winsorized" = list(bool_winsorize = TRUE),
  "no covariates, winsorized" = list(covariates = character(0),
                                     bool_winsorize = TRUE)
)
for(name in colnames(covariate_mat)){
  version_list[[paste0("without ", name)]] <- list(
    covariates = setdiff(colnames(covariate_mat), name)
  )
  version_list[[paste0("only ", name)]] <- list(covariates = name)
}

# The halves -------------------------------------------------------------------

.run_half <- function(split, half){
  idx <- which(rownames(pseudobulk_mat) %in%
                 donor_df[[id_var]][donor_df[[split]] == half])
  sub_count_mat <- pseudobulk_mat[idx, , drop = FALSE]
  # The same library for every gene of a donor.
  library_sum_mat <- matrix(as.numeric(donor_total_vec[idx]) /
                              mean(donor_total_vec[idx] / num_cells_vec[idx]),
                            nrow = length(idx), ncol = ncol(sub_count_mat))
  res_list <- lapply(version_list, function(version){
    name_vec <- if(is.null(version$covariates)) colnames(covariate_mat) else
      version$covariates
    .donor_wls_version(count_mat = sub_count_mat,
                       library_sum_mat = library_sum_mat,
                       arm_vec = arm_vec[idx],
                       covariate_mat = covariate_mat[idx, name_vec,
                                                     drop = FALSE],
                       bool_offset = !identical(version$bool_offset, FALSE),
                       bool_weights = !identical(version$bool_weights, FALSE),
                       bool_winsorize = isTRUE(version$bool_winsorize))
  })
  res_list[["raw (donor sums)"]] <- list(
    lfc_vec = .raw_version(sub_count_mat, library_sum_mat, arm_vec[idx]),
    se_vec = rep(NA_real_, ncol(sub_count_mat))
  )
  res_list
}

# The stored `donor_wls@total` of step 3, to check the reconstruction.
.stored <- function(split, half, method){
  rules_env <- new.env()
  load(file.path(out_dir, paste0(cohort, "_", split, "_half", half,
                                 "_rules.RData")), envir = rules_env)
  list(lfc_vec = rules_env$lfc_mat[colnames(pseudobulk_mat), method],
       status_vec = rules_env$gene_df$status[match(colnames(pseudobulk_mat),
                                                   rules_env$gene_df$gene)])
}

.cor_ok <- function(x_vec, y_vec, idx){
  idx <- idx[is.finite(x_vec[idx]) & is.finite(y_vec[idx])]
  stats::cor(x_vec[idx], y_vec[idx])
}

out_df <- do.call(rbind, lapply(c("split1", "split2"), function(split){
  res_1 <- .run_half(split, 1)
  res_2 <- .run_half(split, 2)
  stored_1 <- .stored(split, 1, "donor_wls@total")
  stored_2 <- .stored(split, 2, "donor_wls@total")
  # The rule's centring shifts every gene of a half by one constant, so the
  # reconstruction must equal the stored rule up to a constant.
  gap_vec <- res_1[["as the rule"]]$lfc_vec - stored_1$lfc_vec
  print(paste0(split, ": reconstruction against the stored donor_wls@total, ",
               "SD of the gap ", signif(stats::sd(gap_vec), 3)))
  stopifnot(stats::sd(gap_vec) < 1e-6)
  bool_capped <- stored_1$status_vec != "estimated" |
    stored_2$status_vec != "estimated"
  set_list <- list(all = seq_along(zero_fraction_vec),
                   zeros_below_80 = which(zero_fraction_vec < 0.8),
                   zeros_80_to_95 = which(zero_fraction_vec >= 0.8 &
                                            zero_fraction_vec < 0.95),
                   zeros_above_95 = which(zero_fraction_vec >= 0.95),
                   rate_estimated = which(!bool_capped),
                   rate_capped = which(bool_capped))
  res_1[["raw@total (stored)"]] <- .stored(split, 1, "raw@total")
  res_2[["raw@total (stored)"]] <- .stored(split, 2, "raw@total")
  res_1[["deseq2 (stored)"]] <- .stored(split, 1, "deseq2")
  res_2[["deseq2 (stored)"]] <- .stored(split, 2, "deseq2")
  do.call(rbind, lapply(names(res_1), function(version){
    do.call(rbind, lapply(names(set_list), function(gene_set){
      idx <- set_list[[gene_set]]
      data.frame(cohort = cohort, split = split, version = version,
                 gene_set = gene_set, num_genes = length(idx),
                 pearson = .cor_ok(res_1[[version]]$lfc_vec,
                                   res_2[[version]]$lfc_vec, idx),
                 sd_lfc = stats::sd(c(res_1[[version]]$lfc_vec[idx],
                                      res_2[[version]]$lfc_vec[idx]),
                                    na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }))
}))
out_df$pearson <- signif(out_df$pearson, 4)
out_df$sd_lfc <- signif(out_df$sd_lfc, 4)
utils::write.csv(out_df, file = file.path(csv_dir,
                                          paste0("donor_wls_diagnosis_",
                                                 cohort, ".csv")),
                 row.names = FALSE)
wide_df <- stats::reshape(out_df[, c("split", "version", "gene_set",
                                     "pearson")],
                          idvar = c("split", "version"), timevar = "gene_set",
                          direction = "wide")
print(wide_df, digits = 2, row.names = FALSE)
print("Done! :)")
