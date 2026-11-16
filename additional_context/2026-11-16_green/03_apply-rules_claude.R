# Green split-half study, step 3: every reporting rule on one fitted half
# Drafted by Claude for Kevin Z. Lin, 2026-11-16

library(Seurat)
library(DESeq2)
library(edgeR)

rm(list = ls())

# Applies every logFC/SE reporting rule of the two autoresearch loops, the
# library-size variants of the rules that read a library, and pseudobulk
# DESeq2 to one fitted half, and stores each rule's per-gene log2 fold
# change and SE. The rule files are the ones the loops left, unedited, and
# the rule list is that of the beyond-score study
# (../2026-11-12_beyond-score/02_apply-rules_claude.R). Nothing is compared
# here; step 4 reads these files.
#
# Run from the repository root:
#   Rscript additional_context/2026-11-16_green/03_apply-rules_claude.R \
#     split1_half1
#
# Writes <out_dir>/green_<label>_rules.RData.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context", "2026-11-16_green")
stopifnot(dir.exists(study_dir))
source(file.path(study_dir, "paths_claude.R"))
loop1_dir <- file.path(repo_dir, "autoresearch_lfcse")
loop2_dir <- file.path(repo_dir, "autoresearch_lfcse2")
# `.load_esvd2()` and `.build_view()`; the file only defines functions.
source(file.path(loop2_dir, "harness", "generator_claude.R"))
.load_esvd2(repo_dir)

fit_label <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(fit_label %in% fit_label_vec)

# The rules --------------------------------------------------------------------

# One file per rule; each defines `report_lfc_se(view)`.
rule_file_vec <- c(
  stored = file.path(loop2_dir, "candidates", "stored_claude.R"),
  null_rescaled = file.path(loop2_dir, "candidates", "null_rescaled_claude.R"),
  deshrunk = file.path(loop2_dir, "candidates", "deshrunk_claude.R"),
  raw = file.path(loop2_dir, "candidates", "raw_claude.R"),
  raw_null = file.path(loop2_dir, "candidates", "raw_null_claude.R"),
  donor_wls = file.path(loop1_dir, "notes", "RULE_nov11_claude.R"),
  gap_trend = file.path(loop2_dir, "candidates", "gap_trend_claude.R"),
  embed_prior = file.path(loop2_dir, "notes", "RULE_nov11_marginal_claude.R"),
  state_adjusted = file.path(loop2_dir, "notes", "RULE_nov11_direct_claude.R")
)
stopifnot(all(file.exists(rule_file_vec)))
rule_list <- lapply(rule_file_vec, function(rule_file){
  rule_env <- new.env(parent = globalenv())
  sys.source(rule_file, envir = rule_env)
  rule_env$report_lfc_se
})

# Which library each rule is run with, beside the fitted (floored) library.
# Every covariate of this design is constant within a donor, so the
# beyond-score study's `_cov` variants (cell-level covariates replaced by
# their donor means) coincide with the plain ones and are not run.
# `raw` is run with the two libraries the recommendation turns on, so that
# it and `donor_wls` can be compared on one library (added 2026-11-16).
variant_list <- list(
  raw = c("esvd_nofloor", "total"),
  donor_wls = c("esvd_nofloor", "esvd_depth", "total", "sf", "tmm"),
  gap_trend = c("esvd_nofloor", "total"),
  raw_null = c("esvd_nofloor", "total"),
  embed_prior = c("esvd_nofloor", "total"),
  state_adjusted = c("esvd_nofloor", "total")
)

# The fitted half --------------------------------------------------------------

load(file.path(out_dir, paste0("green_", fit_label, "_esvd.RData")))
stopifnot(inherits(esvd_obj, "eSVD"))
load(data_file)

# The package's own rebuild of a diet object's counts and covariates, which
# checks both against summaries recorded at the fit.
rebuilt <- eSVD2:::.get_counts_and_covariates(input_obj = esvd_obj,
                                              seurat_obj = seurat_obj)
rm(seurat_obj)
gene_vec <- colnames(rebuilt$dat)
fit <- esvd_obj[[esvd_obj[["latest_Fit"]]]]
covariates <- rebuilt$covariates
stopifnot(identical(rownames(covariates), rownames(fit$x_mat)),
          identical(colnames(covariates), colnames(fit$z_mat)))

# The rules and `.build_view()` name the case-control column `CC_1`.
cc_variable <- paste0(case_control_var, "_", case_control_levels[2])
stopifnot(cc_variable %in% colnames(covariates))
z_mat <- fit$z_mat[gene_vec, , drop = FALSE]
colnames(covariates)[colnames(covariates) == cc_variable] <- "CC_1"
colnames(z_mat)[colnames(z_mat) == cc_variable] <- "CC_1"

pvalue_list <- esvd_obj$pvalue_list
stored_df <- data.frame(
  gene = gene_vec,
  log2fc = as.numeric(esvd_obj$log2fc_vec[gene_vec]),
  log2fc_se = as.numeric(esvd_obj$log2fc_se_vec[gene_vec]),
  teststat = as.numeric(esvd_obj$teststat_vec[gene_vec]),
  df = as.numeric(pvalue_list$df_vec[gene_vec]),
  gaussian_z = as.numeric(pvalue_list$gaussian_teststat[gene_vec]),
  pvalue = 10^(-as.numeric(pvalue_list$log10pvalue[gene_vec])),
  fdr = as.numeric(pvalue_list$fdr_vec[gene_vec]),
  stringsAsFactors = FALSE
)
attr(stored_df, "null_mean") <- pvalue_list$null_mean
attr(stored_df, "null_sd") <- pvalue_list$null_sd

dat <- as.matrix(rebuilt$dat)
rm(rebuilt)
slim <- list(alpha_max = 2 * max(dat), # what `eSVD()` passes to the posterior
             case_control_vec = as.numeric(covariates[, "CC_1"]),
             covariates = covariates,
             dat = dat,
             individual_vec = esvd_obj$individual,
             nuisance_fit_vec = fit$nuisance_vec[gene_vec],
             nuisance_status_vec = fit$nuisance_status[gene_vec],
             stored_df = stored_df,
             x_mat = fit$x_mat,
             y_mat = fit$y_mat[gene_vec, , drop = FALSE],
             z_mat = z_mat)
rm(dat)
view <- .build_view(slim)
rm(slim)
individual_chr <- as.character(view$individual_vec)

# The view's posterior must be the one the fit's statistics came from: the
# arms' means of the donors' mean posterior reproduce the stored ones.
donor_mean_mat <- rowsum(view$posterior_mean_mat, group = individual_chr) /
  as.numeric(table(individual_chr)[sort(unique(individual_chr))])
donor_cc_vec <- view$case_control_vec[match(rownames(donor_mean_mat),
                                            individual_chr)]
case_mean_vec <- colMeans(donor_mean_mat[donor_cc_vec == 1, , drop = FALSE])
control_mean_vec <- colMeans(donor_mean_mat[donor_cc_vec == 0, , drop = FALSE])
posterior_gap <- max(abs(case_mean_vec / esvd_obj$case_mean[gene_vec] - 1),
                     abs(control_mean_vec / esvd_obj$control_mean[gene_vec] - 1))
print(paste0("Largest relative gap between the rebuilt and the stored arm ",
             "means: ", signif(posterior_gap, 3)))
stopifnot(posterior_gap < 1e-6)

# Library variants -------------------------------------------------------------

# A cells x genes library from a per-cell vector.
.expand_library <- function(library_vec, dat){
  library_mat <- matrix(library_vec / mean(library_vec), nrow = nrow(dat),
                        ncol = ncol(dat))
  dimnames(library_mat) <- dimnames(dat)
  library_mat
}

# One library at a time (each is a cells x genes matrix of doubles).
#   esvd_nofloor  exp of the fit's linear predictor over every covariate but
#                 the arm, without the posterior's floor at 0.1
#   esvd_depth    the same from the intercept and Log_UMI alone
#   total         the cell's summed counts
#   sf            the cell's summed counts, rescaled within each donor so
#                 that the donor's sum is DESeq2's median-of-ratios size
#                 factor of the pseudobulk
#   tmm           the cell's summed counts times edgeR's TMM factor of its
#                 donor's pseudobulk
.make_library <- function(view, variant){
  dat <- view$dat
  individual_chr <- as.character(view$individual_vec)
  if(variant %in% c("esvd_nofloor", "esvd_depth")){
    name_vec <- setdiff(colnames(view$covariates), "CC_1")
    if(variant == "esvd_depth") name_vec <- c("Intercept", "Log_UMI")
    library_mat <- exp(tcrossprod(view$covariates[, name_vec],
                                  view$z_mat[, name_vec]))
    dimnames(library_mat) <- dimnames(dat)
    return(library_mat)
  }
  total_vec <- rowSums(dat)
  if(variant == "total") return(.expand_library(total_vec, dat))
  pseudobulk_mat <- t(rowsum(dat, group = individual_chr))
  if(variant == "sf"){
    donor_total_vec <- colSums(pseudobulk_mat)
    size_factor_vec <- tryCatch(
      DESeq2::estimateSizeFactorsForMatrix(pseudobulk_mat),
      error = function(e){
        DESeq2::estimateSizeFactorsForMatrix(pseudobulk_mat, type = "poscounts")
      }
    )
    return(.expand_library(total_vec * size_factor_vec[individual_chr] /
                             donor_total_vec[individual_chr], dat))
  }
  stopifnot(variant == "tmm")
  tmm_vec <- edgeR::calcNormFactors(pseudobulk_mat, method = "TMM")
  names(tmm_vec) <- colnames(pseudobulk_mat)
  .expand_library(total_vec * tmm_vec[individual_chr], dat)
}

# DESeq2 -----------------------------------------------------------------------

# DESeq2 on pseudobulk counts, one column per donor. The design holds every
# covariate of the fit but the depth, and the arm last. Default `DESeq()`
# and `results()`.
.run_deseq2 <- function(view){
  individual_chr <- as.character(view$individual_vec)
  pseudobulk_mat <- t(rowsum(view$dat, group = individual_chr))
  first_idx <- match(colnames(pseudobulk_mat), individual_chr)
  name_vec <- setdiff(colnames(view$covariates),
                      c("Intercept", "Log_UMI", "CC_1"))
  col_df <- as.data.frame(view$covariates[first_idx, c(name_vec, "CC_1"),
                                          drop = FALSE])
  rownames(col_df) <- colnames(pseudobulk_mat)
  design_formula <- stats::as.formula(paste("~", paste(c(name_vec, "CC_1"),
                                                       collapse = " + ")))
  storage.mode(pseudobulk_mat) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_mat,
                                        colData = col_df,
                                        design = design_formula)
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  res <- DESeq2::results(dds, name = "CC_1")
  idx_vec <- match(colnames(view$dat), rownames(res))
  data.frame(gene = colnames(view$dat),
             log2fc = res$log2FoldChange[idx_vec],
             se = res$lfcSE[idx_vec],
             pvalue = res$pvalue[idx_vec],
             padj = res$padj[idx_vec],
             base_mean = res$baseMean[idx_vec],
             stringsAsFactors = FALSE)
}

# Running the rules ------------------------------------------------------------

.run_rule <- function(rule_fn, view){
  gene_vec <- colnames(view$dat)
  time_start <- proc.time()[["elapsed"]]
  value <- tryCatch(rule_fn(view), error = function(e){e})
  seconds <- proc.time()[["elapsed"]] - time_start
  if(inherits(value, "error") || !is.data.frame(value) ||
     !identical(as.character(value$gene), gene_vec)){
    error <- if(inherits(value, "error")) conditionMessage(value) else
      "malformed return"
    return(list(log2fc = rep(NA_real_, length(gene_vec)),
                se = rep(NA_real_, length(gene_vec)),
                error = error, seconds = seconds))
  }
  list(log2fc = as.numeric(value$log2fc), se = as.numeric(value$se),
       error = "", seconds = seconds)
}

result_list <- list()
for(rule in names(rule_list)){
  # No rule draws random numbers; the seed makes that harmless if one ever
  # does.
  set.seed(10)
  result_list[[rule]] <- .run_rule(rule_list[[rule]], view)
  print(paste0(rule, ": ", round(result_list[[rule]]$seconds), " s ",
               result_list[[rule]]$error))
}

# The fitted library's floor, before the library is exchanged.
variant_vec <- unique(unlist(variant_list))
floor_fraction_vec <- NULL
library_floored_mat <- view$library_mat
for(variant in variant_vec){
  view$library_mat <- .make_library(view = view, variant = variant)
  if(variant == "esvd_nofloor"){
    floor_fraction_vec <- colMeans(view$library_mat < 0.1)
  }
  for(rule in names(variant_list)){
    if(!variant %in% variant_list[[rule]]) next
    label <- paste0(rule, "@", variant)
    set.seed(10)
    result_list[[label]] <- .run_rule(rule_list[[rule]], view)
    print(paste0(label, ": ", round(result_list[[label]]$seconds), " s ",
                 result_list[[label]]$error))
  }
}
view$library_mat <- library_floored_mat
rm(library_floored_mat)

deseq2_df <- .run_deseq2(view)
result_list[["deseq2"]] <- list(log2fc = deseq2_df$log2fc, se = deseq2_df$se,
                                error = "", seconds = NA_real_)

# Order: each rule followed by its variants, DESeq2 last.
method_vec <- unlist(lapply(names(rule_list), function(rule){
  c(rule, paste0(rule, "@", variant_list[[rule]])[
    seq_along(variant_list[[rule]])])
}))
method_vec <- c(method_vec, "deseq2")
stopifnot(setequal(method_vec, names(result_list)))
result_list <- result_list[method_vec]

lfc_mat <- sapply(result_list, function(res){res$log2fc})
se_mat <- sapply(result_list, function(res){res$se})
dimnames(lfc_mat) <- list(gene_vec, method_vec)
dimnames(se_mat) <- list(gene_vec, method_vec)
error_vec <- sapply(result_list, function(res){res$error})
seconds_vec <- sapply(result_list, function(res){res$seconds})

# Per-gene and per-donor records -----------------------------------------------

gene_df <- data.frame(
  gene = gene_vec,
  zero_fraction = colMeans(view$dat == 0),
  mean_count = colMeans(view$dat),
  # Share of the gene's cells whose fitted library is below the posterior's
  # floor of 0.1.
  floor_fraction = floor_fraction_vec,
  status = as.character(view$nuisance_status_vec),
  esvd_pvalue = view$stored_df$pvalue,
  esvd_fdr = view$stored_df$fdr,
  deseq2_pvalue = deseq2_df$pvalue,
  deseq2_padj = deseq2_df$padj,
  deseq2_base_mean = deseq2_df$base_mean,
  stringsAsFactors = FALSE
)
donor_cell_vec <- table(individual_chr)
donor_summary_df <- data.frame(
  donor_id = names(donor_cell_vec),
  num_cells = as.numeric(donor_cell_vec),
  case_control = view$case_control_vec[match(names(donor_cell_vec),
                                             individual_chr)],
  stringsAsFactors = FALSE
)
null_mean <- attr(view$stored_df, "null_mean")
null_sd <- attr(view$stored_df, "null_sd")
num_cells <- nrow(view$dat)
fit_seconds <- seconds

date_of_run <- Sys.time()
session_info <- utils::sessionInfo()
note <- paste0("Every reporting rule of the two autoresearch loops, the ",
               "library variants and pseudobulk DESeq2 on the ", fit_label,
               " fit. lfc_mat and se_mat are genes x methods (log2); gene_df ",
               "has the gene-level records and the two tests' p-values; ",
               "error_vec is empty where a rule ran.")
save(lfc_mat, se_mat, gene_df, donor_summary_df, error_vec, seconds_vec,
     null_mean, null_sd, num_cells, fit_label, fit_seconds, posterior_gap,
     date_of_run, session_info, note,
     file = file.path(out_dir, paste0("green_", fit_label, "_rules.RData")))

print(error_vec[error_vec != ""])
print("Done! :)")
