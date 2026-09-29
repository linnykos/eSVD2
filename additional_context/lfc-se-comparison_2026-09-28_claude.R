# eSVD2's log2 fold change and its SE, beside DESeq2, dreamlet and NEBULA
# Drafted by Claude for Kevin Z. Lin, 2026-09-28
#
# Supplementary evidence for eSVD2 1.1.0 (`compute_log_fold_change`). It is
# NOT a unit test and does not ship: it needs three Bioconductor / CRAN
# packages that are not in `Suggests`. The shipped tests are T-LFC-01 to
# T-LFC-18 in tests/testthat/test_compute_log_fold_change_claude.R.
#
# Run from the package root:
#   Rscript additional_context/lfc-se-comparison_2026-09-28_claude.R
#
# One simulated cohort from `generate_null()`: 20 individuals (10 case, 10
# control), 50 cells each, 100 genes. Genes 1-10 carry a planted effect of
# +/- 1 on the natural-log scale. Everything is reported on the log2 scale;
# NEBULA's natural-log output is divided by ln 2.

library(DESeq2)
library(dreamlet)
library(nebula)
library(SingleCellExperiment)

rm(list = ls())

devtools::load_all(".", quiet = TRUE)

# Simulating the cohort --------------------------------------------------------

set.seed(10)
sim <- generate_null(cell_per_person = 50,
                     num_genes = 100,
                     num_individuals = 20)
count_mat <- as.matrix(sim$obs_mat)
# No underscore in the gene names: several of the packages below rewrite them.
colnames(count_mat) <- paste0("gene", seq_len(ncol(count_mat)))
gene_vec <- colnames(count_mat)
covariates <- sim$covariates
individual_vec <- as.character(sim$metadata_individual)
cc_vec <- covariates[, "CC"]

truth_vec <- log(
  apply(exp(sim$nat_mat[cc_vec == 1, , drop = FALSE]), 2, mean) /
    apply(exp(sim$nat_mat[cc_vec == 0, , drop = FALSE]), 2, mean)
) / log(2)
names(truth_vec) <- gene_vec

# eSVD2, in the order `eSVD()` runs it -----------------------------------------

esvd_obj <- suppressWarnings(
  initialize_esvd(dat = count_mat,
                  covariates = covariates,
                  metadata_individual = factor(individual_vec),
                  bool_intercept = TRUE,
                  case_control_variable = "CC",
                  k = 5,
                  lambda = 0.1,
                  metadata_case_control = cc_vec,
                  verbose = 0)
)
esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                               fit_name = "fit_Init",
                                               omitted_variables = "Log_UMI")
esvd_obj <- suppressWarnings(
  opt_esvd(input_obj = esvd_obj,
           l2pen = 0.1,
           max_iter = 50,
           offset_variables = setdiff(colnames(esvd_obj$covariates), "CC"),
           tol = 1e-6,
           fit_name = "fit_First",
           fit_previous = "fit_Init",
           verbose = 0)
)
esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                               fit_name = "fit_First",
                                               omitted_variables = "Log_UMI")
esvd_obj <- suppressWarnings(
  opt_esvd(input_obj = esvd_obj,
           l2pen = 0.1,
           max_iter = 50,
           offset_variables = NULL,
           tol = 1e-6,
           fit_name = "fit_Second",
           fit_previous = "fit_First",
           verbose = 0)
)
esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                               fit_name = "fit_Second",
                                               omitted_variables = NULL)
esvd_obj <- suppressWarnings(estimate_nuisance(input_obj = esvd_obj,
                                               verbose = 0))
esvd_obj <- compute_posterior(input_obj = esvd_obj,
                              alpha_max = 2 * max(count_mat))
esvd_obj <- compute_test_statistic(input_obj = esvd_obj)
esvd_obj <- compute_pvalue(input_obj = esvd_obj)
esvd_df <- report_results(esvd_obj)
nuisance_vec <- esvd_obj[[esvd_obj[["latest_Fit"]]]]$nuisance_vec

# The three comparators --------------------------------------------------------

## The per-individual metadata, shared by the two pseudobulk methods.
individual_df <- unique(data.frame(Individual = individual_vec,
                                   CC = cc_vec,
                                   Sex = covariates[, "Sex"],
                                   Age = covariates[, "Age"]))
rownames(individual_df) <- individual_df$Individual
individual_df$CC <- factor(individual_df$CC, levels = c(0, 1))
individual_df$Sex <- factor(individual_df$Sex)

## DESeq2 on pseudobulk counts (one column per individual).
pseudobulk_mat <- t(rowsum(count_mat, group = individual_vec))
pseudobulk_mat <- pseudobulk_mat[, rownames(individual_df)]
stopifnot(all(colnames(pseudobulk_mat) == rownames(individual_df)))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_mat,
                                      colData = individual_df,
                                      design = ~ Sex + Age + CC)
dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
deseq_res <- DESeq2::results(dds, name = "CC_1_vs_0")
deseq_df <- data.frame(logFC = deseq_res$log2FoldChange,
                       logFC_se = deseq_res$lfcSE,
                       row.names = rownames(deseq_res))

## dreamlet: pseudobulk, voom weights, a linear model per gene. Its table has
## no SE column; the SE is logFC / t.
cell_df <- data.frame(Individual = individual_vec,
                      CC = factor(cc_vec, levels = c(0, 1)),
                      Sex = factor(covariates[, "Sex"]),
                      Age = covariates[, "Age"],
                      cell_type = "all",
                      row.names = rownames(count_mat))
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = Matrix::Matrix(t(count_mat), sparse = TRUE)),
  colData = cell_df
)
pseudobulk_sce <- dreamlet::aggregateToPseudoBulk(sce,
                                                  assay = "counts",
                                                  cluster_id = "cell_type",
                                                  sample_id = "Individual",
                                                  verbose = FALSE)
processed <- suppressWarnings(
  dreamlet::processAssays(pseudobulk_sce,
                          formula = ~ Sex + Age + CC,
                          min.cells = 1,
                          min.count = 1,
                          min.samples = 2,
                          min.prop = 0,
                          quiet = TRUE)
)
dreamlet_fit <- suppressWarnings(
  dreamlet::dreamlet(processed, formula = ~ Sex + Age + CC, quiet = TRUE)
)
dreamlet_tab <- as.data.frame(
  dreamlet::topTable(dreamlet_fit, coef = "CC1", number = Inf)
)
dreamlet_df <- data.frame(logFC = dreamlet_tab$logFC,
                          logFC_se = dreamlet_tab$logFC / dreamlet_tab$t,
                          row.names = dreamlet_tab$ID)

## NEBULA: a negative binomial mixed model on the cells, with a random effect
## per individual. Natural-log output.
nebula_data <- nebula::group_cell(count = t(count_mat),
                                  id = individual_vec,
                                  pred = stats::model.matrix(~ CC + Sex + Age,
                                                             data = cell_df),
                                  offset = apply(count_mat, 1, sum))
if(is.null(nebula_data)){
  # `group_cell()` returns NULL when the cells are already grouped.
  nebula_data <- list(count = t(count_mat),
                      id = individual_vec,
                      pred = stats::model.matrix(~ CC + Sex + Age,
                                                 data = cell_df),
                      offset = apply(count_mat, 1, sum))
}
nebula_fit <- suppressWarnings(
  nebula::nebula(count = nebula_data$count,
                 id = nebula_data$id,
                 pred = nebula_data$pred,
                 offset = nebula_data$offset,
                 model = "NBGMM",
                 verbose = FALSE)
)
nebula_df <- data.frame(logFC = nebula_fit$summary$logFC_CC1 / log(2),
                        logFC_se = nebula_fit$summary$se_CC1 / log(2),
                        row.names = nebula_fit$summary$gene)

# Putting them side by side ----------------------------------------------------

method_list <- list(eSVD2 = esvd_df[, c("logFC", "logFC_se")],
                    DESeq2 = deseq_df,
                    dreamlet = dreamlet_df,
                    NEBULA = nebula_df)
shared_vec <- gene_vec
for(method_name in names(method_list)){
  method_df <- method_list[[method_name]]
  keep_vec <- rownames(method_df)[is.finite(method_df$logFC) &
                                    is.finite(method_df$logFC_se)]
  shared_vec <- intersect(shared_vec, keep_vec)
}
print(paste0("Genes with a finite estimate and SE from all four methods: ",
             length(shared_vec), " of ", length(gene_vec)))

diverged_vec <- nuisance_vec[shared_vec] > 1e4
print(paste0("Of those, genes whose eSVD2 nuisance estimate exceeds 1e4: ",
             sum(diverged_vec)))

summarize_method <- function(method_df, gene_subset){
  logfc_vec <- method_df[gene_subset, "logFC"]
  se_vec <- method_df[gene_subset, "logFC_se"]
  z_vec <- (logfc_vec - truth_vec[gene_subset]) / se_vec

  c(cor_with_truth = stats::cor(logfc_vec, truth_vec[gene_subset]),
    rmse = sqrt(mean((logfc_vec - truth_vec[gene_subset])^2)),
    median_se = stats::median(se_vec),
    calibration = sqrt(mean((logfc_vec - truth_vec[gene_subset])^2) /
                         mean(se_vec^2)),
    coverage_2se = mean(abs(z_vec) <= 2))
}

subset_list <- list(
  "all shared genes" = shared_vec,
  "eSVD2 nuisance not diverged" = shared_vec[!diverged_vec],
  "eSVD2 nuisance diverged" = shared_vec[diverged_vec]
)
for(subset_name in names(subset_list)){
  gene_subset <- subset_list[[subset_name]]
  if(length(gene_subset) < 3) next
  print(paste0("---- ", subset_name, " (", length(gene_subset), " genes)"))
  summary_mat <- t(sapply(method_list, summarize_method,
                          gene_subset = gene_subset))
  print(round(summary_mat, 3))
}

print("---- Spearman correlation of the log2 fold changes (all shared genes)")
logfc_mat <- sapply(method_list, function(method_df){
  method_df[shared_vec, "logFC"]
})
print(round(stats::cor(logfc_mat, method = "spearman"), 3))

print("---- Spearman correlation of the SEs (all shared genes)")
se_mat <- sapply(method_list, function(method_df){
  method_df[shared_vec, "logFC_se"]
})
print(round(stats::cor(se_mat, method = "spearman"), 3))

print("---- Median over genes of (eSVD2 SE) / (other SE), not-diverged genes")
keep_vec <- shared_vec[!diverged_vec]
ratio_vec <- sapply(setdiff(names(method_list), "eSVD2"), function(method_name){
  stats::median(method_list$eSVD2[keep_vec, "logFC_se"] /
                  method_list[[method_name]][keep_vec, "logFC_se"])
})
print(round(ratio_vec, 3))

print(utils::sessionInfo())
