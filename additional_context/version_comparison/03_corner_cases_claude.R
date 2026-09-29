# Corner-case behavior of ONE version of eSVD2
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Each case perturbs one small simulated cohort in a way NEWS.md says the
# devel version handles differently, runs the full step-by-step pipeline, and
# records what happened: finished silently, finished with warnings, or stopped
# (and at which stage), plus the number of non-finite outputs and one
# case-specific number (`metric_value`).
#
# Run from the package root, once per version:
#   Rscript additional_context/version_comparison/03_corner_cases_claude.R master
#   Rscript additional_context/version_comparison/03_corner_cases_claude.R devel
#
# Writes output/<label>/corner_cases.csv.

rm(list = ls())

label <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(label %in% c("master", "devel"))

cmp_dir <- file.path("additional_context", "version_comparison")
out_dir <- file.path(cmp_dir, "output", label)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lib_dir <- file.path(cmp_dir, "lib", label)
library(eSVD2, lib.loc = lib_dir)
print(paste0("Loaded eSVD2 ", utils::packageVersion("eSVD2", lib.loc = lib_dir),
             " for `", label, "`"))

source(file.path(cmp_dir, "helpers_claude.R"))

# Helpers ----------------------------------------------------------------------

.format_sim <- function(sim){
  numeric_var_vec <- intersect("Age", colnames(sim$covariate_df))
  if(length(numeric_var_vec) == 0) numeric_var_vec <- NULL
  eSVD2::format_covariates(dat = sim$dat,
                           covariate_df = sim$covariate_df,
                           rescale_numeric_variables = numeric_var_vec)
}

# Keeps the cells in `cell_idx` of a simulated cohort.
.subset_cells <- function(sim,
                          cell_idx){
  sim$dat <- sim$dat[cell_idx, , drop = FALSE]
  sim$covariate_df <- sim$covariate_df[cell_idx, , drop = FALSE]
  sim$individual_vec <- droplevels(sim$individual_vec[cell_idx])
  sim$cc_vec <- sim$cc_vec[cell_idx]
  sim
}

# Runs one case and returns its row of the results table. `covariates` is
# formatted from `sim` unless supplied; `metric_fun` maps the gene table to the
# case-specific number.
.run_case <- function(case_name,
                      sim,
                      description,
                      covariates = NULL,
                      metric_fun = NULL,
                      metric_name = ""){
  if(is.null(covariates)){
    covariates <- tryCatch(.format_sim(sim), error = function(e){e})
  }
  if(inherits(covariates, "error")){
    res <- list(error_message = conditionMessage(covariates),
                esvd_obj = NULL,
                stage = "format_covariates",
                warning_vec = character(0))
  } else {
    res <- .run_pipeline(dat = sim$dat,
                         covariates = covariates,
                         cc_var = "CC_1",
                         individual_vec = sim$individual_vec)
  }

  gene_df <- NULL
  if(is.null(res$error_message)){
    gene_df <- .extract_gene_df(esvd_obj = res$esvd_obj,
                                gene_vec = colnames(sim$dat))
  }
  metric_value <- NA
  if(!is.null(gene_df) && !is.null(metric_fun)) metric_value <- metric_fun(gene_df)

  outcome <- "silent"
  if(length(res$warning_vec) > 0) outcome <- "warning"
  if(!is.null(res$error_message)) outcome <- "error"
  message_val <- ifelse(is.null(res$error_message), "",
                        paste0(res$stage, ": ", res$error_message))
  method_val <- NA
  if(!is.null(res$esvd_obj$pvalue_list$method)){
    method_val <- res$esvd_obj$pvalue_list$method
  }

  print(paste0(case_name, ": ", outcome))
  list(row_df = data.frame(case = case_name,
                           description = description,
                           outcome = outcome,
                           error = message_val,
                           warnings = paste0(unique(res$warning_vec),
                                             collapse = " || "),
                           num_nonfinite_teststat = ifelse(is.null(gene_df), NA,
                                                           sum(!is.finite(gene_df$teststat))),
                           num_nonfinite_p = ifelse(is.null(gene_df), NA,
                                                    sum(!is.finite(gene_df$log10p))),
                           num_sig = ifelse(is.null(gene_df), NA,
                                            sum(gene_df$fdr < 0.05, na.rm = TRUE)),
                           method = method_val,
                           null_mean = ifelse(is.null(res$esvd_obj), NA,
                                              as.numeric(res$esvd_obj$pvalue_list$null_mean)),
                           null_sd = ifelse(is.null(res$esvd_obj), NA,
                                            as.numeric(res$esvd_obj$pvalue_list$null_sd)),
                           metric_name = metric_name,
                           metric_value = metric_value),
       gene_df = gene_df)
}

# The cases --------------------------------------------------------------------

# 10 individuals (5 / 5) x 30 cells, 100 genes, 10 of them DE.
base_sim <- .simulate_cohort(num_cells_per_individual = 30,
                             num_de = 10,
                             num_genes = 100,
                             num_individuals = 10,
                             seed_number = 7)
case_list <- list()

## The reference run, which several cases compare against.
reference_res <- .run_case("reference", base_sim,
                           description = "Unperturbed cohort (reference)")
case_list$reference <- reference_res$row_df
reference_teststat_vec <- reference_res$gene_df$teststat

## Determinism: the same data with a different ambient RNG state.
set.seed(1)
tmp1 <- .run_case("determinism", base_sim, description = "")
set.seed(2)
tmp2 <- .run_case("determinism", base_sim,
                  description = "Same data, fitted after set.seed(1) and after set.seed(2)",
                  metric_fun = function(gene_df){
                    max(abs(gene_df$teststat - tmp1$gene_df$teststat))
                  },
                  metric_name = "max |diff| in Welch statistic between the two fits")
case_list$determinism <- tmp2$row_df

## Sparse input, which is how a Seurat object hands over its counts.
sim <- base_sim
sim$dat <- methods::as(methods::as(sim$dat, "dMatrix"), "CsparseMatrix")
case_list$sparse_input <- .run_case(
  "sparse_input", sim,
  description = "Same counts as a dgCMatrix",
  metric_fun = function(gene_df){
    max(abs(gene_df$teststat - reference_teststat_vec))
  },
  metric_name = "max |diff| in Welch statistic vs dense reference"
)$row_df

## An NA in a sparse count matrix.
sim$dat@x[5] <- NA
case_list$sparse_na <- .run_case(
  "sparse_na", sim,
  description = "One NA among the nonzero entries of a dgCMatrix"
)$row_df

## Two genes with no counts at all.
sim <- base_sim
sim$dat[, c(3, 50)] <- 0
case_list$all_zero_gene <- .run_case(
  "all_zero_gene", sim,
  description = "Two genes are all-zero"
)$row_df

## A covariate that duplicates another.
covariates <- .format_sim(base_sim)
covariates <- cbind(covariates, Age_copy = covariates[, "Age"])
case_list$collinear_covariate <- .run_case(
  "collinear_covariate", base_sim,
  description = "A covariate column duplicated (Age twice)",
  covariates = covariates
)$row_df

## Sex perfectly confounded with case-control status.
sim <- base_sim
sim$covariate_df$Sex <- factor(ifelse(sim$cc_vec == 1, "M", "F"))
case_list$confounded_design <- .run_case(
  "confounded_design", sim,
  description = "Sex identical to case-control status"
)$row_df

## Few genes, where locfdr's spline fit fails.
sim <- base_sim
sim$dat <- sim$dat[, 1:30, drop = FALSE]
case_list$few_genes <- .run_case(
  "few_genes", sim,
  description = "Only 30 genes"
)$row_df

## Very strong up-regulation of nearly Poisson genes, the case for which the
## Welch statistic is largest.
effect_vec <- c(rep(4, 5), rep(0, 95))
sim <- .simulate_cohort(cc_effect_vec = effect_vec,
                        num_cells_per_individual = 30,
                        num_genes = 100,
                        num_individuals = 10,
                        nuisance_range = c(500, 1000),
                        seed_number = 8)
case_list$extreme_de <- .run_case(
  "extreme_de", sim,
  description = "5 genes up-regulated by exp(4), nearly Poisson",
  metric_fun = function(gene_df){max(gene_df$gaussian_teststat)},
  metric_name = "largest Gaussian statistic"
)$row_df

## One case individual only.
cell_idx <- which(base_sim$cc_vec == 0 |
                    base_sim$individual_vec == "indiv6")
case_list$one_case_individual <- .run_case(
  "one_case_individual", .subset_cells(base_sim, cell_idx),
  description = "The case arm has a single individual"
)$row_df

## One individual contributes only two cells.
drop_idx <- which(base_sim$individual_vec == "indiv1")[-(1:2)]
case_list$two_cell_individual <- .run_case(
  "two_cell_individual", .subset_cells(base_sim, -drop_idx),
  description = "One control individual has only 2 cells"
)$row_df

## Half of one individual's cells labelled as case.
sim <- base_sim
idx_vec <- which(sim$individual_vec == "indiv1")[1:15]
sim$cc_vec[idx_vec] <- 1
sim$covariate_df$CC[idx_vec] <- "1"
case_list$individual_in_both_arms <- .run_case(
  "individual_in_both_arms", sim,
  description = "Half of one individual's cells labelled case, half control"
)$row_df

corner_df <- do.call(rbind, case_list)
corner_df$version <- label
utils::write.csv(corner_df, file.path(out_dir, "corner_cases.csv"),
                 row.names = FALSE)
