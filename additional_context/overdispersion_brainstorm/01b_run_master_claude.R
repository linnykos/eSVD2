# Runs master (3d5f7bf) on the four regimes added for the brainstorm
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# additional_context/version_comparison ran master on its own six regimes
# only. This runs it on weak_de, low_count, wide_rate and trend, so that
# "master" in the report is master's own output in every regime and not the
# legacy cap applied to devel's estimate. It attaches the master library, so it
# must run in its own R process.
#
# Run from the package root, after 00_simulate_extra_claude.R:
#   Rscript additional_context/overdispersion_brainstorm/01b_run_master_claude.R
#
# Writes output/master_extra_genes.csv.

library(parallel)

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
work_dir <- file.path("additional_context", "overdispersion_brainstorm")
data_dir <- file.path(work_dir, "data")
out_dir <- file.path(work_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(eSVD2, lib.loc = file.path(cmp_dir, "lib", "master"))
print(paste0("Loaded eSVD2 ",
             utils::packageVersion("eSVD2",
                                   lib.loc = file.path(cmp_dir, "lib",
                                                       "master"))))
source(file.path(cmp_dir, "helpers_claude.R"))

# Fitting ----------------------------------------------------------------------

stem_vec <- sub("\\.rds$", "", list.files(data_dir, pattern = "\\.rds$"))

res_list <- parallel::mclapply(stem_vec, function(stem){
  sim <- readRDS(file.path(data_dir, paste0(stem, ".rds")))
  covariates <- eSVD2::format_covariates(dat = sim$dat,
                                         covariate_df = sim$covariate_df,
                                         rescale_numeric_variables = "Age")
  res <- .run_pipeline(dat = sim$dat,
                       covariates = covariates,
                       cc_var = "CC_1",
                       individual_vec = sim$individual_vec)
  stopifnot(is.null(res$error_message))
  cbind(data.frame(regime = sub("__rep.*$", "", stem),
                   rep = as.integer(sub("^.*__rep", "", stem))),
        .extract_gene_df(esvd_obj = res$esvd_obj, gene_vec = sim$gene_vec))
}, mc.cores = 8)
stopifnot(!any(sapply(res_list, function(x){inherits(x, "try-error")})))

gene_df <- do.call(rbind, res_list)
num_col_vec <- sapply(gene_df, is.numeric)
gene_df[num_col_vec] <- lapply(gene_df[num_col_vec], signif, digits = 6)
utils::write.csv(gene_df, file.path(out_dir, "master_extra_genes.csv"),
                 row.names = FALSE)
print(paste0("Ran master on ", length(stem_vec), " data sets"))
