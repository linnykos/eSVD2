# Fits every simulated cohort once with devel and caches the fit
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# The candidates of 02_run_candidates_claude.R differ only in the nuisance
# rate, so the expensive part (initialize_esvd through estimate_nuisance) is
# done once here. Reuses the data and the devel library that
# additional_context/version_comparison/run_all_claude.sh built, and adds the
# regimes of 00_simulate_extra_claude.R.
#
# Run from the package root:
#   Rscript additional_context/overdispersion_brainstorm/01_cache_fits_claude.R

library(parallel)

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
out_dir <- file.path("additional_context", "overdispersion_brainstorm")
data_dir_vec <- c(file.path(cmp_dir, "data"), file.path(out_dir, "data"))
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

library(eSVD2, lib.loc = file.path(cmp_dir, "lib", "devel_1.1.0"))
source(file.path(cmp_dir, "helpers_claude.R"))

# Fitting ----------------------------------------------------------------------

file_vec <- list.files(data_dir_vec, pattern = "\\.rds$", full.names = TRUE)
stem_vec <- sub("\\.rds$", "", basename(file_vec))
stopifnot(anyDuplicated(stem_vec) == 0)

res_list <- parallel::mclapply(seq_along(stem_vec), function(i){
  stem <- stem_vec[i]
  sim <- readRDS(file_vec[i])
  if(is.null(sim$covariates)){
    numeric_var_vec <- intersect("Age", colnames(sim$covariate_df))
    if(length(numeric_var_vec) == 0) numeric_var_vec <- NULL
    covariates <- eSVD2::format_covariates(dat = sim$dat,
                                           covariate_df = sim$covariate_df,
                                           rescale_numeric_variables = numeric_var_vec)
    cc_var <- "CC_1"
  } else {
    covariates <- sim$covariates
    cc_var <- "CC"
  }

  # The whole pipeline, so that the cached object is exactly what devel
  # produced; the candidates overwrite its last three stages.
  res <- .run_pipeline(dat = sim$dat,
                       covariates = covariates,
                       cc_var = cc_var,
                       individual_vec = sim$individual_vec)
  stopifnot(is.null(res$error_message))

  saveRDS(list(cc_var = cc_var,
               esvd_obj = res$esvd_obj,
               gene_vec = sim$gene_vec,
               is_de_vec = sim$is_de_vec,
               nuisance_true_vec = sim$nuisance_true_vec,
               truth_log2fc_vec = sim$truth_log2fc_vec),
          file = file.path(cache_dir, paste0(stem, ".rds")))
  stem
}, mc.cores = 8)

stopifnot(all(unlist(res_list) == stem_vec))
print(paste0("Cached ", length(stem_vec), " fits in ", cache_dir))
