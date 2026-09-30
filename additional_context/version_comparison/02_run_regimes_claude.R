# Fits every simulated regime with ONE version of eSVD2
# Drafted by Claude for Kevin Z. Lin, 2026-09-29; updated 2026-09-29 for 1.2.0
#
# The same script serves every version: the step-by-step API (initialize_esvd
# through compute_pvalue) has the same function names and arguments at
# 3d5f7bf and at 1.2.0, except for `cap_multiplier`, which only 1.2.0 has.
#
# Run from the package root, in this order:
#   Rscript additional_context/version_comparison/02_run_regimes_claude.R master
#   Rscript additional_context/version_comparison/02_run_regimes_claude.R devel
#   Rscript additional_context/version_comparison/02_run_regimes_claude.R devel_nocap
#   Rscript additional_context/version_comparison/02_run_regimes_claude.R devel_swap
#
# `devel` is 1.2.0 at its default `cap_multiplier = 10`. `devel_nocap` is 1.2.0
# with `cap_multiplier = Inf`, whose rates are those of 1.1.0 (the uncapped
# estimate this comparison reported on first); it is kept as a reference arm.
#
# `devel_swap` is an ablation, not a version: the devel code, with the
# per-gene nuisance rate replaced by the value master estimated on the same
# data set, right after estimate_nuisance(). Where devel_swap reproduces master,
# the difference between the versions is the nuisance estimate; where it
# reproduces devel, the difference lies elsewhere. It needs master's output.
#
# Writes output/<label>/genes.csv (one row per gene per data set) and
# output/<label>/runs.csv (one row per data set).

rm(list = ls())

label <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(label %in% c("master", "devel", "devel_nocap", "devel_swap"))
lib_label <- ifelse(label == "master", "master", "devel")

cmp_dir <- file.path("additional_context", "version_comparison")
data_dir <- file.path(cmp_dir, "data")
out_dir <- file.path(cmp_dir, "output", label)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lib_dir <- file.path(cmp_dir, "lib", lib_label)
library(eSVD2, lib.loc = lib_dir)
print(paste0("Loaded eSVD2 ", utils::packageVersion("eSVD2", lib.loc = lib_dir),
             " for `", label, "`"))

# NULL leaves the argument out of the call, which master requires.
cap_multiplier <- NULL
if(label == "devel_nocap") cap_multiplier <- Inf

if(label == "devel_swap"){
  master_gene_df <- utils::read.csv(file.path(cmp_dir, "output", "master",
                                              "genes.csv"))
}

source(file.path(cmp_dir, "helpers_claude.R"))

# Fitting every data set -------------------------------------------------------

stem_vec <- sub("\\.rds$", "", list.files(data_dir, pattern = "\\.rds$"))
run_list <- vector("list", length(stem_vec))
gene_list <- vector("list", length(stem_vec))
names(run_list) <- stem_vec
names(gene_list) <- stem_vec

for(stem in stem_vec){
  regime <- sub("__rep.*$", "", stem)
  rep_idx <- as.integer(sub("^.*__rep", "", stem))
  sim <- readRDS(file.path(data_dir, paste0(stem, ".rds")))

  # generate_null() hands over formatted covariates; every other regime hands
  # over raw metadata, which each version formats with its own
  # format_covariates().
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

  nuisance_override_vec <- NULL
  if(label == "devel_swap"){
    idx_vec <- which(master_gene_df$regime == regime &
                       master_gene_df$rep == rep_idx)
    stopifnot(length(idx_vec) == length(sim$gene_vec),
              all(master_gene_df$gene[idx_vec] == sim$gene_vec))
    nuisance_override_vec <- master_gene_df$nuisance[idx_vec]
  }

  start_time <- Sys.time()
  res <- .run_pipeline(dat = sim$dat,
                       covariates = covariates,
                       cc_var = cc_var,
                       individual_vec = sim$individual_vec,
                       cap_multiplier = cap_multiplier,
                       nuisance_override_vec = nuisance_override_vec)
  elapsed_val <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if(!is.null(res$error_message)){
    print(paste0(stem, ": ERROR in ", res$stage, ": ", res$error_message))
    run_list[[stem]] <- data.frame(regime = regime,
                                   rep = rep_idx,
                                   status = "error",
                                   message = paste0(res$stage, ": ",
                                                    res$error_message),
                                   seconds = elapsed_val,
                                   num_warnings = length(res$warning_vec),
                                   warnings = paste0(unique(res$warning_vec),
                                                     collapse = " || "),
                                   null_mean = NA,
                                   null_sd = NA,
                                   method = NA,
                                   num_nonfinite_z = NA)
    next()
  }

  gene_df <- .extract_gene_df(esvd_obj = res$esvd_obj,
                              gene_vec = sim$gene_vec)
  gene_list[[stem]] <- cbind(data.frame(regime = regime, rep = rep_idx),
                             gene_df)

  # master's multtest() does not report which empirical-null estimator ran.
  method_val <- res$esvd_obj$pvalue_list$method
  if(is.null(method_val)) method_val <- NA
  run_list[[stem]] <- data.frame(regime = regime,
                                 rep = rep_idx,
                                 status = "ok",
                                 message = "",
                                 seconds = elapsed_val,
                                 num_warnings = length(res$warning_vec),
                                 warnings = paste0(unique(res$warning_vec),
                                                   collapse = " || "),
                                 null_mean = as.numeric(res$esvd_obj$pvalue_list$null_mean),
                                 null_sd = as.numeric(res$esvd_obj$pvalue_list$null_sd),
                                 method = method_val,
                                 num_nonfinite_z = sum(!is.finite(gene_df$gaussian_teststat)))
  print(paste0(stem, ": ", round(elapsed_val, 1), " s, ",
               length(res$warning_vec), " warnings"))
}

run_df <- do.call(rbind, run_list)
run_df$version <- label
gene_df <- do.call(rbind, gene_list)
gene_df$version <- label
utils::write.csv(run_df, file.path(out_dir, "runs.csv"), row.names = FALSE)
# Six significant digits keeps the three tracked gene tables near 2 MB each.
num_col_vec <- sapply(gene_df, is.numeric)
gene_df[num_col_vec] <- lapply(gene_df[num_col_vec], signif, digits = 6)
utils::write.csv(gene_df, file.path(out_dir, "genes.csv"), row.names = FALSE)
print(paste0("Total: ", round(sum(run_df$seconds), 1), " s"))
