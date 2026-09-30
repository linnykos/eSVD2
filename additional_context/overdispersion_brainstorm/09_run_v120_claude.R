# Runs the shipped cap (eSVD2 1.2.0) end to end on every brainstorm regime
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# The dry-runs 02-08 prototyped the cap by overwriting the 1.1.0 rate with
# pmin(MLE, 10 * median_i s_ji) after estimate_nuisance(). 1.2.0 implements
# it inside estimate_nuisance(), with three differences from the prototype:
#   1. a boundary gene (D_j <= 0) is set to the cap itself, whatever value the
#      optimizer stopped at;
#   2. the floor is min_val * median_i s_ji, not an absolute min_val;
#   3. the cap is computed from the library the package itself builds.
# This script fits every data set with 1.2.0 at its default cap_multiplier = 10,
# with no override, and puts its results beside the prototype's `cap_10s`, so
# the report can say whether the shipped cap is the one it evaluated.
#
# Run from the package root, after 02_run_candidates_claude.R:
#   Rscript additional_context/overdispersion_brainstorm/09_run_v120_claude.R
#
# Writes output/report_v120_by_rep.csv (one row per data set) and
# output/report_v120_agreement.csv (one row per regime).

library(parallel)

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
work_dir <- file.path("additional_context", "overdispersion_brainstorm")
out_dir <- file.path(work_dir, "output")

lib_dir <- file.path(cmp_dir, "lib", "devel")
library(eSVD2, lib.loc = lib_dir)
stopifnot(utils::packageVersion("eSVD2", lib.loc = lib_dir) >= "1.2.0")
print(paste0("Loaded eSVD2 ", utils::packageVersion("eSVD2", lib.loc = lib_dir)))
source(file.path(cmp_dir, "helpers_claude.R"))

# Fitting ----------------------------------------------------------------------

path_vec <- c(list.files(file.path(cmp_dir, "data"), pattern = "\\.rds$",
                         full.names = TRUE),
              list.files(file.path(work_dir, "data"), pattern = "\\.rds$",
                         full.names = TRUE))

res_list <- parallel::mclapply(path_vec, function(path){
  stem <- sub("\\.rds$", "", basename(path))
  sim <- readRDS(path)

  # generate_null() hands over formatted covariates; every other regime hands
  # over raw metadata (as in version_comparison/02_run_regimes_claude.R).
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

  res <- .run_pipeline(dat = sim$dat,
                       covariates = covariates,
                       cc_var = cc_var,
                       individual_vec = sim$individual_vec)
  stopifnot(is.null(res$error_message))
  gene_df <- .extract_gene_df(esvd_obj = res$esvd_obj, gene_vec = sim$gene_vec)
  mle_vec <- res$esvd_obj[[res$esvd_obj[["latest_Fit"]]]]$nuisance_mle_vec
  cbind(data.frame(regime = sub("__rep.*$", "", stem),
                   rep = as.integer(sub("^.*__rep", "", stem)),
                   is_de = sim$is_de_vec),
        gene_df,
        nuisance_mle = as.numeric(mle_vec))
}, mc.cores = 8)
stopifnot(!any(sapply(res_list, function(x){inherits(x, "try-error")})))
v120_df <- do.call(rbind, res_list)

# Discoveries and cap counts per data set --------------------------------------

by_rep_list <- lapply(split(v120_df, list(v120_df$regime, v120_df$rep),
                            drop = TRUE), function(x){
  data.frame(regime = x$regime[1],
             rep = x$rep[1],
             fp = sum(x$fdr < 0.05 & !x$is_de),
             tp = sum(x$fdr < 0.05 & x$is_de),
             num_de = sum(x$is_de),
             num_capped = sum(x$nuisance_status %in% c("capped", "boundary")),
             num_boundary = sum(x$nuisance_status == "boundary"),
             num_failed = sum(x$nuisance_status == "failed"),
             num_genes = nrow(x))
})
by_rep_df <- do.call(rbind, by_rep_list)
by_rep_df <- by_rep_df[order(by_rep_df$regime, by_rep_df$rep), ]
utils::write.csv(by_rep_df, file.path(out_dir, "report_v120_by_rep.csv"),
                 row.names = FALSE)

# Agreement with the prototype's cap_10s ---------------------------------------

# output/genes.csv is gitignored (made by 02_run_candidates_claude.R); the
# agreement table is skipped, with a message, when it is absent.
proto_path <- file.path(out_dir, "genes.csv")
if(file.exists(proto_path)){
  proto_df <- utils::read.csv(proto_path)
  proto_df <- proto_df[proto_df$candidate == "cap_10s",
                       c("regime", "rep", "gene", "nuisance", "fdr",
                         "gaussian_teststat")]
  colnames(proto_df)[4:6] <- c("nuisance_proto", "fdr_proto",
                               "gaussian_teststat_proto")
  both_df <- merge(v120_df, proto_df, by = c("regime", "rep", "gene"))
  stopifnot(nrow(both_df) == nrow(v120_df))

  agreement_list <- lapply(split(both_df, both_df$regime), function(x){
    ratio_vec <- x$nuisance / x$nuisance_proto
    data.frame(regime = x$regime[1],
               num_genes = nrow(x),
               frac_rate_within_1pct = mean(abs(ratio_vec - 1) < 0.01),
               max_rate_ratio = max(ratio_vec),
               min_rate_ratio = min(ratio_vec),
               max_abs_diff_teststat = max(abs(x$gaussian_teststat -
                                                 x$gaussian_teststat_proto)),
               num_fdr_calls_differ = sum((x$fdr < 0.05) !=
                                            (x$fdr_proto < 0.05)))
  })
  agreement_df <- do.call(rbind, agreement_list)
  num_col_vec <- sapply(agreement_df, is.numeric)
  agreement_df[num_col_vec] <- lapply(agreement_df[num_col_vec], signif,
                                      digits = 4)
  utils::write.csv(agreement_df,
                   file.path(out_dir, "report_v120_agreement.csv"),
                   row.names = FALSE)
  print(agreement_df)
} else {
  print("output/genes.csv not found; skipping the agreement table")
}
print(paste0("Ran 1.2.0 on ", length(path_vec), " data sets"))
