rm(list = ls())

# Real-data checks on Tati's four cohorts (Writeup50)
# Drafted by Claude for Kevin Z. Lin, 2026-11-07
#
# Reads the Writeup50 per-cohort CSVs (one row per gene, the four methods'
# logFC / SE / p-value) and the empirical-null table, plus Tati's diet eSVD
# objects for the nuisance status and the Gaussianized statistics. No refit.
# Writes to output/:
#   01_sign-agreement-recentered.csv   eSVD-DE sign agreement under logFC
#                                      recentering, per ordered pair
#   01_sign-agreement-cross-method.csv same genes (A's DEGs under method M_A),
#                                      sign in B under every method M_B
#   01_lfc-cross-method.csv            within-cohort agreement of the four
#                                      methods' logFC and SE
#   01_null-by-status.csv              the empirical null fitted within
#                                      nuisance-status strata
#   01_deg-by-status.csv               eSVD-DE's DEGs split by status, with
#                                      their sign in the other cohorts
# Runs on Kevin's laptop in about a minute.

library(locfdr)

exploration_dir <- getwd()
source(file.path(exploration_dir, "helpers_claude.R"))
.load_esvd2_devel(exploration_dir)

# Paths (laptop) ---------------------------------------------------------------

was2code_dir <- file.path("~", "Library", "CloudStorage", "Dropbox",
                          "Collaboration-and-People", "archive", "tati")
csv_dir <- file.path(was2code_dir, "git", "Was2CoDE_analysis", "csv", "kevin",
                     "Writeup50")
tati_out_dir <- file.path(was2code_dir, "out", "tati")
output_dir <- file.path(exploration_dir, "output")

dataset_vec <- c("prater", "gabitto", "sun", "green")
method_vec <- c("deseq2", "dreamlet", "nebula", "esvd")
writeup_vec <- c(prater = "Writeup44",
                 gabitto = "Writeup45",
                 sun = "Writeup46",
                 green = "Writeup47")
fdr_threshold <- 0.05
# The Rmd's convergence caps (NEBULA's non-converged genes).
lfc_cap <- 10
se_cap <- 100

# Reading ----------------------------------------------------------------------

csv_list <- lapply(dataset_vec, function(dataset){
  utils::read.csv(file.path(csv_dir, paste0("Writeup50_", dataset,
                                            "_de-results.csv")),
                  stringsAsFactors = FALSE)
})
names(csv_list) <- dataset_vec
gene_vec <- csv_list[[1]]$gene
stopifnot(all(sapply(csv_list, function(df){identical(df$gene, gene_vec)})))
null_df <- utils::read.csv(file.path(csv_dir,
                                     "Writeup50_esvd-empirical-null.csv"),
                           stringsAsFactors = FALSE)

# The Rmd's method table: converged flag and the DEG flag (BH within the
# converged genes).
.method_table <- function(csv_df, method){
  lfc_vec <- csv_df[[paste0(method, "_logFC")]]
  se_vec <- csv_df[[paste0(method, "_se")]]
  pvalue_vec <- csv_df[[paste0(method, "_pvalue")]]
  bool_converged_vec <- is.finite(lfc_vec) & is.finite(se_vec) &
    is.finite(pvalue_vec) & abs(lfc_vec) <= lfc_cap & se_vec <= se_cap
  padj_vec <- rep(NA_real_, length(lfc_vec))
  padj_vec[bool_converged_vec] <- stats::p.adjust(
    pvalue_vec[bool_converged_vec], method = "BH")
  data.frame(gene = csv_df$gene,
             lfc = lfc_vec,
             se = se_vec,
             pvalue = pvalue_vec,
             bool_converged = bool_converged_vec,
             bool_deg = !is.na(padj_vec) & padj_vec <= fdr_threshold,
             stringsAsFactors = FALSE)
}
stat_list <- lapply(dataset_vec, function(dataset){
  tmp <- lapply(method_vec, function(method){
    .method_table(csv_list[[dataset]], method)
  })
  names(tmp) <- method_vec
  tmp
})
names(stat_list) <- dataset_vec

# 1. eSVD-DE sign agreement under logFC recentering -----------------------------

# The discovery set is A's eSVD-DE DEGs (stored p-values) within the genes
# converged in both cohorts; only B's logFC is recentred, since the sign of a
# DEG in A does not change under a shift that is small against its effect.
# Variants: the empirical-null shift `- null_mean * se` (option A), the
# median over the universe, and the mean over the central 90% of genes.
.recentre <- function(lfc_vec, se_vec, bool_universe, variant, null_mean){
  if(variant == "stored") return(lfc_vec)
  if(variant == "null_shift") return(lfc_vec - null_mean * se_vec)
  if(variant == "median") return(lfc_vec - stats::median(lfc_vec[bool_universe]))
  if(variant == "trimmed_mean"){
    return(lfc_vec - mean(lfc_vec[bool_universe], trim = 0.05))
  }
  stop("unknown variant")
}
variant_vec <- c("stored", "null_shift", "median", "trimmed_mean")

recentred_row_list <- list()
for(cohort_a in dataset_vec){
  for(cohort_b in setdiff(dataset_vec, cohort_a)){
    a_df <- stat_list[[cohort_a]]$esvd
    b_df <- stat_list[[cohort_b]]$esvd
    bool_universe <- a_df$bool_converged & b_df$bool_converged
    idx_discovery <- which(bool_universe & a_df$bool_deg)
    null_mean_b <- null_df$null_mean[null_df$dataset == cohort_b]
    for(variant in variant_vec){
      lfc_b_vec <- .recentre(b_df$lfc, b_df$se, bool_universe, variant,
                             null_mean_b)
      recentred_row_list[[length(recentred_row_list) + 1]] <- data.frame(
        cohort_a = cohort_a,
        cohort_b = cohort_b,
        variant = variant,
        num_discovery = length(idx_discovery),
        shift_b = stats::median(b_df$lfc[bool_universe]) -
          stats::median(lfc_b_vec[bool_universe]),
        frac_same_sign = mean(sign(a_df$lfc[idx_discovery]) ==
                                sign(lfc_b_vec[idx_discovery])),
        stringsAsFactors = FALSE)
    }
  }
}
recentred_df <- do.call(rbind, recentred_row_list)
utils::write.csv(recentred_df,
                 file = file.path(output_dir,
                                  "01_sign-agreement-recentered.csv"),
                 row.names = FALSE)
print("Sign agreement of eSVD-DE under recentering (fraction same sign):")
print(stats::reshape(recentred_df[, c("cohort_a", "cohort_b", "variant",
                                      "frac_same_sign")],
                     idvar = c("cohort_a", "cohort_b"),
                     timevar = "variant", direction = "wide"),
      digits = 3)

# 2. Same genes, every method's sign in B ---------------------------------------

# Separates "eSVD-DE's DEGs in A replicate less" from "eSVD-DE's logFC in B is
# noisier": for the DEGs of A under method M_A (universe: converged under M_A
# in A and under M_B in B), the fraction whose sign in B under method M_B
# equals their sign in A under M_A.
cross_row_list <- list()
for(cohort_a in dataset_vec){
  for(cohort_b in setdiff(dataset_vec, cohort_a)){
    for(method_a in method_vec){
      a_df <- stat_list[[cohort_a]][[method_a]]
      for(method_b in method_vec){
        b_df <- stat_list[[cohort_b]][[method_b]]
        bool_universe <- a_df$bool_converged & b_df$bool_converged
        idx_discovery <- which(bool_universe & a_df$bool_deg)
        cross_row_list[[length(cross_row_list) + 1]] <- data.frame(
          cohort_a = cohort_a,
          cohort_b = cohort_b,
          method_a = method_a,
          method_b = method_b,
          num_discovery = length(idx_discovery),
          frac_same_sign = if(length(idx_discovery) > 0){
            mean(sign(a_df$lfc[idx_discovery]) ==
                   sign(b_df$lfc[idx_discovery]))
          } else NA_real_,
          # How far B's estimates of these genes are from zero in SE units
          # under M_B: a sign is reliable only when |lfc| is a few SE.
          median_abs_z_b = if(length(idx_discovery) > 0){
            stats::median(abs(b_df$lfc[idx_discovery] /
                                b_df$se[idx_discovery]))
          } else NA_real_,
          stringsAsFactors = FALSE)
      }
    }
  }
}
cross_df <- do.call(rbind, cross_row_list)
utils::write.csv(cross_df,
                 file = file.path(output_dir,
                                  "01_sign-agreement-cross-method.csv"),
                 row.names = FALSE)
print("Pairs into Prater: sign in Prater under each method, by discovery method")
print(cross_df[cross_df$cohort_b == "prater",
               c("cohort_a", "method_a", "method_b", "num_discovery",
                 "frac_same_sign")],
      digits = 3)

# 3. Within-cohort agreement of the four methods ---------------------------------

lfc_row_list <- list()
for(dataset in dataset_vec){
  bool_all <- Reduce(`&`, lapply(stat_list[[dataset]], function(df){
    df$bool_converged
  }))
  for(method_a in method_vec){
    for(method_b in method_vec){
      a_df <- stat_list[[dataset]][[method_a]]
      b_df <- stat_list[[dataset]][[method_b]]
      lfc_row_list[[length(lfc_row_list) + 1]] <- data.frame(
        dataset = dataset,
        method_a = method_a,
        method_b = method_b,
        num_genes = sum(bool_all),
        spearman_lfc = stats::cor(a_df$lfc[bool_all], b_df$lfc[bool_all],
                                  method = "spearman"),
        spearman_se = stats::cor(a_df$se[bool_all], b_df$se[bool_all],
                                 method = "spearman"),
        median_se_ratio = stats::median(a_df$se[bool_all] /
                                          b_df$se[bool_all]),
        frac_same_sign_all = mean(sign(a_df$lfc[bool_all]) ==
                                    sign(b_df$lfc[bool_all])),
        stringsAsFactors = FALSE)
    }
  }
}
lfc_df <- do.call(rbind, lfc_row_list)
utils::write.csv(lfc_df,
                 file = file.path(output_dir, "01_lfc-cross-method.csv"),
                 row.names = FALSE)
print("Spearman of logFC with DESeq2, within each cohort:")
print(lfc_df[lfc_df$method_b == "deseq2" & lfc_df$method_a != "deseq2",
             c("dataset", "method_a", "spearman_lfc", "spearman_se",
               "median_se_ratio")],
      digits = 3)

# 4. The empirical null within nuisance-status strata ---------------------------

# If the stored SE's miscalibration differs between capped and estimated
# genes, one global `null_sd` cannot be right for both. The diet objects
# carry the Gaussianized statistics and the status.
.load_esvd <- function(dataset){
  writeup <- writeup_vec[[dataset]]
  env <- new.env()
  load(file.path(tati_out_dir, writeup,
                 paste0(writeup, "_", dataset, "_esvd.RData")),
       envir = env)
  env[[paste0(dataset, "_esvd")]]$original
}

.fit_null <- function(z_vec){
  res <- tryCatch({
    fit <- locfdr::locfdr(z_vec, plot = 0)
    c(fit$fp0["mlest", "delta"], fit$fp0["mlest", "sigma"])
  }, warning = function(w){c(NA, NA)}, error = function(e){c(NA, NA)})
  # The package's own fallback, so a stratum locfdr refuses still gets a
  # (cruder) estimate.
  if(anyNA(res)){
    simple <- eSVD2:::.multtest_simple(z_vec, observed_quantile = c(0.05, 0.95))
    res <- c(simple$null_mean, simple$null_sd)
    method <- "simple"
  } else {
    method <- "locfdr"
  }
  list(method = method, null_mean = unname(res[1]), null_sd = unname(res[2]))
}

status_row_list <- list()
deg_row_list <- list()
for(dataset in dataset_vec){
  esvd_obj <- .load_esvd(dataset)
  z_vec <- esvd_obj$pvalue_list$gaussian_teststat
  status_vec <- as.character(esvd_obj[[esvd_obj$latest_Fit]]$nuisance_status)
  stopifnot(length(z_vec) == length(status_vec))
  stratum_list <- list(
    all = rep(TRUE, length(z_vec)),
    estimated = status_vec == "estimated",
    capped_or_boundary = status_vec %in% c("capped", "boundary")
  )
  for(stratum in names(stratum_list)){
    bool_vec <- stratum_list[[stratum]] & is.finite(z_vec)
    fit <- .fit_null(z_vec[bool_vec])
    status_row_list[[length(status_row_list) + 1]] <- data.frame(
      dataset = dataset,
      stratum = stratum,
      num_genes = sum(bool_vec),
      null_method = fit$method,
      null_mean = fit$null_mean,
      null_sd = fit$null_sd,
      sd_z = stats::sd(z_vec[bool_vec]),
      median_abs_z = stats::median(abs(z_vec[bool_vec])),
      median_se = stats::median(esvd_obj$log2fc_se_vec[bool_vec]),
      stringsAsFactors = FALSE)
  }
  # 5. Which genes eSVD-DE calls, by status, and whether their sign holds
  # in the other cohorts under DESeq2 (the method with the most DEGs there).
  esvd_df <- stat_list[[dataset]]$esvd
  stopifnot(identical(names(z_vec), esvd_df$gene))
  for(stratum in c("estimated", "capped_or_boundary")){
    bool_stratum <- stratum_list[[stratum]] & esvd_df$bool_converged
    bool_deg <- bool_stratum & esvd_df$bool_deg
    for(cohort_b in setdiff(dataset_vec, dataset)){
      b_df <- stat_list[[cohort_b]]$deseq2
      idx_vec <- which(bool_deg & b_df$bool_converged)
      deg_row_list[[length(deg_row_list) + 1]] <- data.frame(
        dataset = dataset,
        stratum = stratum,
        cohort_b = cohort_b,
        num_genes = sum(bool_stratum),
        num_deg = sum(bool_deg),
        frac_deg = sum(bool_deg) / sum(bool_stratum),
        median_abs_lfc_deg = stats::median(abs(esvd_df$lfc[bool_deg])),
        median_se_deg = stats::median(esvd_df$se[bool_deg]),
        frac_same_sign_deseq2_b = mean(sign(esvd_df$lfc[idx_vec]) ==
                                         sign(b_df$lfc[idx_vec])),
        stringsAsFactors = FALSE)
    }
  }
  rm(esvd_obj)
}
status_df <- do.call(rbind, status_row_list)
deg_df <- do.call(rbind, deg_row_list)
utils::write.csv(deg_df,
                 file = file.path(output_dir, "01_deg-by-status.csv"),
                 row.names = FALSE)
print("eSVD-DE DEGs by nuisance status, and their sign in B under DESeq2:")
print(deg_df, digits = 3)
utils::write.csv(status_df,
                 file = file.path(output_dir, "01_null-by-status.csv"),
                 row.names = FALSE)
print("Empirical null by nuisance status:")
print(status_df, digits = 3)

print("Done! :)")
