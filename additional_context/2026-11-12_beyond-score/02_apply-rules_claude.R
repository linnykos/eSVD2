# Beyond-score study, step 2: every rule on every cohort, gene by gene
# Drafted by Claude for Kevin Z. Lin, 2026-11-12

rm(list = ls())

library(DESeq2)
library(edgeR)

# Applies every reporting rule the two autoresearch loops considered, the
# library-size variants of the rules that read a library, and pseudobulk
# DESeq2, to every cohort of one bank, and stores the per-gene log2 fold
# change and SE of each. Nothing is scored here; step 3 reads these files.
#
# Run from the repository root:
#   VECLIB_MAXIMUM_THREADS=1 Rscript additional_context/2026-11-12_beyond-score/02_apply-rules_claude.R <bank3|bank2|bank1> [--resume]
# Cores: `ESVD2_CORES`, default 5.
#
# Writes <ESVD2_BANK_DIR>/lfcse-beyond-score/<bank>/<cohort>.rds, outside
# the repository, and logs to output/apply_<bank>.log beside this script.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context",
                       "2026-11-12_beyond-score")
stopifnot(dir.exists(study_dir))
loop1_dir <- file.path(repo_dir, "autoresearch_lfcse")
loop2_dir <- file.path(repo_dir, "autoresearch_lfcse2")
source(file.path(loop2_dir, "harness", "generator_claude.R"))
.load_esvd2(repo_dir)

arg_vec <- commandArgs(trailingOnly = TRUE)
bank <- arg_vec[1]
stopifnot(bank %in% c("bank1", "bank2", "bank3"))
bool_resume <- "--resume" %in% arg_vec
num_cores <- as.integer(Sys.getenv("ESVD2_CORES", unset = "5"))

root_dir <- Sys.getenv("ESVD2_BANK_DIR", unset = "")
if(root_dir == "") root_dir <- tools::R_user_dir("eSVD2", which = "cache")
bank_dir <- file.path(root_dir, c(bank1 = "lfcse-bank-1",
                                  bank2 = "lfcse-bank-2",
                                  bank3 = "lfcse-bank-3")[[bank]])
eval_dir <- file.path(root_dir, "lfcse-beyond-score", bank)
dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- file.path(study_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(output_dir, paste0("apply_", bank, ".log"))

.append_log <- function(text){
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", text, "\n"),
      file = log_file, append = TRUE)
}

# The rules --------------------------------------------------------------------

# One file per rule; each defines `report_lfc_se(view)`. The first five and
# `gap_trend` are the reference rows of both loops, `donor_wls` the first
# loop's kept rule, `embed_prior` and `state_adjusted` the second loop's.
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

# Which library variants each rule is run with, beside the fitted library.
# The posterior-based rules (`stored`, `null_rescaled`, `deshrunk`) read no
# library that can be exchanged.
variant_list <- list(
  donor_wls = c("esvd_nofloor", "esvd_depth", "total", "total_cov", "sf_cov",
                "tmm_cov", "oracle"),
  gap_trend = c("esvd_nofloor", "total_cov"),
  raw_null = c("esvd_nofloor", "total_cov"),
  embed_prior = c("esvd_nofloor", "total_cov"),
  state_adjusted = c("esvd_nofloor", "total_cov")
)

# Library variants -------------------------------------------------------------

# Columns of the design that vary within a donor, the depth aside.
.cell_level_columns <- function(covariates, individual_vec){
  individual_chr <- as.character(individual_vec)
  bool_vec <- apply(covariates, 2, function(col_vec){
    any(tapply(col_vec, individual_chr, function(x){diff(range(x))}) > 0)
  })
  setdiff(colnames(covariates)[bool_vec], c("Intercept", "Log_UMI", "CC_1"))
}

# The design with each cell-level covariate replaced by its donor mean, so
# that a rule's donor-level regression adjusts for it, as a pseudobulk
# analysis would.
.donor_mean_covariates <- function(covariates, individual_vec){
  individual_chr <- as.character(individual_vec)
  for(name in .cell_level_columns(covariates, individual_vec)){
    covariates[, name] <- stats::ave(covariates[, name], individual_chr)
  }
  covariates
}

# A cells x genes library from a per-cell vector.
.expand_library <- function(library_vec, dat){
  library_mat <- matrix(library_vec / mean(library_vec), nrow = nrow(dat),
                        ncol = ncol(dat))
  dimnames(library_mat) <- dimnames(dat)
  library_mat
}

# Every library the variants need, each cells x genes.
#   esvd_nofloor  exp of the fit's linear predictor over every covariate but
#                 the arm, without the posterior's floor at 0.1
#   esvd_depth    the same from the intercept and Log_UMI alone
#   total         the cell's summed counts
#   sf            the cell's summed counts, rescaled within each donor so
#                 that the donor's sum is DESeq2's median-of-ratios size
#                 factor of the pseudobulk
#   tmm           the cell's summed counts times edgeR's TMM factor of its
#                 donor's pseudobulk
#   oracle        the generator's library times the true effect of the
#                 cell-level covariates (third bank only)
.make_libraries <- function(view, cohort_truth, population){
  dat <- view$dat
  individual_chr <- as.character(view$individual_vec)
  library_name_vec <- setdiff(colnames(view$covariates), "CC_1")
  esvd_nofloor_mat <- exp(tcrossprod(view$covariates[, library_name_vec],
                                     view$z_mat[, library_name_vec]))
  esvd_depth_mat <- exp(tcrossprod(view$covariates[, c("Intercept", "Log_UMI")],
                                   view$z_mat[, c("Intercept", "Log_UMI")]))
  dimnames(esvd_nofloor_mat) <- dimnames(dat)
  dimnames(esvd_depth_mat) <- dimnames(dat)

  total_vec <- rowSums(dat)
  pseudobulk_mat <- t(rowsum(dat, group = individual_chr))
  donor_total_vec <- colSums(pseudobulk_mat)
  size_factor_vec <- tryCatch(
    DESeq2::estimateSizeFactorsForMatrix(pseudobulk_mat),
    error = function(e){
      DESeq2::estimateSizeFactorsForMatrix(pseudobulk_mat, type = "poscounts")
    }
  )
  tmm_vec <- edgeR::calcNormFactors(pseudobulk_mat, method = "TMM")
  names(tmm_vec) <- colnames(pseudobulk_mat)
  sf_cell_vec <- total_vec * size_factor_vec[individual_chr] /
    donor_total_vec[individual_chr]
  tmm_cell_vec <- total_vec * tmm_vec[individual_chr]

  library_list <- list(esvd_nofloor = esvd_nofloor_mat,
                       esvd_depth = esvd_depth_mat,
                       total = .expand_library(total_vec, dat),
                       sf = .expand_library(sf_cell_vec, dat),
                       tmm = .expand_library(tmm_cell_vec, dat))
  if(!is.null(cohort_truth$library_true_vec)){
    oracle_mat <- cohort_truth$library_true_vec *
      exp(outer(cohort_truth$pct_vec, as.numeric(population$pct_effect_vec)) +
            outer(cohort_truth$phase_vec,
                  as.numeric(population$phase_effect_vec)))
    dimnames(oracle_mat) <- dimnames(dat)
    library_list$oracle <- oracle_mat
  }
  library_list
}

# The view a variant hands to a rule.
.variant_view <- function(view, variant, library_list, covariates_donor){
  if(variant %in% c("esvd_nofloor", "esvd_depth", "total", "oracle")){
    view$library_mat <- library_list[[variant]]
    return(view)
  }
  view$library_mat <- library_list[[sub("_cov$", "", variant)]]
  view$covariates <- covariates_donor
  view
}

# DESeq2 -----------------------------------------------------------------------

# DESeq2 on pseudobulk counts, one column per donor. The design holds every
# covariate of the fit but the depth, cell-level ones as their donor means,
# and the arm last. Default `DESeq()` and `results()`.
.run_deseq2 <- function(view, covariates_donor){
  individual_chr <- as.character(view$individual_vec)
  pseudobulk_mat <- t(rowsum(view$dat, group = individual_chr))
  first_idx <- match(colnames(pseudobulk_mat), individual_chr)
  name_vec <- setdiff(colnames(covariates_donor),
                      c("Intercept", "Log_UMI", "CC_1"))
  col_df <- as.data.frame(covariates_donor[first_idx, c(name_vec, "CC_1"),
                                           drop = FALSE])
  rownames(col_df) <- colnames(pseudobulk_mat)
  design_formula <- stats::as.formula(paste("~", paste(c(name_vec, "CC_1"),
                                                       collapse = " + ")))
  storage.mode(pseudobulk_mat) <- "integer"
  gene_vec <- colnames(view$dat)
  res <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_mat,
                                          colData = col_df,
                                          design = design_formula)
    dds <- suppressMessages(suppressWarnings(DESeq2::DESeq(dds,
                                                           quiet = TRUE)))
    DESeq2::results(dds, name = "CC_1")
  }, error = function(e){e})
  if(inherits(res, "error")){
    return(list(df = data.frame(gene = gene_vec, log2fc = NA_real_,
                                se = NA_real_, pvalue = NA_real_,
                                padj = NA_real_, stringsAsFactors = FALSE),
                error = conditionMessage(res)))
  }
  idx_vec <- match(gene_vec, rownames(res))
  list(df = data.frame(gene = gene_vec,
                       log2fc = res$log2FoldChange[idx_vec],
                       se = res$lfcSE[idx_vec],
                       pvalue = res$pvalue[idx_vec],
                       padj = res$padj[idx_vec],
                       stringsAsFactors = FALSE),
       error = "")
}

# The cohorts of a bank --------------------------------------------------------

# One row per cohort: where its public cohort and its truth are, and its
# factor labels. Banks 1 and 2 contribute their test split only, the split
# no rule was tuned on.
.cohort_table <- function(bank){
  if(bank == "bank3"){
    manifest_df <- utils::read.csv(file.path(bank_dir, "manifest.csv"),
                                   stringsAsFactors = FALSE)
    manifest_df <- manifest_df[manifest_df$status == "ok", ]
    cohort_df <- do.call(rbind, lapply(c("a", "b"), function(replicate){
      data.frame(cohort = paste0(manifest_df$id, "_", replicate),
                 id = manifest_df$id,
                 replicate = replicate,
                 file_public = file.path("examples",
                                         paste0(manifest_df$id, "_", replicate,
                                                "_public.rds")),
                 file_truth = manifest_df$file_truth,
                 num_cells = manifest_df[[paste0("num_cells_", replicate)]],
                 stringsAsFactors = FALSE)
    }))
    return(cohort_df)
  }
  manifest_df <- utils::read.csv(file.path(bank_dir, "manifest_test.csv"),
                                 stringsAsFactors = FALSE)
  manifest_df <- manifest_df[manifest_df$status == "ok", ]
  if(bank == "bank2"){
    return(data.frame(cohort = manifest_df$id,
                      id = manifest_df$id,
                      replicate = "a",
                      file_public = manifest_df$file_public,
                      file_truth = manifest_df$file_truth,
                      num_cells = manifest_df$num_cells,
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(c("a", "b"), function(replicate){
    data.frame(cohort = paste0(manifest_df$id, "_", replicate),
               id = manifest_df$id,
               replicate = replicate,
               file_public = manifest_df$file_public,
               file_truth = manifest_df$file_truth,
               num_cells = manifest_df[[paste0("num_cells_", replicate)]],
               stringsAsFactors = FALSE)
  }))
}

# The slim cohort, both truths, and the generator's pieces of one cohort.
.load_cohort <- function(bank, cohort_row){
  public <- readRDS(file.path(bank_dir, cohort_row$file_public))
  truth <- readRDS(file.path(bank_dir, cohort_row$file_truth))
  replicate <- cohort_row$replicate
  if(bank == "bank3"){
    return(list(slim = public,
                truth_marginal_vec = as.numeric(truth$truth_marginal_vec),
                truth_direct_vec = as.numeric(truth$truth_direct_vec),
                cohort_truth = truth$cohort_list[[replicate]],
                population = truth$population,
                meta = truth$meta))
  }
  if(bank == "bank2"){
    return(list(slim = public,
                truth_marginal_vec = as.numeric(truth$truth_log2fc_vec),
                truth_direct_vec = as.numeric(truth$population$cc_effect_vec) /
                  log(2),
                cohort_truth = list(),
                population = truth$population,
                meta = truth$meta))
  }
  # The first bank's two cohorts of a pair have different true effects and
  # no cell states, so the two truths coincide.
  truth_vec <- as.numeric(truth[[replicate]]$truth_log2fc_vec)
  list(slim = public[[replicate]],
       truth_marginal_vec = truth_vec,
       truth_direct_vec = truth_vec,
       cohort_truth = list(),
       population = truth$population,
       meta = truth$meta)
}

# One cohort -------------------------------------------------------------------

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

.run_cohort <- function(cohort_idx){
  cohort_row <- cohort_df[cohort_idx, , drop = FALSE]
  out_file <- file.path(eval_dir, paste0(cohort_row$cohort, ".rds"))
  if(bool_resume && file.exists(out_file)) return("kept")
  time_start <- proc.time()[["elapsed"]]
  loaded <- .load_cohort(bank = bank, cohort_row = cohort_row)
  view <- .build_view(loaded$slim)
  gene_vec <- colnames(view$dat)
  library_list <- .make_libraries(view = view,
                                  cohort_truth = loaded$cohort_truth,
                                  population = loaded$population)
  covariates_donor <- .donor_mean_covariates(view$covariates,
                                             view$individual_vec)
  cell_level_vec <- .cell_level_columns(view$covariates, view$individual_vec)

  result_list <- list()
  for(rule in names(rule_list)){
    # No rule draws random numbers; the seed makes that harmless if one
    # ever does.
    set.seed(10)
    result_list[[rule]] <- .run_rule(rule_list[[rule]], view)
    for(variant in variant_list[[rule]]){
      if(variant == "oracle" && is.null(library_list$oracle)) next
      set.seed(10)
      result_list[[paste0(rule, "@", variant)]] <- .run_rule(
        rule_list[[rule]],
        .variant_view(view = view, variant = variant,
                      library_list = library_list,
                      covariates_donor = covariates_donor)
      )
    }
  }
  deseq2_res <- .run_deseq2(view = view, covariates_donor = covariates_donor)
  result_list[["deseq2"]] <- list(log2fc = deseq2_res$df$log2fc,
                                  se = deseq2_res$df$se,
                                  error = deseq2_res$error,
                                  seconds = NA_real_)
  # DESeq2 without the cell-level covariates in its design: what a
  # pseudobulk analysis gives when it leaves them out. The same fit where
  # the cohort has none.
  deseq2_nocell_res <- deseq2_res
  if(length(cell_level_vec) > 0){
    keep_vec <- setdiff(colnames(covariates_donor), cell_level_vec)
    deseq2_nocell_res <- .run_deseq2(
      view = view,
      covariates_donor = covariates_donor[, keep_vec, drop = FALSE]
    )
  }
  result_list[["deseq2_nocell"]] <- list(log2fc = deseq2_nocell_res$df$log2fc,
                                         se = deseq2_nocell_res$df$se,
                                         error = deseq2_nocell_res$error,
                                         seconds = NA_real_)

  method_vec <- names(result_list)
  lfc_mat <- sapply(result_list, function(res){res$log2fc})
  se_mat <- sapply(result_list, function(res){res$se})
  dimnames(lfc_mat) <- list(gene_vec, method_vec)
  dimnames(se_mat) <- list(gene_vec, method_vec)

  individual_chr <- as.character(view$individual_vec)
  gene_df <- data.frame(
    gene = gene_vec,
    truth_marginal = loaded$truth_marginal_vec,
    truth_direct = loaded$truth_direct_vec,
    zero_fraction = colMeans(view$dat == 0),
    mean_count = colMeans(view$dat),
    # Share of the gene's cells whose fitted library sits at the
    # posterior's floor of 0.1.
    floor_fraction = colMeans(library_list$esvd_nofloor < 0.1),
    status = as.character(view$nuisance_status_vec),
    esvd_pvalue = view$stored_df$pvalue,
    esvd_fdr = if(is.null(view$stored_df$fdr)){
      stats::p.adjust(view$stored_df$pvalue, method = "BH")
    } else {
      view$stored_df$fdr
    },
    deseq2_pvalue = deseq2_res$df$pvalue,
    deseq2_padj = deseq2_res$df$padj,
    deseq2_nocell_padj = deseq2_nocell_res$df$padj,
    stringsAsFactors = FALSE
  )
  out <- list(cohort = cohort_row$cohort,
              id = cohort_row$id,
              replicate = cohort_row$replicate,
              bank = bank,
              meta = loaded$meta,
              gene_df = gene_df,
              lfc_mat = lfc_mat,
              se_mat = se_mat,
              error_vec = sapply(result_list, function(res){res$error}),
              seconds_vec = sapply(result_list, function(res){res$seconds}),
              cell_level_vec = cell_level_vec,
              num_cells = nrow(view$dat),
              num_donors = length(unique(individual_chr)),
              null_mean = attr(view$stored_df, "null_mean"),
              null_sd = attr(view$stored_df, "null_sd"))
  tmp_file <- paste0(out_file, ".tmp")
  saveRDS(out, file = tmp_file)
  file.rename(tmp_file, out_file)
  .append_log(paste0(cohort_row$cohort, " done in ",
                     round(proc.time()[["elapsed"]] - time_start), " s (",
                     cohort_idx, " of ", nrow(cohort_df), ")"))
  "done"
}

# Running ----------------------------------------------------------------------

cohort_df <- .cohort_table(bank)
# The largest cohorts first, so that the workers finish together.
cohort_df <- cohort_df[order(-cohort_df$num_cells, cohort_df$cohort), ]
rownames(cohort_df) <- NULL
only_val <- Sys.getenv("ESVD2_ONLY", unset = "")
if(only_val != "") cohort_df <- cohort_df[cohort_df$cohort %in% only_val, ]
.append_log(paste0("apply started: ", bank, ", ", nrow(cohort_df),
                   " cohorts, ", num_cores, " cores, resume = ", bool_resume))

status_list <- parallel::mclapply(seq_len(nrow(cohort_df)), function(idx){
  tryCatch(.run_cohort(idx), error = function(e){
    .append_log(paste0(cohort_df$cohort[idx], " FAILED: ", conditionMessage(e)))
    paste0("failed: ", conditionMessage(e))
  })
}, mc.cores = num_cores, mc.preschedule = FALSE)
status_vec <- sapply(status_list, function(x){
  if(is.character(x)) x[1] else "worker died"
})
print(table(sub(":.*", "", status_vec)))
.append_log(paste0("apply finished: ", sum(status_vec %in% c("done", "kept")),
                   " of ", nrow(cohort_df), " cohorts"))
stopifnot(all(status_vec %in% c("done", "kept")))
print("Done! :)")
