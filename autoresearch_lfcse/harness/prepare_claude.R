rm(list = ls())

# The bank builder of the lfcse loop
# Drafted by Claude for Kevin Z. Lin, 2026-11-11
#
# Builds the bank `lfcse-bank-1`: 144 cells x 2 splits, each example a pair
# of cohorts generated and fitted end to end by eSVD2 1.2.0, so that a
# scoring refits nothing. Design: notes/MEMO_claude.md, section 3.
#
# Run from the repository root:
#   VECLIB_MAXIMUM_THREADS=1 Rscript autoresearch_lfcse/harness/prepare_claude.R [--resume]
# Without `--resume` the bank folder is emptied first; with it, an example
# whose two files exist is kept and only its checksums are recomputed.
# Cores: `ESVD2_BUILD_CORES`, default 5 (Kevin's limit for the build,
# 2026-11-11). About 45 minutes on 5 cores.
#
# Writes, under the bank folder (`.bank_dir()`, outside the repository):
#   <split>/<id>_public.rds   the two slim fitted cohorts; nothing of the truth
#   <split>/<id>_truth.rds    the truth, the DESeq2 reference, the labels
#   manifest_<split>.csv      one row per example: files, checksums, factor
#                             labels, seed, realized values
#   BANK_INFO.md              generator version, date, design
# and logs to autoresearch_lfcse/output/prepare.log.

library(DESeq2)

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
loop_dir <- file.path(repo_dir, "autoresearch_lfcse")
stopifnot(dir.exists(file.path(loop_dir, "harness")))
source(file.path(loop_dir, "harness", "generator_claude.R"))
.load_esvd2(repo_dir)

bool_resume <- "--resume" %in% commandArgs(trailingOnly = TRUE)
num_cores <- as.integer(Sys.getenv("ESVD2_BUILD_CORES", unset = "5"))
stopifnot(!is.na(num_cores), num_cores >= 1)
bank_dir <- .bank_dir()
output_dir <- file.path(loop_dir, "output")
log_file <- file.path(output_dir, "prepare.log")
split_vec <- c("dev", "test")

.append_log <- function(text){
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", text, "\n"),
      file = log_file, append = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if(!bool_resume && dir.exists(bank_dir)){
  unlink(bank_dir, recursive = TRUE)
}
for(split in split_vec){
  dir.create(file.path(bank_dir, split), recursive = TRUE,
             showWarnings = FALSE)
}
.append_log(paste0("prepare started: ", .bank_version(), ", ", num_cores,
                   " cores, resume = ", bool_resume, ", bank at ", bank_dir))

# The tasks --------------------------------------------------------------------

design_df <- .bank_design()
task_df <- merge(data.frame(split = split_vec, stringsAsFactors = FALSE),
                 design_df, by = NULL)
task_df$id <- .example_id(cell_idx = task_df$cell_idx, split = task_df$split)
task_df$seed <- mapply(.example_seed, cell_idx = task_df$cell_idx,
                       split = task_df$split)
# The largest examples first, so that the workers finish together.
size_vec <- task_df$donors * ifelse(task_df$cells == "equal", 50, 190)
task_df <- task_df[order(-size_vec, task_df$id), ]
rownames(task_df) <- NULL
factor_vec <- c("donors", "cells", "rates", "effects", "covariates",
                "donor_sd")

.run_task <- function(task_idx){
  task <- task_df[task_idx, , drop = FALSE]
  public_file <- file.path(task$split, paste0(task$id, "_public.rds"))
  truth_file <- file.path(task$split, paste0(task$id, "_truth.rds"))
  row_df <- task[, c("id", "split", "cell_idx", factor_vec, "seed")]
  row_df$file_public <- public_file
  row_df$file_truth <- truth_file

  # A file is written under a temporary name and renamed, so a file under
  # its final name is complete; on resume both must still read and carry
  # the task's seed, otherwise the example is rebuilt.
  bool_have <- FALSE
  if(bool_resume &&
     all(file.exists(file.path(bank_dir, c(public_file, truth_file))))){
    truth <- tryCatch(readRDS(file.path(bank_dir, truth_file)),
                      error = function(e){NULL})
    public <- tryCatch(readRDS(file.path(bank_dir, public_file)),
                       error = function(e){NULL})
    bool_have <- !is.null(truth) && !is.null(public) &&
      identical(truth$meta$seed, task$seed)
    rm(public)
  }
  if(bool_have){
    info_df <- truth$info_df
  } else {
    time_start <- proc.time()[["elapsed"]]
    example <- tryCatch(.make_example(design_row = task, split = task$split),
                        error = function(e){
                          list(error_message = conditionMessage(e))
                        })
    if(!is.null(example$error_message)){
      .append_log(paste0(task$id, " FAILED: ", example$error_message))
      row_df$status <- "failed"
      row_df$error <- example$error_message
      return(row_df)
    }
    info_df <- example$info_df
    truth <- c(example$truth,
               list(info_df = info_df,
                    meta = as.list(task[, c("id", "split", "cell_idx",
                                            factor_vec, "seed")])))
    for(pair in list(list(example$public, public_file),
                     list(truth, truth_file))){
      tmp_file <- file.path(bank_dir, paste0(pair[[2]], ".tmp"))
      saveRDS(pair[[1]], file = tmp_file)
      file.rename(tmp_file, file.path(bank_dir, pair[[2]]))
    }
    .append_log(paste0(task$id, " done in ",
                       round(proc.time()[["elapsed"]] - time_start), " s (",
                       task_idx, " of ", nrow(task_df), ")"))
  }
  row_df$md5_public <- unname(tools::md5sum(file.path(bank_dir, public_file)))
  row_df$md5_truth <- unname(tools::md5sum(file.path(bank_dir, truth_file)))
  row_df$status <- "ok"
  row_df$error <- ""
  cbind(row_df, info_df)
}

# Building ---------------------------------------------------------------------

time_start <- proc.time()[["elapsed"]]
row_list <- parallel::mclapply(seq_len(nrow(task_df)), .run_task,
                               mc.cores = num_cores, mc.preschedule = FALSE)
build_minutes <- (proc.time()[["elapsed"]] - time_start) / 60

# A worker that died returns a try-error or NULL instead of a row.
bool_row <- sapply(row_list, is.data.frame)
if(!all(bool_row)){
  .append_log(paste0(sum(!bool_row), " tasks returned no row: ",
                     paste(task_df$id[!bool_row], collapse = ", ")))
  stop(sum(!bool_row), " tasks returned no row; see ", log_file,
       " and rerun with --resume")
}
column_vec <- unique(unlist(lapply(row_list, colnames)))
row_list <- lapply(row_list, function(row_df){
  row_df[setdiff(column_vec, colnames(row_df))] <- NA
  row_df[, column_vec, drop = FALSE]
})
manifest_df <- do.call(rbind, row_list)
manifest_df <- manifest_df[order(manifest_df$split, manifest_df$cell_idx), ]

for(split in split_vec){
  utils::write.csv(manifest_df[manifest_df$split == split, ],
                   file = file.path(bank_dir,
                                    paste0("manifest_", split, ".csv")),
                   row.names = FALSE)
}

# BANK_INFO.md -----------------------------------------------------------------

num_ok <- sum(manifest_df$status == "ok")
info_vec <- c(
  paste0("# ", .bank_version()),
  "",
  paste0("Built ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by ",
         "`autoresearch_lfcse/harness/prepare_claude.R` with eSVD2 ",
         as.character(utils::packageVersion("eSVD2")), " and R ",
         as.character(getRversion()), ", on ", num_cores, " cores in ",
         round(build_minutes), " minutes."),
  "",
  paste0("Generator version `", .bank_version(), "` ",
         "(`harness/generator_claude.R`). One example is a pair of cohorts ",
         "(A, B) of 600 genes sharing intercepts, loadings, rates and ",
         "covariate effects, with case-control effects correlated at 0.8 ",
         "and a size ratio of 2, each fitted end to end at rank 3 with the ",
         "nuisance cap at 10. Factors: donors (10, 20, 40), cells per donor ",
         "(50 equal; unequal log-normal, median 150, 30 to 600), rates ",
         "(2 to 8; 0.5 to 200), effects (sparse, sparse_up, dense), ",
         "covariates (none; donor sex and age), donor SD (",
         paste(.donor_sd_levels(), collapse = ", "), "): 144 cells, one ",
         "example per cell per split."),
  "",
  paste0("Examples: ", num_ok, " of ", nrow(manifest_df), " built (",
         sum(manifest_df$split == "dev" & manifest_df$status == "ok"),
         " dev, ",
         sum(manifest_df$split == "test" & manifest_df$status == "ok"),
         " test)."),
  "",
  paste0("The public files hold no cells x genes matrix of doubles; the ",
         "scorer rebuilds the library and the posterior from the stored ",
         "fit with `.build_view()`.")
)
writeLines(info_vec, con = file.path(bank_dir, "BANK_INFO.md"))

.append_log(paste0("prepare finished: ", num_ok, " of ", nrow(manifest_df),
                   " examples in ", round(build_minutes), " minutes"))
print(table(manifest_df$split, manifest_df$status))
view_diff <- max(manifest_df$view_diff_a, manifest_df$view_diff_b,
                 na.rm = TRUE)
print(paste0("Largest gap between a rebuilt view and its fit's posterior: ",
             view_diff))
# The scorer's view must be the fit's own posterior.
stopifnot(view_diff < 1e-8)
print("Done! :)")
