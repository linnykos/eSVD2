# Beyond-score study, step 1: build the bank `lfcse-bank-3`
# Drafted by Claude for Kevin Z. Lin, 2026-11-12

rm(list = ls())

# Builds 216 examples, each one population and two replicate cohorts
# generated and fitted end to end by eSVD2 1.2.0 (432 fits).
#
# Run from the repository root:
#   VECLIB_MAXIMUM_THREADS=1 Rscript additional_context/2026-11-12_beyond-score/01_build-bank_claude.R [--resume]
# Without `--resume` the bank folder is emptied first; with it, an example
# whose files exist and carry the task's seed is kept.
# Cores: `ESVD2_BUILD_CORES`, default 5 (Kevin's limit, 2026-11-12).
#
# Writes, under the bank folder (`.bank_dir()`, outside the repository):
#   examples/<id>_<a|b>_public.rds   the slim fitted cohort
#   examples/<id>_truth.rds          the population, both cohorts' truth
#   manifest.csv                     one row per example
# and logs to output/build.log beside this script.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context",
                       "2026-11-12_beyond-score")
stopifnot(dir.exists(study_dir))
source(file.path(repo_dir, "autoresearch_lfcse2", "harness",
                 "generator_claude.R"))
source(file.path(study_dir, "generator_claude.R"))
.load_esvd2(repo_dir)

bool_resume <- "--resume" %in% commandArgs(trailingOnly = TRUE)
num_cores <- as.integer(Sys.getenv("ESVD2_BUILD_CORES", unset = "5"))
stopifnot(!is.na(num_cores), num_cores >= 1)
bank_dir <- .bank_dir()
output_dir <- file.path(study_dir, "output")
log_file <- file.path(output_dir, "build.log")

.append_log <- function(text){
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", text, "\n"),
      file = log_file, append = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if(!bool_resume && dir.exists(bank_dir)) unlink(bank_dir, recursive = TRUE)
dir.create(file.path(bank_dir, "examples"), recursive = TRUE,
           showWarnings = FALSE)
.append_log(paste0("build started: ", .bank_version(), ", ", num_cores,
                   " cores, resume = ", bool_resume, ", bank at ", bank_dir))

# The tasks --------------------------------------------------------------------

task_df <- .bank_design()
task_df$id <- .example_id(cell_idx = task_df$cell_idx)
task_df$seed <- .example_seed(cell_idx = task_df$cell_idx)
# The largest examples first, so that the workers finish together.
size_vec <- task_df$donors * ifelse(task_df$cells == "equal", 50, 190)
task_df <- task_df[order(-size_vec, task_df$id), ]
rownames(task_df) <- NULL
factor_vec <- .factor_names()
replicate_vec <- .replicate_names()

.run_task <- function(task_idx){
  task <- task_df[task_idx, , drop = FALSE]
  public_file_vec <- file.path("examples", paste0(task$id, "_", replicate_vec,
                                                  "_public.rds"))
  truth_file <- file.path("examples", paste0(task$id, "_truth.rds"))
  row_df <- task[, c("id", "cell_idx", factor_vec, "seed")]
  row_df$file_truth <- truth_file

  bool_have <- FALSE
  if(bool_resume &&
     all(file.exists(file.path(bank_dir, c(public_file_vec, truth_file))))){
    truth <- tryCatch(readRDS(file.path(bank_dir, truth_file)),
                      error = function(e){NULL})
    bool_have <- !is.null(truth) && identical(truth$meta$seed, task$seed)
  }
  if(bool_have){
    info_df <- truth$info_df
  } else {
    time_start <- proc.time()[["elapsed"]]
    example <- tryCatch(.make_example(design_row = task),
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
                    meta = as.list(task[, c("id", "cell_idx", factor_vec,
                                            "seed")])))
    # Written under a temporary name and renamed, so that a file under its
    # final name is complete. The truth goes last: its presence marks the
    # example as built.
    file_vec <- c(public_file_vec, truth_file)
    object_list <- c(unname(example$public_list), list(truth))
    for(i in seq_along(file_vec)){
      tmp_file <- file.path(bank_dir, paste0(file_vec[i], ".tmp"))
      saveRDS(object_list[[i]], file = tmp_file)
      file.rename(tmp_file, file.path(bank_dir, file_vec[i]))
    }
    .append_log(paste0(task$id, " done in ",
                       round(proc.time()[["elapsed"]] - time_start), " s (",
                       task_idx, " of ", nrow(task_df), ")"))
  }
  row_df$status <- "ok"
  row_df$error <- ""
  cbind(row_df, info_df)
}

# Building ---------------------------------------------------------------------

time_start <- proc.time()[["elapsed"]]
row_list <- parallel::mclapply(seq_len(nrow(task_df)), .run_task,
                               mc.cores = num_cores, mc.preschedule = FALSE)
build_minutes <- (proc.time()[["elapsed"]] - time_start) / 60

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
manifest_df <- manifest_df[order(manifest_df$cell_idx), ]
utils::write.csv(manifest_df, file = file.path(bank_dir, "manifest.csv"),
                 row.names = FALSE)

.append_log(paste0("build finished: ", sum(manifest_df$status == "ok"),
                   " of ", nrow(manifest_df), " examples in ",
                   round(build_minutes), " minutes"))
print(table(manifest_df$status))
view_diff <- max(manifest_df$view_diff_a, manifest_df$view_diff_b,
                 na.rm = TRUE)
print(paste0("Largest gap between a rebuilt view and its fit's posterior: ",
             view_diff))
stopifnot(view_diff < 1e-8)
print("Done! :)")
