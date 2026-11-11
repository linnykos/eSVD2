rm(list = ls())

# The scorer of the lfcse loop
# Drafted by Claude for Kevin Z. Lin, 2026-11-11
#
# FIXED during a run: this file is the same judge for every commit. It
# scores `report_lfc_se()` of the trial file on every example of one split
# of the bank `lfcse-bank-1` and prints the summary block that
# harness/record_result.sh reads. Contract: notes/MEMO_claude.md, sections
# 1, 2 and 5.
#
# Run from the repository root:
#   Rscript autoresearch_lfcse/harness/evaluate_claude.R [dev|test] > autoresearch_lfcse/output/run.log 2>&1
#
# Environment:
#   ESVD2_BANK_DIR   root of the bank (default: R's user cache for eSVD2)
#   ESVD2_CORES      examples scored in parallel (default 8, the memo's)
#   ESVD2_BUDGET     wall-clock budget of a scoring, seconds (default 600)
#   ESVD2_TRIAL      another file to score in place of the trial file, for
#                    the reference rows of BANK.md; its table goes to
#                    output/diagnostics/ and no run file is touched
#
# Writes output/last_scores.csv (dev) or output/test_scores.csv (test), one
# row per example, and on test appends a row to output/test_runs.tsv.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
loop_dir <- file.path(repo_dir, "autoresearch_lfcse")
stopifnot(dir.exists(file.path(loop_dir, "harness")))
source(file.path(loop_dir, "harness", "generator_claude.R"))

arg_vec <- commandArgs(trailingOnly = TRUE)
split <- if(length(arg_vec) >= 1) arg_vec[1] else "dev"
stopifnot(split %in% c("dev", "test"))
num_cores <- as.integer(Sys.getenv("ESVD2_CORES", unset = "8"))
budget_seconds <- as.numeric(Sys.getenv("ESVD2_BUDGET", unset = "600"))
# Seconds one cohort's call of the trial function may take (memo, section 2).
cohort_seconds <- 30
# The validity gate: the pooled slope of the reported logFC on the truth
# over weak-effect genes must lie in this interval (memo, section 5).
gate_range <- c(0.9, 1.1)

default_trial_file <- file.path(loop_dir, "report_rule_claude.R")
trial_file <- Sys.getenv("ESVD2_TRIAL", unset = "")
bool_reference <- trial_file != ""
if(!bool_reference) trial_file <- default_trial_file
output_dir <- file.path(loop_dir, "output")
dir.create(file.path(output_dir, "diagnostics"), recursive = TRUE,
           showWarnings = FALSE)
bank_dir <- .bank_dir()
time_start <- proc.time()[["elapsed"]]

commit_val <- tryCatch(system2("git", c("rev-parse", "--short=7", "HEAD"),
                               stdout = TRUE, stderr = FALSE),
                       error = function(e){"unknown"})
trial_md5 <- if(file.exists(trial_file)){
  unname(tools::md5sum(trial_file))
} else {
  "missing"
}

# A failed scoring prints the short block and stops.
.fail <- function(status, error_vec){
  cat("---\n")
  cat("status:         ", status, "\n", sep = "")
  for(error in utils::head(unique(error_vec), 5)){
    cat("error:          ", gsub("\n", " ", error), "\n", sep = "")
  }
  cat("seconds:        ",
      round(proc.time()[["elapsed"]] - time_start, 1), "\n", sep = "")
  cat("commit:         ", commit_val, "\n", sep = "")
  cat("trial_md5:      ", trial_md5, "\n", sep = "")
  quit(save = "no", status = 0)
}

# The bank ---------------------------------------------------------------------

manifest_file <- file.path(bank_dir, paste0("manifest_", split, ".csv"))
if(!file.exists(manifest_file)){
  .fail("crash", paste0("no manifest at ", manifest_file,
                        "; run harness/prepare_claude.R"))
}
manifest_df <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
if(any(manifest_df$status != "ok")){
  .fail("crash", paste0(sum(manifest_df$status != "ok"),
                        " examples of the manifest were not built"))
}
md5_public_vec <- unname(tools::md5sum(file.path(bank_dir,
                                                 manifest_df$file_public)))
md5_truth_vec <- unname(tools::md5sum(file.path(bank_dir,
                                                manifest_df$file_truth)))
bool_changed <- is.na(md5_public_vec) | is.na(md5_truth_vec) |
  md5_public_vec != manifest_df$md5_public |
  md5_truth_vec != manifest_df$md5_truth
if(any(bool_changed)){
  .fail("crash", paste0("the bank differs from its manifest (",
                        sum(bool_changed), " examples, first ",
                        manifest_df$id[which(bool_changed)[1]],
                        "); rebuild it and re-baseline"))
}

# The tripwire -----------------------------------------------------------------

# The trial function sees the public view and nothing else. The file is
# parsed to tokens (comments excluded) and refused when a symbol or a
# string constant names a way to read or write files, to evaluate text, to
# reach another environment or the call stack, or anything of the bank,
# the truth, the harness or its outputs. A tripwire, not a proof.
.forbidden_symbol_vec <- function(){
  c("readRDS", "saveRDS", "load", "save", "save.image", "read.csv",
    "read.table", "read.delim", "readLines", "writeLines", "readline",
    "scan", "file", "url", "gzfile", "bzfile", "xzfile", "pipe",
    "textConnection", "rawConnection", "socketConnection", "unz",
    "write", "write.csv", "write.table", "sink", "dput", "dget", "dump",
    "unserialize", "serialize", "readBin", "writeBin", "readChar",
    "writeChar", "list.files", "list.dirs", "dir", "file.exists",
    "file.path", "normalizePath", "path.expand", "getwd", "setwd",
    "tempfile", "tempdir", "unlink", "file.remove", "file.copy",
    "source", "sys.source", "eval", "evalq", "eval.parent", "parse",
    "str2lang", "str2expression", "body", "formals", "args",
    "Sys.getenv", "Sys.setenv", "Sys.unsetenv", "Sys.glob", "Sys.time",
    "sys.call", "sys.calls", "sys.function", "sys.frame", "sys.frames",
    "sys.parent", "sys.parents", "sys.on.exit", "sys.status",
    "parent.frame", "parent.env", "environment", "globalenv", "emptyenv",
    "topenv", "baseenv", "as.environment", "new.env", "local", "attach",
    "with", "within", "pos.to.env", "search", "searchpaths",
    "get", "get0", "mget", "assign", "exists", "rm", "delayedAssign",
    "makeActiveBinding", "lockBinding", "unlockBinding", "match.fun",
    "do.call", "Recall", "match.arg", "match.call",
    "getFromNamespace", "asNamespace", "getNamespace", "loadNamespace",
    "requireNamespace", "library", "require", "attachNamespace",
    "getExportedValue", "getAnywhere", "getS3method", "getOption",
    "options", "system", "system2", "shell", "Sys.which", "commandArgs",
    "quit", "q", "on.exit", "reg.finalizer", "trace", "untrace",
    "setTimeLimit", "setSessionTimeLimit", "Sys.sleep", "set.seed",
    "withCallingHandlers", "signalCondition",
    "mclapply", "mcparallel", "mccollect", "makeCluster", "parLapply",
    "dyn.load", ".Call", ".External", ".C", ".Fortran", ".Internal",
    ".Primitive", ".GlobalEnv", ".BaseNamespaceEnv", "Rscript")
}

# Substrings refused inside any symbol or string constant.
.forbidden_word_vec <- function(){
  c("truth", "manifest", "bank", "harness", "evaluate", "prepare",
    "generator", "deseq", "edger", "limma", "nebula", "dreamlet", "apeglm",
    "ashr", "best_scores", "last_scores", "test_scores", "results.tsv",
    "run.log", "meta", "gene_class", "population", "donor_df")
}

# Substrings refused inside a string constant only.
.forbidden_string_vec <- function(){
  c("/", "\\\\", "~", "\\.rds", "\\.csv", "\\.rdata", "\\.rda", "\\.tsv",
    "\\.txt", "\\.r$", "esvd2_", "::")
}

.allowed_package_vec <- function(){
  c("base", "stats", "Matrix", "locfdr", "eSVD2")
}

# Returns character(0) when the file passes, otherwise the refusals.
.tripwire <- function(file){
  expr_list <- tryCatch(parse(file = file, keep.source = TRUE),
                        error = function(e){e})
  if(inherits(expr_list, "error")){
    return(paste0("the trial file does not parse: ",
                  conditionMessage(expr_list)))
  }
  token_df <- utils::getParseData(expr_list)
  token_df <- token_df[token_df$terminal & token_df$token != "COMMENT", ]
  refusal_vec <- character(0)
  refuse <- function(text, line){
    refusal_vec <<- c(refusal_vec, paste0("line ", line, ": ", text))
  }

  bool_name <- token_df$token %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL",
                                     "SYMBOL_FORMALS", "SYMBOL_SUB",
                                     "SLOT", "SPECIAL")
  bool_string <- token_df$token == "STR_CONST"
  for(i in which(bool_name | bool_string)){
    text <- token_df$text[i]
    bare <- gsub("^[\"'`]|[\"'`]$", "", text)
    if(bare %in% .forbidden_symbol_vec()){
      refuse(paste0("`", bare, "` is not available to the trial file"),
             token_df$line1[i])
    }
    for(word in .forbidden_word_vec()){
      if(grepl(word, bare, ignore.case = TRUE)){
        refuse(paste0("`", bare, "` names `", word, "`"), token_df$line1[i])
      }
    }
    if(bool_string[i]){
      for(word in .forbidden_string_vec()){
        if(grepl(word, bare, ignore.case = TRUE)){
          refuse(paste0("the string ", text, " contains `", word, "`"),
                 token_df$line1[i])
        }
      }
    }
  }
  for(i in which(token_df$token %in% c("NS_GET_INT", "LEFT_ASSIGN") &
                 token_df$text %in% c(":::", "<<-"))){
    refuse(paste0("`", token_df$text[i], "` is not available"),
           token_df$line1[i])
  }
  for(i in which(token_df$token == "RIGHT_ASSIGN")){
    refuse("right assignment is not available", token_df$line1[i])
  }
  for(i in which(token_df$token == "SYMBOL_PACKAGE")){
    if(!token_df$text[i] %in% .allowed_package_vec()){
      refuse(paste0("package `", token_df$text[i], "` is not allowed"),
             token_df$line1[i])
    }
  }

  # The top level holds function definitions and named constants only: any
  # other expression would run in the loader's frame.
  .is_constant <- function(expr){
    if(is.numeric(expr) || is.character(expr) || is.logical(expr)){
      return(TRUE)
    }
    is.call(expr) && identical(expr[[1]], quote(`-`)) && length(expr) == 2 &&
      is.numeric(expr[[2]])
  }
  for(expr in expr_list){
    bool_ok <- is.call(expr) && identical(expr[[1]], quote(`<-`)) &&
      is.name(expr[[2]]) &&
      ((is.call(expr[[3]]) && identical(expr[[3]][[1]], quote(`function`))) ||
         .is_constant(expr[[3]]))
    if(!bool_ok){
      refuse(paste0("a top-level expression that is neither a function ",
                    "definition nor a named constant: ",
                    substr(paste(deparse(expr), collapse = " "), 1, 60)), "?")
    }
  }
  refusal_vec
}

refusal_vec <- if(file.exists(trial_file)){
  .tripwire(trial_file)
} else {
  paste0("no trial file at ", trial_file)
}
if(length(refusal_vec) > 0) .fail("crash", paste0("tripwire: ", refusal_vec))

.load_esvd2(repo_dir)

# The token scan alone is defeated by a name built at run time
# (`lapply(x, paste0("read", "RDS"))`), so the trial file is loaded into an
# environment whose only ancestor is a curated copy of base R rooted at the
# empty environment: every base function except those that read or write
# files, evaluate text, or reach an environment, the call stack, the system
# or a namespace, plus a `::` that serves the allowed packages. A name the
# trial file constructs is looked up along this chain and finds nothing
# that is not in it.
.denied_name_pattern <- function(){
  paste0("^(read|write|file|dir|list\\.files|list\\.dirs|unlink|Sys\\.|",
         "sys\\.|dyn|eval|parse|source|load|save|get|mget|assign|attach|",
         "detach|library|require|system|socket|url|open|close|pipe|fifo|",
         "gzfile|bzfile|xzfile|unz|gzcon|seek|trace|untrace|browser|debug|",
         "undebug|options|setwd|getwd|normalizePath|path\\.|basename|",
         "dirname|temp|download|sink|dput|dget|dump|serialize|unserialize|",
         "scan|shell|nsl|mem|gc|R\\.|Rprof|socket|make|lock|unlock|",
         "register|set\\.seed|RNG|truncate|flush|isatty|stdin|stdout|",
         "stderr|invokeRestart|computeRestarts|bindingIs|substitute|quote|",
         "bquote|sQuote|match\\.fun|match\\.call|Recall|reg\\.|",
         "interactive|commandArgs|quit|q$|on\\.exit|body|formals|args|",
         "delayedAssign|rapply|emptyenv|globalenv|topenv|baseenv|local|",
         "with|within|exists|rm$|remove|search|pos\\.to|date|format\\.Date|",
         "Vectorize|kronecker|Negate|Position|Find|lazyLoad|setHook|",
         "packageEvent|autoload|cat$|curl|find\\.package|serverSocket|",
         "addTaskCallback|removeTaskCallback|taskCallback|callCC)",
         "|[Ee]nv|[Nn]amespace|[Cc]onnection|[Bb]inding|[Ff]rame$|",
         "^\\.(Call|External|C$|Fortran|Internal|Primitive|Device|",
         "Library|libPaths|Platform|Options|GlobalEnv|BaseNamespaceEnv|",
         "First|Last|onLoad|onAttach|packages|dynLibs|getNamespace|",
         "userHooksEnv|S3|S4|knownS3|mergeImportMethods|doTrace|",
         "handleSimpleError|signalSimpleWarning|tryResumeInterrupt|",
         "Script|standard_regexps|rowNamesDF|TAOCP|traceback|__)")
}

.make_trial_parent <- function(){
  base_env <- baseenv()
  parent_env <- new.env(parent = emptyenv())
  # Names the pattern would deny by accident and a rule needs.
  keep_vec <- c("function", "seq_len", "seq_along", "data.frame", "seq",
                "sequence", "seq.int", "nlevels", "levels", "make.names",
                "make.unique", "isTRUE", "isFALSE", "format", "formatC",
                "as.data.frame", "is.data.frame")
  for(name in ls(base_env, all.names = TRUE)){
    value <- base_env[[name]]
    bool_denied <- grepl(.denied_name_pattern(), name) ||
      name %in% .forbidden_symbol_vec()
    bool_keep <- name %in% keep_vec ||
      (!bool_denied && (is.function(value) ||
                          name %in% c("pi", "T", "F", "letters", "LETTERS",
                                      ".Machine")))
    if(bool_keep){
      base::assign(name, value, envir = parent_env)
    }
  }
  # `tryCatch` and `try` let a rule fall back; the clock is on the
  # scorer's side, so catching a condition gains nothing else.
  for(name in c("tryCatch", "try", "suppressWarnings", "stop", "warning")){
    base::assign(name, base_env[[name]], envir = parent_env)
  }
  allowed_vec <- .allowed_package_vec()
  denied_pattern <- .denied_name_pattern()
  base::assign("::", function(pkg, name){
    pkg <- as.character(substitute(pkg))
    name <- as.character(substitute(name))
    if(!pkg %in% allowed_vec){
      stop("package `", pkg, "` is not available to the trial file")
    }
    if(pkg == "base"){
      if(!exists(name, envir = parent_env, inherits = FALSE)){
        stop("`base::", name, "` is not available to the trial file")
      }
      return(parent_env[[name]])
    }
    if(grepl(denied_pattern, name)){
      stop("`", pkg, "::", name, "` is not available to the trial file")
    }
    getExportedValue(pkg, name)
  }, envir = parent_env)
  if(exists(":::", envir = parent_env, inherits = FALSE)){
    rm(list = ":::", envir = parent_env)
  }
  parent_env
}
trial_env <- new.env(parent = .make_trial_parent())
load_res <- tryCatch(sys.source(trial_file, envir = trial_env),
                     error = function(e){e})
if(inherits(load_res, "error")){
  .fail("crash", paste0("the trial file does not load: ",
                        conditionMessage(load_res)))
}
defined_vec <- ls(trial_env, all.names = TRUE)
shadow_vec <- defined_vec[sapply(defined_vec, function(name){
  any(sapply(c("base", "stats", "Matrix", "locfdr", "eSVD2", "parallel",
               "utils", "tools", "methods"),
             function(pkg){exists(name, envir = asNamespace(pkg),
                                  inherits = FALSE)}))
})]
if(length(shadow_vec) > 0){
  .fail("crash", paste0("the trial file redefines `", shadow_vec,
                        "`, a function of base R or of a loaded package"))
}
if(!is.function(trial_env$report_lfc_se)){
  .fail("crash", "the trial file does not define `report_lfc_se()`")
}

# Calling the trial function ---------------------------------------------------

# The call runs in a forked child whose global environment is emptied
# first, so the trial function cannot reach this script's objects; the
# child returns the value and nothing else. The clock is on this side: a
# child still running at the limit is killed.
.call_trial <- function(view, time_limit, seed_number){
  trial_fn <- trial_env$report_lfc_se
  # Evaluated here: the child's global environment is emptied before the
  # promises would otherwise be forced.
  force(view)
  force(seed_number)
  # The child's random stream is fixed by the example and the cohort, so a
  # rule that draws random numbers scores the same on every scoring; the
  # trial file itself cannot seed.
  job <- parallel::mcparallel({
    rm(list = ls(envir = globalenv(), all.names = TRUE), envir = globalenv())
    set.seed(seed_number)
    tryCatch(trial_fn(view),
             error = function(e){
               structure(list(message = conditionMessage(e)),
                         class = "trial_error")
             })
  }, silent = TRUE, mc.set.seed = FALSE)
  call_start <- proc.time()[["elapsed"]]
  res <- parallel::mccollect(job, wait = FALSE, timeout = time_limit)
  if(is.null(res) && proc.time()[["elapsed"]] - call_start < time_limit - 1){
    return(structure(list(message = "the child process died"),
                     class = "trial_error"))
  }
  if(is.null(res)){
    tools::pskill(job$pid, signal = tools::SIGKILL)
    parallel::mccollect(job, wait = TRUE)
    return(structure(list(message = paste0("no value within ", time_limit,
                                           " seconds")),
                     class = "trial_timeout"))
  }
  if(is.null(res[[1]])){
    return(structure(list(message = "the child process returned nothing"),
                     class = "trial_error"))
  }
  res[[1]]
}

# The contract of the returned value. Returns "" or the reason.
.check_return <- function(value, gene_vec){
  if(inherits(value, "trial_error")) return(value$message)
  if(inherits(value, "try-error")) return(as.character(value))
  if(!is.data.frame(value)) return("the returned value is not a data frame")
  if(!all(c("gene", "log2fc", "se") %in% colnames(value))){
    return("the returned data frame lacks `gene`, `log2fc` or `se`")
  }
  if(nrow(value) != length(gene_vec)){
    return(paste0("returned ", nrow(value), " rows for ", length(gene_vec),
                  " genes"))
  }
  if(!identical(as.character(value$gene), gene_vec)){
    return("`gene` is not `colnames(view$dat)` in order")
  }
  if(!is.numeric(value$log2fc) || !is.numeric(value$se)){
    return("`log2fc` or `se` is not numeric")
  }
  if(any(!is.finite(value$log2fc)) || any(!is.finite(value$se))){
    return(paste0(sum(!is.finite(value$log2fc) | !is.finite(value$se)),
                  " genes have a non-finite `log2fc` or `se`"))
  }
  if(any(value$se <= 0)) return("`se` is not positive for every gene")
  ""
}

# Scoring ----------------------------------------------------------------------

# The Gaussian negative log-likelihood of the truth, per gene (nats).
.nll_vec <- function(lfc_vec, se_vec, truth_vec){
  0.5 * log(2 * pi * se_vec^2) + (lfc_vec - truth_vec)^2 / (2 * se_vec^2)
}

# The Figure 4 moments as Writeup50 computes them (copied from the SE
# exploration's helpers): S_A^2 = mean(lfc_A^2) - mean(SE_A^2),
# rho_A = S_A^2 / mean(lfc_A^2), k = S_A / S_B,
# r = mean(lfc_A lfc_B) / (S_A S_B), on the genes whose SE is at or below
# the 95% quantile of each cohort's SEs.
.figure4_errors <- function(lfc_a_vec, lfc_b_vec, se_a_vec, se_b_vec,
                            truth_a_vec, truth_b_vec){
  bool_vec <- se_a_vec <= stats::quantile(se_a_vec, 0.95) &
    se_b_vec <= stats::quantile(se_b_vec, 0.95)
  lfc_a_vec <- lfc_a_vec[bool_vec]
  lfc_b_vec <- lfc_b_vec[bool_vec]
  signal_ms_a <- mean(lfc_a_vec^2) - mean(se_a_vec[bool_vec]^2)
  signal_ms_b <- mean(lfc_b_vec^2) - mean(se_b_vec[bool_vec]^2)
  truth_ms_a <- mean(truth_a_vec[bool_vec]^2)
  truth_ms_b <- mean(truth_b_vec[bool_vec]^2)
  # The realized reliability of A (the truth's share of the mean squared
  # logFC), and the error of its estimate from the reported SEs.
  rho_truth <- truth_ms_a / mean(lfc_a_vec^2)
  rho_err <- abs(signal_ms_a - truth_ms_a) / mean(lfc_a_vec^2)
  if(signal_ms_a <= 0 || signal_ms_b <= 0){
    return(c(k_err = NA_real_, r_err = NA_real_, rho_err = rho_err,
             rho_truth = rho_truth))
  }
  k_val <- sqrt(signal_ms_a / signal_ms_b)
  r_val <- mean(lfc_a_vec * lfc_b_vec) / sqrt(signal_ms_a * signal_ms_b)
  truth_k <- sqrt(truth_ms_a / truth_ms_b)
  truth_r <- mean(truth_a_vec[bool_vec] * truth_b_vec[bool_vec]) /
    sqrt(truth_ms_a * truth_ms_b)
  c(k_err = abs(k_val - truth_k), r_err = abs(r_val - truth_r),
    rho_err = rho_err, rho_truth = rho_truth)
}

# Per-cohort diagnostics from the returned triple and the truth.
.score_cohort <- function(value, truth, status_vec){
  lfc_vec <- value$log2fc
  se_vec <- value$se
  truth_vec <- as.numeric(truth$truth_log2fc_vec)
  class_vec <- truth$gene_class_vec
  nll_vec <- .nll_vec(lfc_vec, se_vec, truth_vec)
  bool_null <- class_vec == "null"
  bool_weak <- class_vec == "weak"
  z_vec <- lfc_vec / se_vec
  padj_vec <- stats::p.adjust(2 * stats::pnorm(-abs(z_vec)), method = "BH")

  deseq2_df <- truth$deseq2_df
  bool_deseq2 <- is.finite(deseq2_df$deseq2_lfc) &
    is.finite(deseq2_df$deseq2_se) & deseq2_df$deseq2_se > 0
  deseq2_nll_vec <- .nll_vec(deseq2_df$deseq2_lfc[bool_deseq2],
                             deseq2_df$deseq2_se[bool_deseq2],
                             truth_vec[bool_deseq2])
  bool_capped <- status_vec != "estimated"
  mean_or_na <- function(vec){if(length(vec) > 0) mean(vec) else NA_real_}

  c(score = mean(nll_vec),
    mse = mean((lfc_vec - truth_vec)^2),
    coverage_95 = mean(abs(lfc_vec - truth_vec) <=
                         stats::qnorm(0.975) * se_vec),
    z_sd_null = if(sum(bool_null) > 1) stats::sd(z_vec[bool_null]) else
      NA_real_,
    bias_null = if(any(bool_null)) stats::median(lfc_vec[bool_null]) else
      NA_real_,
    weak_cross = sum(lfc_vec[bool_weak] * truth_vec[bool_weak]),
    weak_square = sum(truth_vec[bool_weak]^2),
    wald_false = if(any(bool_null)) 100 * mean(padj_vec[bool_null] <= 0.05)
    else NA_real_,
    vs_deseq2 = mean(nll_vec[bool_deseq2]) - mean(deseq2_nll_vec),
    score_capped = mean_or_na(nll_vec[bool_capped]),
    score_estimated = mean_or_na(nll_vec[!bool_capped]),
    median_se = stats::median(se_vec))
}

factor_vec <- c("donors", "cells", "rates", "effects", "covariates",
                "donor_sd")
deadline <- time_start + budget_seconds

# One example: both cohorts' views go to the trial function before the
# truth file is read.
.score_example <- function(example_idx){
  row_df <- manifest_df[example_idx, c("id", "split", "cell_idx",
                                       factor_vec)]
  example_start <- proc.time()[["elapsed"]]
  fail_row <- function(status, error){
    row_df$status <- status
    row_df$error <- error
    row_df
  }
  if(example_start > deadline){
    return(fail_row("timeout", "the scoring's budget ran out"))
  }
  public <- readRDS(file.path(bank_dir, manifest_df$file_public[example_idx]))
  value_list <- list()
  status_list <- list()
  for(cohort in c("a", "b")){
    view <- .build_view(public[[cohort]])
    value <- .call_trial(view = view, time_limit = cohort_seconds,
                         seed_number = manifest_df$seed[example_idx] +
                           if(cohort == "a") 1 else 2)
    if(inherits(value, "trial_timeout")){
      return(fail_row("timeout", paste0("cohort ", cohort, ": ",
                                        value$message)))
    }
    reason <- .check_return(value = value, gene_vec = colnames(view$dat))
    if(reason != ""){
      return(fail_row("crash", paste0("cohort ", cohort, ": ", reason)))
    }
    value_list[[cohort]] <- value
    status_list[[cohort]] <- as.character(view$nuisance_status_vec)
    rm(view)
  }

  truth <- readRDS(file.path(bank_dir, manifest_df$file_truth[example_idx]))
  score_a_vec <- .score_cohort(value = value_list$a, truth = truth$a,
                               status_vec = status_list$a)
  score_b_vec <- .score_cohort(value = value_list$b, truth = truth$b,
                               status_vec = status_list$b)
  mean_vec <- colMeans(rbind(score_a_vec, score_b_vec), na.rm = TRUE)
  mean_vec[is.nan(mean_vec)] <- NA_real_
  # The two sums of the gate's slope are pooled, not averaged.
  mean_vec[c("weak_cross", "weak_square")] <-
    (score_a_vec + score_b_vec)[c("weak_cross", "weak_square")]
  figure4_vec <- .figure4_errors(
    lfc_a_vec = value_list$a$log2fc, lfc_b_vec = value_list$b$log2fc,
    se_a_vec = value_list$a$se, se_b_vec = value_list$b$se,
    truth_a_vec = as.numeric(truth$a$truth_log2fc_vec),
    truth_b_vec = as.numeric(truth$b$truth_log2fc_vec)
  )
  row_df$status <- "ok"
  row_df$error <- ""
  row_df <- cbind(row_df, as.data.frame(as.list(mean_vec)),
                  as.data.frame(as.list(figure4_vec)))
  row_df$score_a <- score_a_vec[["score"]]
  row_df$score_b <- score_b_vec[["score"]]
  row_df$seconds <- proc.time()[["elapsed"]] - example_start
  row_df
}

row_list <- parallel::mclapply(seq_len(nrow(manifest_df)), .score_example,
                               mc.cores = num_cores, mc.preschedule = FALSE)
bool_row <- sapply(row_list, is.data.frame)
if(!all(bool_row)){
  error_vec <- sapply(row_list[!bool_row], function(x){
    if(is.null(x)) "a worker died" else as.character(x)[1]
  })
  .fail("crash", paste0("scorer worker: ", error_vec))
}
column_vec <- unique(unlist(lapply(row_list, colnames)))
row_list <- lapply(row_list, function(row_df){
  row_df[setdiff(column_vec, colnames(row_df))] <- NA
  row_df[, column_vec, drop = FALSE]
})
score_df <- do.call(rbind, row_list)
score_df$bank <- .bank_version()

# The per-example table --------------------------------------------------------

trial_name <- sub("\\.R$", "", basename(trial_file))
table_file <- if(bool_reference){
  file.path(output_dir, "diagnostics",
            paste0("reference_", trial_name, "_", split, "_scores.csv"))
} else if(split == "dev"){
  file.path(output_dir, "last_scores.csv")
} else {
  file.path(output_dir, "test_scores.csv")
}
utils::write.csv(score_df, file = table_file, row.names = FALSE)

if(any(score_df$status == "timeout")){
  .fail("timeout", score_df$error[score_df$status == "timeout"])
}
if(any(score_df$status != "ok")){
  .fail("crash", paste0(score_df$id[score_df$status != "ok"], ": ",
                        score_df$error[score_df$status != "ok"]))
}

# The summary block ------------------------------------------------------------

num_examples <- nrow(score_df)
se_of_mean <- function(vec){stats::sd(vec) / sqrt(length(vec))}
slope_weak <- sum(score_df$weak_cross) / sum(score_df$weak_square)
bool_gate <- slope_weak >= gate_range[1] && slope_weak <= gate_range[2]

delta_vs_best <- NA_real_
paired_se <- NA_real_
best_file <- file.path(output_dir, "best_scores.csv")
best_note <- ""
if(split == "dev" && !bool_reference && file.exists(best_file)){
  best_df <- utils::read.csv(best_file, stringsAsFactors = FALSE)
  bool_same <- setequal(best_df$id, score_df$id) &&
    all(best_df$bank == .bank_version()) && all(best_df$status == "ok")
  if(bool_same){
    diff_vec <- score_df$score - best_df$score[match(score_df$id, best_df$id)]
    delta_vs_best <- mean(diff_vec)
    paired_se <- se_of_mean(diff_vec)
  } else {
    best_note <- "best_scores.csv is from another bank or split; ignored"
  }
}

fmt <- function(x, digits = 4){
  if(is.na(x)) "NA" else formatC(x, format = "f", digits = digits)
}
line <- function(key, value){
  cat(formatC(paste0(key, ":"), width = -28), value, "\n", sep = "")
}
seconds_val <- proc.time()[["elapsed"]] - time_start

if(best_note != "") print(best_note)
cat("---\n")
line("score", fmt(mean(score_df$score)))
line("score_se", fmt(se_of_mean(score_df$score)))
line("rmse", fmt(sqrt(mean(score_df$mse))))
line("coverage_95", fmt(mean(score_df$coverage_95)))
line("z_sd_null", fmt(mean(score_df$z_sd_null, na.rm = TRUE)))
line("bias_null", fmt(mean(score_df$bias_null, na.rm = TRUE)))
line("slope_weak", fmt(slope_weak))
line("wald_false", fmt(mean(score_df$wald_false, na.rm = TRUE), 3))
line("k_err", fmt(mean(score_df$k_err, na.rm = TRUE)))
line("r_err", fmt(mean(score_df$r_err, na.rm = TRUE)))
line("figure4_valid", fmt(mean(!is.na(score_df$k_err)), 3))
line("rho_err", fmt(mean(score_df$rho_err)))
line("score_vs_deseq2", fmt(mean(score_df$vs_deseq2, na.rm = TRUE)))
for(factor_name in factor_vec){
  for(level in sort(unique(score_df[[factor_name]]))){
    line(paste0("score_", factor_name, "_", level),
         fmt(mean(score_df$score[score_df[[factor_name]] == level])))
  }
}
line("score_status_capped", fmt(mean(score_df$score_capped, na.rm = TRUE)))
line("score_status_estimated",
     fmt(mean(score_df$score_estimated, na.rm = TRUE)))
line("delta_vs_best", fmt(delta_vs_best))
line("paired_se", fmt(paired_se))
line("gate", paste0(if(bool_gate) "pass" else "fail", " (slope_weak ",
                    fmt(slope_weak), " in [", gate_range[1], ", ",
                    gate_range[2], "])"))
line("seconds", fmt(seconds_val, 1))
line("n_examples", num_examples)
line("bank", .bank_version())
line("budget", budget_seconds)
line("commit", commit_val)
line("trial_md5", trial_md5)
line("status", if(bool_gate) "ok" else "gate_failed")

if(split == "test" && !bool_reference){
  runs_file <- file.path(output_dir, "test_runs.tsv")
  if(!file.exists(runs_file)){
    cat("time\tcommit\ttrial_md5\tbank\tscore\tscore_se\tn_examples\n",
        file = runs_file)
  }
  cat(paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), commit_val,
            trial_md5, .bank_version(), fmt(mean(score_df$score)),
            fmt(se_of_mean(score_df$score)), num_examples, sep = "\t"),
      "\n", sep = "", file = runs_file, append = TRUE)
}
