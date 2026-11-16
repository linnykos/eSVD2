# Gabitto split-half study, step 2: eSVD-DE on one half of one split
# Drafted by Claude for Kevin Z. Lin, 2026-11-15

library(Seurat)

rm(list = ls())

# Fits eSVD-DE (eSVD2 1.2.0, the settings of Tati's Writeup45 fit of the
# whole cohort) on the donors of one half of one split. One process per
# fit, so the four fits run side by side.
#
# Run from the repository root:
#   VECLIB_MAXIMUM_THREADS=2 Rscript \
#     additional_context/2026-11-15_gabitto/02_fit-esvd_claude.R split1_half1
#
# Writes <out_dir>/gabitto_<label>_esvd.RData (the diet eSVD object) and,
# while it runs, <out_dir>/gabitto_<label>_esvd_tmp.RData.

# Paths and settings -----------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context", "2026-11-15_gabitto")
stopifnot(dir.exists(study_dir))
source(file.path(study_dir, "paths_claude.R"))

fit_label <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(fit_label %in% fit_label_vec)
split <- sub("_half[12]$", "", fit_label)
half <- as.integer(sub("^split[12]_half", "", fit_label))

# eSVD2 1.2.0 is the working tree's build; the system library has 1.0.1.2.
lib_dir <- file.path(repo_dir, "additional_context", "version_comparison",
                     "lib", "devel")
.libPaths(c(normalizePath(lib_dir), .libPaths()))
library(eSVD2)
stopifnot(as.character(utils::packageVersion("eSVD2")) == "1.2.0")

# The half's cells -------------------------------------------------------------

load(data_file)
load(file.path(out_dir, "gabitto_splits.RData"))
donor_vec <- donor_df$donor_id[donor_df[[split]] == half]
keep_vec <- colnames(seurat_obj)[as.character(seurat_obj$donor_id) %in%
                                   donor_vec]
seurat_obj <- subset(seurat_obj, cells = keep_vec)
print(paste0(fit_label, ": ", length(donor_vec), " donors, ",
             ncol(seurat_obj), " cells, ", nrow(seurat_obj), " genes"))
print(table(unique(seurat_obj@meta.data[, c("donor_id",
                                            "CERAD_binary")])$CERAD_binary))

# The fit ----------------------------------------------------------------------

# `eSVD_helper()` seeds nothing itself; Tati's script set this seed before
# its call.
set.seed(10)
time_start <- Sys.time()
esvd_obj <- eSVD2::eSVD_helper(
  batch_var_prefix = NULL,
  case_control_levels = case_control_levels,
  case_control_var = case_control_var,
  categorical_vars = categorical_vars,
  id_var = id_var,
  numerical_vars = numerical_vars,
  seurat_obj = seurat_obj,
  verbose = 2,
  intermediate_save = file.path(out_dir, paste0("gabitto_", fit_label,
                                                "_esvd_tmp.RData"))
)
seconds <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
stopifnot(inherits(esvd_obj, "eSVD"))
print(paste0("Fitted in ", round(seconds), " s"))

date_of_run <- Sys.time()
session_info <- utils::sessionInfo()
note <- paste0("eSVD2::eSVD_helper() of eSVD2 1.2.0 on the ", fit_label,
               " donors of Writeup43_gabitto_cleaned_only-counts.RData, ",
               "case-control CERAD_binary, adjusted for sex, APOE4status, ",
               "PMI, Ageatdeath (Tati's Writeup45 settings; every other ",
               "argument at its default, so k = 30 and the diet object). ",
               "donor_vec is the half's donors.")
save(esvd_obj, donor_vec, fit_label, seconds, date_of_run, session_info, note,
     file = file.path(out_dir, paste0("gabitto_", fit_label, "_esvd.RData")))

print("Done! :)")
