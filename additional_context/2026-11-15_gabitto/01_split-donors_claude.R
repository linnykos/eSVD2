# Gabitto split-half study, step 1: two donor-disjoint halves, twice
# Drafted by Claude for Kevin Z. Lin, 2026-11-15

library(Seurat)

rm(list = ls())

# Splits the 80 Gabitto donors into two halves, stratified by CERAD arm x sex
# x APOE4 status, once per seed. No donor is in both halves of a split. The
# two splits use different seeds, so which donor goes where differs.
#
# Run from the repository root:
#   Rscript additional_context/2026-11-15_gabitto/01_split-donors_claude.R
#
# Writes <out_dir>/gabitto_splits.RData and csv/split_donors.csv,
# csv/split_balance.csv beside this script.

# Paths ------------------------------------------------------------------------

repo_dir <- getwd()
study_dir <- file.path(repo_dir, "additional_context", "2026-11-15_gabitto")
stopifnot(dir.exists(study_dir))
source(file.path(study_dir, "paths_claude.R"))
csv_dir <- file.path(study_dir, "csv")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Settings ---------------------------------------------------------------------

seed_vec <- c(split1 = 10, split2 = 20)
strata_vars <- c("CERAD_binary", "sex", "APOE4status")

# The donors -------------------------------------------------------------------

load(data_file)
metadata_df <- seurat_obj@meta.data
donor_vars <- c("donor_id", strata_vars, "PMI", "Ageatdeath")
donor_df <- unique(metadata_df[, donor_vars])
stopifnot(anyDuplicated(donor_df$donor_id) == 0)
donor_df$donor_id <- as.character(donor_df$donor_id)
donor_df$num_cells <- as.numeric(table(metadata_df$donor_id)[donor_df$donor_id])
donor_df <- donor_df[order(donor_df$donor_id), ]
rownames(donor_df) <- NULL
donor_df$stratum <- apply(donor_df[, strata_vars], 1, paste0, collapse = "_")

# Splitting --------------------------------------------------------------------

# Within a stratum the donors are permuted and dealt alternately to the two
# halves. A stratum with an odd count has one donor left over; which half
# starts the deal alternates over the odd strata (in a random order), so
# the halves' totals differ by at most one donor per arm's worth of odd
# strata rather than every leftover landing in half 1.
.split_donors <- function(donor_df, seed_number = 10){
  set.seed(seed_number)
  half_vec <- rep(NA_integer_, nrow(donor_df))
  stratum_vec <- sample(unique(donor_df$stratum))
  first_half <- sample(2, 1)
  for(stratum in stratum_vec){
    idx_vec <- which(donor_df$stratum == stratum)
    # `sample(x)` on a length-1 vector permutes seq_len(x), so permute the
    # positions rather than the values.
    idx_vec <- idx_vec[sample(length(idx_vec))]
    deal_vec <- rep(c(first_half, 3 - first_half), length.out = length(idx_vec))
    half_vec[idx_vec] <- deal_vec
    if(length(idx_vec) %% 2 == 1) first_half <- 3 - first_half
  }
  stopifnot(!anyNA(half_vec))
  half_vec
}

for(split in names(seed_vec)){
  donor_df[[split]] <- .split_donors(donor_df, seed_number = seed_vec[[split]])
}
# The two splits must be different partitions, in either labelling.
stopifnot(!all(donor_df$split1 == donor_df$split2),
          !all(donor_df$split1 == 3 - donor_df$split2))

# Balance ----------------------------------------------------------------------

balance_df <- do.call(rbind, lapply(names(seed_vec), function(split){
  do.call(rbind, lapply(seq_len(2), function(half){
    do.call(rbind, lapply(c("Control", "Case"), function(arm){
      sub_df <- donor_df[donor_df[[split]] == half &
                           donor_df$CERAD_binary == arm, ]
      data.frame(split = split,
                 half = half,
                 arm = arm,
                 num_donors = nrow(sub_df),
                 num_cells = sum(sub_df$num_cells),
                 num_female = sum(sub_df$sex == "female"),
                 num_apoe4 = sum(sub_df$APOE4status == "Y"),
                 mean_pmi = round(mean(sub_df$PMI), 2),
                 mean_age = round(mean(sub_df$Ageatdeath), 2),
                 stringsAsFactors = FALSE)
    }))
  }))
}))
print(balance_df)
print(table(split1 = donor_df$split1, split2 = donor_df$split2))

utils::write.csv(donor_df, file = file.path(csv_dir, "split_donors.csv"),
                 row.names = FALSE)
utils::write.csv(balance_df, file = file.path(csv_dir, "split_balance.csv"),
                 row.names = FALSE)

date_of_run <- Sys.time()
session_info <- utils::sessionInfo()
note <- paste("Two donor-disjoint half splits of the 80 Gabitto donors,",
              "stratified by CERAD_binary x sex x APOE4status; seeds 10",
              "(split1) and 20 (split2). donor_df has one row per donor with",
              "its half (1 or 2) under each split.")
save(donor_df, balance_df, seed_vec, date_of_run, session_info, note,
     file = file.path(out_dir, "gabitto_splits.RData"))

print("Done! :)")
