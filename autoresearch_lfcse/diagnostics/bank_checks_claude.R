rm(list = ls())

# Realized values of the bank against the memo's targets
# Drafted by Claude for Kevin Z. Lin, 2026-11-11
#
# Reads the two manifests of `lfcse-bank-1` (and, when present, the
# baseline's per-example tables) and prints the tables that
# notes/BANK.md reports: examples built, fit time, warnings, the share of
# capped genes per level of `rates`, the cells per donor, the sex imbalance
# between the arms, the population k and r, and the realized reliability
# per level of `effects` under the baseline.
#
# Run from the repository root:
#   Rscript autoresearch_lfcse/diagnostics/bank_checks_claude.R

repo_dir <- getwd()
loop_dir <- file.path(repo_dir, "autoresearch_lfcse")
source(file.path(loop_dir, "harness", "generator_claude.R"))
bank_dir <- .bank_dir()

manifest_df <- do.call(rbind, lapply(c("dev", "test"), function(split){
  utils::read.csv(file.path(bank_dir, paste0("manifest_", split, ".csv")),
                  stringsAsFactors = FALSE)
}))
print("Examples by split and status:")
print(table(manifest_df$split, manifest_df$status))

both <- function(name){
  c(manifest_df[[paste0(name, "_a")]], manifest_df[[paste0(name, "_b")]])
}
label <- function(name){rep(manifest_df[[name]], times = 2)}
summarize <- function(value_vec, group_vec){
  round(do.call(rbind, tapply(value_vec, group_vec, function(vec){
    c(min = min(vec), median = stats::median(vec), mean = mean(vec),
      max = max(vec))
  })), 3)
}

print("Fit seconds per cohort, by donors and cells:")
print(summarize(both("fit_seconds"), paste(label("donors"), label("cells"))))
print(paste0("Total fit time: ", round(sum(both("fit_seconds")) / 3600, 2),
             " core-hours over ", 2 * nrow(manifest_df), " cohorts"))
print(paste0("Cohorts with at least one warning: ",
             sum(both("num_warnings") > 0), "; cohorts with a gene that ",
             "was all zero: ", sum(both("num_zero_genes") > 0)))
print(paste0("Largest gap between a rebuilt view and its fit's posterior: ",
             max(both("view_diff"))))

print("Share of genes capped, by rates:")
print(summarize(both("frac_capped"), label("rates")))
print("Share of genes capped, by rates and cells and donors:")
print(summarize(both("frac_capped"),
                paste(label("rates"), label("cells"), label("donors"))))
print("Share of genes at the boundary, by rates:")
print(summarize(both("frac_boundary"), label("rates")))

print("Cells per cohort, by donors and cells:")
print(summarize(both("num_cells"), paste(label("donors"), label("cells"))))
bool_unequal <- label("cells") == "unequal"
print(paste0("Unequal cells per donor: smallest donor ",
             min(both("min_cells")[bool_unequal]), ", largest ",
             max(both("max_cells")[bool_unequal])))

print("Share of sex = 1 among cases and among controls, by donors:")
print(summarize(both("sex_case"), label("donors")))
print(summarize(both("sex_control"), label("donors")))

print("Population k and r (uncentred moments of the truth), by effects:")
print(summarize(manifest_df$truth_k, manifest_df$effects))
print(summarize(manifest_df$truth_r, manifest_df$effects))

for(split in c("dev", "test")){
  score_file <- file.path(loop_dir, "output",
                          paste0("baseline_", split, "_scores.csv"))
  if(!file.exists(score_file)) next
  score_df <- utils::read.csv(score_file, stringsAsFactors = FALSE)
  print(paste0("Realized reliability of cohort A under the baseline (",
               split, "), by effects:"))
  print(summarize(score_df$rho_truth, score_df$effects))
  print(paste0("... and by effects and donors (", split, "):"))
  print(summarize(score_df$rho_truth,
                  paste(score_df$effects, score_df$donors)))
}
print("Done! :)")
