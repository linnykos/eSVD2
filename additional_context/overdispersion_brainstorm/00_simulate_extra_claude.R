# Four extra regimes for the overdispersion brainstorm
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# The six regimes of additional_context/version_comparison all have 15 genes
# with a large effect and mean counts of 1 to 4.5, so every candidate finds
# every DE gene and they cannot be told apart on power. These three use the
# same generator with a weaker signal, lower counts, and a wider range of
# rates. A fourth, `trend`, ties each gene's rate to its expression
# (helpers_simulate_claude.R), which no other regime does.
#
# Run from the package root:
#   Rscript additional_context/overdispersion_brainstorm/00_simulate_extra_claude.R

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
data_dir <- file.path("additional_context", "overdispersion_brainstorm", "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(cmp_dir, "helpers_claude.R"))
source(file.path("additional_context", "overdispersion_brainstorm",
                 "helpers_simulate_claude.R"))

# The regimes ------------------------------------------------------------------

num_genes <- 300
num_reps <- 10
weak_effect_vec <- c(rep(c(0.3, -0.3), length.out = 30),
                     rep(0, num_genes - 30))

regime_spec_list <- list(
  weak_de = list(cc_effect_vec = weak_effect_vec),
  low_count = list(cc_effect_vec = 2 * weak_effect_vec,
                   gene_intercept_range = c(-2.5, -0.5)),
  wide_rate = list(cc_effect_vec = weak_effect_vec,
                   nuisance_range = c(0.5, 200))
)

for(rep_idx in seq_len(num_reps)){
  for(regime in names(regime_spec_list)){
    # Offset by 6 so that no seed repeats one of the version comparison.
    seed_number <- 1000 * rep_idx + 6 + match(regime, names(regime_spec_list))
    arg_list <- c(regime_spec_list[[regime]],
                  list(num_cells_per_individual = 30,
                       num_genes = num_genes,
                       seed_number = seed_number))
    res <- do.call(.simulate_cohort, arg_list)
    file_stem <- sprintf("%s__rep%02d", regime, rep_idx)
    saveRDS(res, file = file.path(data_dir, paste0(file_stem, ".rds")))
  }
}

for(rep_idx in seq_len(num_reps)){
  res <- .simulate_cohort_trend(num_genes = num_genes,
                                seed_number = 1000 * rep_idx + 30)
  saveRDS(res, file = file.path(data_dir,
                                sprintf("trend__rep%02d.rds", rep_idx)))
}
print(paste0("Wrote ", num_reps * (length(regime_spec_list) + 1),
             " data sets to ", data_dir))
