# Simulated cohorts for the master (3d5f7bf) vs devel (1.1.0) comparison
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Writes one .rds per regime into additional_context/version_comparison/data/,
# so both package versions are fitted on byte-identical inputs. Needs no
# version of eSVD2 except for the `generate_null` regime, which uses the devel
# library's `generate_null()` (the function is unchanged since 3d5f7bf).
#
# Run from the package root:
#   Rscript additional_context/version_comparison/01_simulate_data_claude.R

rm(list = ls())

cmp_dir <- file.path("additional_context", "version_comparison")
data_dir <- file.path(cmp_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(cmp_dir, "helpers_claude.R"))

# The regimes ------------------------------------------------------------------

# 20 individuals (10 case / 10 control) x 30 cells, 300 genes: large enough
# that locfdr usually converges, small enough that one fit takes a few seconds.
library(eSVD2, lib.loc = file.path(cmp_dir, "lib", "devel"))

num_genes <- 300
num_reps <- 10
num_cells_per_individual <- 30
strong_effect_vec <- c(rep(c(3, -3), length.out = 15),
                       rep(0, num_genes - 15))

regime_spec_list <- list(
  baseline = list(),
  null = list(num_de = 0),
  strong_de = list(cc_effect_vec = strong_effect_vec),
  near_poisson = list(nuisance_range = c(50, 1000)),
  minimal_design = list(bool_covariates = FALSE)
)

truth_list <- list()
for(rep_idx in seq_len(num_reps)){
  for(regime in c(names(regime_spec_list), "generate_null")){
    seed_number <- 1000 * rep_idx + match(regime, c(names(regime_spec_list),
                                                    "generate_null"))
    if(regime == "generate_null"){
      res <- .simulate_generate_null(num_cells_per_individual = num_cells_per_individual,
                                     num_genes = num_genes,
                                     num_individuals = 20,
                                     seed_number = seed_number)
    } else {
      arg_list <- c(regime_spec_list[[regime]],
                    list(num_cells_per_individual = num_cells_per_individual,
                         num_genes = num_genes,
                         seed_number = seed_number))
      res <- do.call(.simulate_cohort, arg_list)
    }

    file_stem <- sprintf("%s__rep%02d", regime, rep_idx)
    saveRDS(res, file = file.path(data_dir, paste0(file_stem, ".rds")))

    nuisance_vec <- res$nuisance_true_vec
    if(is.null(nuisance_vec)) nuisance_vec <- rep(NA, length(res$gene_vec))
    truth_list[[file_stem]] <- data.frame(regime = regime,
                                          rep = rep_idx,
                                          gene = res$gene_vec,
                                          is_de = res$is_de_vec,
                                          truth_log2fc = as.numeric(res$truth_log2fc_vec),
                                          nuisance_true = as.numeric(nuisance_vec),
                                          mean_count = as.numeric(colMeans(res$dat)))
  }
}

truth_df <- do.call(rbind, truth_list)
out_dir <- file.path(cmp_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(truth_df, file.path(out_dir, "truth.csv"), row.names = FALSE)
print(paste0("Wrote ", length(truth_list), " data sets to ", data_dir))
