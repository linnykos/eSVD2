#!/bin/bash
# Reruns every dry-run behind additional_context/OVERDISPERSION_BRAINSTORM.md
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Needs additional_context/version_comparison/run_all_claude.sh to have been
# run first: it builds lib/devel, the six original regimes, and master's
# output. Run from the package root:
#   bash additional_context/overdispersion_brainstorm/run_all_claude.sh
#
# About 4 minutes on 8 cores.

set -e
work_dir=additional_context/overdispersion_brainstorm
mkdir -p "$work_dir/output/logs"

for script in 00_simulate_extra 01_cache_fits 01b_run_master 02_run_candidates \
              03_summarize 04_downstream 05_cap_sweep 06_more_cells \
              07_deseq2_variants 08_report_tables; do
  echo "== $script"
  Rscript "$work_dir/${script}_claude.R" > "$work_dir/output/logs/$script.txt" 2>&1
done
# rmarkdown needs pandoc: either on the PATH, or through RSTUDIO_PANDOC set to
# the folder that holds RStudio's copy.
echo "== knitting the report"
Rscript -e "rmarkdown::render('$work_dir/overdispersion_cap_claude.Rmd', quiet = TRUE)"
echo "Done. Report: $work_dir/overdispersion_cap_claude.html. Tables are in $work_dir/output/summary_*.csv and output/logs/."
