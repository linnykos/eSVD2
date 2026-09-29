#!/usr/bin/env bash
# Runs the whole master (3d5f7bf) vs devel (1.1.0) comparison end to end.
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Run from the package root:
#   bash additional_context/version_comparison/run_all_claude.sh
#
# About 6 minutes on an Apple-silicon laptop (master and devel fitted in
# parallel). The Rmd reads only output/, so re-knitting needs none of this.
set -euo pipefail

cmp_dir="additional_context/version_comparison"
log_dir="${cmp_dir}/output/logs"
mkdir -p "${log_dir}"
SECONDS=0

bash "${cmp_dir}/00_install_versions_claude.sh" > "${log_dir}/00_install.txt" 2>&1
Rscript "${cmp_dir}/01_simulate_data_claude.R" > "${log_dir}/01_simulate.txt" 2>&1

Rscript "${cmp_dir}/02_run_regimes_claude.R" master > "${log_dir}/02_master.txt" 2>&1 &
Rscript "${cmp_dir}/02_run_regimes_claude.R" devel > "${log_dir}/02_devel.txt" 2>&1 &
wait
# The ablation reads master's nuisance estimates, so it runs after master.
Rscript "${cmp_dir}/02_run_regimes_claude.R" devel_swap > "${log_dir}/02_devel_swap.txt" 2>&1 &
Rscript "${cmp_dir}/03_corner_cases_claude.R" master > "${log_dir}/03_master.txt" 2>&1
Rscript "${cmp_dir}/03_corner_cases_claude.R" devel > "${log_dir}/03_devel.txt" 2>&1
wait

echo "total_seconds,${SECONDS}" > "${cmp_dir}/output/wall_time.csv"
echo "Finished in ${SECONDS} s"
