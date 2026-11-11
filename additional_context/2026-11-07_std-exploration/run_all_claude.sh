#!/bin/bash
# Runs the 2026-11-07 standard-error exploration end to end.
# Usage, from this folder:  bash run_all_claude.sh
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p cache output
# rmarkdown needs pandoc >= 2.8; RStudio's copy is used when none is on PATH.
if ! command -v pandoc > /dev/null; then
  export PATH="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64:$PATH"
fi
Rscript 01_real-data-checks_claude.R > output/01_log.txt 2>&1
Rscript 02_replicate-cohorts_claude.R > output/02_log.txt 2>&1
Rscript 03_figure4-pairs_claude.R > output/03_log.txt 2>&1
Rscript -e 'rmarkdown::render("summary_claude.Rmd", quiet = TRUE)'
echo "Done! :)"
