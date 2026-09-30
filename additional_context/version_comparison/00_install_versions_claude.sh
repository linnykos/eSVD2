#!/usr/bin/env bash
# Installs the eSVD2 versions being compared into private libraries, so none
# touches the system library and each can be loaded in its own R process.
# Drafted by Claude for Kevin Z. Lin, 2026-09-29; updated 2026-09-29 for 1.2.0.
#
#   master      : git commit 3d5f7bf (version 1.0.1.07, the GitHub `master`)
#   devel_1.1.0 : git commit d49e402 (version 1.1.0, the uncapped nuisance
#                 rate). The dry-runs of ../overdispersion_brainstorm/ were
#                 built on it and load it from here.
#   devel       : the current working tree (version 1.2.0, the capped rate)
#
# Run from the package root:
#   bash additional_context/version_comparison/00_install_versions_claude.sh
set -euo pipefail

cmp_dir="additional_context/version_comparison"
lib_master="${cmp_dir}/lib/master"
lib_devel_110="${cmp_dir}/lib/devel_1.1.0"
lib_devel="${cmp_dir}/lib/devel"
src_master="$(mktemp -d)/eSVD2"
src_devel_110="$(mktemp -d)/eSVD2"

mkdir -p "${lib_master}" "${lib_devel_110}" "${lib_devel}" \
  "${src_master}" "${src_devel_110}"

# `git archive` exports the tree of the commit without touching the checkout.
git archive 3d5f7bf | tar -x -C "${src_master}"
R CMD INSTALL --no-test-load --library="${lib_master}" "${src_master}"

git archive d49e402 | tar -x -C "${src_devel_110}"
R CMD INSTALL --library="${lib_devel_110}" "${src_devel_110}"

# The working tree carries the comparison folder itself; R CMD INSTALL honors
# .Rbuildignore only when building, so install from a built tarball instead.
build_dir="$(mktemp -d)"
R CMD build --no-build-vignettes --no-manual . && mv eSVD2_*.tar.gz "${build_dir}/"
R CMD INSTALL --library="${lib_devel}" "${build_dir}"/eSVD2_*.tar.gz

for lib_label in master devel_1.1.0 devel; do
  Rscript -e "cat('${lib_label}:', as.character(packageVersion('eSVD2', lib.loc = '${cmp_dir}/lib/${lib_label}')), '\n')"
done
