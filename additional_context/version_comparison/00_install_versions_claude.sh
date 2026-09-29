#!/usr/bin/env bash
# Installs the two eSVD2 versions being compared into private libraries, so
# neither touches the system library and both can be loaded in separate R
# processes. Drafted by Claude for Kevin Z. Lin, 2026-09-29.
#
#   master : git commit 3d5f7bf (version 1.0.1.07, the GitHub `master`)
#   devel  : the current working tree (version 1.1.0)
#
# Run from the package root:
#   bash additional_context/version_comparison/00_install_versions_claude.sh
set -euo pipefail

cmp_dir="additional_context/version_comparison"
lib_master="${cmp_dir}/lib/master"
lib_devel="${cmp_dir}/lib/devel"
src_master="$(mktemp -d)/eSVD2"

mkdir -p "${lib_master}" "${lib_devel}" "${src_master}"

# `git archive` exports the tree of the commit without touching the checkout.
git archive 3d5f7bf | tar -x -C "${src_master}"
R CMD INSTALL --no-test-load --library="${lib_master}" "${src_master}"

# The working tree carries the comparison folder itself; R CMD INSTALL honors
# .Rbuildignore only when building, so install from a built tarball instead.
build_dir="$(mktemp -d)"
R CMD build --no-build-vignettes --no-manual . && mv eSVD2_*.tar.gz "${build_dir}/"
R CMD INSTALL --library="${lib_devel}" "${build_dir}"/eSVD2_*.tar.gz

Rscript -e "cat('master:', as.character(packageVersion('eSVD2', lib.loc = '${lib_master}')), '\n')"
Rscript -e "cat('devel: ', as.character(packageVersion('eSVD2', lib.loc = '${lib_devel}')), '\n')"
