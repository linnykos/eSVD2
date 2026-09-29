# CLAUDE_kevin.md — Kevin's Context

> **Current state only.** Every narrative section below is updated *in place* each session — overwrite, don't append. The External Locations table is keyed by (location, machine): replace the matching row when a path changes, add a row only for a new pair, and delete rows that no longer exist. The append-only dated log lives in `HISTORY_kevin.md`.
>
> **Owned by Kevin.** Only Kevin's session writes this file and `HISTORY_kevin.md`. Other collaborators may read it for context but must not edit it.

## About Kevin
- Role in project: package author and maintainer (`aut`, `cre`); first author of the eSVD-DE paper
- Background: statistics / biostatistics, single-cell methods; Department of Biostatistics, University of Washington
- Email: kzlin@uw.edu

## External Locations (per-machine paths)
Resolves the location names declared in the master `CLAUDE.md` → *External Locations* to real paths **on Kevin's machines**. One row per (location, machine) pair.

| Location name | Machine | Path | Notes |
|---|---|---|---|
| `EXAMPLES_REPO` | — | *(not recorded)* | Clone of `linnykos/eSVD2_examples`; add a row when it is checked out somewhere |
| `PAPER_DATA` | — | *(not recorded)* | Public datasets (GSE136831, GSE135893, Smillie, Velmeshev); add a row when downloaded |
| `WAS2CODE_REPO` | personal laptop (macOS), Dropbox | `~/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/tati/git/Was2CODE` | Under `archive/`, so treat as frozen. `R/esvd_helper.R` (60 lines) is the file being imported into `eSVD2`; copy it in rather than depending on this path |

Name the machine specifically enough that another collaborator can tell whether it is reachable to them. **A missing row means unknown; only an explicit *(not present)* row means known-absent.**

## Project Status (as of 2026-09-28)

**Goal: get `eSVD2` onto CRAN.** Correctness first; efficiency is explicitly out
of scope for now.

**Version is now `1.1.0`** (session 12). It adds the log2 fold change and its
standard error to the ordinary output: `log2fc_vec`, `log2fc_se_vec`,
`case_var`, `control_var` on the `eSVD` object from both test paths, a
`logFC_se` column in `report_results()`, and the exported
`compute_log_fold_change()`. The statistic is proposal 2 of the Was2CoDE
wiki page `code-esvd2.md`. **Session 12's work is uncommitted on `devel`**
(sessions 1-11 are committed, last commit `a95e533`); Kevin vets the tests,
then commits.

**`R CMD check --as-cran` on the 1.1.0 tarball: `Status: 2 NOTEs`**
(`New submission`, and a missing HTML Tidy on this machine), the same two as
before, with the vignette built and tarball 1.0 MB. Suite under
`NOT_CRAN=true`: **996 pass / 0 fail / 0 warnings / 0 skip** in 41 s across
30 files, up from 688. Three of the new tests (T-LFC-15, -16, -17) are
`skip_on_cran()`, so a bare `R CMD check` does not run them.

**Read `additional_context/UNIT_TEST_PLAN.md` section 2.18 first** for the
new feature: the 18 tests (T-LFC-01 to -18) with their oracles, the ten
breakages they were shown to catch, the five findings, the three corrections
made to the tests after their first run, and Q-LFC-1 to -3.
`additional_context/CRAN_READINESS.md` section 0 is still the session-10 text
and does not mention 1.1.0.

**What still stands between the package and a submission** is decisions,
not code: the three Q-LFC questions; `Authors@R` for Yixuan Qiu and Kathryn
Roeder; whether to keep the misspelled argument names; confirming the ASD
tutorials as pkgdown articles; confirming the rank-deficiency refusal; then a
Windows/Linux check (`devtools::check_win_devel()`, rhub) and one sanitizer
run.

## Key Methodological Details

- **The log2 SE is `(1/ln 2) * sqrt(case_var/(n1*case_mean^2) +
  control_var/(n0*control_mean^2))`**, `n` counting individuals, the
  variances being the mixture variances the Welch statistic divides by (no
  Bessel). `.compute_log2_fold_change()` is the only place the formula lives.
- **The within term is most of that SE** (median 98% of `log2fc_se^2` on
  F-SMALL, 95% on `generate_null()`), so the SE is 4 to 8 times the SD of a
  fit-fixed bootstrap over individuals. Against cohorts that are *refitted*
  it is about right: calibration ratio 0.77, coverage of +/- 2 SE 0.96, on
  genes with an ordinary nuisance estimate. Refitting moves the fold change
  far more than resampling individuals does.
- **Where the nuisance estimate diverges the SE is anti-conservative.** About
  12% of genes on the toy cohorts get a rate near 1e7 (true rates 0.1 to 10);
  the posterior collapses onto the fit, the within term vanishes, coverage
  falls to 0.36, and the Welch statistic is in the tens. This is Q10's
  downstream face. The diverged genes also inflate the geometric mean that
  `bool_stabilize_underdispersion` divides every gene's rate by.
- **The depth adjustment shifts every fold change.** `Log_UMI` is the log of
  the observed total, which moves with the DE genes: on F-SMALL five planted
  genes raise a case cell's total by 2^0.23 and the 35 null genes come out at
  a median of -0.21. A test's truth must be the generator's `nat_mat` ratio,
  and a 3-SE tolerance hides this bias.
- **`generate_null()`'s "null_large_var" genes are not null for a ratio of
  means.** The arms share a log-scale mean but differ in within-individual SD
  (0.1 vs 0.75), so the ratio of arithmetic means is exp(0.276), log2FC 0.40.
  Only the "null_interleaved" genes (even positions from 12) have truth 0.
  The Welch statistic is a contrast of arithmetic means too, so T-PROP-06,
  which treats genes 11 onward as null, is counting non-null genes as null
  (not yet investigated).
- **Against DESeq2 / dreamlet / NEBULA** on one simulated cohort (10 vs 10
  individuals, 100 genes): fold changes agree at Spearman 0.97 to 0.98; the
  eSVD2 SE is a median 1.8 to 2.1 times theirs and its rank correlation with
  theirs is 0.69 to 0.71 (theirs with each other: 0.97 to 0.99).
- **`compute_test_statistic.default()` warns on a non-positive arm mean**
  (fold change `NA`). Mean-zero Gaussian matrices, which several older tests
  use, trigger it; T-TSTAT-06 now asserts it.
- **`expect_equal(NULL, NULL)` passes**, so a test of a field that does not
  exist yet is green. The new test file guards every section-A call with a
  length check for this reason.
- **`.combine_two_named_lists()` never overwrites an existing entry**
  (T-UTIL-03 pins it), so `param` goes stale when any stage is rerun with
  different settings. The two test functions now overwrite their own two
  entries (`test_case_individuals`, `test_control_individuals`), because the
  SE divides by their lengths. Every other stage's `param` entries still go
  stale on a rerun.
- **`report_results()` reads `logFC` from `log2fc_vec`** and recomputes it
  from the arm means only for an object built before 1.1.0.
- **`nuisance_vec` is the Gamma *rate* `β = 1/γ`**, the reciprocal of the paper's
  overdispersion `γ_j`. Larger `nuisance_vec` = *less* overdispersion. Stated
  in `?eSVD`, both `?estimate_nuisance` methods and `?compute_posterior`.
- **The pipeline is numerically chaotic, and this matters for tests.** A
  4e-16 change in `z_mat` (switching the reparameterization from `lm()` to
  an algebraically identical QR solve) flipped whether `locfdr` converges on
  the 18-gene `F-TINY` fixture, turning eleven silent gene-status tests into
  warning tests. `.svd_start_vector()` makes identical inputs give identical
  output; it does not make nearby inputs give nearby output. Any test that
  compares two runs must go through the same code path or use a tolerance
  set with this in mind. The tests muffle the `locfdr` fallback warning via
  `.muffle_locfdr_fallback()` in `helper-fixtures.R`.
- **`locfdr` needs roughly 40+ genes.** It converges on the 40-gene example
  chain and fails on 18; `multtest()` warns whenever it falls back.
- **Collinear covariates are now refused twice**, by name: at
  `initialize_esvd()` (QR rank of `covariates`) and at
  `reparameterization_esvd_covariates()`. Before, `lm()` returned `NA`
  coefficients that silently entered `z_mat` and killed the second
  `opt_esvd` with "missing value where TRUE/FALSE needed". Consequence:
  `format_covariates(variables_enumerate_all = ...)` output can no longer be
  fed to `initialize_esvd(bool_intercept = TRUE)`.
- **The reparameterization is plain least squares on the design matrix**
  (`qr.coef` / `qr.fitted`), so covariate names with spaces or parentheses
  survive; `as.data.frame()` used to `make.names()` them.
- **`.initialize_coefficient()`'s GLM fallback** (zero or one covariate left
  after intercept and offsets) now carries the library-size offset. With
  `Intercept + Log_UMI + CC_1` only, it used to fit intercepts without depth
  adjustment.
- **Line-search failures are counted in C++ and warned about once in R.**
  `constr_newton` returns `linesearch_failed`; `opt_x`/`opt_yz` attach
  `num_linesearch_failed` as an attribute; `opt_esvd.default` sums and warns.
  No `Rcpp::warning()` remains in a C++ frame. Never fires on the fixtures.
- **`gamma_rate`'s bracket fix has a large downstream consequence** (session
  8; unchanged): Poisson-like genes get enormous rates, so
  `bool_stabilize_underdispersion` rescales every gene by `10^-mean(log10)`.
  Kevin's call (Q10 in the readiness doc).
- **Efron's truncated MLE** is `(delta0, log sigma0, p0)` under `L-BFGS-B`
  with `p0 ∈ [1e-4, 1]`; `.multtest_simple` is moment-corrected for
  truncation. `.multtest_locfdr` still catches warnings (Q7).
- **Bessel: no correction** (Kevin, 2026-08-29).
- **`.t_to_gaussian()`** is the one place `Ẑ = Φ⁻¹(F_df(T̂))` is computed.
- **`bool_diet = TRUE` keeps the final fit** (Q5 of session 9; unchanged).
- **Three levels share `.which_all_zero()`**; `.reinsert_genes()` pads by
  name. Seurat rewrites `gene_1` to `gene-1`; the examples rename first.
- **`print()` behind `verbose` stays** (readiness §5.3 resolved as no
  change): CRAN's reviewer boilerplate explicitly accepts `if(verbose)
  cat()`, and Kevin's style guide mandates `print(paste0())`.
- The public API is the 17 `export()` lines in `NAMESPACE`; `opt_x`,
  `opt_yz`, `data_loader`, `esvd_family` are documented `@keywords internal`.
- `devtools::document()` (roxygen2 8.1.0 here) rewrites `RoxygenNote:
  7.3.3` to `Config/roxygen2/version: 8.1.0`; reverted in sessions 9 and 10
  (Q11 in the readiness doc).
- `R CMD build` with the vignette needs pandoc; RStudio's copy at
  `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64`
  prepended to `PATH` works. Run with `_R_CHECK_FORCE_SUGGESTS_=false`
  no longer needed.

## Open Questions / Next Steps

**Decisions for Kevin, newest first.** Q-LFC-1 to -3 are in
`UNIT_TEST_PLAN.md` section 2.18; Q1-Q12 are in `CRAN_READINESS.md` section
0.3.

1. **Q-LFC-1**: the SE is unreliable where the nuisance estimate diverges.
   Document only (done), repair the divergence upstream (Q10), or flag such
   genes in `report_results()`.
2. **Q-LFC-2**: keep the warning from `compute_test_statistic.default()` on a
   non-positive arm mean, or return `NA` silently from the matrix method.
3. **Q-LFC-3**: the names `log2fc_vec` / `log2fc_se_vec` and `logFC_se`.
4. Is the mixture variance (per-cell posterior variance undivided by cells per
   individual) the intended variance for a reported SE? It is what the test
   uses and what the wiki's proposal 2 specifies; the alternative is the
   variance of the individual means alone.
5. `Authors@R`: add Yixuan Qiu and Kathryn Roeder as `aut`? One spelling of
   Kevin's name across `DESCRIPTION` / `LICENSE`.
6. Rename `bool_library_includes_interept` / `library_multipler` now or in a
   later minor version (recommend later, with a deprecation shim).
7. Confirm ASD tutorials as `vignettes/articles/` pkgdown articles.
8. Confirm refusing rank-deficient covariates at `initialize_esvd()` (and the
   `variables_enumerate_all` consequence).
9. Keep the aggregated line-search warning, or downgrade to `verbose`.
10. `multtest()`'s fallback warning on < ~40 genes: keep as is?
11. `generate_null()` gene names `gene_1` vs Seurat's `gene-1`.
12. Regenerated 400 x 60 legacy fixture vs porting the legacy tests to the
    built fixtures.
13. `param` going stale on a rerun of `opt_esvd`, `estimate_nuisance` or
    `compute_posterior` (see Key Methodological Details): fix in
    `.combine_two_named_lists()`, or leave?
14. T-PROP-06 counts `generate_null()`'s "null_large_var" genes as null;
    they are not. Investigate whether the test still means what it says.
15. Carried over: nuisance blow-up (now with a measured consequence, item 1);
    `.multtest_locfdr` catching warnings; `bool_diet` keeping the fit; the
    non-determinism amplification.

**Ready to start, blocked on nothing:**

16. Vet `tests/testthat/test_compute_log_fold_change_claude.R`, then commit
    session 12 on `devel`.
17. Windows / Linux checks (`devtools::check_win_devel()`,
    `rhub::rhub_check()`) and one ASan/UBSan run.
18. The nine suggested tests in `CRAN_READINESS.md` section 0.4, once the
    decisions they depend on are made.
19. Rebuild the pkgdown site (`docs/`) so it lists `compute_log_fold_change`.
20. Optional: `.compute_df()` could now read the stored `case_var` /
    `control_var` instead of recomputing them; `.split_individuals_by_arm()`
    refactor (four copies of the case/control derivation).
