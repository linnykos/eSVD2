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

## Project Status (as of 2026-09-06)

**Goal: get `eSVD2` onto CRAN.** Correctness first; efficiency is explicitly out
of scope for now.

Ten package sessions in (session 11 on 2026-09-06 changed nothing here; it
wrote the vetting lessons into the shared `vet-r-package` skill). **The package now passes `R CMD check --as-cran` with the
vignette built: `Status: 2 NOTEs`** (`New submission`, and a missing HTML
Tidy on this machine), tarball **1.0 MB**, suite **687 pass / 0 fail /
0 warnings / 1 skip** in 20 s. Session 9's commit `b36376c` is the last
commit; **this session's work is uncommitted** — 45 modified/renamed files
plus `NEWS.md`, `cran-comments.md`, two new test files. Review, then commit
on `devel`.

**Read `additional_context/CRAN_READINESS.md` §0 first.** It is rewritten
each session and holds: what changed (§0.1, nine new correctness fixes N1–N9
plus packaging), the state of every audit item (§0.2), **twelve questions
for Kevin** (§0.3) and nine suggested tests (§0.4). The §1–§5 audit below it
is the frozen 2026-08-27 text. `TEST_RUN_REPORT.md` Part 7 and
`UNIT_TEST_PLAN.md` are unchanged.

**What still stands between the package and a submission** is decisions,
not code: `Authors@R` for Yixuan Qiu and Kathryn Roeder; version `1.0.2` vs
`1.1.0`; whether to keep the misspelled argument names; confirming the ASD
tutorials as pkgdown articles; confirming the new rank-deficiency refusal;
then a Windows/Linux check (`devtools::check_win_devel()`, rhub) and one
sanitizer run.

## Key Methodological Details

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

**Decisions for Kevin — the full list with context is `CRAN_READINESS.md`
§0.3 (Q1–Q12).** In order of consequence:

1. `Authors@R`: add Yixuan Qiu and Kathryn Roeder as `aut`? One spelling of
   Kevin's name across `DESCRIPTION` / `LICENSE`.
2. Version: `1.0.2` (set) or `1.1.0`.
3. Rename `bool_library_includes_interept` / `library_multipler` now or in a
   later minor version (recommend later, with a deprecation shim).
4. Confirm ASD tutorials as `vignettes/articles/` pkgdown articles.
5. Confirm refusing rank-deficient covariates at `initialize_esvd()` (and the
   `variables_enumerate_all` consequence).
6. Keep the aggregated line-search warning, or downgrade to `verbose`.
7. `multtest()`'s fallback warning on < ~40 genes: keep as is?
8. `generate_null()` gene names `gene_1` vs Seurat's `gene-1`.
9. Regenerated 400 × 60 legacy fixture vs porting the legacy tests to the
   built fixtures.
10. Carried over: nuisance blow-up; `.multtest_locfdr` catching warnings;
    `bool_diet` keeping the fit; the non-determinism amplification.

**Ready to start, blocked on nothing:**

11. Commit this session's work on `devel` after review.
12. Windows / Linux checks (`devtools::check_win_devel()`,
    `rhub::rhub_check()`) and one ASan/UBSan run.
13. The nine suggested tests in `CRAN_READINESS.md` §0.4, once the decisions
    they depend on are made.
14. Optional: `.split_individuals_by_arm()` refactor (four copies of the
    case/control derivation) and §1.6 (`.compute_df()` recomputes).
