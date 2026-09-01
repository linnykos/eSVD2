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

## Project Status (as of 2026-09-01)

**Goal: get `eSVD2` onto CRAN.** Correctness first; efficiency is explicitly out
of scope for now.

Nine sessions in. **The code is now fixed to the test suite.** Kevin reviewed
the suite (2026-09-01) and gave the go-ahead; this session changed 19 `R/`
files, 4 `src/` files, `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `README.md`, and
added `R/eSVD_helper_claude.R` (`filter_cohort()`, `eSVD_helper()`, the
`gene_status` reinsertion). **Suite: 644 pass / 0 fail / 0 skip** (was 522 /
37 / 28). `R CMD check --as-cran` on the built tarball: **1 WARNING, 3 NOTEs**
(was 5 WARNINGs, 3 NOTEs); the one WARNING (`MASS` in tests) was fixed after
that run and the tests pass under the check itself. **Nothing is committed** —
`git status` shows 41 modified files plus 4 new; review then commit on `devel`.

The working documents:

- **`additional_context/TEST_RUN_REPORT.md`** — now has **Part 7** (this
  session): every fix by defect, the five tests and one fixture that had to
  change and why, and the non-determinism finding. Read this first.
- **`additional_context/UNIT_TEST_PLAN.md`** — the spec the code was fixed to.
  Every decision in it (§9) is implemented; §2.16 / §2.17 are built exactly to
  the five-step order.
- **`additional_context/CRAN_READINESS.md`** — the audit. **Now further behind**:
  §1.1, §1.2, §1.5–§1.8, §2.1–§2.5, §2.9 and the §7 landmines are done; §2.6,
  §2.7, §2.8 (tarball is **14.8 MB**) and §5 are not. Needs a pass.
- **`additional_context/RMPFR_REPORT.md`** — settled; `Rmpfr` is gone.

**What still blocks CRAN**, in order: (1) the tarball is 14.8 MB against the
5 MB limit — `tests/assets/synthetic_data.RData` (9.7 MB) and
`initialize_esvd1.rda` (1.3 MB) feed the twelve legacy tests; (2) the two
`asd*.Rmd` vignettes (undeclared `sparseMatrixStats`, multi-GB downloads);
(3) `DESCRIPTION` polish (§2.6: version, title, `Authors@R`); (4) `print()` →
`message()` package-wide (§5.3).

**Code review of the diff (`/code-review high`, 2026-09-01)** verified ten
findings, tabulated in `TEST_RUN_REPORT.md` §7.4. Six were fixed in-session.
**Four are open and are the first items below.**

## Key Methodological Details

- **`nuisance_vec` is the Gamma *rate* `β = 1/γ`**, the reciprocal of the paper's
  overdispersion `γ_j`. Larger `nuisance_vec` = *less* overdispersion. Now
  stated in `?eSVD`; the other roxygen blocks still need it (§5.4).
- **The pipeline was not deterministic and is numerically sensitive.** `irlba`
  starts from a random vector; two identical `eSVD()` runs differed by 2e-9 at
  initialization and by **up to 8** in the Welch statistics of the strongest
  genes (6.70 vs 7.41; 13.99 vs 14.29), while the rest agreed to 1e-5. Fixed
  for reproducibility with `.svd_start_vector()` (golden-ratio sequence through
  `qnorm`, no RNG touched) for both `irlba` and `RSpectra`; identical runs now
  agree exactly. **The amplification is real and unexplained** — most likely
  a flat direction between the latent factors and the free covariate
  coefficients in the second fit. T-OPT-03's "deterministic" only tested
  `opt_esvd` from a fixed start.
- **`gamma_rate`'s bracket fix has a large downstream consequence.** With the
  `max(s)` cap gone (session 8, Kevin-approved), Poisson-like genes get the
  true MLE, which is enormous: on `.small_esvd_obj()` (true rates 2.3–7.6)
  the estimates run **18 to 6e7, 10 of 20 genes above 1e4,
  `mean(log10) = 4.55`**, so `bool_stabilize_underdispersion = TRUE` fires and
  rescales **every** gene's nuisance by `10^-4.55`. Before the fix the cap at
  ≈266 hid this. On pure Poisson data with the correct mean, `gamma_rate`
  throws Boost's "no root" (caught; falls to `log_gamma_rate`'s `exp(10)`
  clamp). That half the fixture's genes look under-dispersed given the fitted
  mean is itself odd — the fit's `mean_mat` may be off. **Kevin's call.**
- **Efron's truncated MLE was wrong in the model, not the optimizer.**
  `theta` was optimized as a free parameter, dropping `p0 <= 1`, the constraint
  that ties the window count to the null mass; without it `sigma0` is nearly
  unidentified on a 90% window (2.57 at 200 genes). Now `(delta0, log sigma0,
  p0)` under `L-BFGS-B` with `p0 ∈ [1e-4, 1]`: 0.967 at 200, 0.983 at 1000
  (locfdr 0.987). `.multtest_simple` is moment-corrected for truncation
  (1.001 vs raw 0.79). Both verified across 50 seeds; one non-convergence in
  50 at N = 1000, reported through `convergence` and treated as a failure.
- **`.multtest_locfdr` still catches warnings**, and T-MT-03 pins that a
  200-gene run falls back. locfdr's routine "f(z) misfit" warning fires on
  large heavy-tailed gene sets while its `mlest` null is still fine, so real
  datasets may be moved off `locfdr` (this was already true before this
  session; what is new is that `multtest()` now **warns** every time).
- **Bessel: no correction** (Kevin, 2026-08-29); `teststat_eSVD =
  teststat_Welch * sqrt(n/(n-1))` is the contract (T-TSTAT-01a).
- **`.t_to_gaussian()`** is the one place `Ẑ = Φ⁻¹(F_df(T̂))` is computed, on
  the log scale mirrored through zero, used by both pipelines.
- **`bool_diet = TRUE` now keeps the final fit** (`fit_Second`: `x_mat`,
  `y_mat`, `z_mat`, `nuisance_vec`) and drops only `dat`, `covariates`,
  `fit_Init`, `fit_First`. It used to drop all three fits, leaving
  `latest_Fit` dangling. T-GS-11 and the §2.16 padding spec both presume the
  fit survives. One-line revert if Kevin disagrees.
- **Three levels share one predicate**, `.which_all_zero()` in `R/utils.R`
  (`colSums(na.rm = TRUE) == 0`): `eSVD_helper` labels, `eSVD` and
  `initialize_esvd` refuse. `.reinsert_genes()` pads by **name**, not
  position, so it does not depend on Seurat's feature order.
- **Seurat rewrites `gene_1` to `gene-1`.** The fixture's genes are now
  `gene1`; any name-based comparison between Seurat-derived output and a
  matrix built outside Seurat has to avoid underscores.
- `exportPattern` is gone; the public API is exactly the 17 `export()` lines
  in `NAMESPACE`. Tests reach internals through the namespace, and pass under
  `R CMD check`.
- The old `tests/assets/*.RData` fixtures are still used by the twelve legacy
  test files and are what makes the tarball 14.8 MB.

## Open Questions / Next Steps

**From the code review, open (ordered by consequence):**

1. **The nuisance blow-up** (Key Methodological Details, third bullet). Is a
   rate of 6e7 for a gene with true rate 5 a fit problem or a data problem,
   and should `bool_stabilize_underdispersion` really rescale all genes by
   `10^-4.55`? Options: cap the MLE at a documented maximum (e.g. `exp(10)`,
   matching `log_gamma_rate`); or make the stabilization robust (median, or a
   trimmed mean, of `log10`). **Needs Kevin.**
2. **`.multtest_locfdr` catching warnings.** Catch only errors, or only the
   specific "CM estimation failed" warnings, so a large real dataset stays on
   `locfdr`? T-MT-03 would need to change. **Needs Kevin.**
3. **Four copies of the case/control-individual derivation** (`eSVD`,
   `compute_test_statistic.eSVD`, `.compute_df`, `compute_test_per_gene`) with
   differing error messages. Refactor into one `.split_individuals_by_arm()`.
   Mechanical; do it.
4. **Helper memory**: `eSVD_helper` holds a transposed copy of the counts
   across the whole `eSVD()` call, and `.reinsert_genes` re-allocates every
   element even when nothing was removed. Efficiency, so deferred by policy;
   the early return is a two-line fix.

**Decisions for Kevin, from this session:**

5. Keep or revert **`bool_diet` keeping the final fit** (above).
6. The **non-determinism amplification** — investigate before submission, or
   accept and document?
7. `LICENSE` stub says `YEAR: 2026`, `COPYRIGHT HOLDER: Kevin Z. Lin` (the
   readiness doc suggested 2024; `Authors@R` says "Kevin Z Lin"). Pick one
   name form for all three places.

**Blocking the package's shape (unchanged):**

8. `Authors@R` roles for Yixuan Qiu and Kathryn Roeder (§2.6).
9. `asd.Rmd` / `asd-preprocess.Rmd`: vignettes or pkgdown articles? (§2.7)
   They are what puts `sparseMatrixStats` back in the check output.
10. The 14.8 MB tarball: shrink or replace `tests/assets/` and port the
    twelve legacy tests to the built fixtures (UNIT_TEST_PLAN.md §1).

**Ready to start, blocked on nothing:**

11. Commit this session's work on `devel` (after Kevin's review).
12. `CRAN_READINESS.md` pass to absorb sessions 3–9.
13. `print()`/`cat()` → `message()` package-wide (§5.3), then the §6 verbose
    tests become `expect_message()`.
14. State the rate/scale inversion in every nuisance roxygen block (§5.4).
15. Re-run `R CMD check --as-cran` after the `MASS` removal and quote it.
