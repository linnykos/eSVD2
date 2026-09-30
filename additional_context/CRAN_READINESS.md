# eSVD2 → CRAN: correctness, tests, and packaging

**Status:** working document, started 2026-08-27. Focus is **correctness**, not
efficiency, per the current goal. Efficiency items are listed only where they are
also a correctness or CRAN-policy hazard.

**How to read this.** Every claim is tagged:

- **[verified]** — I ran it in R 4.5.1 in this session and the output is quoted.
- **[inspection]** — read from the source; the reasoning is given, but not executed.
- **[policy]** — a CRAN rule, not a bug.

Items are ordered by how much they can hurt: silent wrong numbers first, then
things that block submission, then things a reviewer will make you fix, then test
debt. **§1–§5 below are the original audit of 2026-08-27 and are kept as
written; the current state of every item is in §0, directly below**, which is
rewritten each session. §0 also carries the open questions and the suggested
tests.

**Companion document.** The proposed test suite lives in
`additional_context/UNIT_TEST_PLAN.md` — ~200 tests with IDs, oracles, and the
open questions that block them. §3 below is now a pointer to it. Every test
citation in this document (`T-PVAL-01`, `T-CPP-FAM-03`, …) refers to that file.

Equation numbers refer to Lin, Qiu & Roeder (2024), indexed in
`additional_context/summary.md` as `lin2024esvd`.

**Parameterization warning that recurs throughout.** The paper's overdispersion
`γ_j` is the Gamma **scale**; the code's `nuisance_vec` is the Gamma **rate**
`β = 1/γ_j` (hence `src/gamma_rate.cpp`). So paper `μ/γ` ↔ code
`mean_mat * nuisance_vec`, and paper `ℓ + 1/γ` ↔ code `library_mat + nuisance_vec`.
The code is internally consistent; the *documentation* never states the inversion,
which is a documentation defect (§5.4), not a numerical one.

---

## 0. Status as of 2026-09-29 (version 1.2.0, commit `1a0a536`) — read this first

**Bottom line.** Both commands ran on HEAD of `devel`, with RStudio's pandoc
prepended to `PATH`, on macOS 14.2.1 / arm64 with R 4.5.1, on 2026-09-29:
`R CMD build` (vignette included), then `R CMD check --as-cran
eSVD2_1.2.0.tar.gz`. The result is **`Status: 2 NOTEs`, 0 ERRORs, 0
WARNINGs**. The two NOTEs are `New submission`, which is unavoidable, and
`'tidy' doesn't look like recent enough HTML Tidy`, which comes from this
machine and not the package.

- **Tests**: `[ FAIL 0 | WARN 0 | SKIP 20 | PASS 1608 ]` under the check (re-run after session 20's fixes, still 2 NOTEs). The
  20 skips are `skip_on_cran()` and source-tree tests.
- **Examples**: all run; the slowest is `eSVD`, at 0.54 s.
- **Vignette**: rebuilds in 65 s.
- **Tarball**: 1,556,652 bytes.
- **Full suite under `NOT_CRAN=true`**: 2110 pass / 0 fail / 0 warnings /
  0 skip, across 34 files and 298 `test_that` blocks.

Two NOTEs is the same result as session 10. The third NOTE that sessions 16 and
17 saw (`unable to verify current time`) came from the check's clock service
and did not recur. **What stands between the package and a submission is the
list of decisions in §0.3, a Windows/Linux check and one sanitizer run. No
known defect stands in the way.**

### 0.1 What changed since session 10 (2026-09-02)

The package went from 1.0.2 to 1.2.0 in two steps. `NEWS.md` has both in full.

**1.1.0 (commit `d49e402`, session 12): the log2 fold change and its SE.**
The fit carries `log2fc_vec` and `log2fc_se_vec`, the arm variances, and the
new exported `compute_log_fold_change()`. `report_results()` gets a
`logFC_se` column. The SE is the delta-method SE of `log2(case_mean /
control_mean)`, with `n` counting individuals and the mixture variances the
Welch statistic already divides by. Tests: `UNIT_TEST_PLAN.md` §2.18. Three
questions remain (Q-LFC-1 to -3, in §0.3).

**The uncapped nuisance rate (sessions 13 and 14).** Sessions 13 and 14
compared master (`3d5f7bf`) with 1.1.0 on simulated cohorts
(`version_comparison/`) and brainstormed fixes (`OVERDISPERSION_BRAINSTORM.md`).
The comparison found that the §4.3 bracket fix of session 8 had a large
downstream cost:

- A gene whose counts are no more variable than Poisson around the fit has
  no finite maximum-likelihood rate, and 1.1.0 gave it about 1e7.
- Its posterior then followed the fit, and its statistic was inflated.
- Under the null, false discoveries at FDR 0.05 rose from 0.2 to 2.7 per 300
  genes, and to 63 in a near-Poisson regime.

Master had been protected by an accidental cap at the gene's largest library
size, at which 78% to 100% of genes sat. This was Q10 of session 10.

**1.2.0 (commit `1a0a536`, sessions 15 to 17): the cap.** Kevin chose Idea 1
of the brainstorm, which gives the rule
`nuisance_vec = max(min(MLE, c · m_j), min_val · m_j)`, with `m_j` the gene's
median library size and `c = cap_multiplier = 10`. Around it came:

- a per-gene `nuisance_status` (`estimated` / `capped` / `boundary` /
  `failed`), also a column of `report_results()`;
- the rate before the cap, stored as `nuisance_mle_vec`;
- `recompute_pvalue()`, which redoes the test at another cap without
  refitting, and for a `bool_diet` object rebuilds the counts and design from
  the Seurat object;
- `plot_nuisance()` and `plot_fitted_vs_observed()`, with `ggplot2` and
  `ggrepel` in `Suggests`.

The decisions and the nuances of the implementation are in the "Decision"
section of `OVERDISPERSION_BRAINSTORM.md`, and the tests in
`UNIT_TEST_PLAN.md` §2.19 to §2.21. Kevin vetted the code before it was
committed. The rerun of `version_comparison/version_comparison_claude.Rmd`
compares master against 1.2.0 on the same simulated cohorts.

**Correctness fixes found on the way** (each pinned by a test):

- `estimate_nuisance()` and `compute_posterior()` now overwrite their `param`
  entries on a rerun. Before, `.combine_two_named_lists()` kept the first
  call's settings, so a redo of the posterior used stale settings (T-CAP-08,
  T-REDO-12).
- The two test functions refresh the recorded individuals, which the SE
  divides by (T-LFC-19).
- The rule that picks the library columns existed three times. It is now one
  function, `.nuisance_library_idx()`, pinned against hand-written column
  sets before the refactor (T-CAP-10, T-POST-13, T-PGENE-01).

### 0.2 Where each audit item stands

Unchanged since session 10 except where marked **(new)**.

| Item | State |
|---|---|
| §1.1 `qnorm(pt())` saturation | done (session 9), `.t_to_gaussian()` |
| §1.2 `verbose = 2` throws | done |
| §1.3 formula-based reparameterization | done (session 10) |
| §1.4 `scale()` type change / doc direction | done (session 10) |
| §1.5 pipeline drift, weak equivalence test | done. **(new)** The library-column rule is now shared (T-PGENE-01) |
| §1.6 `.compute_df()` recomputes group stats | open, low priority. It could now read the stored `case_var` / `control_var` of 1.1.0 |
| §1.7 `sigma0` unconstrained | done |
| §1.8 small items | done |
| §2.1 `sparseMatrixStats` | done |
| §2.2 `Rmpfr` | done |
| §2.3 `exportPattern` | done. **(new)** The public API is now 21 `export()` lines, including `compute_log_fold_change`, `recompute_pvalue`, `plot_nuisance` and `plot_fitted_vs_observed` |
| §2.4 `LICENSE` stub | done |
| §2.5 compiler `override` | done |
| §2.6 `DESCRIPTION` | done except `Authors@R` (Q1). **(new)** Version is 1.2.0; `ggplot2 (>= 3.4.0)` and `ggrepel` added to `Suggests` |
| §2.7 vignettes | done (articles; Q4). **(new)** `eSVD2.Rmd` now draws the two diagnostic plots |
| §2.8 tarball | done: 1.56 MB, up from 1.0 MB because the vignette draws five figures |
| §2.9 `MASS` | done |
| §4.1 pointer null checks | not done; cosmetic |
| §4.2 `Rcpp::warning` in C++ frame | done |
| §4.3 `gamma_rate` bracket | done (session 8). **(new)** Its downstream consequence, the diverging rate, is resolved by the cap of 1.2.0, applied in R so that `gamma_rate` stays an MLE |
| §4.4 valgrind / sanitizers | **not done**; recommended once before submission |
| §5.1 to §5.6 | as in session 10 (§5.5 partly done) |

### 0.3 Questions for Kevin (blocking or affecting submission)

**Resolved since session 10.**

- The old Q2 (version number) is resolved: the version is 1.2.0.
- The old Q10's nuisance blow-up is resolved by the cap.
- The three 1.2.0 decisions were confirmed on 2026-09-29:
  - a boundary gene is set to the cap itself;
  - `bool_diet` drops `covariates`;
  - `min_val` is in the units of the cap.

**Open, newest first.**

1. **One real data set under 1.2.0** (Idea 11 of the brainstorm). How far the
   cap moves the published results is unmeasured. It needs a `PAPER_DATA`
   path in `CLAUDE_kevin.md`.
2. **Q-LFC-1.** Under the cap the ±2 SE interval covers the true log2 fold
   change for at least 99% of genes of every status in five of six simulated
   regimes. The measurement is from the version comparison rerun, 10 cohorts
   each. The exception is `strong_de` (0.72 to 0.85), from the `Log_UMI`
   depth bias. Without the cap, diverged genes had 0 to 0.47. Remaining:
   reword the caveat in `?compute_log_fold_change`, which still describes the
   uncapped case.
3. **Q-LFC-2.** `compute_test_statistic.default()` warns on a non-positive
   arm mean. Keep the warning, or return `NA` silently?
4. **Q-LFC-3.** The names `log2fc_vec` / `log2fc_se_vec`, and `logFC_se` in
   `report_results()`.
5. **`plot_fitted_vs_observed()` calls `set.seed(seed_number)`** with a
   default of 10, the house convention. This resets the caller's random
   stream when the plot is drawn inside a loop. Keep it, default to `NULL`, or
   restore the stream on exit?
6. **`param` still goes stale** on a rerun of `initialize_esvd()`,
   `opt_esvd()` or the reparameterization; the later stages now overwrite
   their entries. Fix it for every stage in `.combine_two_named_lists()`, or
   leave it?
7. **T-PROP-06** counts `generate_null()`'s "null_large_var" genes as null.
   For a ratio of arithmetic means their truth is log2FC 0.40. Investigate
   whether the test still means what it says.
8. **`Authors@R`** (was Q1): add Yixuan Qiu and Kathryn Roeder as `aut`? Also
   use one spelling of Kevin's name across `DESCRIPTION` and `LICENSE`.
9. **Misspelled argument names** (was Q3): `bool_library_includes_interept`
   and `library_multipler`. The recommendation is to rename them in a later
   minor version, with a deprecation shim.
10. **Unchanged from session 10:**
    - Q4 (ASD tutorials as pkgdown articles);
    - Q5 (refusing rank-deficient covariates);
    - Q6 (the aggregated line-search warning);
    - Q7 (`multtest()`'s fallback warning under about 40 genes);
    - Q8 (`generate_null()` gene names);
    - Q9 (the regenerated legacy fixture);
    - Q11 (`RoxygenNote` against `Config/roxygen2/version`);
    - carried over from Q10: whether `.multtest_locfdr` should catch warnings,
      `bool_diet` keeping the final fit, and the non-determinism
      amplification.
11. **Before submitting:**
    - run the check on Windows and Linux (`devtools::check_win_devel()`,
      `rhub::rhub_check()`);
    - run once with sanitizers (§4.4);
    - rebuild the pkgdown site (`docs/`), which does not yet list the four
      functions added in 1.1.0 and 1.2.0.

### 0.4 Suggested new tests (from working through the code)

The session-10 list, minus what 1.2.0 made moot, plus two new entries.

| ID | What it would pin | Why it is missing |
|---|---|---|
| T-PG-06 | `compute_test_per_gene()` vs the matrix path with `bool_adjust_covariates = TRUE` **and** `pseudocount = 1` | T-PGENE-01 covers the library settings, not these two branches |
| T-CN-04 | A run in which the line search actually fails, asserting one warning whose count matches `sum(attr(...))` | Needs a reproducible failing case |
| T-INIT-11 | `eSVD()` on a Seurat object whose count layer is a dense `matrix` | Needs a dense-layer Seurat fixture |
| T-ESVD-13 | Non-syntactic metadata levels through the whole `eSVD()` path | Seurat's own handling is untested here |
| T-FMT-10 | `format_covariates(variables_enumerate_all = ...)` then `initialize_esvd()` | Blocked on Q5 |
| T-NUIS-06 | `bool_use_log = TRUE` agrees with `FALSE` on genes with status `estimated` | Now well-defined: the status says which genes are interior |
| T-MT-05 | `multtest()` with more than about 500 genes stays on `locfdr` | Depends on Q7 |
| T-SVD-05 | `.svd_start_vector()` is identical across two R sessions | Needs `callr` |
| T-LFC-20 **(new)** | T-LFC-17's coverage assertion restricted to capped and boundary genes, now that the cap keeps them from collapsing | Measured in simulation (Q-LFC-1); not yet a test |
| T-POST-14 **(new)** | A fixture on which `alpha_max` changes a statistic | `alpha_max` changed nothing on F-TINY at 1, 5, 86 or 1000; whether it ever binds is unknown |

### 0.5 What session 10 changed (kept for its reasoning)

*Session-10 bottom line, superseded by §0 above:* `R CMD build` (vignette included) + `R CMD check --as-cran`
on the tarball: **`Status: 2 NOTEs`, 0 ERRORs, 0 WARNINGs** (macOS 14.2.1 / arm64, R 4.5.1, 2026-09-02). The two NOTEs are `New submission` (unavoidable) and the environmental `'tidy' doesn't look like recent enough HTML Tidy` (this machine, not the package). Examples all run, slowest 0.55 s; vignette rebuilds in 64 s; under the check the suite reports `[ FAIL 0 | WARN 0 | SKIP 3 | PASS 685 ]` (the two extra skips read the source tree, which is absent in a tarball). Test suite: **687 pass / 0 fail / 0 warnings /
1 skip** (`T-PROP-06`, `skip_on_cran`), 20 s. Tarball **1.0 MB** (was 13.2 MB).
What remains before submission is a short list of decisions for Kevin
(§0.3), not code.


**Test-suite warnings (Kevin's question).** A warning that escapes a test is
not a CRAN failure — `R CMD check` only fails tests on errors — but each one
is noise that hides a real one, and the check log prints them. Both were
tests that produced the *intended* warning without asserting it:
`T-REP-09` (rank-deficiency warning from `.identification`) now uses
`expect_warning()`, and the `estimate_nuisance.default` extreme-`mu` test
asserts the "nuisance estimation failed for 1 of 10 gene(s)" warning. A
third source appeared mid-session and is worth knowing about: on the
20-gene `F-TINY` fixture `multtest()` legitimately falls back from
`locfdr` and warns, and **whether `locfdr` converges on 18 genes is a
numerical coin flip** — a change to the reparameterization that altered
`z_mat` by 4e-16 flipped eleven gene-status tests from silent to warning.
The test helpers `.helper_run()` / `.esvd_run()` now muffle exactly that
warning (`.muffle_locfdr_fallback()` in `helper-fixtures.R`); every other
warning still surfaces. This is the same sensitivity recorded in
`TEST_RUN_REPORT.md` §7.3, seen at a smaller scale.

**Correctness (new this session, none previously listed):**

| # | What | Where |
|---|---|---|
| N1 | **The single-covariate GLM fallback in `.initialize_coefficient()` dropped the library-size offset.** `glmnet` needs two predictors, so a design with only `Intercept + Log_UMI + CC_1` (i.e. `eSVD()` with no other covariates) went through `stats::glm()` *without* `offset = Log_UMI`, initializing every gene's intercept without adjusting for depth (intercept 1.13 vs the correct −1.6 on a Poisson(3) toy). Also `colnames(df) <- "y"` on a two-column frame left the predictor unnamed. Both fixed; oracle is `stats::glm(y ~ x + offset(log_umi))`. Tests `T-INIT-09`, `T-INIT-10`. | `R/initialization.R` |
| N2 | **§1.3 was still open and is now fixed.** `reparameterization_esvd_covariates()` regressed each factor with `stats::lm(x ~ ., data = as.data.frame(...))`, which (a) `make.names()`-mangles covariate names with spaces or parentheses into names that are not columns of `z_mat`, and (b) returns `NA` for aliased columns, which propagated into `z_mat` and killed the *second* `opt_esvd` with "missing value where TRUE/FALSE needed". Replaced by `qr()` / `qr.coef()` / `qr.fitted()` on the design matrix (numerically identical: max diff 4e-16 on `F-TINY`), with a rank check that names the aliased columns. Tests `T-REP-05`, `T-REP-06`. | `R/reparameterization.R` |
| N3 | **`initialize_esvd()` now refuses a rank-deficient `covariates`** by name. Previously glmnet's ridge silently fit it and the failure surfaced two stages later (N2). Consequence: the legacy test `initialize_esvd works` passed all 20 individual indicators and had to drop them, as `eSVD()` and the fixture script already do. **See Q5.** | `R/initialization.R`, `tests/testthat/test_initialization.R` |
| N4 | **`eSVD()` errored on a numeric or character categorical variable** (`droplevels()` on a non-factor). The single-level test now uses `unique()`; the variable is converted to a factor two lines later anyway. Found by the new `eSVD()` example. Test `T-ESVD-12`. | `R/eSVD.R` |
| N5 | **§1.4's `scale()` type change was still in the code** in both `compute_posterior.default` and `compute_test_per_gene` (T-POST-07 passed by coercion). Replaced by `10^(log10(v) - mean(log10(v)))`; `nuisance_vec` stays a named numeric vector. Roxygen for `bool_stabilize_underdispersion` now says what the code does and in which direction. | `R/posterior.R`, `R/compute_test_per_gene.R` |
| N6 | **§4.2 fixed.** `Rcpp::warning()` no longer fires inside the C++ line search. `constr_newton` returns `linesearch_failed`; `opt_x`/`opt_yz` count failures and attach `num_linesearch_failed` as an attribute; `opt_esvd.default` raises **one** aggregated R warning. `ObjectiveYZ` gained the five missing `override`s. Tests `T-CN-01..03`. **See Q6.** | `src/constrained_newton.cpp`, `src/optimization.cpp`, `R/optimization.R`, `R/utils.R` |
| N7 | `opt_esvd.default` refuses `NA` in `x_init`/`y_init`/`z_init`/`covariates` and stops with a named message when the objective becomes non-finite, instead of R's generic error from the convergence test. Test `T-OPT-06`. | `R/optimization.R` |
| N8 | `eSVD()` computed `alpha_max` as `2 * max(dat@x)`, which errors if the Seurat count layer is dense. Now `2 * max(dat)`. | `R/eSVD.R` |
| N9 | 48 `T`/`F` literals → `TRUE`/`FALSE`; six scalar `&`/`|` → `&&`/`||`; `.opt_esvd_setup_z_mat()` returns explicitly; stray no-op expression removed from `generate_null()`; empty `if` block removed from `compute_test_per_gene()`. | package-wide |

**Packaging:**

- `DESCRIPTION`: `Version: 1.0.2`; Title in title case without the package
  name; Description no longer starts with "This package" and cites
  `Lin, Qiu and Roeder (2024) <doi:10.1186/s12859-024-05724-7>`; `URL`
  gains the GitHub repo, `BugReports` added. `Suggests` trimmed to what is
  used: `knitr`, `numDeriv`, `rmarkdown` (**was missing** — the vignette
  is `rmarkdown::html_vignette`), `Seurat`, `SeuratObject`, `testthat`,
  `withr`. Dropped `devtools`, `EnhancedVolcano` (Bioconductor),
  `npregfast` (unused).
- **Tarball 13.2 MB → 1.0 MB.** `tests/assets/synthetic_data.RData`
  regenerated at 20 individuals × 20 cells × 60 genes with `compress = "xz"`
  (432 KB; was 2000 × 150 and 9.7 MB); it keeps only the objects a test
  reads. `tests/assets/data_generation.R` rewritten to match (the fixture
  is what that script produces under `devtools::load_all()`).
  `initialize_esvd1.rda` (1.3 MB) was loaded by no test and is deleted.
  **See Q10.**
- **§2.7 resolved as option (b):** `asd.Rmd`, `asd-preprocess.Rmd` and the
  three PNGs moved to `vignettes/articles/` (the `usethis::use_article()`
  layout: pkgdown renders them under the same `articles/` URLs, `R CMD
  build` excludes them via `^vignettes/articles$`). Their vignette-index
  blocks are removed, the PNGs are included by relative path, and
  `devtools::session_info()` → `sessionInfo()`. **See Q4.**
- `eSVD2.Rmd`: `library(Seurat)` replaced by a `requireNamespace()` guard on
  the one UMAP chunk; `T`/`F` fixed; `sessionInfo()`.
- **`@examples` on every exported function** (21 Rd files). The chain is
  `generate_null(15 cells × 8 individuals × 40 genes)` → each step; the whole
  chain runs in about 0.5 s and `locfdr` converges at 40 genes. The three
  Seurat-dependent examples (`eSVD`, `eSVD_helper`, `filter_cohort`) are
  guarded by `requireNamespace("SeuratObject")`; the first two are also in
  `\donttest{}` (they still run under `--as-cran`).
- `@returns` for `opt_x`, `opt_yz`, `data_loader`, `esvd_family`,
  `print.esvd_data_loader` (was the placeholder "description"); all five
  marked `@keywords internal`. `.reparameterize` is `@noRd`. §5.4: the
  rate/scale inversion is now stated in `estimate_nuisance.eSVD`,
  `estimate_nuisance.default`, `compute_posterior.default` and `eSVD`.
  Typos (`outputed`, `expoential`, `flexbility`, "defficiency" in a warning)
  fixed; `sfari_df` `@source` pointed at the Velmeshev paper and now points
  at SFARI; `gandal_df`'s example carried a personal `/Users/kevinlin/...`
  path.
- `NEWS.md` (ships) and `cran-comments.md` (`.Rbuildignore`d) created.
  README: version, installation, dependency paragraph; `compatibile`.
- Junk removed: `inst/include/.Rapp.history` (0 bytes, would have shipped),
  `vignettes/.DS_Store`, the empty `tests/testthat/_snaps/`.

**Deliberately not changed:**

- **§5.3 `print()` → `message()`: resolved as *no change*.** Every
  `print()`/`cat()` in `R/` is behind `verbose`. CRAN's own reviewer
  boilerplate ("Instead of print()/cat() rather use message()/warning() or
  **if(verbose) cat(..)**") explicitly accepts verbose-gated output, and
  Kevin's R style guide mandates `print(paste0(...))` for progress. The one
  ungated `message()` (in `report_results` when results are missing) is
  fine as is.
- `Config/testthat/edition: 3` and removing `context()` from the twelve
  legacy files: not done. Edition 3 changes `expect_equal()` semantics and
  would have to be re-verified across 687 tests for no CRAN benefit.
- The misspelled argument names `bool_library_includes_interept` and
  `library_multipler` (see Q3).

---

## 1. Correctness defects (silent wrong numbers)

### 1.1 `qnorm(pt(t, df))` saturates to `+Inf` for up-regulated genes, silently disabling the empirical null — **[verified]**

`R/compute_pvalue.R:78` and `R/compute_test_per_gene.R:~300`:

```r
gaussian_teststat <- stats::qnorm(stats::pt(teststat_vec, df = df_vec))
```

This implements `Ẑ_j = Φ⁻¹(F_df(T̂_j))` (paper, "Performing multiple testing
correction"). In double precision the composition is **asymmetric**: the lower
tail stays finite far longer than the upper tail, because `pt` underflows to `0`
gracefully but saturates at exactly `1`.

Measured at `df = 18`:

| `t` | `pt(t, df)` | `qnorm(pt(t, df))` |
|---|---|---|
| −40 | 2.43e−19 | −8.915 |
| −10 | 4.47e−09 | −5.750 |
| +10 | 1.000000 | +5.750 |
| **+40** | **1** | **`Inf`** |

The numerically stable form gives the correct finite value:
`qnorm(pt(-40, 18, log.p = TRUE), log.p = TRUE)` = `-8.915293`.

**The damage is not the `Inf` itself — it is what happens downstream.** I traced
the whole chain:

1. One `Inf` enters `multtest()`.
2. `.multtest_locfdr()` calls `locfdr::locfdr(teststat_vec, plot = 0)`, which
   **errors**: `'to' must be a finite number` — [verified].
3. That error is swallowed by the `tryCatch` at `R/multtest.R:56–63`, which
   returns `NA`.
4. `multtest()` falls through to `.multtest_truncatedGauss()`, which computes
   `mean`/`sd` of a vector containing `Inf` → `NA`.
5. It falls through again to `.multtest_simple()` — the method whose own comment
   reads *"this is a purposefully overly simple method, meant to use when
   everything else has failed. Realistically, it is not a good estimator."*

So a single strongly up-regulated gene silently downgrades the entire empirical
null — the paper's headline mechanism for controlling Type-1 error — to an
estimator the authors themselves disclaim. **And the user cannot tell.**
`multtest()` returns `method` in its list, but `compute_pvalue()` (and
`compute_test_per_gene()`) keep only `fdr_vec`, `null_mean`, `null_sd` and
**discard `method`**, so `pvalue_list` carries no record of which estimator ran.

Note this is exactly the regime the method is designed to find: real DE genes in
a cohort with 15+ individuals per arm routinely produce |t| well past 10.

**Fix.**
1. Use the log-scale composition:
   `stats::qnorm(stats::pt(teststat_vec, df = df_vec, log.p = TRUE), log.p = TRUE)`.
   This is a drop-in replacement — it is exact for the finite cases and finite
   where the current form is `Inf`.
2. Store `method` in `pvalue_list`, and `warning()` when the fallback chain
   drops below `locfdr`. A silent degradation of the null model is not
   acceptable in a published method.
3. Guard `multtest()` at entry: `stopifnot(all(is.finite(teststat_vec)))` after
   the fix, so a future regression surfaces loudly.

**Test (new, `test_compute_pvalue.R`):** feed a test-statistic vector containing
`t = 40, df = 18`; assert `gaussian_teststat` is finite, assert
`pvalue_list$method == "locfdr"`, and assert monotonicity of
`log10pvalue` in `|teststat|`. See `UNIT_TEST_PLAN.md` `T-PVAL-01`, `T-PVAL-02`,
`T-PVAL-03`, `T-MT-03`, `T-MT-04`.

---

### 1.2 `opt_esvd.default(verbose = 2)` throws — **[verified]**

`R/optimization.R:230–231`:

```r
print("Residual of loss: ", resid)
print("Threshold for termination: ", thresh)
```

`print()` has no `...`-of-values contract; the second argument binds to
`print.default`'s `digits`. With `resid` a small double this errors:

```
Error: invalid printing digits 0
```

[verified]. So the most verbose diagnostic mode of the package's core optimizer
crashes instead of printing. It has clearly never been executed — which is itself
the finding: `verbose` branches are untested throughout (§3.6).

**Fix.** `message("Residual of loss: ", resid)` etc. Prefer `message()` over
`print()`/`cat()` package-wide anyway (§5.3).

---

### 1.3 `reparameterization_esvd_covariates()` breaks on non-syntactic covariate names — **[inspection]**

`R/reparameterization.R:~120`:

```r
tmp_df <- cbind(x_mat[,ell], covariate_mat2)
colnames(tmp_df)[1] <- "x"
tmp_df <- as.data.frame(tmp_df)
lm_res <- stats::lm(x ~ . , data = tmp_df)
coef_vec <- stats::coef(lm_res)
names(coef_vec)[1] <- "Intercept"
...
z_mat[j, names(coef_vec)] <- z_mat[j, names(coef_vec)] + coef_vec * y_mat[j, ell]
```

Two failure modes, both reachable from ordinary user input:

- **Name mangling.** `as.data.frame()` applies `make.names()`. A covariate column
  produced by `format_covariates()` is `paste0(var, "_", lvl)`, and `lvl` comes
  from the user's factor levels — e.g. a diagnosis level `"ASD (severe)"` yields
  the column `Diagnosis_ASD (severe)`, which `make.names()` rewrites to
  `Diagnosis_ASD..severe.`. `names(coef_vec)` then carries the mangled name,
  which is not a column of `z_mat` → `subscript out of bounds`. Note
  `format_covariates()` performs **no** name sanitization, so nothing upstream
  prevents this.
- **Rank deficiency.** If `covariate_mat2` is collinear (the paper explicitly
  warns that including individual one-hot vectors makes `C` collinear), `lm`
  drops aliased terms and `coef()` returns `NA` for them. Those `NA`s propagate
  straight into `z_mat`, and from there into `μ̂` and every posterior. There is no
  check.

**Fix.** Do the projection in linear algebra rather than via a formula — this is
a least-squares fit of `x_mat[,ell]` on `cbind(1, covariate_mat2)`, so use
`qr.solve()` / `qr.coef()` on the design matrix directly and keep the original
column names. Guard with `stopifnot(all(is.finite(coef_vec)))` and check the QR
rank, erroring with a message that names the collinear covariates.

**Test:** run `reparameterization_esvd_covariates()` with (a) a factor level
containing a space and parentheses, (b) a deliberately collinear covariate pair;
assert an informative error rather than `NA`-contamination or a subscript error.
See `UNIT_TEST_PLAN.md` `T-REP-05`, `T-REP-06`, `T-FMT-09`.

---

### 1.4 `bool_stabilize_underdispersion` — code and documentation disagree, and `scale()` silently changes the type — **[verified for the type change; inspection for the doc mismatch]**

`R/posterior.R:628–630`:

```r
if(bool_stabilize_underdispersion & mean(log10(nuisance_vec)) > 0) {
  nuisance_vec <- 10^(scale(log10(nuisance_vec), center = T, scale = F))
}
```

- **Doc mismatch.** The roxygen says the rescaling happens when *"the global mean
  over-dispersion is less than 1"*, but the guard fires when
  `mean(log10(nuisance_vec)) > 0`, i.e. when the geometric mean of `nuisance_vec`
  is **greater** than 1. Given `nuisance_vec` is the *rate* `1/γ` (see the
  warning at the top), a large value means *low* overdispersion — so the code may
  well be right and the prose wrong. Either way one of them must change, and a
  reader cannot currently tell which.
- **Type change.** `scale()` returns an `n × 1` **matrix** and moves `names()` to
  `rownames()` — [verified]: for `c(a=10, b=100, c=1000)` the result has
  `class = c("matrix","array")`, `names(out)` is `NULL`, `rownames(out)` is
  `a b c`. Downstream `sweep(..., STATS = nuisance_vec)` and `nuisance_vec[j]`
  still work by coercion, so this is currently latent rather than active — but it
  means `nuisance_vec` has two different types depending on a data-dependent
  branch, which is exactly the kind of thing that becomes a bug on the next edit.
- `&` should be `&&` here (scalar guard).

**Fix.** `nuisance_vec <- 10^(log10(nuisance_vec) - mean(log10(nuisance_vec)))`,
preserving names and type. Reconcile the roxygen with the intended direction and
say which way "over/under-dispersion" runs in this parameterization.

**Test:** assert `is.numeric(nuisance_vec) && !is.matrix(nuisance_vec)` and that
`names()` survive, on both branches of the guard. See `UNIT_TEST_PLAN.md`
`T-POST-07`; the direction question is `T-POST-08` / `Q-POST-1`.

---

### 1.5 The two pipelines can drift, and the one test comparing them can't detect it — **[inspection]**

`compute_test_per_gene()` is documented as reproducing
`compute_posterior()` + `compute_test_statistic()` + `compute_pvalue()` without
allocating `n × p` matrices. Two problems:

- **Defaults already differ.** `compute_posterior.eSVD()` has
  `library_min = 0.1`; `compute_test_per_gene()` has `library_min = 1e-2`.
  Called with defaults, the two paths give different answers. (`eSVD()` passes
  `library_min` explicitly, so its two branches agree — but a user calling either
  function directly gets a silent discrepancy.)
- **The equivalence test is too weak.** `test_compute_test_per_gene.R:45`:

  ```r
  expect_true(abs(sum(res1$teststat_vec - res2$teststat_vec)) <= 1e-3)
  ```

  This is `abs(sum(differences))`, so per-gene errors of opposite sign cancel.
  A path that got every gene wrong by ±5 in alternating signs would pass.

**Fix.** Align the defaults (pick 0.1, matching `compute_posterior`, and say so).
Replace the assertion with `expect_equal(res1$teststat_vec, res2$teststat_vec,
tolerance = 1e-8)`, which compares element-wise. Extend it to `case_mean`,
`control_mean`, `pvalue_list$df_vec` and `pvalue_list$gaussian_teststat`.

**This is the single highest-value test in the package** — it pins one
implementation against the other across the whole posterior/test/p-value stack.

---

### 1.6 `.compute_df()` recomputes group statistics instead of reusing them — **[inspection]**

`compute_pvalue()` calls `.compute_df()`, which re-derives
`case_individuals`, the averaging matrix, and both group means and variances —
the exact quantities `compute_test_statistic()` already computed and threw away.
Beyond the wasted work, the two copies of the code can diverge under editing, and
the df would then be computed from different group statistics than the test
statistic it is paired with. It is also the reason `compute_pvalue()` requires
`input_obj$dat` to still exist, which conflicts with `bool_diet = TRUE` in
`eSVD()`.

**Fix.** Have `compute_test_statistic()` store `case_gaussian_var`,
`control_gaussian_var`, `n1`, `n2` on the object; make `.compute_df()` a pure
function of those four. Keep `.compute_df()` exported-as-internal for the tests.

---

### 1.7 `.multtest_truncatedGauss()` optimizes `sigma0` unconstrained — **[inspection]**

`R/multtest.R:~120`. The objective takes `param_vec = c(delta0, sigma0, theta)`
and bounds only `theta` (`if(theta > 0.99 | theta < 0.01) return(Inf)`).
`sigma0` is free, and Nelder-Mead is free to step it negative; `stats::dnorm`
and `stats::pnorm` with `sd < 0` return `NaN` (with a warning) — [verified for
the `NaN` behaviour, `pnorm(0, 0, -1)` is `NaN`]. I did **not** reproduce a run
where Nelder-Mead actually went negative, so treat this as a latent hazard rather
than a demonstrated failure.

**Fix.** Optimize `log(sigma0)`, or return `Inf` for `sigma0 <= 0` alongside the
existing `theta` guard. Also add `optim(..., control = list(maxit = ...))` and
check `optim_res$convergence`, which is currently ignored.

---

### 1.8 Small correctness/robustness items — **[inspection]**

| Where | Issue | Fix |
|---|---|---|
| `R/data_management.R` (last `else`) | `stopifnot("what_obj is not found")` — a character is not a logical. It *does* error, with the confusing message `"what_obj is not found" is not TRUE` — [verified] | `stop("what_obj not recognized: ", what_obj)` |
| `R/fisher_test.R` | `stats::dhyper(..., log = F)` uses `F` for `FALSE` | Use `FALSE`; `T`/`F` are rebindable variables. Applies package-wide — `F`/`T` appear as defaults in `initialize_esvd`, `estimate_nuisance.*`, `compute_posterior.*`, `.svd_safe` callers |
| `R/reparameterization.R` `.identification()` | `sqrt(eigen_sym$values)` with no non-negativity guard; the function already `warning()`s about rank deficiency but proceeds | Clamp at 0 or error; a `NaN` here silently poisons `x_mat`/`y_mat` |
| `R/reparameterization.R` `.irlba_custom`, `.rpsectra_custom` | `check_stability & K > 5` — `&` on scalars | `&&` |
| `R/optimization.R` | If `covariates = NULL`, `z_mat` stays `NULL` and is returned as such; works only because `.opt_esvd_format_matrices()` short-circuits on `is.null(covariates)` | Add an explicit test for the no-covariate path (§3.3) |
| `R/optimization_helper.R` `.opt_esvd_setup_z_mat()` | Returns the value of an `if`/`else` whose branches are assignments — correct, but by accident | Return `z_mat` explicitly |
| `R/multtest.R`, `R/compute_pvalue.R` | `Rmpfr::pnorm()` on plain doubles is **identical** to `stats::pnorm()` — [verified]: `Rmpfr::pnorm(-40, log.p=TRUE)` returns a `numeric`, exactly equal to `stats::pnorm(-40, log.p=TRUE)` = −804.6084 | See §2.2 — this dependency buys nothing and costs a system library |

---

## 2. CRAN blockers (submission will be rejected)

**Hard blockers, in the order I would fix them** (the subsection numbering below
is not the priority order):

| # | Blocker | Effort |
|---|---|---|
| §2.8 | Tarball is **13.1 MB**; CRAN's limit is **5 MB** | small — delete unused fixture objects |
| §2.1 | `sparseMatrixStats` (Bioconductor) in `Imports` | small — ~4 lines of `Matrix` |
| §2.4 | `LICENSE` is not a valid DCF stub | trivial — 2 lines |
| §2.5 | Install-time compiler warnings | trivial — add `override` ×5 |
| §2.9 | `MASS` used in tests, not declared | trivial — one `Suggests` entry |
| §2.3 | 14 undocumented exported objects | medium — needs `eSVD()` docs |
| §2.2 | `Rmpfr` dependency buys nothing | small |
| §2.6 | `DESCRIPTION` field wording/version | small |
| §2.7 | Vignettes need Bioconductor packages | medium — decide vignette vs. article |

### 2.1 `sparseMatrixStats` is a Bioconductor package in `Imports` — **[verified]**

`available.packages()` against CRAN: `sparseMatrixStats` is **not on CRAN**
(neither is `EnhancedVolcano`, in `Suggests`). CRAN packages may only depend on
packages available from CRAN or a repository declared in
`Additional_repositories`, and a Bioconductor package in **`Imports`** is a hard
stop — CRAN's check machines will not install it. **[policy]**

It is used in exactly one place, `.compute_matrix_sd()` in
`R/reparameterization.R`:

```r
if(inherits(x = mat, what = c('dgCMatrix', 'dgTMatrix'))) {
  sd_vec <- sparseMatrixStats::colSds(mat)
} else {
  sd_vec <- matrixStats::colSds(mat)
}
```

...and `.compute_matrix_sd()` is only ever reached with `sd_vec = NULL` from
`.svd_safe()` callers (`.initialize_residuals`, `.factorize_matrix`), so the
branch is **dead code in practice**.

**Options, in order of preference:**
1. Drop the dependency and compute sparse column SDs from `Matrix` primitives:
   `colMeans` and `colMeans(x^2)` give `sd = sqrt(n/(n-1) * (E[x²] − E[x]²))`.
   ~4 lines, no dependency, and it removes a Bioconductor edge entirely.
2. Move to `Suggests` + `Additional_repositories: https://bioconductor.org/packages/release/bioc`
   and guard with `requireNamespace()`.
3. Submit to Bioconductor instead of CRAN (the `biocViews:` field in
   `DESCRIPTION` suggests this was once the plan).

Whichever is chosen, decide it before anything else in this document — it changes
what the package *is*.

### 2.2 `Rmpfr` is a heavy system dependency that provides no benefit — **[verified]**

`Rmpfr` requires the GMP and MPFR C libraries. The README already warns:
*"We have noted that `Rmpfr` is sometimes tricky to install due to its required
C++ libraries."* `DESCRIPTION` declares **no `SystemRequirements:` field**, which
is itself a check finding. **[policy]**

And it does nothing here. `Rmpfr::pnorm` is used only on plain doubles
(`R/multtest.R`, `R/compute_pvalue.R`, `R/compute_test_per_gene.R`). Its methods
dispatch to `stats::pnorm` for non-`mpfr` input, so — [verified] —
`Rmpfr::pnorm(-40, mean=0, sd=1, log.p=TRUE)` returns a base `numeric`,
`identical()` to `stats::pnorm(-40, log.p=TRUE)`. No extra precision is obtained
because no argument is ever an `mpfr` object.

**Fix.** Replace all `Rmpfr::pnorm` with `stats::pnorm` and drop `Rmpfr` from
`Imports`. This removes a system-library dependency, several classes of install
failure, and one line of README apology, at zero numerical cost. (If arbitrary
precision *is* genuinely wanted in the far tail, it has to be done deliberately —
`Rmpfr::mpfr(x, precBits = ...)` first — but §1.1's `log.p` fix already recovers
the range that matters.)

### 2.3 `exportPattern("^[[:alpha:]]+")` exports 14 undocumented objects — **[verified]**

`R/zzz.R` carries `@exportPattern "^[[:alpha:]]+"`, so the NAMESPACE exports
every non-dot-prefixed object in the package. Cross-referencing the 89 functions
defined under `R/` against the 33 files in `man/`, these are exported with **no
documentation**:

```
eSVD                     data_loader_description   test_data_loader
objfn_all_r              objfn_Xi_r                objfn_YZj_r
grad_Xi_r                grad_YZj_r
hessian_Xi_r             hessian_YZj_r
feas_Xi_r                feas_YZj_r
gamma_rate               log_gamma_rate
```

`R CMD check` reports this as a **WARNING** (`Undocumented code objects`), which
blocks submission. **[policy]**

Two distinct problems here:

- **The raw Rcpp bindings are public API.** `objfn_Xi_r`, `hessian_YZj_r`,
  `test_data_loader`, etc. take external pointers created by `data_loader()` and
  `esvd_family()`. If a user calls one with a pointer restored from a saved
  `.rds`, the `Rcpp::XPtr` is null and dereferencing it segfaults (§4.1). These
  should not be exported at all.
- **`eSVD()` — the package's main entry point — is exported by accident.** It has
  no roxygen block, no `.Rd`, and no test. The one function a new user is most
  likely to call is the least documented thing in the package. This is the single
  biggest documentation gap.

**Fix.** Delete `@exportPattern` from `R/zzz.R`; add explicit `@export` to the
intended public API; write a full roxygen block for `eSVD()` including
`@examples`. Keep `@useDynLib` / `@importFrom(Rcpp, evalCpp)`.

### 2.4 The `LICENSE` file is not a valid license stub — **[verified]**

`R CMD check` reports:

```
* checking DESCRIPTION meta-information ... NOTE
License stub is invalid DCF.
```

`DESCRIPTION` declares `License: MIT + file LICENSE`. For that declaration, R
requires `LICENSE` to be a **two-line DCF stub**, not the license text — the MIT
text itself lives in R's shared license database. The current file is the full
20-line MIT text beginning `Copyright (c) 2016 C. Kevin Lin`.

**Fix.** Replace the entire file with exactly:

```
YEAR: 2024
COPYRIGHT HOLDER: Kevin Z. Lin
```

Two side notes worth resolving at the same time: the year is **2016**, and the
name is **"C. Kevin Lin"**, which matches neither the current year nor the
`Authors@R` entry ("Kevin Z Lin"). Pick one form of the name and use it in
`DESCRIPTION`, `LICENSE`, and the README.

### 2.5 Compiler warnings flagged as significant — **[verified]**

```
* checking whether package ‘eSVD2’ can be installed ... WARNING
Found the following significant warnings:
  ./data_loader.h:106:25: warning: 'operator++' overrides a member function but is not marked 'override' [-Winconsistent-missing-override]
  optimization.cpp:73:19: warning: 'grad' ...      [-Winconsistent-missing-override]
  optimization.cpp:81:19: warning: 'hessian' ...   [-Winconsistent-missing-override]
  optimization.cpp:89:10: warning: 'direction' ... [-Winconsistent-missing-override]
  optimization.cpp:97:10: warning: 'feas' ...      [-Winconsistent-missing-override]
```

CRAN treats "significant warnings" at install time as a submission blocker.

**Fix.** Add the `override` specifier to the five virtual-function overrides
named above (they implement `VecIterator::operator++` and the four
`Objective` pure virtuals declared in `src/objective.h`). This is mechanical and
also protects against a signature drifting out of sync with the base class — a
class of bug that `override` exists to catch and that would otherwise silently
create an overload instead of an override.

### 2.6 `DESCRIPTION` cleanup — **[policy]**

| Field | Problem | Fix |
|---|---|---|
| `Version: 1.0.1.07` | Four components with a leading zero in the last; CRAN wants `major.minor.patch[.dev]` and reads `07` oddly | `1.0.2` |
| `Title` | Sentence-case, too long, repeats the package name ("eSVD2 for performing the eSVD-DE for…") | Title Case, ≤65 chars, no package name, no "R package" |
| `Description` | Does not cite the method paper in CRAN's required form | Add `Lin, Qiu and Roeder (2024) <doi:10.1186/s12859-024-05724-7>` |
| `Depends: R (>= 3.5.0), Rcpp` | `Rcpp` in `Depends` attaches it to the user's search path | Move `Rcpp` to `Imports` |
| `biocViews:` | Bioconductor-only field; CRAN flags it as non-standard | Remove (unless §2.1 resolves toward Bioconductor) |
| `SystemRequirements:` | Absent, though `Rmpfr` needs GMP/MPFR | Moot once §2.2 removes `Rmpfr`; otherwise declare it |
| `Suggests: devtools` | Never used in tests or vignettes (only in `tests/assets/data_generation.R`, which is not run at check time) | Drop |
| `Authors@R` | Yixuan Qiu and Kathryn Roeder are credited in the paper's author contributions ("KZL and YQ coded eSVD-DE in R and C++") but are not in `Authors@R` | Add with appropriate roles; also add `Roeder` as `aut` if intended. **Needs Kevin's decision.** |
| `LazyData: true` | Only valid with a `data/` dir — fine here, but confirm the four `.rda` are compressed (`tools::checkRdaFiles`) | Verify |

### 2.7 Vignettes load `Suggests` packages unconditionally — **[inspection]**

`vignettes/asd.Rmd:34` and `asd-preprocess.Rmd:34` are **evaluated** chunks that
run `library(Seurat)`, `library(SeuratObject)`, `library(EnhancedVolcano)`,
`library(sparseMatrixStats)`. `EnhancedVolcano` and `sparseMatrixStats` are
Bioconductor-only (§2.1), so on CRAN's machines these vignettes fail to build.
The heavy analysis chunks are correctly `eval = FALSE`, and both vignettes
require external downloads (a UCSC `rawMatrix.zip`, a Dropbox `.RData`), so they
can never really run at check time. **[policy]**

**Fix.** Either (a) guard the setup chunk:

```r
knitr::opts_chunk$set(eval = requireNamespace("Seurat", quietly = TRUE) &&
                             requireNamespace("EnhancedVolcano", quietly = TRUE))
```

and use `requireNamespace()` rather than `library()`, or (b) move `asd.Rmd` and
`asd-preprocess.Rmd` out of `vignettes/` into `pkgdown/articles/` so they ship on
the website but not in the tarball. **(b) is cleaner** — they are tutorials
requiring multi-GB downloads, not vignettes. `eSVD2.Rmd` is self-contained
(the author states both simulations complete in under 2 minutes) and should stay,
but must be timed against CRAN's 10-minute total check budget.

### 2.8 The tarball is 13.1 MB — CRAN's limit is 5 MB — **[verified]**

`R CMD check --as-cran` reports `Size of tarball: 13111968 bytes`. CRAN rejects
source packages over **5 MB** without a prior arrangement. This is the single
largest blocker and, happily, the easiest to fix.

The cause is `tests/assets/synthetic_data.RData` (10.1 MB on disk). It holds
eleven objects, of which most are never used by any test:

| Object | Class | In-memory size | Needed? |
|---|---|---|---|
| `eSVD_obj` | `eSVD` | 8.7 MB | a **fully fitted** object — the tests refit anyway |
| `dat` | `dgCMatrix` | 2.9 MB | **yes** |
| `library_mat` | matrix 2000×150 | 2.5 MB | generation by-product |
| `gamma_mat` | matrix 2000×150 | 2.3 MB | generation by-product |
| `nat_mat_nolib` | matrix 2000×150 | 2.3 MB | generation by-product |
| `covariates` | matrix 2000×23 | 0.5 MB | **yes** |
| `metadata` | data.frame | 0.1 MB | **yes** |
| `session_info` | `session_info` | 36 KB | no — provenance, belongs in the script |
| `nuisance_vec`, `true_cc_status`, `date_of_run` | — | ~3 KB | small, keep |

`tests/assets/initialize_esvd1.rda` adds another 1.4 MB, and the three vignette
PNGs (`asd-umap.png`, `volcano.png`, `volcano_hk.png`) add ~1.5 MB.

**Fix, in order of impact:**
1. Regenerate `synthetic_data.RData` keeping only `dat`, `covariates`,
   `metadata`, `nuisance_vec`, `true_cc_status` — and save with
   `save(..., compress = "xz")`. This alone should drop ~15 MB of in-memory
   objects. Update `tests/assets/data_generation.R` to match so provenance stays
   honest.
2. Shrink the fixture itself. 2000 cells × 150 genes is far larger than any unit
   test needs; 200 cells × 40 genes over 8 individuals exercises every code path
   and runs faster. CRAN also caps total check time at 10 minutes.
3. Move the vignette PNGs to the pkgdown site, or downsample them. If §2.7
   resolves by moving `asd*.Rmd` to `pkgdown/articles/`, the PNGs go with them.
4. Re-verify with `R CMD build` and check the reported tarball size.

Note this interacts with §2.7: the PNGs exist only to display results the
`eval = FALSE` chunks cannot produce at build time.

### 2.9 `MASS` is used in the test suite but not declared — **[verified]**

```
* checking for unstated dependencies in ‘tests’ ... WARNING
'::' or ':::' import not declared from: ‘MASS’
```

`tests/testthat/test_reparameterization.R` calls `MASS::mvrnorm()` at 12 sites
(lines 44, 45, 75, 76, 101, 102, 120, 121, 142, 143, 161, 162). It passes locally
only because `MASS` ships with R and happens to be installed.

**Fix.** Add `MASS` to `Suggests`. Alternatively drop the dependency —
`mvrnorm(n, mu, Sigma)` here is only ever used with `diag(5)`, `2*diag(5)` and
one `toeplitz(5:1)`, all of which are reproducible with
`matrix(rnorm(n*5), n, 5) %*% chol(Sigma)`. Dropping it is slightly preferable
since it removes a fixture dependency from the test suite entirely.

(Unrelated but adjacent: `src/gamma_rate.cpp:77` mentions `library(MASS)` inside
a comment block of example R code. That is inert and needs no action.)

---

## 3. Missing unit tests → see `UNIT_TEST_PLAN.md`

**This section has moved.** The full proposed test suite now lives in
`additional_context/UNIT_TEST_PLAN.md`, which supersedes what used to be here:
everything §3 listed is absorbed there, expanded with the C++ backend, error
messages, fixtures, and an explicit list of the questions that block individual
tests. Each test carries an ID (`T-<AREA>-nn`) so it can be referred to in
review, plus an **oracle** field saying where its expected answer comes from —
which is the field that separates a real correctness test from a snapshot of
current behaviour.

**Read `UNIT_TEST_PLAN.md` before writing any test.** It is a proposal awaiting
review, not an agreed plan: ~200 `test_that` blocks against the current ~40, and
16 open questions that need Kevin's answer before the corresponding tests can be
written.

What it adds that this section did not have:

- **§1 Fixtures.** The single largest piece of work. The current
  `synthetic_data.RData` cannot support the suite — too large (§2.8), stores a
  fully-fitted object so tests read back their own output, and has only one
  shape, so no degenerate case is reachable. Proposes `F-TINY` / `F-SMALL` /
  `F-DEGEN` / `F-NULL` / `F-DERIV`.
- **§3 C++ tests** — data loader (dense int / dense double / sparse as each
  other's oracle, the `Flag::na` and all-zero-column paths), `numDeriv` gradient
  and Hessian checks for all 7 families, constrained-Newton behaviour, `opt_x` /
  `opt_yz`, `gamma_rate` vs `log_gamma_rate`, and external-pointer hygiene.
- **§5 Error-message tests** — 28 of them, one per malformed-input path.
- **§7 A regression table** mapping every defect in §1, §2 and §4 of *this*
  document to the test that must fail before its fix and pass after. That table
  is the acceptance criterion for §6 below.
- **§9 Sixteen open questions**, each blocking at least one test. Four of them
  (`Q-POST-1`, `Q-TSTAT-1`, `Q-REP-1`, `Q-GAM-1`) are about intended behaviour
  and only Kevin can answer them.

Four defects were found while writing that plan which are **not** listed in §1 of
this document and probably should be:

1. **Four of seven families are unusable at their documented defaults.**
   `opt_esvd.default`'s default `nuisance_vec = rep(NA, ncol(input_obj))`, but
   `gaussian`, `curved_gaussian`, `neg_binom` and `neg_binom2` all consume
   `gamma` — so `objfn_all_r` returns `NA` and `opt_esvd.default(family =
   "gaussian")` dies with `Error: missing value where TRUE/FALSE needed`, naming
   neither the family nor the parameter — **[verified]**. (Test `T-OPT-05`.)
2. **`format_covariates()` drops the *first* factor level; its roxygen says the
   last** — **[verified]**: `factor(c("a","a","b","b","c","c"))` yields columns
   `g_b`, `g_c`. (Test `T-FMT-02`.)
3. **`format_covariates()` rescales only the variables named in
   `rescale_numeric_variables`; its roxygen says it "rescales all the numerical
   variables"** — [inspection]. (Test `T-FMT-06`.)
4. **`data_loader()` returns a null external pointer without erroring** for an S4
   that is not `dgCMatrix` (e.g. a `dgeMatrix` — plausible user input) or a dense
   matrix that is neither integer nor numeric. The `Rcpp::stop("unsupported
   matrix type")` is only reached for a non-S4 non-matrix. The failure surfaces
   later as `Error: external pointer is not valid` — **[verified]**. (Tests
   `T-CPP-LOAD-05`, `T-VAL-22`.)

And one claim in §4.1 below is **overstated**: a serialized-and-restored `XPtr`
does *not* segfault under the current Rcpp. `saveRDS`/`readRDS` round-tripping
`esvd_family()` or `data_loader()` output and then calling `feas_Xi_r()` or
`objfn_all_r()` gives a clean R error, `Error: external pointer is not valid`,
from Rcpp's checked `XPtr(SEXP)` constructor — **[verified]**. The guard is still
worth adding (the message is unhelpful, and we would then own the guarantee
rather than inheriting it from a dependency), but this is not a segfault risk.
The live pointer risk is different and is `T-CPP-PTR-04`: `DenseDataLoader` holds
an `Eigen::Ref` to the **R matrix's own memory**, so a loader outliving its R
matrix dangles.


## 4. C++ backend (`src/`)

Restricting to correctness; the paper already flags optimization speed as future
work and that is explicitly out of scope right now.

### 4.1 External pointers are not null-checked — **[inspection]**

`src/esvd_family.cpp` returns `Rcpp::XPtr<Distribution>(distr, true)`, and
`src/data_loader.cpp` similarly. Consumers do
`Rcpp::XPtr<Distribution> distr = family["internal"];` and dereference without
checking. An `XPtr` that has been serialized and restored (`saveRDS`/`readRDS`,
`save`/`load`, or a `parallel` worker) has a **null** address, and dereferencing
it **segfaults R**.

This is reachable today because §2.3 exports `objfn_Xi_r` and friends, and
because `eSVD()` has an `intermediate_save` argument that `save()`s the object
mid-pipeline. The saved object does not currently *contain* a family/loader
pointer, so I have not demonstrated a crash — but the exported API makes one easy
to construct, and CRAN treats segfaults as a hard rejection.

**Fix.** A guard at every entry point:

```cpp
if (R_ExternalPtrAddr(family["internal"]) == NULL)
    Rcpp::stop("eSVD2: invalid family pointer (was this object restored from disk?)");
```

**Test.** `saveRDS()` an `esvd_family()` result, `readRDS()` it, call
`objfn_all_r()`, and `expect_error()` — not a crash. See `UNIT_TEST_PLAN.md`
`T-CPP-PTR-01..04`. **Note the correction recorded in §3:** this already gives a
clean `Error: external pointer is not valid` rather than a segfault
— **[verified]** — so the guard is about the *message*, not about crash safety.
The live pointer risk is `T-CPP-PTR-04` instead.

### 4.2 `Rcpp::warning()` inside a C++ frame holding Eigen objects — **[inspection]**

`src/constrained_newton.cpp:52`:

```cpp
Rcpp::warning("line search failed, returning the initial x");
```

`Rcpp::warning` routes to `Rf_warning`, which under `options(warn = 2)` becomes
an error and `longjmp`s — past the destructors of the `MatrixXd`/`VectorXd` and
`NumericVector` objects live on that stack. This is the classic R-longjmp/C++-RAII
hazard. I have not reproduced a leak or crash, so this is a code-review concern
rather than a demonstrated defect; but line-search failure is a real, reachable
condition, and CRAN reviewers do look for this pattern.

**Fix.** Set a status flag and return it to R, letting the R wrapper issue the
warning. Alternatively use `Rcpp::Rcerr` for the diagnostic.

Related: line-search failure is currently *only* signalled by this warning. The
returned `List` carries `step = 0.0` — so `opt_x`/`opt_yz` can silently return an
unchanged iterate. `opt_esvd.default`'s convergence test then sees no change and
`break`s, reporting convergence. **A failed optimization is currently
indistinguishable from a converged one.** Surfacing the status flag fixes both
problems at once.

### 4.3 `gamma_rate()`'s upper-bound search — **[inspection]**

`src/gamma_rate.cpp`, the bracket search:

```cpp
for(int i = 0; i < max_try; i++) {
    const double new_ub = gamma * ub;
    std::pair<double,double> new_dvals = deriv(new_ub);
    if(new_dvals.second <= 0.0) break;   // keeps the OLD ub
    ub = new_ub;
}
```

The comment above it says to shrink *"until `[l(ub)]''>0` and `[l(gamma*ub)]''<0`"*,
and on `break` the old `ub` is retained — so the intent may well be exactly this.
But if the loop exhausts `max_try` without ever satisfying the condition, `ub` is
used anyway with no signal, and `newton_raphson_iterate` is handed a bracket that
may not bracket a root; Boost then returns a bound rather than a root. The
function returns that value as if it were an MLE. `.nuisance_in_sequence()` in R
accepts any finite positive result without further checking.

**Fix.** Return a convergence status alongside the estimate, and have
`.nuisance_in_sequence()` fall through to `log_gamma_rate` on non-convergence
(it already has that fallback — it just never learns it is needed). Add a comment
resolving the bracket-direction ambiguity so the next reader does not have to
re-derive it.

**Test.** `test_gamma_rate.R` currently checks one well-behaved configuration.
Add: all-zero counts for a gene; a gene with a single non-zero count; very large
`mu` with tiny `s`; and assert `gamma_rate` and `exp(log_gamma_rate)` agree to
~1e-4 across a grid (they estimate the same quantity by different routes, so they
are each other's oracle — a second free equivalence test, like §1.5). See
`UNIT_TEST_PLAN.md` `T-CPP-GAM-01..08`; note the file's `T-CPP-GAM-02` proposes
reusing the R reference implementation already sitting in the comment block at
`src/gamma_rate.cpp:58–75` as an independent oracle, and `T-CPP-GAM-08` flags
that a `mu` shorter than `x` is an out-of-bounds read rather than an error.

> **Afterword, 2026-09-29.** The bracket was fixed in session 8, so the
> search can now return a rate above the library size. That exposed the
> likelihood's boundary: for a gene that is no more variable than Poisson
> around the fit, no finite rate maximizes it, and the rate ran to about 1e7.
> This inflated the test (§0.1). Version 1.2.0 caps the rate in R
> (`.apply_nuisance_cap()`) and leaves `gamma_rate` an MLE. See
> `OVERDISPERSION_BRAINSTORM.md`, "Decision".

### 4.4 Packaging hygiene — **[policy]**

- `src/*.o`, `src/*.so` were sitting in the working tree. Now excluded by
  `.Rbuildignore` (and `R CMD build` reports `cleaning src`), and by
  `.gitignore`. **Done this session.**
- No `src/Makevars`. Not required — `LinkingTo: Rcpp, RcppEigen, BH` is
  sufficient and R defaults to C++17 — but confirm no compiler warnings under
  `-Wall -pedantic`, which CRAN uses.
- Symbol registration is present and correct (`R_init_eSVD2` with
  `R_useDynDynamicSymbols(dll, FALSE)`); no action needed.
- Run the check with `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` and, for the C++,
  under valgrind or at minimum `-fsanitize=address,undefined` once before
  submission. CRAN runs both.

---

## 5. Documentation and packaging

### 5.1 `eSVD()` needs a full roxygen block
See §2.3. It is the entry point; it currently has none. Include `@examples`
using `generate_null()` so the example runs at check time without downloads.

### 5.2 `@return` is missing on several documented functions — **[policy]**
`opt_x`, `opt_yz` (`R/optimization_helper.R`), `data_loader`, `esvd_family` have
roxygen blocks with no `@return`. `R CMD check` flags missing `\value{}` as a
WARNING. Note `data_loader`'s existing `@returns description` is a placeholder.

### 5.3 Replace `print()` with `message()` package-wide — **[policy]**
`print()` and `cat()` for progress appear throughout (`eSVD.R`,
`initialization.R`, `nuisance.R`, `compute_test_per_gene.R`,
`compute_test_statistic.R`, `reparameterization.R`). CRAN requires diagnostics on
`stderr` and suppressible — `message()`. This also fixes §1.2. `cat('*')` progress
markers in `.initialize_coefficient()` and `estimate_nuisance.default()` should go
behind `message()` too.

### 5.4 Document the rate/scale inversion
Every roxygen mention of "nuisance" or "over-dispersion" should state that
`nuisance_vec` is the Gamma **rate** `β = 1/γ`, the reciprocal of the paper's
`γ_j`, and that **larger** means **less** overdispersed. Right now a reader
comparing `?compute_posterior` to Eq. 14 concludes the code is wrong.

Also state the three modelling assumptions the paper lists (Gamma–Poisson counts;
covariate effects removable in a GLM; DE = difference in *means* between
individuals) somewhere a user will find them — the `eSVD()` help page.

### 5.5 testthat modernization
- `context()` is deprecated in testthat 3e; remove from all 12 files.
- Add `Config/testthat/edition: 3` to `DESCRIPTION`.
- Replace `load("../assets/...")` with `testthat::test_path()`; better, move
  fixtures to `tests/testthat/fixtures/`. The current form depends on the working
  directory and breaks outside `test_check()`.
- `tests/testthat/_snaps/` is empty and was dropped by `R CMD build`
  ("Removed empty directory"). Either use snapshot tests or delete the directory.
- `tests/assets/data_generation.R` is the provenance script for the fixtures —
  good practice, keep it, but it uses `devtools::session_info()`, which is why
  `devtools` is in `Suggests` (§2.6). It is not run at check time, so the
  dependency can go.

### 5.6 README
- The "Installation" section should change from `devtools::install_github()` to
  `install.packages("eSVD2")` once accepted.
- Drop the `Rmpfr` install caveat once §2.2 lands.
- It says "version `1.0.0` as of June 16, 2024" while `DESCRIPTION` says
  `1.0.1.07` — reconcile.

---

## 6. Suggested order of work

Each step is independently verifiable, and later steps depend on earlier ones.

0. **Clear the mechanical blockers in one sitting** — they are independent of
   every design question and each takes minutes: shrink the test fixture (§2.8,
   the 5 MB limit), rewrite `LICENSE` as a 2-line stub (§2.4), add `override`
   ×5 (§2.5), add `MASS` to `Suggests` (§2.9). Re-run the check and confirm the
   tarball is under 5 MB before doing anything else.
1. **Decide CRAN vs. Bioconductor** (§2.1). Everything else follows from this.
   Assuming CRAN: remove `sparseMatrixStats`, remove `Rmpfr` (§2.2), remove
   `biocViews`, move `Rcpp` to `Imports`.
2. **Fix §1.1** (the `log.p` transform) and add the test. This is the one defect
   that changes published-style results, and it is a two-line fix.
3. **Fix §1.2, §1.4, §1.8** — small, local, each with a test.
4. **Kill `exportPattern`** (§2.3); write `eSVD()`'s documentation (§5.1) and
   the missing `@return`s (§5.2). Re-run `devtools::document()`.
5. **Strengthen the two-pipeline equivalence test** (§1.5) and align the
   `library_min` defaults. This becomes the safety net for everything after.
6. **Add the numerical-gradient tests** (§3.2). These make the C++ safe to touch.
7. **Fix §4.1, §4.2, §4.3** — the C++ status-reporting work, now guarded by (6).
8. **Fix §1.3, §1.6, §1.7** — the larger R refactors, now guarded by (5).
9. **Fill the coverage gaps** — now enumerated in `UNIT_TEST_PLAN.md` §2, §4,
   §5 and §6. Add `covr::package_coverage()` to CI and set a floor (C-09).
   `UNIT_TEST_PLAN.md` §10 gives a test-first ordering that differs from this
   one: the harness and fixtures land first, then the two *free* equivalence
   tests (matrix-vs-per-gene, `gamma_rate`-vs-`log_gamma_rate`) as a safety net,
   then the gradient tests, and only then the fixes above.
10. **Vignettes** (§2.7), README (§5.6), `DESCRIPTION` polish (§2.6),
    `cran-comments.md`.
11. **`R CMD check --as-cran`** clean on macOS + Linux + Windows
    (`devtools::check_win_devel()`, `rhub::check_for_cran()`), then submit.

---

## Appendix: `R CMD check --as-cran` results

### 2026-09-29 (version 1.2.0, commit `1a0a536`) — current

R 4.5.1, aarch64-apple-darwin20, macOS Sonoma 14.2.1, Apple clang 15.0.0 for
the package's C++. RStudio's pandoc prepended to `PATH`. Built and checked
outside the repository.

```
R CMD build eSVD2                       # vignette built; tarball 1,556,652 bytes
R CMD check --as-cran eSVD2_1.2.0.tar.gz
```

**Result: `Status: 2 NOTEs`.**

```
* checking CRAN incoming feasibility ... [3s/13s] NOTE
New submission
* checking tests ...
  Running 'testthat.R' [11s/11s]    [ FAIL 0 | WARN 0 | SKIP 20 | PASS 1512 ]
 [11s/11s] OK
* checking re-building of vignette outputs ... [64s/65s] OK
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
Skipping checking math rendering: package 'V8' unavailable
```

Everything else is `OK`, including `examples with --run-donttest` (slowest
`eSVD` 0.54 s, then `plot_fitted_vs_observed` 0.45 s). The run was without
`_R_CHECK_FORCE_SUGGESTS_`, so every `Suggests` package (now including
`ggplot2` and `ggrepel`) is installed here.

### 2026-09-02 (session 10)

R 4.5.1, aarch64-apple-darwin20, macOS Sonoma 14.2.1, Apple clang 15.0.0.
pandoc from RStudio's bundle (`/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64`) prepended to `PATH`.

```
R CMD build eSVD2                       # vignette built; tarball 1,020,379 bytes
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran eSVD2_1.0.2.tar.gz
```

**Result: `Status: 2 NOTEs`.**

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kevin Z Lin <kzlin@uw.edu>'
New submission
* checking tests ... Running 'testthat.R' OK      [ FAIL 0 | WARN 0 | SKIP 3 | PASS 685 ]
* checking re-building of vignette outputs ... [63s/64s] OK
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
```

Everything else `OK`, including `examples with --run-donttest` (slowest
example `eSVD`, 0.55 s), `installed package size`, `compiled code`,
`Rd \usage sections`, `for unstated dependencies in vignettes`. The
`_R_CHECK_FORCE_SUGGESTS_=false` flag is now moot: every `Suggests`
package is on CRAN and installed here.

### 2026-08-27 — baseline

Recorded **2026-08-27** on R 4.5.1, aarch64-apple-darwin20, macOS Sonoma 14.2.1,
Apple clang 15.0.0. Command:

```
R CMD build --no-build-vignettes eSVD2
R CMD check --as-cran --no-manual --no-vignettes --no-build-vignettes eSVD2_1.0.1.07.tar.gz
```

**Result: `Status: 5 WARNINGs, 3 NOTEs`.**

| Check | Result | Covered in |
|---|---|---|
| CRAN incoming feasibility | NOTE | §2.4, §2.6, §2.8 |
| whether package can be installed | WARNING | §2.5 |
| DESCRIPTION meta-information | NOTE (`License stub is invalid DCF`) | §2.4 |
| top-level files | NOTE (`pandoc` not installed) | environmental — ignore |
| missing documentation entries | WARNING (14 objects) | §2.3 |
| unstated dependencies in `tests` | WARNING (`MASS`) | §2.9 |
| files in `vignettes` / package vignettes | WARNING ×2 | see caveat below |
| **tests** | **OK** — the existing suite passes | — |
| R code for possible problems | OK | — |
| dependencies in R code | OK | — |
| compiled code | OK | — |

The incoming-feasibility NOTE reads, in full:

```
New submission
Version contains leading zeroes (1.0.1.07)
License components with restrictions and base license permitting such: MIT + file LICENSE
Package has a VignetteBuilder field but no prebuilt vignette index.
The Title field starts with the package name.
The Title field should be in title case.
The Description field should not start with the package name, 'This package' or similar.
Size of tarball: 13111968 bytes
```

**Three caveats — this is a floor, not a ceiling.**

1. **The two vignette WARNINGs and the "no prebuilt vignette index" line are
   artifacts of my `--no-build-vignettes` / `--no-vignettes` flags**, not
   independent defects. A real submission builds vignettes, at which point §2.7
   (Bioconductor packages loaded unconditionally) becomes the live failure
   instead. Do not treat these two as fixed when they disappear.
2. **The dependency policy violations did not fire.** `sparseMatrixStats` and
   `EnhancedVolcano` are installed on this machine and I passed
   `_R_CHECK_FORCE_SUGGESTS_=false`, so `checking package dependencies` reported
   `OK`. On CRAN's machines it will not. §2.1 stands regardless of this output.
3. **Static analysis was clean** (`checking R code for possible problems ... OK`,
   `checking dependencies in R code ... OK`, `checking compiled code ... OK`).
   Every defect in §1 is *semantic* — `R CMD check` cannot see any of them. That
   is the argument for §3: the check being green will never mean the numbers are
   right.

Re-run this command after each step of §6 and update the table above.
