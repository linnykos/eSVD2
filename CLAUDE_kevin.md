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
| `OVERDISPERSION_WIKI` | personal laptop (macOS), Dropbox | `~/Library/CloudStorage/Dropbox/Collaboration-and-People/amywatt/git/overdispersion_wiki` | Read only from this project. Start at `wiki/index.md`; pages are `wiki/pages/<slug>.md` |
| `WAS2CODE_REPO` | personal laptop (macOS), Dropbox | `~/Library/CloudStorage/Dropbox/Collaboration-and-People/archive/tati/git/Was2CODE` | Under `archive/`, so treat as frozen. `R/esvd_helper.R` (60 lines) is the file being imported into `eSVD2`; copy it in rather than depending on this path |

Name the machine specifically enough that another collaborator can tell whether it is reachable to them. **A missing row means unknown; only an explicit *(not present)* row means known-absent.**

## Project Status (as of 2026-09-30)

**Goal: get `eSVD2` onto CRAN.** Correctness first; efficiency is explicitly out
of scope for now.

**Version 1.2.0 (the cap on the nuisance rate) is vetted and committed as
`1a0a536` on `devel`.** Session 18 brought every document under
`additional_context/` up to 1.2.0, and session 19 corrected about 35
misleading comments, roxygen blocks and two message strings in `R/` and
`src/` (no behaviour change; `man/` regenerated). Session 20 changed
`estimate_nuisance.eSVD`'s default to `bool_covariates_as_library = TRUE` and
made both posterior paths and `eSVD()` refuse the same nonsensical settings
through `.check_posterior_args()`. All three are committed as `b0dd8f2`.

**No file in `R/` or `tests/testthat/` carries the `_claude` suffix any
more** (session 21, uncommitted). Six drafted test files had the name of an
older test file, so each was appended to it; the suite is now 28 files. The
scripts and reports under `additional_context/` keep the suffix.

**`R CMD check --as-cran` after the rename (2026-09-30): `Status: 3 NOTEs`**
(`New submission`, HTML Tidy on this machine, and "unable to verify current
time", which is the check failing to reach a time server), tests
`[ FAIL 0 | WARN 0 | SKIP 20 | PASS 1608 ]`. The suite under `NOT_CRAN=true`
gives **2110 pass / 0 fail / 0 warnings / 0 skip**: 298 `test_that` blocks in
28 files.

**What stands between the package and a submission**:

- one real data set under the cap;
- rewording `?compute_log_fold_change`'s caveat;
- Q-LFC-2 and Q-LFC-3;
- `Authors@R`;
- the misspelled argument names;
- confirming the ASD articles and the rank-deficiency refusal;
- a Windows/Linux check, one sanitizer run, and a pkgdown rebuild.

`CRAN_READINESS.md` §0.3 is the list.

## Key Methodological Details

- **The cap is applied in R, in one place, `.apply_nuisance_cap()`**, used by
  `estimate_nuisance()` and by `recompute_pvalue()`; `gamma_rate` stays a
  maximum-likelihood estimator. Status precedence: `failed`, then
  `boundary` (`D_j <= 0`), then `capped`, then `estimated`.
  `nuisance_num_capped` counts every gene the cap replaced, boundary genes
  included.
- **The rule for which covariate columns form the library size lives in one
  place, `.nuisance_library_idx()`** (in `R/nuisance_cap.R`), used by
  `estimate_nuisance.eSVD()`, `compute_posterior.default()`,
  `compute_test_per_gene()` and `plot_fitted_vs_observed()` since session
  17. The eight column sets it must give on F-TINY are written out by hand
  in `.library_column_oracle()` (`helper-fixtures.R`), which T-CAP-10 and
  T-POST-13 check against.
- **`.estimate_nuisance_matrix()` walks the genes once**: the per-gene
  `vapply` returns the rate and the boundary statistic `D_j` together, and
  `.compute_boundary_statistic()` is a per-gene helper (`x_vec, mu_vec,
  s_vec`). T-CAP-11 pins dense against sparse counts through that loop.
- **A boundary gene is set to the cap, not to `min(optimizer's value,
  cap)` (Kevin confirmed, 2026-09-29).** The optimizer's value is arbitrary
  (about 1e7 on one route, exactly `exp(10)` on the log route) and is below
  the cap once the median library size is above about 2200, as with
  `bool_library_includes_interept = FALSE`. `nuisance_mle_vec` keeps the
  optimizer's value.
- **`min_val` is in the units of the cap (Kevin, 2026-09-29): the floor is
  `min_val * median_i s_ji`, so no gene is above its cap** and on the
  unit-free scale of `plot_nuisance()` every rate lies in `[min_val, c]`.
  `.check_min_val()` refuses `min_val >= cap_multiplier` in both
  `estimate_nuisance()` methods and in `recompute_pvalue()` (which reads
  the recorded `nuisance_min_val`). A failed gene gets the floor. The
  earlier absolute floor held a gene with a degenerate fit (median library
  size about 1e-8) thousands of times above its cap; `plot_nuisance()` no
  longer needs to count such genes.
- **`bool_diet = TRUE` drops `covariates` (Kevin, 2026-09-29: the diet is
  to be as small as possible)**, so `recompute_pvalue()` and
  `plot_fitted_vs_observed()` keep rebuilding the design from the Seurat
  object and checking it against the stored summaries.
- **The posterior's settings are validated in one place,
  `.check_posterior_args()` (`R/posterior.R`)**, called by
  `compute_posterior.default`, `compute_test_per_gene` and, on entry,
  `eSVD()`. `nuisance_lower_quantile = NULL` means no floor on both paths.
  `eSVD()` has no `nuisance_lower_quantile` argument and passes `NULL` to the
  check.
- **The cap is tighter than `c` in the posterior for any gene whose library
  is below `library_min = 0.1`**, because `compute_posterior()` floors the
  library there.
- **`recompute_pvalue()` replays `param`**, so every stage it repeats must
  overwrite its own `param` entries on a rerun. `estimate_nuisance()`,
  `compute_posterior()`, `compute_test_statistic()` and
  `compute_test_per_gene()` now do; `initialize_esvd()`, `opt_esvd()` and
  the reparameterization still go through `.combine_two_named_lists()`
  alone and go stale.
- **A diet object is redone from the Seurat object**, subset to the fitted
  cells and analyzed genes *before* the covariates are built (`Log_UMI` and
  the rescaled numeric covariates depend on which cells and genes are
  present), through `.prepare_esvd_covariates()`, the block `eSVD()` itself
  calls. The rebuilt matrices are recognized by summaries stored at the fit
  (gene means, column sums, and sums weighted by `cos(sqrt(2) * i)` so that
  an exchange between cells is noticed). Only an object built by `eSVD()` or
  `eSVD_helper()` of 1.2.0 can be rebuilt; a stage-by-stage object must
  still carry its `dat` and `covariates`.
- **The SD in `plot_fitted_vs_observed()` is the model's marginal one**,
  `sqrt(m * (1 + s / beta))` with the capped `nuisance_vec` and the library
  of `estimate_nuisance()`, not the posterior's rescaled rate. On counts
  drawn from the model at their true rates 1 to 2% of pairs are beyond 3 SD
  (a count is right-skewed); a rate ten times too large gives about 5%, ten
  times too small under 0.1%.
- **Counts below the fit are almost never marked**, because the fit minus
  3 SD is negative for most pairs; the red points of the plot are counts
  above the fit.
- **`alpha_max` changed no statistic on F-TINY** at 1, 5, 86 or 1000, and
  `library_min = 2` none either (every library size is above 2.4). A test
  of the posterior's settings on that fixture must use `library_min >= 20`
  or `pseudocount`.

- **Master vs 1.2.0 (version comparison rerun, session 18).**
  - With the cap, false discoveries are back at master's level, within the
    noise of the counts (SE about 0.5): null 0.2 vs 0.2 per 300 genes, and
    `near_poisson` 1.8 vs 2.6 (type-I 0.065 vs 0.067).
  - Welch r with master is at least 0.986, and `devel_swap` still reproduces
    master.
  - 95% of genes are within twofold of the true rate, against master's 59%.
  - The one visible cost: `generate_null()` finds 8.1 of 10 DE genes, against
    master's 8.8. The uncapped 1.1.0 finds 8.2, so the loss comes from
    estimating the rate, not from the cap.
  - Under the cap, ±2 SE covers the true log2FC for at least 99% of genes of
    every status in five of six regimes; `strong_de` is at 0.72 to 0.85, from
    the `Log_UMI` depth bias.
- **1.2.0 run end to end equals the brainstorm's prototype `cap_10s`**: the
  same FDR calls in all 100 data sets, and rates within 1% except two
  `generate_null` genes that the old absolute floor had lifted
  (`09_run_v120_claude.R`).
- **`cap_multiplier = Inf` reproduces 1.1.0's statistics to every stored
  digit**, so `devel_nocap` in the comparison stands in for 1.1.0.
- **Master is not bit-reproducible on `generate_null()`**: two runs on
  identical data differ by up to 0.006 in the Welch statistic, from its
  random SVD start.
- **Libraries under `version_comparison/lib/`**:
  - `master` (3d5f7bf);
  - `devel_1.1.0` (d49e402), which the brainstorm's scripts 01 to 07 load so
    that their cached 1.1.0 fits stay valid;
  - `devel` (the working tree, 1.2.0).
- **In the bullets below that compare "master" and "devel", "devel" is
  version 1.1.0, the uncapped estimate.** They are the evidence for the
  cap.
- **Master (3d5f7bf) vs devel differ on ordinary data only through the
  nuisance estimate.** logFC agrees at r >= 0.998 in all six simulated
  regimes. Overwriting devel's `nuisance_vec` with master's right after
  `estimate_nuisance()` reproduces master's Welch statistics to within 0.006
  on every data set; no other change (QR reparameterization, GLM-fallback
  offset, SVD start, new empirical-null estimators) reaches the statistic.
- **Master's rate is `min(MLE, max library size of the gene)`**, and 78% to
  100% of genes sit at that cap in the model-generated regimes (39% under
  `generate_null()`), which is why its rates have Spearman about 0 with the
  truth. `pmin(devel rate, max_i s_ji)` reproduces master's rates to a
  median of 0.1% and its Welch statistics to within 0.025.
- **Devel's MLE is accurate away from the boundary; the "two- to threefold
  overestimate" was a units mismatch.** The rate shares units with the
  fitted library size, which contains the gene intercept; the generator's
  library size averages 1. Divided by the gene's median fitted library size,
  the estimate is 4% to 16% above the truth (39% at low counts), Spearman
  0.77. The comparison report now makes the comparison in the fit's units.
- **A rate of about 1e7 is the boundary of the likelihood, not a solver
  failure.** Around the Poisson limit the log-likelihood is
  `l_Poisson + D / (2 * beta)` with `D = sum_i [(A_i - m_i)^2 - A_i] / mu_i`,
  `m = mu * s`; a gene is at the boundary when `D <= 0` (agrees with devel's
  output for 11,997 of 12,000 genes). Fractions: 0.5% ordinary, 2.5% at mean
  counts of 0.1 to 0.6, 12% with a trend, 16% when true rates span 0.5 to
  200, 48% near-Poisson.
- **The true rates miscalibrate the test only for genes that are nearly
  Poisson in truth.** With true rates of 0.5 to 200 the oracle gives 5.7
  false discoveries per 300 genes, 55 of 57 from genes whose true rate is
  above 30 times the library size. With true rates below about 50 (the
  `trend` regime) the oracle is the best choice (1.5 false, 24.3 of 30 true).
  The null SD of the Gaussianized statistic is about 0.3 and grows with the
  rate (0.69 in the top fifth of devel's rates).
- **`bool_stabilize_underdispersion` discards the absolute scale of the
  rates whenever their geometric mean exceeds 1**: every eSVD-model regime
  under devel, but not `generate_null()` (0.7), nor `low_count` (0.5) and
  `generate_null()` (0.4) under master.
- **When the rate follows expression, the cap beats master on the test as
  well as on the estimate.** In the `trend` regime (DESeq2's
  `alpha0 + alpha1 / m`, alpha0 0.1, alpha1 0.02): cap 1.0 false and 21.7
  true discoveries, master 2.3 and 20.7; paired differences -1.3 (SE 0.6)
  and +1.0 (SE 0.45). Master's rate is flat in expression where the truth
  falls twentyfold. The cap's cost is on lowly expressed genes whose true
  rate is near 20: 3.0 true discoveries of 10.7 against the oracle's 5.9.
- **Shrinkage toward a fitted trend gives the best rates when a trend
  exists** (97% within twofold, 24.3 true discoveries, 1.6 false), run at
  600 cells only.
- **A cap `rate <= c * median_i s_ji` is flat in false discoveries for c
  from 1 to 10 and rises from c = 20.** At c = 10 about 11% of genes are
  capped in ordinary data. A cap at the technical-noise level of the
  literature (c of 50 to 500 here) is too loose: 1.2 to 2.9 false
  discoveries under the null.
- **Bounds built from the gene's likelihood weaken as cells are added; a cap
  does not.** At 600 cells the empirical-Bayes MAP rate acts as a soft cap
  near 10 times the library size and matches the hard cap. At 3000 cells
  with true rates spanning 0.5 to 200 it leaves 6.2 false discoveries per
  300 genes against the cap's 3.0 (profile lower bound 5.7, oracle 10.1).
- **As an estimate, both the cap at 10 and the shrinkage beat master
  wherever the rate is identifiable**: 95% and 99% of genes within twofold
  of the truth against master's 60%, Spearman 0.77 against 0.05. In
  near-Poisson data no candidate is within a factor of 8 or ranks the genes.
  Shrinkage followed by the cap equals the cap at 3000 cells.
- **DESeq2's trend centre and `trigamma((m - p)/2)` sampling variance change
  nothing in these simulations**, which have no trend by construction. The
  trigamma value is 0.003 at 600 cells against a measured 0.10 to 2.1.
- **With 3000 cells and true rates of 2 to 8 no gene is at the boundary**
  and the uncapped MLE is calibrated; the ordinary-regime inflation is a
  small-sample effect, the wide-rate inflation is not.
- **A diverged rate is almost always a false discovery**: 74-100% of null
  genes with a diverged rate reach FDR < 0.05, against under 1% of the
  others. Per 300 null genes, devel gives 2.7 false discoveries under the
  global null against master's 0.2. In the near-Poisson regime devel gives
  63 against 2.6, with type-I error at p < 0.05 of 0.28 against 0.07.
- **The t-to-z log-scale fix matters only because of the `gamma_rate`
  change**: master's largest Welch statistic in simulation was 26; devel's
  is 53, and 10 of devel's statistics would be `Inf` under master's
  `qnorm(pt())`.
- **Master called all 30 of 30 genes significant, silently**, on a 30-gene
  cohort (its unconstrained fallback null estimator). Devel called none.
- **NEWS's "zeroes NA counts in sparse matrices" never takes effect via
  `format_covariates()`**: `Log_UMI = log(rowSums(dat))` is already `NA`,
  and `initialize_esvd()` stops on the covariates first.

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
- **Where the nuisance estimate diverges the SE is anti-conservative**
  (measured without a cap). About 12% of genes on the toy cohorts get a rate
  near 1e7 (true rates 0.1 to 10); the posterior collapses onto the fit, the
  within term vanishes, coverage falls to 0.36, and the Welch statistic is
  in the tens. The diverged genes also inflate the geometric mean that
  `bool_stabilize_underdispersion` divides every gene's rate by. Under the
  cap the coverage is restored (see the master vs 1.2.0 bullet).
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
  SE divides by their lengths.
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
  8): Poisson-like genes get enormous rates, so
  `bool_stabilize_underdispersion` rescales every gene by `10^-mean(log10)`.
  The cap of 1.2.0 is Kevin's answer to it (Q10 in the readiness doc).
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
- The public API is the 21 `export()` lines in `NAMESPACE`; `opt_x`,
  `opt_yz`, `data_loader`, `esvd_family` are documented `@keywords internal`.
- `devtools::document()` (roxygen2 8.1.0 here) rewrites `RoxygenNote:
  7.3.3` to `Config/roxygen2/version: 8.1.0`; reverted in sessions 9 and 10
  (Q11 in the readiness doc).
- `R CMD build` with the vignette needs pandoc; RStudio's copy at
  `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64`
  prepended to `PATH` works. Run with `_R_CHECK_FORCE_SUGGESTS_=false`
  no longer needed.

## Open Questions / Next Steps

**Decisions for Kevin, newest first.** Q-LFC-2 and Q-LFC-3 are in
`UNIT_TEST_PLAN.md` §2.18; the full list is `CRAN_READINESS.md` §0.3.

0. **Review and commit the rename of session 21.** Also decide whether the
   26 `*_claude.*` files under `additional_context/` lose the suffix too
   (their scripts call each other by name, and `.gitignore` names
   `01_cache_fits_claude.R`), and whether the `_claude` convention in the
   master `CLAUDE.md` still holds for new files.
0f. **Run one real data set under 1.2.0** (brainstorm Idea 11). It needs a
    `PAPER_DATA` path. This is the only check the default `c = 10` has not had.
0e. `plot_fitted_vs_observed()` resets the caller's random stream
    (`seed_number = 10`). Keep it, default to `NULL`, or restore the stream on
    exit?
1. **Q-LFC-1**: coverage under the cap has been measured (see Key details).
   What remains is to reword the caveat in `?compute_log_fold_change`, which
   still describes the uncapped case.
2. **Q-LFC-2**: keep the warning from `compute_test_statistic.default()` on a
   non-positive arm mean, or return `NA` silently from the matrix method?
3. **Q-LFC-3**: the names `log2fc_vec` / `log2fc_se_vec` and `logFC_se`.
4. Is the mixture variance (per-cell posterior variance, not divided by the
   cells per individual) the intended variance for a reported SE?
5. `Authors@R`: add Yixuan Qiu and Kathryn Roeder as `aut`? Also use one
   spelling of Kevin's name.
6. Rename `bool_library_includes_interept` / `library_multipler` now or
   later? The recommendation is later, with a deprecation shim.
7. Confirm the ASD tutorials as `vignettes/articles/` pkgdown articles.
8. Confirm refusing rank-deficient covariates at `initialize_esvd()`.
9. Keep the aggregated line-search warning, or downgrade it to `verbose`?
10. `multtest()`'s fallback warning on fewer than about 40 genes: keep it as
    is?
11. `generate_null()` gene names: `gene_1` vs Seurat's `gene-1`.
12. The regenerated 400 x 60 legacy fixture, or port the legacy tests to the
    built fixtures?
13. `param` still goes stale on a rerun of `initialize_esvd`, `opt_esvd` or
    the reparameterization. Fix it in `.combine_two_named_lists()`?
14. T-PROP-06 counts `generate_null()`'s "null_large_var" genes as null; they
    are not.
15. Carried over: `.multtest_locfdr` catching warnings; `bool_diet` keeping
    the fit; the amplification of non-determinism.

**Ready to start, blocked on nothing:**

16. Vet `tests/testthat/test_compute_log_fold_change.R`, which was
    committed in `d49e402` without a recorded vet.
16b. Make `sparse_na` work as NEWS says, or reword the NEWS item.
17. Windows / Linux checks (`devtools::check_win_devel()`,
    `rhub::rhub_check()`) and one ASan/UBSan run.
18. The suggested tests in `CRAN_READINESS.md` §0.4, including T-LFC-20 and
    T-POST-14 (whether `alpha_max` ever binds).
19. Rebuild the pkgdown site (`docs/`) so it lists the four new functions.
20. Optional: `.compute_df()` could read the stored `case_var` /
    `control_var`; the `.split_individuals_by_arm()` refactor.
