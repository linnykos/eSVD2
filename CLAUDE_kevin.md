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

## Project Status (as of 2026-08-29)

**Goal: get `eSVD2` onto CRAN.** Correctness first; efficiency is explicitly out
of scope for now.

Seven sessions in. **Nothing in the package's `R/` or `src/` has been changed
yet** — by design. The test suite is now WRITTEN AND RUN against the unmodified
code: 15 new files, **388 pass / 87 fail / 29 skip**, with the existing 79-test
suite untouched and still green. Results and all findings are in
**`additional_context/TEST_RUN_REPORT.md`**; the `Rmpfr` decision is settled in
**`additional_context/RMPFR_REPORT.md`**. The two working documents are:

- **`additional_context/CRAN_READINESS.md`** — the audit. *What is wrong*, with
  every claim tagged `[verified]` / `[inspection]` / `[policy]`. **Now visibly
  behind**: §2.2 (`Rmpfr`) and §2.1 (`sparseMatrixStats`) are both settled more
  strongly than it states, and four findings from sessions 3–4 are not in it at
  all. Needs a pass.
- **`additional_context/UNIT_TEST_PLAN.md`** — the proposed test suite, ~236
  `test_that` blocks against the current ~40, each with an ID (`T-<AREA>-nn`) and
  an explicit **oracle**. **Reviewed four times (2026-08-29), and every question
  in it is now closed.** §9.1–§9.4 record the 30 resolved decisions with their
  consequences; §9.5 lists three defaults to confirm while implementing, none
  blocking. Three sections were added across the review passes: **§2.15** (does
  `Rmpfr` earn its place), **§2.16** (the `gene_status` output), **§2.17** (cohort
  filtering, imported from `WAS2CODE_REPO` and respecified). **The plan is ready
  to implement on Kevin's word.**
- **`additional_context/TEST_RUN_REPORT.md`** — the first run. Five findings in
  neither prior document, the predicted failures confirmed, and what passed.
- **`additional_context/RMPFR_REPORT.md`** — `Rmpfr` is not needed; drop it from
  `DESCRIPTION` outright rather than moving it to `Suggests`.

**Session-1 infrastructure (done):** `.gitignore` / `.Rbuildignore` rewritten so
the tarball contains only the package; large-file pre-commit hook installed and
`core.hooksPath` set; `.RData` untracked; `additional_context/summary.md` written
with the paper's equations mapped to the functions implementing them.

**Baseline check recorded**: `R CMD check --as-cran` gives **5 WARNINGs,
3 NOTEs** on R 4.5.1 / macOS. The existing test suite **passes** — which is the
point: `checking R code for possible problems` reports OK while every defect in
§1 of the readiness doc is live.

**The findings that matter most:**

1. **`qnorm(pt(t, df))` returns `+Inf` for `t = 40, df = 18`** — verified. One
   strongly up-regulated gene makes `locfdr::locfdr()` error, the `tryCatch` in
   `multtest()` swallows it, and the empirical null silently degrades to
   `.multtest_simple()` — the estimator whose own comment disclaims it. The user
   gets no signal because `compute_pvalue()` discards `multtest()`'s `method`
   field. Two-line fix (`log.p = TRUE`). Tests `T-PVAL-01/02`.
2. **Tarball is 13.1 MB against CRAN's 5 MB limit** — verified. Driven by
   `tests/assets/synthetic_data.RData`, which stores a fitted 8.7 MB `eSVD_obj`
   plus three ~2.3 MB generation by-products no test uses. Fixing this and
   building the fixtures the test plan needs is the same piece of work.
3. **Four of seven families are unusable at their documented defaults** —
   `opt_esvd.default` defaults `nuisance_vec = rep(NA, p)`, but `gaussian`,
   `curved_gaussian`, `neg_binom` and `neg_binom2` all consume `gamma`, so the
   objective is `NA` and the call dies with `missing value where TRUE/FALSE
   needed`. Test `T-OPT-05`. Not in the readiness doc.
4. **All-zero genes misbehave in three places at once** — `glmnet` warns
   `"Convergence for 1th lambda value not reached"` and returns the wrong
   lambda's fit; `gamma_rate` and `log_gamma_rate` disagree 21× (the latter
   saturated at its `-10` clamp); the Welch statistic is `0/0`. Removed by the
   new `gene_status` feature.
5. **Two one-commit landmines found in session 4**, neither in either working
   document. `NAMESPACE` has **no `export(eSVD)`** — the main user-facing
   function is exported only by the blanket `exportPattern("^[[:alpha:]]+")`,
   which `Q-CPP-1` removes, so that change would silently un-export `eSVD()`.
   And `eSVD()` calls `SeuratObject::LayerData()` with `SeuratObject` in
   `Suggests` and **no `requireNamespace()` guard anywhere in `R/`**, which CRAN
   requires.

**Decisions made: target CRAN, not Bioconductor.** `sparseMatrixStats` comes out
(reimplement `colSds` in `Matrix`, vendor the C++ only if the test fails).
`Rmpfr` leaves `DESCRIPTION` entirely — the evidence tests run once as a script
under `additional_context/`, not in `tests/`.

## Key Methodological Details

- **`nuisance_vec` is the Gamma *rate* `β = 1/γ`**, the reciprocal of the paper's
  overdispersion `γ_j`. Larger `nuisance_vec` = *less* overdispersion. The code is
  internally consistent; no roxygen block says this, so comparing
  `?compute_posterior` to Eq. 14 makes the code look wrong when it isn't.
  **Resolved (`Q-POST-1`): the code is right and the prose is wrong**, so
  `bool_stabilize_underdispersion`'s fix is a roxygen rewrite, not a code change.
- **Two independent implementations of the same computation exist**, and each is
  the other's free oracle — no new reference implementation needed:
  `compute_posterior` + `compute_test_statistic` + `compute_pvalue` (matrix path)
  versus `compute_test_per_gene` (fused per-gene path); and `gamma_rate` versus
  `log_gamma_rate` in `src/`. `UNIT_TEST_PLAN.md` §10 puts both at step 2 — the
  cheapest possible safety net for every later change.
- Their `library_min` defaults **disagreed** (0.1 vs. 1e-2). Resolved: **0.1
  everywhere**, so `compute_test_per_gene`'s default changes.
- **Why `gamma_rate` vs `log_gamma_rate` was commented out**: on an all-zero gene
  they return `9.706e-4` and `exp(-10) = 4.5e-5` — the second is its clamp
  boundary, not an estimate. Restricting the grid to non-degenerate inputs should
  make the equivalence test pass.
- **A third free oracle sits unused in a comment**: `src/gamma_rate.cpp` lines
  58–75 carry a complete R implementation of the objective, gradient and Hessian.
  `T-CPP-GAM-02` turns it into an independent check on both C++ routines.
- **`numDeriv` has been in `Suggests` since the beginning and never used.**
  Gradient/Hessian checks for the 7 families are the highest-value C++ tests
  available: an analytic-derivative error currently produces a silently
  mis-converged fit with *no* symptom.
- **A failed optimization is indistinguishable from a converged one.**
  `line_search` warns and returns `step = 0.0` with an unchanged iterate;
  `opt_esvd.default`'s convergence test then sees no change and `break`s,
  reporting convergence. `T-CPP-NEWT-02`.
- **`Rmpfr` cannot help, settled empirically.** Doubles in log space match a
  200-bit MPFR computation to relative ≤7.3e-17 across `z ∈ [-1e4, -5]`;
  `2*Rmpfr::pnorm(-40)` underflows to exactly 0 just as `stats::pnorm` does; and
  `stats::p.adjust` silently coerces an `mpfr` vector to double, so even a
  correctly written MPFR pipeline could not deliver precision to the user.
  **The real information loss is elsewhere**: `report_results()` returns
  `10^(-log10pvalue)`, which underflows to 0 above `log10pvalue = 308`, so every
  strongly-DE gene reports `p = 0` and cannot be ranked — while
  `pvalue_list$log10pvalue` holds the distinction. One-line fix: expose it.
- **`.compute_matrix_sd`'s sparse branch is dead code.** Its only in-package
  caller, `.initialize_residuals()`, passes `sd_vec = NULL` and a dense matrix,
  so the `sparseMatrixStats` removal cannot change any eSVD2 result.
- **Removing all-zero genes will change results for the genes that remain**, and
  that is correct rather than a bug. `.initialize_residuals()` takes the SVD of
  `log1p(A) − ZᵀC`; an all-zero column of `A` is not a zero column of that matrix
  but `−ZᵀC`, so empty genes have been pulling on the factorization of every
  other gene. Needs a release note. `Log_UMI` and `alpha_max` are provably
  unaffected.
- **The call graph is settled (`Q-COH-7`) and it is the plan's most consequential
  decision.** `eSVD_helper()` owns *all* preprocessing and postprocessing —
  donor drop, cohort checks, `gene_status` labelling, temporary gene removal, and
  reinsertion after `eSVD()` returns. `eSVD()` and `initialize_esvd()` **error**
  on anything that reaches them. The house rule is *the wrapper filters, the
  pipeline refuses*; apply it to any further preprocessing rather than
  re-deciding. Note this **supersedes `Q-STATUS-1`'s placement** — `gene_status`
  is not on `eSVD()`'s output.
- **The order inside `eSVD_helper()` is fixed**: (1) drop low-cell donors with a
  warning, (2) check the four cohort minima on the *filtered* object, (3) label
  `gene_status` and remove all-zero genes, (4) call `eSVD()`, (5) reinsert. Both
  ordering constraints are now internal to one function rather than contracts
  between two exported ones, which is far cheaper to keep true. (i) The donor drop
  must precede `gene_status`, or a gene expressed only in a dropped donor is
  marked `analyzed` and sent through with an all-zero count vector (T-COH-07).
  (ii) The cohort checks must follow the drop, or they are evaluated against a
  cohort that no longer exists when the model is fit (T-COH-03).
- **Two things the call graph bought, worth remembering when it is tempting to
  refactor.** `Q-STATUS-2` ("BH and the empirical null on analyzed genes only")
  becomes *structural* rather than conventional, because `eSVD()` never receives
  an all-zero gene and so cannot include one. And `compute_test_per_gene` never
  needs to learn about `gene_status`, because reinsertion happens after `eSVD()`
  returns — T-PG-01..06 are untouched.
- **What it costs: "all zero" is now decided in three places** (the helper labels,
  `eSVD()` errors, `initialize_esvd()` errors). One shared internal predicate
  `.which_all_zero()`, or the helper starts passing genes `eSVD()` rejects and the
  user gets an internal error from a function they never called. T-GS-17.
- **Both filters are `Seurat` object subsets, not matrix subsets** — `obj[, cells]`
  and `obj[genes, ]` — because `eSVD()` takes `seurat_obj`. Feature subsetting
  touches variable features, `scale.data` and reductions; none matter to `eSVD()`,
  but the `counts` layer and `meta.data` must survive intact. T-COH-13.
- **With `bool_check_donors` removed (`Q-COH-4`), thresholds are the only escape
  hatch.** Setting one to `0` disables it — every comparison is `<=` or `<`.
  The hazard is the middle: `min_cells_per_id = 1` or `2` weakens the drop rather
  than disabling it, and a surviving 1-cell donor is back on the ±Inf →
  silent-`locfdr`-degradation path. Recommend a `stopifnot` for off-or-safe.
- **`R CMD check` is structurally blind** to every defect in §1 of the readiness
  doc. A green check will never mean the numbers are right; only the tests will.

## Open Questions / Next Steps

**NEW, from the first test run (2026-08-29). Ordered by value:**

0a. **`gamma_rate()` is capped by the library size** — `src/gamma_rate.cpp:181`
   sets the search bracket's upper bound to `Rcpp::max(s)` and only ever shrinks
   it, so the estimated Gamma rate can never exceed the largest library size.
   Verified: correct to β = 1, then pinned at ≈0.9995 for β = 1.2, 1.5, 3, 10,
   while `log_gamma_rate` is right throughout. `estimate_nuisance` defaults to
   `bool_use_log = FALSE`, so every gene with over-dispersion γ < 1 gets the
   wrong nuisance, silently. Almost certainly why T-CPP-GAM-01 was commented out.
   **Stopgap: default `bool_use_log = TRUE`. Fix: let the bracket grow.**
0b. **Does the test statistic want the Bessel correction? — needs Kevin.**
   `.compute_mixture_gaussian_variance()` returns the population variance (÷n);
   Welch's t uses the sample variance (÷(n−1)). Verified to differ by exactly
   √(n/(n−1)) per arm, so the statistic is inflated 9.5% at 6 donors/arm and the
   p-values are anti-conservative. The mixture genuinely is a population object,
   so computing its population variance is defensible; using it in a statistic
   compared to a t distribution is the step that does not follow. **This changes
   every reported p-value.**
0c. **`.multtest_simple` underestimates the null sd by 21%** on exactly N(0,1)
   data — anti-conservative — and it is the last fallback. `truncated_mle` is 12%
   high. Also: at 200 genes the fallback fires on *ordinary* data, because
   `.multtest_locfdr` catches warnings and `locfdr` warns at modest gene counts.
0d. **Q-SVD-3 answered by the test**: the naive `Matrix` rewrite returns NaN
   (catastrophic cancellation). Ship T-SVD-05c's stable sparse form instead —
   it matches `matrixStats` to 1e-10 and needs no vendoring, no `cph` entry.
0e. `format_covariates()`'s "rescale" divides by the root-mean-square, not the
   sd, so it is not standardization. Documentation fix; the unit-invariance it
   does deliver is worth stating.


**Blocking the package's shape:**

1. ~~CRAN or Bioconductor?~~ **Resolved: CRAN.**
2. **`Authors@R` lists only Kevin.** The paper credits Yixuan Qiu with co-writing
   the R and C++ code and Kathryn Roeder with the method. Decide roles before
   submission. (§2.6)
3. **Should `asd.Rmd` / `asd-preprocess.Rmd` stay vignettes?** They need multi-GB
   external downloads and Bioconductor packages. Recommendation: move to
   `pkgdown/articles/`, keep only `eSVD2.Rmd`. (§2.7)

**No open questions in `UNIT_TEST_PLAN.md`.** §9.5 records three defaults taken
in the text, to confirm while implementing rather than before:

4. `eSVD()` errors on the three **correctness** conditions (all-zero gene,
   ≤2-cell donor, one donor in an arm) but **not** on the four cohort *power*
   minima — a legitimately underpowered 15-cell run must stay possible through
   `eSVD()` directly. (T-COH-11)
5. `min_cells_per_id` must be `0` (off) or `>= 3` (safe), never 1–2, which
   weakens the drop rather than disabling it and puts a 1-cell donor back on the
   §1.1 path. (T-COH-06)
6. The function is named `eSVD_helper`, not the imported file's `esvd_helper`.

**Ready to start, blocked on nothing:**

9. Step 0 of `CRAN_READINESS.md` §6 — the mechanical blockers (LICENSE stub,
   `override` ×5, drop `MASS` from the tests) — is independent of everything
   above. Merge its fixture shrink with `UNIT_TEST_PLAN.md` §1.
10. The two one-commit landmines in finding 5 above: `@export` on `eSVD()`, and
    a `requireNamespace("SeuratObject")` guard.
11. `UNIT_TEST_PLAN.md` §10 gives a 9-step test-first ordering, with §2.17/§2.16
    built as one step and `eSVD()`'s three errors (T-COH-11) written first so a
    helper bug fails loudly. **Implementation is gated on Kevin's go-ahead**,
    which he has now four times said he will give separately.
12. `CRAN_READINESS.md` needs a pass to absorb the session-3/4 findings.
