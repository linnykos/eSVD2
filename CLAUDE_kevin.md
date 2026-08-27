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

Name the machine specifically enough that another collaborator can tell whether it is reachable to them. **A missing row means unknown; only an explicit *(not present)* row means known-absent.**

## Project Status (as of 2026-08-27)

**Goal: get `eSVD2` onto CRAN.** Correctness first; efficiency is explicitly out
of scope for now.

Two sessions in. Nothing in the package's `R/` or `src/` has been changed yet —
by design. The two working documents are:

- **`additional_context/CRAN_READINESS.md`** — the audit. *What is wrong*, with
  every claim tagged `[verified]` / `[inspection]` / `[policy]`. Read before
  touching `DESCRIPTION`, `NAMESPACE`, `src/`, or `tests/`.
- **`additional_context/UNIT_TEST_PLAN.md`** — the proposed test suite. *What we
  would have to assert to know it is right.* ~200 `test_that` blocks against the
  current ~40, each with an ID (`T-<AREA>-nn`) and an explicit **oracle** field.
  §3 of the readiness doc is now a pointer to this file. **Awaiting Kevin's
  review** — 16 open questions in its §9 block individual tests.

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
   verified in session 2, and *not* in the readiness doc's §1.
   `opt_esvd.default` defaults `nuisance_vec = rep(NA, p)`, but `gaussian`,
   `curved_gaussian`, `neg_binom` and `neg_binom2` all consume `gamma`, so the
   objective is `NA` and the call dies with `missing value where TRUE/FALSE
   needed`. Test `T-OPT-05`.

**Decision needed before other work** (`CRAN_READINESS.md` §2.1): `eSVD2`
`Imports` `sparseMatrixStats`, which is Bioconductor-only and a hard CRAN stop.
Either drop it (reachable only through a dead branch — ~4 lines of `Matrix`
replaces it) or submit to Bioconductor. Everything else assumes CRAN.

## Key Methodological Details

- **`nuisance_vec` is the Gamma *rate* `β = 1/γ`**, the reciprocal of the paper's
  overdispersion `γ_j`. Larger `nuisance_vec` = *less* overdispersion. The code is
  internally consistent; no roxygen block says this, so comparing
  `?compute_posterior` to Eq. 14 makes the code look wrong when it isn't.
  Several proposed tests assert a *direction* and their expected sign depends on
  this.
- **Two independent implementations of the same computation exist**, and each is
  the other's free oracle — no new reference implementation needed:
  `compute_posterior` + `compute_test_statistic` + `compute_pvalue` (matrix path)
  versus `compute_test_per_gene` (fused per-gene path); and `gamma_rate` versus
  `log_gamma_rate` in `src/`. `UNIT_TEST_PLAN.md` §10 puts both at step 2,
  before anything else is touched, because they are the cheapest possible safety
  net for every later change.
- Their `library_min` defaults **disagree** (0.1 vs. 1e-2), so the two paths
  return different numbers when called directly with defaults. `eSVD()` passes
  the value explicitly, which is why this has gone unnoticed.
- **A third free oracle is sitting unused in a comment**: `src/gamma_rate.cpp`
  lines 58–75 carry a complete R implementation of the objective, gradient and
  Hessian. Test `T-CPP-GAM-02` turns it into an independent check on both C++
  routines.
- **`numDeriv` has been in `Suggests` since the beginning and has never been
  used.** Gradient/Hessian checks for the 7 families are the highest-value C++
  tests available: an analytic-derivative error currently produces a silently
  mis-converged fit with *no* symptom — the loss still decreases and the
  optimizer still terminates.
- **A failed optimization is currently indistinguishable from a converged one.**
  `line_search` warns and returns `step = 0.0` with an unchanged iterate;
  `opt_esvd.default`'s convergence test then sees no change and `break`s,
  reporting convergence. Test `T-CPP-NEWT-02`.
- **`R CMD check` is structurally blind** to every defect in §1 of the readiness
  doc. A green check will never mean the numbers are right; only the tests will.

## Open Questions / Next Steps

**Blocking the package's shape:**

1. **CRAN or Bioconductor?** Blocks the `sparseMatrixStats` decision and hence
   most of the packaging work. (`CRAN_READINESS.md` §2.1; test `T-SVD-05`,
   question `Q-SVD-1`)
2. **`Authors@R` lists only Kevin.** The paper credits Yixuan Qiu with co-writing
   the R and C++ code and Kathryn Roeder with the method. Decide roles before
   submission. (§2.6)
3. **Should `asd.Rmd` / `asd-preprocess.Rmd` stay vignettes?** They need multi-GB
   external downloads and Bioconductor packages. Recommendation: move to
   `pkgdown/articles/`, keep only `eSVD2.Rmd`. (§2.7)

**Blocking specific tests — these are about *intended behaviour*, so only Kevin
can answer them** (full list of 16 in `UNIT_TEST_PLAN.md` §9):

4. **`Q-POST-1`: which direction should `bool_stabilize_underdispersion` fire?**
   The roxygen and the code disagree; because `nuisance_vec` is the rate, the
   code may be right and the prose wrong. The test's expected value depends
   entirely on the answer. (§1.4)
5. **`Q-TSTAT-1`:** an individual with exactly one cell — finite-but-huge
   statistic, or an error? This is precisely the input that produces the `±Inf`
   that breaks `locfdr`.
6. **`Q-REP-1`:** rank-deficient `.reparameterize` — warn, error, or proceed?
   `opt_esvd` currently swallows the failure in a `tryCatch` and silently keeps
   the unreparameterized matrices, so a rank-deficient fit looks like a good one.
7. **`Q-GAM-1`:** what should `gamma_rate` return for an all-zero gene?
8. **`Q-PROP-2`:** where does the null-calibration test live? It is the paper's
   Type-1-error claim made executable and would have caught finding 1 above — and
   it is also the most likely to make a CRAN check flaky on a machine we don't
   control. Recommendation: a slow/CI-only tier.

**Ready to start, blocked on nothing:**

9. Step 0 of `CRAN_READINESS.md` §6 — the mechanical blockers (LICENSE stub,
   `override` ×5, drop `MASS` from the tests) — is independent of every question
   above. The fixture shrink in that step should be merged with
   `UNIT_TEST_PLAN.md` §1, since regenerating the fixtures is the same work.
10. Once the test list is agreed: `UNIT_TEST_PLAN.md` §10 gives a test-first
    ordering (harness → fixtures → the two free equivalence tests → gradient
    tests → fixes) that deliberately differs from `CRAN_READINESS.md` §6.
