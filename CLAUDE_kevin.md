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

Session 1 established the scaffolding and a full audit. Where things stand:

- **Done this session.** `.gitignore` and `.Rbuildignore` rewritten so the CRAN
  tarball contains only the package (verified: the tarball's top level is now
  just `DESCRIPTION LICENSE NAMESPACE R README.md data inst man src tests
  vignettes`). Large-file pre-commit hook installed and `core.hooksPath` set.
  `.RData` untracked from git. `additional_context/summary.md` written with the
  paper's equations mapped to the functions that implement them.
- **The audit lives in `additional_context/CRAN_READINESS.md`.** That is the
  working document for this effort — read it before touching `DESCRIPTION`,
  `NAMESPACE`, `src/`, or `tests/`. Nothing in it has been fixed yet.
- **Baseline check recorded**: `R CMD check --as-cran` gives **5 WARNINGs,
  3 NOTEs** on R 4.5.1 / macOS. The existing test suite **passes**.

**The two findings that matter most:**

1. **`qnorm(pt(t, df))` returns `+Inf` for `t = 40, df = 18`** — verified. One
   strongly up-regulated gene makes `locfdr::locfdr()` error, which the
   `tryCatch` in `multtest()` swallows, silently degrading the empirical null all
   the way down to `.multtest_simple()` — the estimator whose own comment says it
   "is not a good estimator". The user gets no signal, because `compute_pvalue()`
   discards `multtest()`'s `method` field. This affects real analyses, not edge
   cases. Two-line fix (`log.p = TRUE` on both calls).
2. **Tarball is 13.1 MB against CRAN's 5 MB limit** — verified. Driven by
   `tests/assets/synthetic_data.RData`, which stores a fully-fitted 8.7 MB
   `eSVD_obj` plus three ~2.3 MB generation by-products that no test uses.

**Decision needed before other work** (see `CRAN_READINESS.md` §2.1): `eSVD2`
currently `Imports` `sparseMatrixStats`, which is Bioconductor-only and therefore
a hard CRAN stop. Either drop it (it is reachable only through a dead branch —
~4 lines of `Matrix` replaces it) or submit to Bioconductor instead. Everything
else in the plan assumes CRAN.

## Key Methodological Details

- **`nuisance_vec` is the Gamma *rate* `β = 1/γ`**, the reciprocal of the paper's
  overdispersion `γ_j`. Larger `nuisance_vec` = *less* overdispersion. The code is
  internally consistent; no roxygen block says this, so comparing `?compute_posterior`
  to Eq. 14 makes the code look wrong when it isn't. Documenting this is on the list.
- **Two independent implementations of the same computation exist**:
  `compute_posterior` + `compute_test_statistic` + `compute_pvalue` (matrix path)
  versus `compute_test_per_gene` (fused per-gene path). They are each other's best
  oracle. Same for `gamma_rate` vs. `log_gamma_rate` in `src/`. Both pairings
  should become strict equivalence tests.
- Their `library_min` defaults currently **disagree** (0.1 vs. 1e-2), so the two
  paths return different numbers when called directly with defaults. `eSVD()`
  passes the value explicitly, which is why this has gone unnoticed.
- `R CMD check` is **structurally blind** to every defect in §1 of the readiness
  doc — `checking R code for possible problems` reports OK. A green check will
  never mean the numbers are right; only §3's tests will.

## Open Questions / Next Steps

1. **CRAN or Bioconductor?** Blocks the `sparseMatrixStats` decision and hence
   most of the packaging work. (`CRAN_READINESS.md` §2.1)
2. **`Authors@R` currently lists only Kevin.** The paper credits Yixuan Qiu with
   co-writing the R and C++ code and Kathryn Roeder with the method. Decide roles
   and add them before submission. (§2.6)
3. **Should `asd.Rmd` / `asd-preprocess.Rmd` stay vignettes?** They need
   multi-GB external downloads and Bioconductor packages. Recommendation: move to
   `pkgdown/articles/`, keep only `eSVD2.Rmd` as a real vignette. (§2.7)
4. **Which `bool_stabilize_underdispersion` direction is intended?** The roxygen
   and the code disagree about whether the guard fires on high or low mean
   overdispersion. Only Kevin can say which was meant. (§1.4)
5. Work order for the fixes is in `CRAN_READINESS.md` §6. Step 0 (the mechanical
   blockers: fixture size, LICENSE stub, `override`, `MASS`) is independent of
   every open question above and can start immediately.
