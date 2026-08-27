# HISTORY_kevin.md — Kevin's Session Log

> **Append-only, ascending chronological order** (oldest at top, newest at the bottom). Add each session's dated entry to the END of this file. Never read at session startup — consulted only on demand for deep history. Current project state lives in `CLAUDE_kevin.md`.
>
> **Owned by Kevin.** Only Kevin's session appends to this file. Other collaborators may read it but must not edit or rewrite entries.

---

### 2026-08-27 (Session 1 — CRAN-readiness scaffolding and audit)

**Scaffolding**
- Initialized the project-setup three-file pattern: master `CLAUDE.md`,
  `CLAUDE_kevin.md`, `HISTORY_kevin.md`. Installed `.githooks/pre-commit`
  (blocks staged files ≥ 50 MB) and set `git config core.hooksPath .githooks`.
  Each collaborator must repeat that `git config` line after cloning.
- Read both PDFs in `additional_context/` and wrote `summary.md`. Wrote it to
  double as the method spec — every equation is mapped to the function that
  implements it — so future sessions need not reopen the PDFs to reason about
  the code.

**Ignore files (the requested change)**
- Rewrote `.Rbuildignore` and `.gitignore` with deliberately different scopes:
  git keeps `additional_context/` and the `CLAUDE*.md` / `HISTORY_*.md` files;
  `.Rbuildignore` excludes them from the tarball.
- **`.Rbuildignore` does not support `#` comments** — every line is parsed as a
  regex. My first version had explanatory comments and `R CMD build` died with
  `invalid regular expression '# Entries are Perl-compatible regular expressions
  matched (case-insensitively,'`. Rewrote it comment-free; the explanation now
  lives in `CLAUDE.md` instead. Worth remembering — the failure mode is a
  confusing PCRE error, not a warning about comments.
- Untracked `.RData` from git (`git rm --cached`; the file remains on disk).
- Verified the result: the tarball's top level is now exactly
  `DESCRIPTION LICENSE NAMESPACE R README.md data inst man src tests vignettes`,
  and `src/*.o` / `*.so` / `vignettes/.DS_Store` are all excluded.

**Build environment note**
- `R CMD build` run directly against the Dropbox CloudStorage working directory
  hangs indefinitely at "checking for file DESCRIPTION" — it stalls copying, and
  the 278 MB `.git` plus on-demand file hydration appear to be the cause.
  Workaround that works: `rsync -a --exclude .git --exclude docs --exclude
  .Rproj.user --exclude additional_context` into a local temp dir and build there.
  Do this for every future check run; it is not a package problem.

**Audit → `additional_context/CRAN_READINESS.md`**
- Recorded a real baseline: `R CMD check --as-cran` = **5 WARNINGs, 3 NOTEs**
  (R 4.5.1, aarch64-darwin). The existing testthat suite **passes**.
- Every claim in the document is tagged `[verified]` / `[inspection]` /
  `[policy]` so a later reader knows what was actually executed.

**Findings worth carrying forward (rationale not visible in the code):**

- **The `+Inf` bug is a chain, not a line.** `qnorm(pt(t, df))` saturating at
  `t = 40, df = 18` is only the trigger. What makes it serious is that
  `locfdr::locfdr()` *errors* on `Inf` (`'to' must be a finite number`), that
  error is swallowed by the `tryCatch` in `.multtest_locfdr()`, and `multtest()`
  then silently walks down its fallback chain to `.multtest_simple()`. Because
  `compute_pvalue()` keeps `null_mean`/`null_sd` but discards `method`, there is
  no trace in the output that the empirical null was ever abandoned. I traced and
  verified each link separately. The asymmetry matters too: the *lower* tail
  stays finite to `t = −40`, so the degradation is triggered only by
  up-regulated genes.
- **Resolved: the paper and the code disagree about `γ` only in notation.** Eq.
  14 reads `(μ/γ + A)/(1/γ + ℓ)` while the code computes
  `(A + μ·nuisance)/(ℓ + nuisance)`. These match iff `nuisance == 1/γ`, i.e. the
  code's nuisance is the Gamma *rate* — consistent with the C++ being named
  `gamma_rate`. The code is right; nothing documents the inversion.
- **Open: which way `bool_stabilize_underdispersion` should fire.** Roxygen says
  "when the global mean over-dispersion is less than 1"; the guard is
  `mean(log10(nuisance_vec)) > 0`, i.e. geometric mean > 1. Under the rate
  parameterization above the code may well be correct and the prose stale, but I
  can't tell which was intended. Needs Kevin.
- **`exportPattern("^[[:alpha:]]+")` in `R/zzz.R` is doing more damage than it
  looks.** It exports the raw Rcpp bindings (`objfn_Xi_r`, `hessian_YZj_r`,
  `test_data_loader`, …), which take external pointers — a user passing a
  deserialized pointer segfaults R. It also means `eSVD()`, the top-level
  wrapper, is public with no docs and no tests purely by accident.
- **Two natural test oracles already exist in the codebase** and neither is being
  used as one: `compute_test_per_gene` vs. the three-function matrix path, and
  `gamma_rate` vs. `log_gamma_rate`. The existing equivalence assertion is
  `abs(sum(res1 - res2)) <= 1e-3`, which sums *signed* differences and so cannot
  detect per-gene errors that cancel. Element-wise `expect_equal` instead.
- **`opt_esvd.default(verbose = 2)` has never been run.** `print(x, resid)` binds
  `resid` to `digits` and errors with `invalid printing digits 0`. A dead
  giveaway that verbose branches are wholly untested.
- **Unpredicted findings that only came from actually running the check**: the
  13.1 MB tarball vs. CRAN's 5 MB cap (`synthetic_data.RData` carries a fitted
  8.7 MB `eSVD_obj` and three unused ~2.3 MB matrices), the `LICENSE` file being
  full MIT text where `License: MIT + file LICENSE` requires a 2-line DCF stub,
  undeclared `MASS` in the tests, and five `-Winconsistent-missing-override`
  warnings. None of these were visible by reading the source — worth running the
  check early next time rather than auditing first.
- **Verified that `Rmpfr` buys nothing here.** `Rmpfr::pnorm` on plain doubles
  returns a base `numeric` identical to `stats::pnorm`; no argument in the
  package is ever an `mpfr` object. So the GMP/MPFR system dependency the README
  apologizes for is pure cost.
- Decided **not** to change any package code this session — the request was to
  produce the plan, and the ordering in `CRAN_READINESS.md` §6 depends on the
  unresolved CRAN-vs-Bioconductor question.
