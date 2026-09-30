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

### [2026-08-27] (Session 2 — writing the unit-test plan)

- Wrote `additional_context/UNIT_TEST_PLAN.md` (720 lines): the full proposed
  suite, ~200 `test_that` blocks against the current ~40. Replaced §3 of
  `CRAN_READINESS.md` with a pointer to it and added `T-*` test citations at the
  §1.1 / §1.3 / §1.4 / §4.1 / §4.3 "Test" paragraphs and in §6 step 9.
- **Every test carries an explicit `oracle` field**, tagged `[oracle]` /
  `[invariant]` / `[snapshot]`. This is the field to argue about in review: a
  test whose oracle is "whatever the code returns today" locks in the bug. The
  ID scheme (`T-<AREA>-nn`) exists so Kevin can strike or dispute tests by
  reference without quoting them.
- **Decided the plan is fixture-first, not assertion-first.** The bulk of the
  work is regenerating `tests/assets/`, not writing `expect_*` calls. Proposed
  `F-TINY` (120×20, 6 individuals), `F-SMALL` (400×40, 8), `F-DEGEN` (hand-built
  degenerate inputs, never stored), `F-NULL`, `F-DERIV` (per-family feasible
  points). Fitted objects are built at test time from raw inputs, not stored —
  smaller *and* a stronger test than reading back a stored fit.
- **Decided the work order in the plan deliberately contradicts
  `CRAN_READINESS.md` §6**: harness → fixtures → the two *free* equivalence tests
  → gradient tests → fixes. Rationale: the matrix-vs-per-gene and
  `gamma_rate`-vs-`log_gamma_rate` pairs need no new oracle (each is the other's)
  and together cover the posterior/test/p-value stack plus the nuisance
  estimator, so they are the cheapest possible safety net and must exist before
  anything is touched.

**Four defects found while writing the plan that §1 of `CRAN_READINESS.md` does
not list** (all recorded in the new §3 of that file):

- **Four of seven families are unusable at their documented defaults** —
  verified. `opt_esvd.default` defaults `nuisance_vec = rep(NA, ncol(dat))`, but
  `gaussian`, `curved_gaussian`, `neg_binom`, `neg_binom2` all consume `gamma`,
  so `objfn_all_r` returns `NA` and `opt_esvd.default(family = "gaussian")` dies
  with `missing value where TRUE/FALSE needed` — naming neither the family nor
  the parameter. Only `poisson` and `bernoulli` ignore `gamma` entirely. This is
  why six families have no tests: a naive test would fail immediately and it
  would look like the test's fault.
- **`format_covariates()` drops the *first* factor level; the roxygen says the
  last** — verified with `factor(c("a","a","b","b","c","c"))` → columns `g_b`,
  `g_c`. Second doc/code mismatch in the same function: it rescales only the
  variables named in `rescale_numeric_variables`, while the roxygen says "all
  the numerical variables".
- **`data_loader()` returns a null external pointer without erroring** for an S4
  that is not `dgCMatrix` (a `dgeMatrix` is plausible user input) or a dense
  matrix that is neither integer nor numeric — verified. The
  `Rcpp::stop("unsupported matrix type")` is only reachable for a non-S4
  non-matrix, because the two type tests are `if`/`else if` and the `else`
  catches neither fall-through.

**Correction to my own session-1 claim (`CRAN_READINESS.md` §4.1):** a
serialized-and-restored `XPtr` does **not** segfault under the current Rcpp.
`saveRDS`/`readRDS` round-tripping `esvd_family()` or `data_loader()` output and
then calling `feas_Xi_r()` / `objfn_all_r()` gives a clean
`Error: external pointer is not valid` from Rcpp's checked `XPtr(SEXP)`
constructor — verified in a subprocess. The guard is still worth adding (we
should own that guarantee rather than inherit it from a dependency, and the
message is unhelpful), but it is not a crash risk and should not be prioritized
as one. **The real pointer hazard is different and previously unnoticed:**
`DenseDataLoader` holds an `Eigen::Ref` to the *R matrix's own memory*, so a
loader outliving its R matrix dangles. That is `T-CPP-PTR-04` and it is the one I
would write first.

- Open: 16 questions collected in `UNIT_TEST_PLAN.md` §9, each blocking at least
  one test. Four need Kevin specifically because they are about *intended*
  behaviour, not code: `Q-POST-1` (which direction should
  `bool_stabilize_underdispersion` fire), `Q-TSTAT-1` (one-cell individual:
  finite-but-huge or error), `Q-REP-1` (rank-deficient `.reparameterize`: warn,
  error, or proceed — `opt_esvd` currently swallows it in a `tryCatch` so a
  rank-deficient fit is indistinguishable from a good one), `Q-GAM-1` (what
  should `gamma_rate` return for an all-zero gene).
- Open: `Q-PROP-2` — the null-calibration test (`generate_null()` → p-values
  approximately uniform, KS test) is simultaneously the most valuable test in the
  plan (it is the paper's Type-1-error claim, executable, and would have caught
  §1.1) and the most likely to make a CRAN check flaky on a machine we don't
  control. My recommendation is a slow/CI-only tier rather than the CRAN suite.
- Noted `verbose` probes were run against the **installed** `eSVD2 1.0.1.2`, not
  the working tree's `1.0.1.07`. Every `[verified]` claim in the new file carries
  that caveat; re-confirm against a fresh install before acting on any of them.
- Decided again **not** to change package code — the request was explicitly to
  produce the test list for review before any implementation or fixture work.

### [2026-08-29] (Session 3 — Kevin's review pass folded into the test plan)

- Kevin answered 15 of the 16 §9 questions in-line as `[[KZL: ...]]`. Each answer
  is now folded into the affected test as **✅ RESOLVED**, with the verbatim
  comment kept alongside, and §9 is split into 9.1 (resolved, with consequences)
  and 9.2 (still open).
- Open: `Q-FIX-1` is the one question the pass did not answer. `[[KZL: Yes, this
  is a good test]]` agrees with the *approach* (build the fitted object at test
  time) but not with the *shape* question (is 120 × 20 at `k = 2` degenerate),
  which every §2 tolerance depends on. Recorded my own reading — fine for
  structural assertions, use `F-SMALL` for recovery assertions — so it is no
  longer blocking, but it still needs Kevin.
- Several answers cost more than the test they unblocked, and that is now
  recorded at each test rather than left to be rediscovered:
  - `Q-REP-1` (warn and proceed) is inert unless `opt_esvd.default` stops
    swallowing the warning in its `tryCatch` → new T-REP-09/10.
  - `Q-TSTAT-1`'s ≤2-cell check must live in `eSVD()` only. Putting it in
    `compute_test_statistic()` would delete T-TSTAT-01 and T-DF-01, which are
    built on a one-cell-per-individual fixture *because* that degenerate case is
    where the mixture machinery reduces exactly to Welch's t and `stats::t.test`
    becomes an external oracle — the only independent check on the headline
    statistic.
  - `Q-GAM-2` (no convergence status) leaves `.nuisance_in_sequence()`'s
    `log_gamma_rate` fallback permanently unreachable, since it triggers on a
    non-finite result and a saturated result is finite. Now a deliberate choice
    to be commented in the code, not an oversight.
  - `Q-POST-1` (code right, prose wrong) makes T-POST-08 a *documentation* fix.
- Resolved: `Q-SVD-1` → CRAN, absorb the `sparseMatrixStats` function. Opened
  `Q-SVD-3` because "vendor it" is the more expensive of the two readings:
  `colSds` there is C++, so vendoring means a new `src/` file plus a `cph` entry
  in `Authors@R` and `inst/COPYRIGHTS` (MIT permits the copy but CRAN requires
  the attribution). Found the deciding fact — `.compute_matrix_sd`'s sparse
  branch is **dead code**: its only in-package caller `.initialize_residuals()`
  passes `sd_vec = NULL` and a dense matrix. So no eSVD2 result changes either
  way, and the four-line `Matrix` rewrite wins. Third option noted: delete
  `sd_vec`/`rescale` from `.svd_safe()` entirely.
- **Rmpfr — settled empirically, not by inspection** (new §2.15, T-MPFR-01..08).
  The old T-PVAL-07 only showed the *current calls* are no-ops; Kevin asked for
  the stronger question. Verified on R 4.5.1 / macOS:
  - Double `pnorm(z, log.p=TRUE)` matches a genuine 200-bit MPFR computation to
    **relative ≤7.3e-17** across `z ∈ [-1e4, -5]`. Doubles are bit-exact in log
    space over every regime the package can reach.
  - `2*Rmpfr::pnorm(-40)` is **exactly 0**, same as `stats::pnorm` — so `Rmpfr`
    does *not* prevent the `pvalue_vec` underflow it appears to guard. The zeros
    come from the `log.p = FALSE` parameterization, not from precision.
  - `stats::p.adjust` on an `mpfr` vector silently coerces to double: a 1e-350
    p-value comes back as exactly 0. **Even a correct MPFR pipeline could not
    deliver its precision to the user**, because BH is double-only. This is the
    fact that closes the question.
  - Found a dependency-free permanent oracle: the Mills-ratio asymptotic series
    matches `stats::pnorm(log.p=TRUE)` to 4.4e-13 at `z=-20` and ≤1.2e-16 beyond,
    so T-MPFR-05 can pin double adequacy after `Rmpfr` is gone.
  - Found the one place information *is* lost, and it is not a precision problem:
    `report_results()` returns `10^(-log10pvalue)`, which underflows to 0 for
    `log10pvalue > 308`, so every strongly-DE gene reports `p = 0` and cannot be
    ranked — while `pvalue_list$log10pvalue` holds the distinction the whole
    time. One-line fix: expose `log10pvalue` as a column.
  - Open: `Q-MPFR-1` — ship the Rmpfr-dependent tests (needs `Rmpfr` in
    `Suggests`, the dependency they exist to remove) or run them once as a
    recorded experiment. Recommended the latter.
- **`gene_status` — Kevin's new output, specified before tested** (new §2.16,
  T-GS-01..16). Verified the three places an all-zero gene misbehaves today:
  `glmnet` warns `"Convergence for 1th lambda value not reached after
  maxit=100000"` and returns the previous lambda's fit; `gamma_rate` returns
  9.706e-4 while `log_gamma_rate` returns exactly `-10` (its clamp) — a 21×
  disagreement, and the most likely reason T-CPP-GAM-01 was commented out; the
  Welch statistic is 0/0.
- Two non-obvious points about the status feature that the spec now states out
  loud because both are silent when wrong:
  - **BH and the empirical null must be fit on analyzed genes only.** Including
    status-2 genes inflates `n` in BH and raises every real gene's adjusted
    p-value in proportion to how empty the cell type is. Kevin's "put back at the
    VERY end" implies it; T-GS-09/10 enforce it. → `Q-STATUS-2` asks him to
    confirm.
  - **Removing all-zero genes changes the fit for the retained genes**, and that
    is expected, not a bug. `.initialize_residuals()` takes the SVD of
    `log1p(A) − ZᵀC`; an all-zero column of `A` is not a zero column of that
    matrix but `−ZᵀC`, so empty genes have been actively pulling on the
    factorization of every other gene. This is a behaviour change for the release
    note. Conversely `Log_UMI` and `alpha_max` are provably unaffected (T-GS-03,
    T-GS-04).
- Recommended a `factor` with `levels = c("analyzed","all_zero")` rather than a
  bare integer: `as.integer()` still gives Kevin's 1/2, printed objects stay
  readable, and the enum is expected to grow (`Q-STATUS-4`).
- Open: `Q-STATUS-1` (does `initialize_esvd()` filter or error — recommended
  error, so there is one definition of "all zero"), `Q-STATUS-3` (`fdr_vec` = 1 or
  NA at status-2 positions — recommended 1), `Q-TSTAT-2` (is the ≤2-cell
  threshold an argument, and error vs. drop the individual — recommended
  `min_cells_per_individual = 3` and error).
- Still no package code changed. The plan is now ~225 proposed tests against ~40
  today; implementation is explicitly gated on Kevin's go-ahead.

### [2026-08-29] (Session 4 — Kevin's review pass 2 on the test plan)

- Kevin answered the remaining 8 questions. Six matched my recommendations
  (Q-FIX-1, Q-SVD-3, Q-MPFR-1, Q-STATUS-1/2/3); Q-STATUS-4 froze the status enum
  at two levels; Q-TSTAT-2 was answered by pointing at existing code rather than
  by a decision.
- Resolved: `Q-MPFR-1` "outside of the package" → `Rmpfr` leaves `DESCRIPTION`
  entirely rather than moving to `Suggests`; T-MPFR-01..04 run once as a script
  under `additional_context/` and never enter `tests/`.
- Resolved: `Q-STATUS-4` keeps only status 1/2. Noted that this weakens my own
  argument for a `factor` over an integer — with the enum frozen it is now only a
  readability preference, and said so in the spec rather than leaving the
  original justification standing.
- **New external location recorded**: `WAS2CODE_REPO` (see `CLAUDE.md` and the
  per-machine table here). It is under `Dropbox/.../archive/`, so treat it as
  frozen — the plan says to *copy* `esvd_helper.R` in, not to depend on the path.
- Read `esvd_helper.R` (60 lines). It answers Q-TSTAT-2 and brings **three
  cohort-level filters that were not in the plan at all** — `min_cells = 20`,
  `min_ids = 4`, `min_cells_casecontrol = 20`, each `warning()` + `return(NA)`,
  plus `min_cells_per_id = 3`. Written up as new §2.17 with 10 tests.
- **Open: `Q-COH-3` — Kevin's two answers conflict.** Pass 1 said the ≤2-cell
  check should *stop* `eSVD()` early; `esvd_helper.R` *drops* those donors and
  proceeds. T-ESVD-09 and T-VAL-31 cannot be written until this is settled.
  Recommended drop (it is the behaviour actually used in the Was2CODE analyses),
  but with a `warning()` naming the donors — the file currently drops silently.
- **Open: `Q-COH-2` — `min_ids` does not prevent the failure it appears to.** It
  counts donors pooled across arms, so 4 case + 1 control passes `min_ids > 4`
  and still gives `n2 - 1 = 0` in `.compute_df()`, i.e. T-TSTAT-04's `NaN` Welch
  df reached through the front door. `min_cells_casecontrol` has the same gap: it
  counts *cells* per arm, and one donor can supply all 20. Recommended adding
  `min_ids_per_arm = 2`. T-COH-04 is deliberately written to fail against the
  file as it stands.
- **Open: `Q-COH-4`** — `bool_check_donors = FALSE` disables the low-cell-donor
  drop along with the three power filters. The three cell/donor minima are power
  judgements a user may override; a 1-cell donor is a *correctness* problem
  (zero within-individual variance → ±Inf → `locfdr` degrades silently).
  Recommended splitting the switch.
- **Open: `Q-COH-1`** (bare `NA` return is poor as a package API — recommended a
  classed `eSVD_skipped` sentinel with a `reason`) and **`Q-COH-5`** (shape:
  recommended an exported `filter_cohort()` that `eSVD()` calls internally).
- **Ordering constraint found, now T-COH-07**: the donor drop must run *before*
  `gene_status` is computed. A gene expressed only in a dropped donor becomes
  all-zero after the drop; computing status first marks it `analyzed` and sends
  an all-zero count vector through the pipeline — the exact failure §2.16 exists
  to prevent, reintroduced by ordering alone. This is why Q-COH-5's answer
  matters beyond style.
- Two mechanical findings, no decision needed, neither in either working
  document:
  - **`NAMESPACE` has no `export(eSVD)`.** The package's main user-facing
    function is exported only by the blanket `exportPattern("^[[:alpha:]]+")`,
    which `Q-CPP-1` removes — so that change would silently un-export `eSVD()`.
    One `@export` tag, but it must land in the same commit. Test T-ESVD-11.
  - **`SeuratObject` is used unguarded from `Suggests`.** `eSVD()` calls
    `SeuratObject::LayerData()` at `R/eSVD.R:38` and there is not one
    `requireNamespace()` call anywhere in `R/`. CRAN requires conditional use of
    a suggested package. Test T-VAL-34. Ironically `esvd_helper.R` already models
    the right pattern, for `eSVD2` itself.
- Plan now ~236 proposed tests. Still no package code changed — Kevin has said
  again to hold off implementing.

### [2026-08-29] (Session 5 — Kevin's review pass 3 on the test plan)

- All five `Q-COH` questions answered, including a five-step ordered spec for
  `eSVD_helper()`. Resolved: `Q-COH-2` (add `min_ids_per_arm = 2`), `Q-COH-3`
  (**drop** low-cell donors with a warning — settles the pass-1/pass-2 conflict),
  `Q-COH-4` (remove `bool_check_donors` entirely; thresholds are the only knob),
  `Q-COH-5` (keep `eSVD_helper()` *and* export `filter_cohort()`).
- `Q-COH-1`: Kevin declined the classed `eSVD_skipped` sentinel and kept
  `warning()` + `return(NA)`. Consequence recorded rather than re-argued: the
  warning *text* is now the only channel identifying which of the three
  rejections fired, so T-COH-05 asserts three distinct `regexp` fragments.
- `Q-COH-4`'s answer is better than the split switch I proposed, but makes
  "set the threshold to 0" the only escape hatch. Verified that works — every
  comparison is `<=` or `<`, so `0` disables. Added T-COH-12 to make it a tested
  contract, and T-COH-06 for the hazard it opens: `min_cells_per_id = 1` or `2`
  does not disable the drop, it *weakens* it, and a surviving 1-cell donor is
  back on the ±Inf → silent-`locfdr`-degradation path. Recommended
  `stopifnot(min_cells_per_id == 0 || min_cells_per_id >= 3)`.
- **Open: `Q-COH-6` — the stated step order is the reverse of the imported
  file's.** Kevin's steps run the three count checks (1–3) then drop low-cell
  donors (4); `esvd_helper.R` drops first and checks after. Checking first admits
  exactly the cohorts the checks exist to reject: 6 donors of which 2 have 1 cell
  passes a 6-donor check, then becomes 4 — and if both were controls, the
  `min_ids_per_arm = 2` just added is being checked against a donor count that no
  longer exists by the time the model is fit. Recommended drop-first; T-COH-03 is
  written to that order and will fail if the literal step order is implemented.
- **Open: `Q-COH-7` — the layering.** `Q-COH-5` keeps `eSVD_helper()` a *sibling*
  of `eSVD()` rather than a layer inside it, so `eSVD()` called directly (as the
  vignette does) gets no donor filtering and stays on the §1.1 path. Proposed
  resolving it with a rule that is **already decided for exactly this shape** by
  `Q-STATUS-1`: the wrapper filters, the pipeline refuses. Wrote it as a
  three-level table — `eSVD_helper()` filters donors and inherits gene filtering;
  `eSVD()` errors on low-cell donors and filters genes; `initialize_esvd()` errors
  on all-zero genes.
- Sub-conflict inside `Q-COH-7`: Kevin's step 5 puts the all-zero gene handling in
  `eSVD_helper()`, while `Q-STATUS-1` put it in `eSVD()`. §2.16 is written to
  `Q-STATUS-1` so there is one implementation of `gene_status`; if it flips, every
  §2.16 test re-points at the helper and `eSVD()` called directly gains no
  `gene_status` output.
- `Q-COH-5`'s answer incidentally defused the `F-TINY` sizing worry from session
  4: because the helper is a sibling, every §2 test that calls `eSVD()` is
  unaffected by the five thresholds. Only T-COH-08 (the no-op check) needs a
  cohort sized to clear all five.
- Naming noted for the import: the file defines `esvd_helper`, Kevin writes
  `eSVD_helper`. Recommended `eSVD_helper` to sit beside `eSVD` / `eSVD2`.
- Plan now ~239 proposed tests. Still no package code changed — Kevin has now
  three times said to hold off implementing.

### [2026-08-29] (Session 6 — Kevin's review pass 4; every question now closed)

- `Q-COH-6` resolved: **drop donors first, then check** — matches the imported
  file and my recommendation.
- `Q-COH-7` resolved, and **it supersedes `Q-STATUS-1`'s placement**: all
  filtering, `gene_status` labelling, temporary removal *and* reinsertion live in
  `eSVD_helper()`; `eSVD()` errors on anything that gets through. `gene_status`
  therefore moves out of `eSVD()` entirely. This is a different structure from
  the one I proposed (which had `eSVD()` filtering genes and the helper only
  filtering donors), and it is better — recorded why below.
- **The architecture makes `Q-STATUS-2` structural instead of conventional.**
  Under my proposal, "BH and the empirical null on analyzed genes only" was a
  discipline `compute_pvalue()` had to honour, with T-GS-09/T-GS-10 there to catch
  it being broken. Under Kevin's, `eSVD()` never receives an all-zero gene, so it
  cannot include one. T-GS-09/10 demoted from load-bearing to architecture checks.
- Second thing the new call graph paid for: **`compute_test_per_gene` no longer
  needs to know about `gene_status` at all**, because reinsertion happens after
  `eSVD()` returns. T-PG-01..06 are untouched, which they would not have been.
- The definitive order inside `eSVD_helper()` is now fixed and every §2.16/§2.17
  test is written against it: (1) drop low-cell donors with a warning, (2) check
  the four cohort minima on the *filtered* object, (3) label `gene_status` from
  the donor-filtered counts and remove the all-zero genes, (4) call `eSVD()`,
  (5) reinsert at original positions.
- **New cost of the design, now T-GS-17**: "all zero" is decided in three places
  (helper labels, `eSVD()` errors, `initialize_esvd()` errors). They must share
  one internal predicate `.which_all_zero()` or the helper will start passing
  genes `eSVD()` rejects, and the user gets an internal error from a function they
  never called.
- **New implementation concern, now T-COH-13**: both filters are `Seurat` *object*
  subsets, not matrix subsets — `obj[, keep_cells]` and `obj[keep_genes, ]` —
  because `eSVD()` takes `seurat_obj`. Feature subsetting a Seurat object touches
  variable features, `scale.data` and reductions; none matter to `eSVD()`, but the
  tests must assert the `counts` layer and `meta.data` survive intact.
- Recorded three implement-time defaults rather than opening new questions
  (§9.5): `eSVD()` errors on the three *correctness* conditions only, not the
  four cohort power minima (an underpowered 15-cell run must stay possible
  directly); `min_cells_per_id` must be 0 or ≥3, never 1–2, which weakens the
  drop rather than disabling it; and the function is named `eSVD_helper`.
- Plan now ~242 proposed tests, **0 open questions**. §9.1–9.4 record 30 resolved
  decisions with their consequences. Still no package code changed — Kevin has now
  four times said to hold off implementing.

### [2026-08-29] (Session 7 — wrote and ran the test suite; five new findings)

- Wrote 15 new test files (`test_*_claude.R`) plus `helper-fixtures.R`, the
  `Rmpfr` decision script and two reports. **No `R/`, `src/`, `DESCRIPTION` or
  `NAMESPACE` change**, as instructed. New suite: 388 pass / 87 fail / 29 skip.
  Existing suite unchanged at 79 pass / 0 fail — writing alongside rather than
  editing kept the baseline clean, which is what makes the failures readable.
- **Fixtures are built in code with a fixed seed, not shipped as `.rda`.**
  Deviation from the §1 draft; it removes the fixture bytes from the tarball
  entirely, makes provenance readable, and forces the "fitted object" tests to
  exercise the pipeline. F-TINY deliberately has no all-zero genes — with
  `gene_status` unimplemented they would fail every unrelated test for the same
  reason.
- **Finding 1, serious: `gamma_rate()` is capped by the library size.**
  `src/gamma_rate.cpp:181` sets the bracket's upper bound to `Rcpp::max(s)` and
  the refinement loop only ever shrinks it, so the estimated Gamma rate can
  never exceed the largest library size. Verified: tracks the truth to β = 1,
  then pins at ≈0.9995 for β = 1.2, 1.5, 3, 10, while `log_gamma_rate` is
  correct throughout. Rescaling `s` moves the answer, which proves the cause;
  the R objective in the `gamma_rate.cpp` comment block confirms the likelihood
  at the true β is strictly higher than at the returned value. Since
  `estimate_nuisance` defaults to `bool_use_log = FALSE`, every gene with
  over-dispersion γ < 1 gets the wrong nuisance. **Almost certainly why
  T-CPP-GAM-01 was commented out** — above the cap the two routines cannot agree.
- **Finding 2: the test statistic omits the Bessel correction.**
  `.compute_mixture_gaussian_variance()` returns the *population* variance
  (÷n); Welch's t uses the sample variance (÷(n−1)). On the one-cell fixture the
  two differ by exactly √(n/(n−1)) per arm — verified to 1e-8. The statistic is
  inflated 9.5% at 6 donors/arm, 5.4% at 10, 2.6% at 20, and the p-values are
  anti-conservative while being referred to a t distribution. Open: is this
  intended? The mixture is a population object, so its population variance is
  defensible; using it in a t statistic is the step that does not follow.
- **Finding 3: `.multtest_simple` underestimates the null sd by 21%** on exactly
  N(0,1) data (0.79 vs locfdr's 0.99), i.e. anti-conservative; `truncated_mle`
  is 12% high, i.e. conservative. This quantifies §1.1's silent degradation for
  the first time. Also found that at 200 genes the fallback fires on *ordinary*
  data, because `.multtest_locfdr` catches warnings and `locfdr` warns at modest
  gene counts — returning null mean 0.70, sd 2.57.
- **Finding 4: `format_covariates()`'s "rescale" is not standardization.**
  `scale(x, center = FALSE)` divides by the root-mean-square, not the sd, so a
  covariate with mean 40 / sd 8 comes out with sd 0.2. The property it does
  deliver — invariance to the recorded unit — is real and worth documenting; the
  roxygen word "rescales" is what needs to change.
- **Finding 5: Q-SVD-3 answered by the test, and neither prior option wins.**
  The naive `Matrix` rewrite returns **NaN** on a large-mean/small-variance
  column (catastrophic cancellation). But vendoring C++ is still unnecessary:
  T-SVD-05c ships a stable sparse form — squared deviations over stored
  non-zeros plus `(n − nnz)·m²` — that matches `matrixStats` to 1e-10.
  Recommend shipping that.
- **All 28 analytic gradient/Hessian checks pass across all seven families.**
  The plan's largest block and its highest-value C++ test: the derivatives are
  correct. `numDeriv` has finally earned its `Suggests` entry.
- Q-POST-1 confirmed empirically: the stabilization fires on
  `mean(log10(nuisance_vec)) > 0`, so the code is right and the prose is wrong.
- Not root-caused, and flagged as such rather than claimed: T-PROP-08a (scale
  equivariance) fails by up to 38% even with the clamps loosened. Most likely my
  test — it reuses `nuisance_vec` from the unscaled data instead of re-fitting.
  A conclusive version needs a full re-fit.
- `test_gene_status_claude.R` (15) and `test_cohort_filter_claude.R` (12) skip
  cleanly on `exists("eSVD_helper")`, so the features can be built against them.
  Four of their tests already run; T-COH-11 confirms `eSVD()` has no Q-COH-7
  guard yet, and T-COH-13 caught a real Seurat detail — feature subsetting needs
  `subset(features =)`, not `obj[genes, ]`.
- `Rmpfr` settled and reported in `additional_context/RMPFR_REPORT.md`: doubles
  match a 200-bit MPFR oracle to ≤7.3e-17 in log space; `2*Rmpfr::pnorm(-40)`
  underflows to 0 exactly as `stats::pnorm` does; and `p.adjust` silently
  coerces an mpfr vector to double, so MPFR precision could not reach the user
  even if used correctly. Recommend dropping it from `DESCRIPTION` outright.

### [2026-08-29] (Session 8 — Kevin's answers applied; first code changes)

- **First changes to `R/` and `src/`.** Suite went 388 pass / 87 fail →
  **442 / 36**; existing 79-test suite still green.
- **`src/gamma_rate.cpp` bracket fixed** (Kevin: "let the bound grow"). Four
  changes: grow `ub` while `[l(ub)]' > 0`; shrink `lb` until `[l(lb)]' > 0` so
  the root is bracketed both sides; start Newton from the geometric mean rather
  than a fixed 1.0; raise the iteration cap 10 → 200 (Boost bisects when a
  Newton step leaves the bracket, so easy inputs are unaffected). Removed the
  old retry loop that shrank `lb` on `[l(b*)]'' > 0` — a workaround for the bad
  bracket that pulled the search off the root once the bracket was correct.
  Verified against an independent `stats::optimize` at 20 (max(s), beta)
  combinations: 18/20 missed before (worst 84%), 0/20 after (all ≤6e-5).
- Consequence worth remembering: `gamma_rate` and `log_gamma_rate` now agree to
  1e-12 across beta = 0.1–10, where they differed 20× above beta = 1. That is
  almost certainly why the equivalence half of `test_gamma_rate.R` was
  commented out; it is now writable.
- **`.check_cohort_is_testable()` added** to `R/compute_test_statistic.R` and
  called from `compute_test_statistic.default`, `.compute_df` and
  `compute_test_per_gene` — all three, so the two pipelines cannot diverge on
  the inputs the guard exists for. Errors on (i) an arm with <2 individuals
  (the real numerical failure: df = 0, then NaN everywhere) and (ii) a donor
  with < `min_cells_per_individual` cells, default 3, `0` disables.
- Recorded that (ii) is a POLICY not a numerical necessity: verified a 1-cell
  donor yields perfectly finite statistics even with zero posterior variance.
  Kevin's call, now explicit in the code rather than implied.
- **The external `stats::t.test` oracle survived the new guard.** Rebuilt the
  fixture as 3 identical cells per donor with zero posterior variance — the
  reduction depends on the zero variance and identical cells, not the cell
  count. Still exact to 1e-10. This was the risk I flagged back in plan §2.8.
- **Bessel: Kevin says no correction**, so the population variance is deliberate
  and the code is right. T-DF-01 deleted (wrong oracle); T-TSTAT-01a keeps the
  contract `teststat_eSVD = teststat_Welch * sqrt(n/(n-1))`, exact at n = 3, 5,
  10, 25. My "Finding 2" was not a defect.
- **T-PROP-08 cut.** Reason recorded in the file: multiplying counts by c is not
  sequencing c times as deep. `Poisson(c*l*lambda)` has mean AND variance
  `c*l*mu`; scaling observed counts gives variance `c^2*l*mu`. The rescaled data
  violates the model's mean-variance relation, so nothing is preserved. A real
  Poisson invariance would be thinning.
- **Open, newly exposed: `log_gamma_rate` is now the binding constraint.** Its
  bounds are log(beta) in [-10, 10], so beta <= 22026. On the fixture one gene
  of twenty is essentially Poisson (over-dispersion 2.6e-7, MLE 3.8e6) and the
  log route returns exactly exp(10) while `gamma_rate` now returns the MLE. The
  asymmetry has flipped; the two are not interchangeable at the extremes.
- **Q-NUIS-1's ±30% is about the estimator, not our code.** Measured median
  relative error 0.257 at 400 cells/gene, 0.107 at 2000, 0.053 at 10000 —
  gene-by-gene ±30% needs ~10000 cells, which no fixture affords in 90 s. The
  Gamma rate is weakly identified when over-dispersion is small (flat
  likelihood). T-NUIS-05 split: 05 asserts our code reproduces an independent
  `stats::optimize` (exact, and it is the assertion that would have caught the
  bracket cap); 05b checks the median across genes. Both pass.
- All ten self-audit defects repaired. Two remain unfixed and are noted:
  T-VAL-34 and T-ESVD-11 read package source via `test_path("..","..")` and will
  silently skip under `R CMD check`.
- Corrected an earlier wrong claim: `.reparameterize` **already** warns and
  proceeds on rank-deficient input (Q-REP-1's desired behaviour). The real
  defect is one level up — `reparameterization_esvd_covariates` warns, proceeds,
  produces non-finite values, then dies in `eigen()`.
- `devtools::document()` swapped `RoxygenNote: 7.3.3` for
  `Config/roxygen2/version: 8.1.0` in DESCRIPTION. Reverted — a local-toolchain
  artifact, not a change Kevin asked for. `man/` regeneration was kept.

### [2026-09-01] (Session 9 — the code fixed to the reviewed suite)

- Kevin: "I've checked all the unit tests. Can you go ahead and fix the
  codebase?" Taken as the go-ahead. Suite went 522 pass / 37 fail / 28 skip →
  **644 / 0 / 0**; `R CMD check --as-cran` 5 WARN / 3 NOTE → 1 WARN / 3 NOTE
  (the WARN, `MASS` in tests, fixed after the run). Uncommitted.
- Full fix table is `TEST_RUN_REPORT.md` Part 7. Non-obvious rationale only
  below.
- **Efron's truncated MLE**: the optimizer was fine; the model dropped the
  `p0 <= 1` constraint by treating `theta` as free. Verified the old code's
  2.57 was the true maximum of its objective via an independent
  reimplementation before changing the model. Constrained version verified
  across 50 seeds at N = 200 and 1000.
- **`.multtest_simple`** corrected by the moments of a normal truncated at the
  null quantiles. Exact when the window is all-null; mildly high otherwise
  (1.26 on a 10% mixture at ±4). It is the last fallback, so acceptable.
- **Non-determinism found**: identical `eSVD()` runs differed by up to 8 in
  Welch statistics, from a 2e-9 `irlba` start-vector difference. Fixed the
  reproducibility (`.svd_start_vector()`), NOT the sensitivity; recorded as
  open. This was what made T-COH-08 and T-GS-05/09/10 fail, not the helper.
- **Seurat renames `gene_1` → `gene-1`.** Fixture genes are now `gene1`. This
  was also the real cause of T-COH-13's `subset()` error.
- **Five tests were unsatisfiable as written** and were changed (report §7.2):
  T-PVAL-01 recomputed the naive form inline; T-PVAL-01b contradicted T-MT-04;
  T-MT-02a pinned the buggy values and contradicted T-MT-02; T-MT-08 pinned the
  bias itself; T-MPFR-08 wanted two zero doubles to differ; T-COH-04 indexed
  alphabetically-sorted factor levels. Plus T-SVD-05c shadowed the package
  function with a local copy (code review).
- **`bool_diet = TRUE` keeps `fit_Second`** — a deliberate behaviour change
  because the spec and T-GS-11 presume it; flagged for Kevin.
- **`eSVD()`'s `grep(id_var, ...)`** (unchanged line) dropped any covariate
  column containing `id_var` as a substring; now removes the exact
  `<id_var>_<level>` names. Found by the code review.
- `eSVD_helper` forwards `min_cells_per_individual = min_cells_per_id` unless
  the caller passes it, so `min_cells_per_id = 0` really disables the check.
- `compute_test_per_gene` had the `min_cells_per_individual` argument but
  never used it (session 8 claimed all three call sites); now it does.
- `multtest()` errors (not warns) when no estimator yields a finite null;
  `null_sd == 0` still proceeds with a step-function p-value (T-MT-07).
- Code-review findings NOT fixed, recorded in `CLAUDE_kevin.md` Open
  Questions 1–4: the nuisance blow-up after the bracket fix (10/20 fixture
  genes above 1e4, stabilization rescales all by 10^-4.55), locfdr's misfit
  warning triggering the fallback on real data, the 4× duplicated
  case/control derivation, and the helper's memory duplication.
- Mechanical CRAN items done in passing: `Rmpfr` and `sparseMatrixStats` out
  of `DESCRIPTION` and code (`.sparse_col_sds()` replaces the latter), `Rcpp`
  moved to Imports, `biocViews` removed, `withr` added to Suggests, `LICENSE`
  two-line stub, `override` ×5, `exportPattern` removed with explicit
  `@export` on the 17 public functions, `requireNamespace("SeuratObject")`.
- Build needs pandoc. On the personal laptop (macOS) the shell PATH has none;
  RStudio's copy at
  `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64`
  works when prepended before `R CMD build`.
- `devtools::document()` again swapped `RoxygenNote` for
  `Config/roxygen2/version: 8.1.0`; reverted again.

### 2026-09-02 (Session 10 — test warnings, CRAN polish pass, readiness doc §0)
- Kevin asked for the two test-suite warnings fixed and a CRAN polish pass,
  with questions collected in `CRAN_READINESS.md`. Both warnings were tests
  that produced the intended warning without asserting it (`T-REP-09`,
  `estimate_nuisance.default` extreme-`mu`); now `expect_warning()`.
- `R CMD check --as-cran` with the vignette built: `Status: 2 NOTEs`
  (`New submission`; HTML Tidy missing locally). Tarball 13.2 MB → 1.0 MB.
  Suite 687 / 0 / 0 warnings / 1 skip, 20 s. Nothing committed.
- Found and fixed four correctness defects not in the audit: (N1) the
  single-covariate GLM fallback in `.initialize_coefficient()` dropped the
  library-size offset; (N2) `reparameterization_esvd_covariates()` via
  `lm(x ~ ., as.data.frame())` mangled non-syntactic names and propagated
  `NA` for aliased columns — replaced by a QR solve, numerically identical
  (4e-16); (N4) `eSVD()` errored on a non-factor categorical variable
  (`droplevels`); (N5) §1.4's `scale()` type change was still present in
  both posterior paths. Plus (N3) `initialize_esvd()` now refuses a
  rank-deficient design by name, (N6) §4.2 `Rcpp::warning` replaced by a
  counted flag and one R warning, (N7) NA/non-finite guards in
  `opt_esvd.default`.
- Decision: §5.3 (`print()` → `message()`) resolved as *no change*; CRAN's
  reviewer text accepts `if(verbose) cat()`, and the style guide mandates
  `print(paste0())`.
- Decision (Claude's, reversible, flagged as Q4): ASD tutorials moved to
  `vignettes/articles/` (pkgdown articles, `.Rbuildignore`d); their PNGs
  included by relative path; `EnhancedVolcano` and `devtools` dropped from
  `Suggests`; `rmarkdown` added (was missing).
- Decision (Claude's, flagged as Q9): legacy fixture regenerated at
  20 individuals × 20 cells × 60 genes with xz (432 KB), keeping only what
  tests read; `initialize_esvd1.rda` was loaded by no test and is deleted.
  Legacy `initialize_esvd works` now drops the individual indicators
  (rank-deficient by construction); legacy per-gene equivalence test is
  element-wise with matched `library_min`/`alpha_max`.
- Added `@examples` to all exported functions (40-gene `generate_null` chain,
  ≤ 0.55 s each; `locfdr` converges at 40 genes, fails at 18). Seurat ones
  guarded by `requireNamespace()`, two in `\donttest{}`.
- Empirical finding: a 4e-16 change in `z_mat` flipped `locfdr`'s
  convergence on the 18-gene fixture (eleven gene-status tests went from
  silent to warning). Handled with `.muffle_locfdr_fallback()` in the test
  helpers; recorded as new evidence for the amplification item.
- `devtools::document()` again swapped `RoxygenNote` for
  `Config/roxygen2/version: 8.1.0`; reverted again (Q11).
- Open: Q1–Q12 in `CRAN_READINESS.md` §0.3; suggested tests §0.4.

### 2026-09-06 (Session 11 — folding the vetting history back into the `vet-r-package` skill)
- No package code, tests, or docs changed this session; the eSVD2 tree is
  exactly as session 10 left it (still uncommitted on `devel`).
- Kevin asked for the sessions 1–10 findings to be written back into the
  shared `vet-r-package` skill (`claude_skills` repo, sibling of this
  project's Dropbox folder). Edited in place, uncommitted there: `pitfalls.md`
  173 → 473 lines, plus targeted additions to `SKILL.md`, `plan-template.md`,
  `cran-readiness.md`.
- Decision: entries are written by bug class in the skill's generic
  vocabulary (unit / observation / feature), never as eSVD2 narrative, so the
  file stays a checklist rather than a session log. Policy items (CRAN
  findings) went to `cran-readiness.md`; bug classes to `pitfalls.md`;
  process lessons (check before reading source, oracle field on every test
  bullet, multi-pass review folding, `code-review` on the fix diff, tests
  for unbuilt features skipping on `exists()`) to `SKILL.md` and the plan
  template.
- Checked the new file with two fresh agents on a planted-bug snippet, one
  given the file and one not. The control found the local bugs on its own;
  the file's marginal value was the fixture- and cohort-level classes
  (self-reading fixture, discarded `method`, fixture below `locfdr`'s range,
  filter order, per-arm counting) and the regime axes. The with-file agent
  flagged five entries as ambiguous and five bug classes with no entry
  (user string as regex, dependency return object indexed by position,
  design built from all metadata, dependency default slot, shared downstream
  stage not covered by an equivalence test); all ten fixed before finishing.
- Open: the skill edits are uncommitted in `claude_skills`; Kevin reviews
  and commits there.

### [2026-09-28] (Session 12 — log2 fold change and its standard error; version 1.1.0)
- Kevin asked for the fitted model to report a log fold change and its SE,
  specifically proposal 2 of the Was2CoDE wiki page `code-esvd2.md` (the
  log2-scale delta-method SE), folded into the ordinary output, with tests he
  will vet himself. Not wanted: the model-coefficient SE, and the bootstrap
  as a feature (it is a test oracle only).
- Found at startup: sessions 9-11's work is committed (`a95e533`), the tree
  was clean, and the suite baseline is 688 pass / 0 fail / 0 skip, not the
  687 / 1 skip the state file recorded.
- Decisions, approved by Kevin as a plan-mode plan (D1-D7): the SE uses the
  mixture variance the Welch statistic uses, not the variance of individual
  means alone; names `log2fc_vec` / `log2fc_se_vec` / `logFC_se`; a
  non-positive arm mean gives `NA` plus a warning; an object saved before
  1.1.0 is refused by `compute_log_fold_change()` but still reported by
  `report_results()` with `logFC_se = NA`; version 1.1.0.
- Rationale for storing `case_var` / `control_var` rather than recomputing:
  a `bool_diet = TRUE` object has dropped `dat` and the posterior matrices,
  so nothing is left to recompute from. For the same reason
  `compute_test_per_gene()` now records the arms' individuals in `param`,
  which only the matrix path did.
- Finding: the within term (average per-cell posterior variance, undivided
  by cells) is a median 98% of `log2fc_se^2` on F-SMALL and 95% on
  `generate_null()`; the SE is 7.8 and 4.3 times the SD of a fit-fixed
  bootstrap over individuals.
- Finding: against 30 *refitted* `generate_null()` cohorts the SE is about
  right on genes with an ordinary nuisance estimate (394 gene-cohorts:
  calibration ratio 0.774, coverage of +/- 2 SE 0.959). So the within term is
  what accounts for the fit's uncertainty, even though it is not a
  donor-sampling variance.
- Finding: on genes whose nuisance estimate diverges (56 gene-cohorts, rate
  near 1e7 against true rates of 0.1 to 10) the SE is anti-conservative:
  calibration ratio 1.538, coverage 0.357. On F-SMALL 7 of 40 genes diverge
  and 5 of them miss their truth by more than 3 SE (planted gene 2 by 12.5).
  This is the measured consequence of Q10.
- Finding: the delta method understates the between-individual SD by up to
  20% at a between-individual CV near 1 with 4 individuals per arm; within
  CV <= 0.2 the bootstrap / delta ratio is 0.98 to 1.03.
- Finding: the depth adjustment biases every fold change when DE genes are a
  large share of the counts (F-SMALL: null genes at a median of -0.21).
- Finding: `generate_null()`'s "null_large_var" genes have a true log2 fold
  change of 0.40 for a ratio of arithmetic means; they are null only for a
  difference of mean logs.
- Open, not investigated: the Welch statistic is also a contrast of
  arithmetic means, so those genes are not null for eSVD2's own test either,
  yet T-PROP-06 counts genes 11 onward as null when it checks that p-values
  are uniform. Whether that test passes because the effect is small beside
  the SE or because the empirical null absorbs it has not been checked.
- Three corrections to my own tests after their first run, all disclosed in
  `UNIT_TEST_PLAN.md` section 2.18: T-LFC-15's tolerance was applied outside
  the CV regime it was derived for, and its "SE never below the bootstrap"
  assertion was not a property of the statistic; T-LFC-16 failed on the one
  planted gene with a diverged nuisance and now restricts its 3-SE assertion
  to not-diverged genes; T-LFC-18 was added for the mechanism. A first
  attempt to fix T-LFC-16 by "composition-correcting" the truth was wrong
  (it assumed a `Log_UMI` coefficient of 1; `generate_null()` fits 0.24) and
  was dropped for the generator's own truth.
- T-LFC-02 passed before any implementation existed, because
  `expect_equal(NULL, NULL)` succeeds; every section-A call now goes through
  a helper that checks the fields' lengths.
- Teeth: ten deliberate breakages in a scratch copy, each caught by at least
  one test (table in section 2.18).
- Two existing test files were edited, one more than the plan said:
  `test_report_results.R:26` (the pinned column set) and T-TSTAT-06 in
  `test_compute_test_statistic_claude.R`, whose mean-zero Gaussian fixture
  now triggers the D4 warning and asserts it.
- Comparison script `additional_context/lfc-se-comparison_2026-09-28_claude.R`
  (one cohort, 10 vs 10 individuals, 100 genes): fold changes agree with
  DESeq2 / dreamlet / NEBULA at Spearman 0.965 to 0.978; eSVD2's SE is a
  median 1.81 / 1.87 / 2.05 times theirs; calibration ratio 0.41 for eSVD2
  against 0.87 / 0.94 / 0.80.
- `NEWS.md`: the 1.0.2 section is relabelled "Development version, not
  released", since "First CRAN submission" now belongs to 1.1.0.
- Open: Q-LFC-1 (what to do about the diverged-nuisance SE), Q-LFC-2 (the
  warning from the matrix method), Q-LFC-3 (names).
- `/code-review` on the diff returned nine findings; seven were fixed.
  Fixed: `report_results()` recomputed `logFC` from the arm means, so it
  could show `Inf` beside an `NA` SE; the arm sizes in `param` went stale on
  a rerun because `.combine_two_named_lists()` never overwrites; an invalid
  gene got an SE of `NaN` rather than `NA` (inputs were blanked, not
  outputs); a single `NA` input exempted a gene from validation; T-LFC-15
  and -16 put thresholds on a chaotic fit without `skip_on_cran()`; the ASD
  article read `log2fc_vec`, which a pre-1.1.0 saved object lacks (reverted
  to the expression that works on every object); README still said 1.0.2.
  Not fixed, by decision: the warning from the matrix method (that is
  Q-LFC-2), and `.compute_df()` recomputing the variances (the optional
  refactor).
- T-LFC-08b and T-LFC-19 were added and T-LFC-12 extended for the review's
  findings. Teeth were re-run on the final code with 14 breakages.
- Result: suite 996 pass / 0 fail / 0 warnings / 0 skip under
  `NOT_CRAN=true` (688 before); `R CMD build` + `R CMD check --as-cran` on
  `eSVD2_1.1.0.tar.gz` gives `Status: 2 NOTEs` (`New submission`; HTML Tidy
  not recent enough on this machine), tests OK, examples OK, vignette
  rebuilt OK, tarball 1,044,809 bytes.
- `devtools::document()` again rewrote `RoxygenNote: 7.3.3` to
  `Config/roxygen2/version: 8.1.0`; reverted each time (Q11).
- Open: whether `param` going stale on a rerun should be fixed for every
  stage in `.combine_two_named_lists()`.
- Nothing committed; session 12 is uncommitted on `devel`.

### [2026-09-29] (Session 13 — master (3d5f7bf) vs devel (1.1.0) comparison on simulated data)
- Built `additional_context/version_comparison/`: private libraries for both versions, six simulated regimes x 10 replicates (300 genes, 20 individuals x 30 cells), 12 corner cases, and a knitted report `version_comparison_claude.html`; `run_all_claude.sh` takes 385 s end to end.
- Both versions are driven by one script, because the step-by-step API has identical names and arguments; `eSVD()` was not used since it needs Seurat.
- Folder `output/` rather than `results/`, because the repo's `.gitignore` ignores every `results/`; `lib/` and `data/` are gitignored as regenerable.
- The `devel_swap` ablation (devel code, master's nuisance rates) reproduced master's Welch statistics on all 60 data sets, attributing the whole ordinary-data difference to the `gamma_rate` change.
- Finding: logFC unchanged (r >= 0.998); devel's null false discoveries at FDR 0.05 rise from 0.2 to 2.7 per 300 genes, and to 63 in the near-Poisson regime (type-I 0.28 at p < 0.05).
- Finding: master's nuisance estimates have Spearman about 0 with the truth; devel's about 0.6, but biased upward two- to threefold and divergent on some genes.
- Finding: the corner-case fixes behave as NEWS says (collinear/confounded covariates, all-zero genes, one-individual arms, too-few-cell individuals now refused early by name); master silently called 30/30 genes significant on a 30-gene cohort.
- Finding: NEWS's sparse-NA zeroing is unreachable through `format_covariates()`, which already puts NA in `Log_UMI`.
- Finding: determinism was not a practical issue in master (1e-11 difference across seeds); devel is bit-identical.
- Finding: the logFC SE covers truth at >= nominal for ordinary-rate genes except in `strong_de` (0.84, the Log_UMI depth bias); coverage collapses (0-0.47) for diverged-rate genes.
- Open: what to do about `gamma_rate` (Q10), now the top decision before CRAN.

### [2026-09-29] (Session 14 — brainstorm on the inflated p-values from the uncapped nuisance rate)
- Wrote `additional_context/OVERDISPERSION_BRAINSTORM.md` (eleven ideas, each with a go/no-go result or task) at Kevin's request that it live in `additional_context/`; `brainstorming_kevin.md` holds only a pointer to it.
- Built `additional_context/overdispersion_brainstorm/`: devel is fitted once per data set and cached, and each candidate only overwrites `nuisance_vec` before `compute_posterior()`, so 18 candidates on 90 data sets take about 15 s.
- Added three regimes (`weak_de`, `low_count`, `wide_rate`) because the six of session 13 find every DE gene under every candidate and cannot separate them on power.
- Finding: master's rate is `min(MLE, max_i s_ji)`, and 78% to 100% of genes are at that cap in the model-generated regimes.
- Finding: devel's "two- to threefold overestimate" (session 13) was a units mismatch between the generator's rate and the fit's; unit-free, the excess is 4% to 16%.
- Finding: a diverged rate is the boundary of the likelihood (Poisson at least as likely as any finite rate; Pearson statistic 0.95 to 0.99), matching the 49.9% divergence of `lause-2021` in the overdispersion wiki.
- Finding: the true rates give 5.7 false discoveries per 300 genes when they span 0.5 to 200, so the test is calibrated by bounding the spread of the rates and not by estimating them better.
- Finding: capping at `c * median_i s_ji` is flat in false discoveries for c from 1 to 10; empirical-Bayes shrinkage and a 90% profile-likelihood lower bound perform the same.
- Finding: a Welch variance between individuals only has null SD 1.5 under the theoretical null and loses most power under the empirical null; rescaling the statistic within tenths of the rate is worse than not rescaling.
- Decision: the cap was prototyped in R after `gamma_rate`, not in C++, so that `gamma_rate` stays an MLE and T-CPP-GAM-05 keeps its meaning.
- The first version of the empirical-Bayes prior centred on the median over all genes, which is infinite once half the genes are at the boundary; it now uses the interior genes, which biases the centre low in near-Poisson data.
- Added location `OVERDISPERSION_WIKI` (read only from this project) to the master `CLAUDE.md` and its path to `CLAUDE_kevin.md`.
- `.gitignore` now excludes the brainstorm's `data/`, `cache/` (258 MB), `output/genes.csv` and `output/gene_summaries.csv`; the summary tables and logs stay tracked.
- Open: which route to take (Q10); the real-data check (Idea 11) has not been run because no `PAPER_DATA` path is recorded.
- Nothing in `R/`, `src/` or `tests/` was changed, so no `R CMD check` was run this session.
- Follow-up in the same session, prompted by Kevin asking how the empirical-Bayes idea was implemented: measured what the prototype did, and added `06_more_cells_claude.R` (three regimes at 150 cells per individual).
- Finding: at 600 cells the MAP rate is a soft cap near 10 times the library size (prior variance 0.25 against a sampling variance of 0.096), which is why it matched the hard cap.
- Finding: at 3000 cells with true rates of 0.5 to 200, shrinkage leaves 6.2 false discoveries against the cap's 3.0, so Ideas 3 and 4 were downgraded from "go" to "go only with a cap" in the brainstorm.
- Finding: at 3000 cells with true rates of 2 to 8, no gene is at the boundary and the uncapped MLE is calibrated.
- The prototype differs from DESeq2's design in four ways: no trend in expression (the centre is one unit-free constant), no Cox-Reid term, a Wald-type sampling variance in place of the trigamma formula, and no rule exempting outlying genes.
- Second follow-up, on Kevin's question of why the prototype departs from DESeq2 and his criterion that the estimate must improve on master: added `07_deseq2_variants_claude.R` and section 3.5 of the brainstorm.
- Finding: scored as estimates, the cap at 10 and the shrinkage put 95% and 99% of genes within twofold of the truth against master's 60%; in near-Poisson data nothing does.
- Finding: a fitted trend centre and DESeq2's trigamma sampling variance leave false discoveries unchanged; the simulations have no trend by construction, so whether one exists in real data is open (Idea 11).
- Finding: shrinkage followed by the cap gives exactly the cap's results at 3000 cells in `wide_rate`.
- Open: Kevin's criterion has two parts that can disagree (accuracy of the rate, calibration of the test); the cap satisfies both, the shrinkage alone only the first at 3000 cells.
- Decision (Kevin): the fix is the cap at 10 times the gene's median library size; it is not to be implemented in `R/` until he has read a report explaining it.
- Wrote and knitted `additional_context/overdispersion_brainstorm/overdispersion_cap_claude.Rmd`: takeaway, formulas, what master and devel each do, the change, the comparison with master, five nuances, recommendation; citations are quotations as the overdispersion wiki's footnotes record them, not rechecked against the PDFs.
- The report reads only small tracked tables (`output/report_*.csv`, made by `08_report_tables_claude.R`), because the per-gene tables it would otherwise need are gitignored.
- Added the `trend` regime (`helpers_simulate_claude.R`), a copy of the version comparison's generator with the rate set by DESeq2's parametric trend; the copy was made so that the session-13 helper and its seeds stay untouched.
- Ran master itself on the four added regimes (`01b_run_master_claude.R`), so "master" in the report is never the legacy cap standing in for it; the two agree.
- Finding: in `trend` the cap beats master on false discoveries (1.0 against 2.3) and true discoveries (21.7 against 20.7), and a cap as tight as master's is worse there (2.8 to 3.5 false discoveries at c of 2 and 1).
- Finding: the boundary condition is `D = sum [(A - m)^2 - A] / mu <= 0`, from the first-order expansion of the log-likelihood around Poisson; it matches devel's output for 11,997 of 12,000 genes.
- Correction: the claim that the true rates do not calibrate the test held only for genes with a true rate above about 30 times the library size; it had been drawn from `wide_rate` alone.
- Correction: the claim that `bool_stabilize_underdispersion` rescales in every regime was unchecked and wrong for `generate_null()` (and for `low_count` under master).
- Reference DOIs in the report were checked against Crossref; the one for Dai, Bao and Bao recalled from memory was wrong and was replaced by 10.1016/j.spl.2012.08.017.
- Open: shrinkage toward a fitted trend at 3000 cells and with a cap as backstop; master at 3000 cells; one real data set.

### [2026-09-29] (Session 15 — the cap implemented as version 1.2.0: `cap_multiplier`, `recompute_pvalue()`, the diagnostic plots)
- Kevin read the brainstorm and chose Idea 1, `min(MLE, c * median_i s_ji)` with `c = 10` a user argument, plus a redo of the p-values at another `c` on a frozen fit and diagnostic plots; nothing under `additional_context/` is to be touched until he has vetted the code.
- Decisions from Kevin before the code was written: a `bool_diet = TRUE` object gets its counts back by the user passing the Seurat object again; the plots draw the unit-free rate `beta_j / median_i s_ji` on a log10 axis with the cap at `c`; the per-gene status is a factor on the fit and a column of `report_results()`; master's exact bound `max_i s_ji` is not offered, `c = 1` is documented as close to it.
- The laptop crashed after the code and tests were written and before documentation, the by-hand check, `R CMD check`, the review and this log; the approved plan survived as `~/.claude/plans/validated-prancing-kazoo.md` and the session resumed from it.
- Finding (by looking at the rendered plot, which no test had done): `min_val` is an absolute floor applied after a cap that is in units of the library size, so a gene whose median library size is about 1e-8 has the rate `min_val`, thousands of times above its cap. 2 of 100 genes on `generate_null()` (40 cells per person, 12 individuals, `k = 3`): their `Log_UMI` coefficients are -4.0 and -1.8 against a median of 0.37, that is, the fit of the gene is degenerate.
- Resolved: the rule `max(min(MLE, cap), min_val)` is kept as planned, because `compute_posterior()` floors the library at `library_min = 0.1` and such a gene's posterior does not follow the fit (statistics 0.6 and 1.0); `plot_nuisance()` now counts such genes in its subtitle and both help pages say why.
- Finding (code review): a boundary gene reached the cap only when the optimizer's stopping value was above it, and that value is exactly `exp(10)` on the log route; with a median library size above about 2200 all 30 boundary genes of a simulation kept 22026 and were not counted as capped.
- Decision (Claude, for Kevin to confirm): a boundary gene is set to the cap itself, the reading of `min(MLE, cap)` with an infinite MLE; with `cap_multiplier = Inf` it keeps the optimizer's value as in 1.1.0. This is the one change to the rule of the approved plan.
- Finding (code review): `compute_posterior.eSVD()` recorded its settings with `.combine_two_named_lists()`, so after a second call `param` described the first and the redo at the same cap returned statistics up to 0.99 away; it now overwrites its entries, as the test functions have since session 12.
- Finding (code review): the rebuilt covariates were recognized by column sums, which do not change when a covariate is exchanged between two individuals with as many cells; the redo then ran silently on the wrong design. `eSVD()` now also records sums weighted by `cos(sqrt(2) * i)` over cells, for the covariates and for the counts.
- Open: the reviewer's root cause for the previous item is that `bool_diet` drops the n x r `covariates` while keeping the n x k `x_mat`; keeping `covariates` on a diet object would remove the rebuild of the design and its check altogether.
- A test written to pin the stale-`param` bug was vacuous on its first draft: `alpha_max = 5, library_min = 2` does not move any statistic on F-TINY (every library size is above 2.4, and `alpha_max` from 1 to 1000 changed nothing). Its own precondition caught this; it now uses `library_min = 20, pseudocount = 1`.
- Open: `alpha_max` had no effect on any statistic of F-TINY at 1, 5, 86 or 1000; whether it is ever binding was not investigated.
- Mutation sweep in a scratch copy, 26 deliberate breakages of the new code: 25 turn a test red. The one that does not is equivalent (reversing the gene order consistently on both axes of `plot_fitted_vs_observed()`).
- The sweep found one gap, now closed by T-DIAG-08b: no test saw a count marked from below, because on counts drawn from the model the fit minus 3 SD is almost always negative and the bar's lower end is 0.
- Not done, from the review, as refactors of files this work did not own or as performance: three copies of the rule that picks the library columns (`.nuisance_library_idx()`, `compute_test_per_gene()`, `compute_posterior.default()`), and the second pass over the genes in `.compute_boundary_statistic()`.
- Open: `plot_fitted_vs_observed()` calls `set.seed(seed_number)` with default 10, as the house style prescribes, which resets the caller's stream when the plot is drawn inside a simulation loop.
- Verification: suite under `NOT_CRAN=true` 1868 pass / 0 fail / 0 warnings / 0 skip in 66 s across 33 files; `R CMD check --as-cran` on `eSVD2_1.2.0.tar.gz` gave `Status: 2 NOTEs` (`New submission`, HTML Tidy); tarball 1.55 MB, up from 1.0 MB because the vignette now draws five figures.
- On `generate_null()` (100 genes): 11 genes capped at `c = 10`, 4 of them at the boundary; the largest Welch statistic is 6.7 at `c = 10`, 6.5 at `c = 1` and 30.7 without a cap, with 8, 9 and 12 genes at FDR < 0.05; 1.6% of cell-gene pairs lie beyond 3 SD.
- Test IDs and oracles are in the headers of the three new test files and not in `UNIT_TEST_PLAN.md`, by Kevin's instruction; nothing is committed.

### [2026-09-29] (Session 16 — Kevin's first answers on 1.2.0: the floor `min_val` moves into the units of the cap)
- Kevin asked for an explanation of the boundary-gene rule (given, decision still his), confirmed that a diet object must not keep `covariates` ("as small as possible"; the code already drops them, so the Seurat rebuild and its check stay), asked whether `min_val` was in library-size units (it was not: an absolute `1e-4` applied after the cap) and directed that it should be, and asked what the two deferred refactors were (explained, not requested).
- Resolved: the rate is now `max(min(MLE, c * m_j), min_val * m_j)` with `m_j` the gene's median library size; `.check_min_val()` refuses a floor that is not one positive finite number below `cap_multiplier`, in both `estimate_nuisance()` methods and in `recompute_pvalue()`, which reads the recorded `nuisance_min_val`. The default `1e-4` is unchanged in value but is now a multiple of `m_j`.
- Consequence: no gene can sit above its cap, so the "above the cap" count in the subtitle of `plot_nuisance()` and its help-page paragraph are gone; a failed gene now gets `min_val * m_j` rather than `min_val`, and the warning says so.
- Tests: T-CAP-03b rewritten as a grid over `(cap_multiplier, min_val)` asserting every unit-free rate lies in `[min_val, cap_multiplier]`; new T-CAP-03c (the rule by hand on three genes with medians 0.5, 2, 8), T-CAP-03d (the argument check, nine bad values, both methods) and T-REDO-13 (a redo at a cap at or below the recorded floor is refused); T-CAP-06 and T-DIAG-03b updated to the new floor. Test IDs stay in the file headers, not in `UNIT_TEST_PLAN.md`.
- Not changed, for Kevin to decide: the boundary gene set to the cap itself (item 0a). Under the plan's literal `min(MLE, cap)` such a gene keeps the optimizer's stopping value whenever that value is below its cap, which is exactly `exp(10)` on the log route and about 1e7 on the first, that is, whenever the median library size is above about 2200 or 1e6.
- Not done: the two refactors, which are quality and speed, not correctness. (i) The rule that picks the library columns of `covariates` exists three times: `.nuisance_library_idx()` (used by `estimate_nuisance.eSVD()` and `plot_fitted_vs_observed()`), an inline copy in `compute_posterior.default()` and another in `compute_test_per_gene()`; the older two wrap it in `unique()`, which `which(%in%)` makes irrelevant, so the three agree today but can drift. (ii) `.estimate_nuisance_matrix()` walks the genes twice, once for the MLE and once in `.compute_boundary_statistic()` for `D_j`, extracting each column of `dat` (a `dgCMatrix` on real data) both times; `D_j` could be computed in the first loop.
- Verification: suite under `NOT_CRAN=true` 1892 pass / 0 fail / 0 warnings / 0 skip in 64 s across 33 files (one new expectation first failed on a names mismatch in the test itself, fixed); `R CMD check --as-cran` on the rebuilt `eSVD2_1.2.0.tar.gz` gave `Status: 3 NOTEs`: the two of session 15 plus `unable to verify current time`, a clock-service NOTE of the check environment. `devtools::document()` again rewrote `RoxygenNote` and was reverted. Nothing is committed.

### [2026-09-29] (Session 17 — Kevin confirms the three 1.2.0 decisions; the two refactors, pinned by tests first)
- Kevin confirmed: the boundary gene is set to the cap itself (item 0a closed), the diet object does not keep `covariates` (0c), and the floor `min_val` in the units of the cap as implemented in session 16 (0d). He asked for both deferred refactors, with unit tests written first so the refactor can be shown not to change behaviour.
- Tests written before the refactor, all passing on the old code: T-CAP-10 (`.nuisance_library_idx()` against eight column sets written out by hand, `.library_column_oracle()` in `helper-fixtures.R`; F-TINY's columns are Intercept, Log_UMI, Age, CC_1, Sex_M, so the case-control column sits between library columns), T-POST-13 (`compute_posterior.default()` components equal a posterior built by hand from those column sets, with `alpha_max = NULL`, `library_min = NULL`, `nuisance_lower_quantile = 0` and no stabilization, so the library columns are the whole test), T-PGENE-01 in the new `test_compute_test_per_gene_claude.R` (the per-gene path equals the matrix path on the six admissible settings of the three library booleans), and T-CAP-11 (`.estimate_nuisance_matrix()` returns the same list on dense and sparse counts, both routes).
- Two of my first-draft expectations were wrong about the test, not the code: the helper returns column positions in design order, not sorted, and a fixture defined inside one test file is invisible to another (moved to `helper-fixtures.R`).
- Beside the tests, a golden snapshot outside the repo (`golden.R` in the session scratchpad): 25 outputs of `.estimate_nuisance_matrix()` (dense, sparse, both routes, an NA column, an Inf entry, `cap_multiplier = Inf`), `compute_posterior.default()` (12 settings) and `compute_test_per_gene()` (6 settings), saved before the refactor and bit-identical (`identical()`) after it.
- Refactor (i): `compute_posterior.default()` and `compute_test_per_gene()` call `.nuisance_library_idx()`; the inline copies and their `unique()` are gone. Refactor (ii): `.estimate_nuisance_matrix()` uses one `vapply` over the genes returning `c(rate, D_j)` (a 2 x p matrix even for p = 1), and `.compute_boundary_statistic()` takes one gene's `x_vec, mu_vec, s_vec`; T-CAP-04b was rewritten to the per-gene signature, with the same definition as its oracle.
- Not measured: the speed gain of (ii); it is one fewer column extraction per gene, which matters on a `dgCMatrix`.
- Verification after the refactor: suite under `NOT_CRAN=true` 2014 pass / 0 fail / 0 warnings / 0 skip in 63 s across 34 files (up from 1892 across 33); `R CMD check --as-cran` on the rebuilt `eSVD2_1.2.0.tar.gz` gave `Status: 3 NOTEs`, the same three as session 16. Nothing is committed.

### [2026-09-29] (Session 18 — `additional_context/` brought up to 1.2.0; the comparison and the cap report rerun)
- Kevin asked for every `.md` under `additional_context/` to be updated for 1.2.0 (vetted and committed as `1a0a536`), a "Decision" section in the brainstorm, and both knitted reports rerun against the current code; nothing in `R/`, `src/` or `tests/` changed.
- Decision (Claude): `version_comparison/` gained a `devel_nocap` arm (1.2.0 at `cap_multiplier = Inf`), so the first edition's comparison survives as a reference; its statistics equal the old 1.1.0 run to every stored digit.
- Decision (Claude): the 1.1.0 install moved to `lib/devel_1.1.0` and the brainstorm's scripts 01 to 07 now load it, because reinstalling `lib/devel` as 1.2.0 would silently have changed what they reproduce.
- Finding: master vs 1.2.0 false discoveries are within noise in every regime (null 0.2 vs 0.2, `near_poisson` 1.8 vs 2.6, `baseline` 2.1 vs 1.8); power is lower only on `generate_null()` (8.1 vs 8.8 of 10), which the uncapped arm shares (8.2).
- Finding: under the cap, ±2 SE covers the true log2FC for at least 99% of genes of every status in five of six regimes (`strong_de` 0.72 to 0.85), against 0 to 0.47 for diverged genes without it; this answers the measurement half of Q-LFC-1.
- Finding: `09_run_v120_claude.R` (new) shows 1.2.0 end to end calls every gene of all 100 brainstorm data sets as the prototype `cap_10s` did; the boundary-gene rule never bound, because no median library in these simulations exceeds about 8.
- Finding: master's Welch statistics on `generate_null()` moved by up to 0.006 between two runs on identical data (its random SVD start); it explains the ablation's 0.006 there.
- Correction: section 4.1 of the comparison report now compares rates in the fit's units (item 0b closed), and the brainstorm's "Corrections this implies elsewhere" are marked made.
- Finding: `.gitignore`'s `*cache*` had been ignoring `overdispersion_brainstorm/01_cache_fits_claude.R`, a script `run_all_claude.sh` needs; an exception line was added.
- Open: the comment in `.apply_nuisance_cap()` says the floor can only move a failed gene, but T-CAP-03c lifts an estimated gene with a tiny MLE; not edited, since `R/` was out of scope.
- Verification: `R CMD check --as-cran` on the HEAD tarball gave `Status: 2 NOTEs` (`New submission`, HTML Tidy), tests `[ FAIL 0 | WARN 0 | SKIP 20 | PASS 1512 ]`; suite under `NOT_CRAN=true` 2014 pass / 0 fail / 0 skip in 64 s; both reports knit. Nothing is committed.

### [2026-09-29] (Session 19 — misleading comments corrected in `R/` and `src/`)
- Kevin asked for every potentially misleading code comment to be fixed; three parallel read-only audits covered all of `R/` and `src/`, and each finding was checked against the code before editing.
- Fixed about 35 comments and roxygen blocks in 16 R files and 6 C++ files. The ones a reader was most likely to act on:
  - `initialize_esvd`'s description claimed two GLMs and a deviance-test p-value;
  - `estimate_nuisance`'s `bool_covariates_as_library` text was copied from `bool_adjust_covariates`;
  - `alpha_max` was called a cap on the numerator (it caps the fitted mean);
  - `bool_stabilize_underdispersion` spoke of "under-dispersion";
  - `gamma_rate.cpp` described a search that no longer runs;
  - the Bernoulli start mapping was written reversed;
  - `curved_gaussian`'s nuisance was called the CV (it is mean/sd);
  - `library_multipler`'s variance description matched no family;
  - `opt_esvd`'s `tol` was called a zero threshold.
- Resolved: the floor comment in `.apply_nuisance_cap()`. `estimate_nuisance()` floors before storing `nuisance_mle_vec`, so from both callers the helper's floor changes nothing. A gene lifted by the floor keeps the status `estimated`, and a failed gene's `nuisance_mle_vec` is the floor, not an estimate.
- Two message strings were also changed: `multtest()`'s error no longer blames `qnorm(pt())` saturation, and `opt_esvd`'s line-search warning says "last accepted step". No test matches either.
- Pointers to `CRAN_READINESS.md` in code now name `additional_context/`; the file is not shipped with the package.
- Open: two code inconsistencies the audit found, left unfixed because the request was about comments: the `bool_covariates_as_library` defaults differ between `estimate_nuisance.eSVD` (FALSE) and the posterior (TRUE); and `compute_test_per_gene` lacks `compute_posterior.default`'s refusal of `bool_adjust_covariates` with `bool_covariates_as_library`.
- Verification:
  - `git diff -U0` of `R/` and `src/` shows no changed line outside comments, roxygen and the two strings;
  - `devtools::document()` rewrote `RoxygenNote` again, and that was reverted;
  - the suite under `NOT_CRAN=true` gives 2014 pass / 0 fail / 0 skip;
  - `R CMD check --as-cran` gives `Status: 2 NOTEs` (`New submission`, HTML Tidy).
- Nothing is committed.

### [2026-09-29] (Session 20 — the library default aligned; the posterior paths refuse the same settings)
- At Kevin's direction, `estimate_nuisance.eSVD()` now defaults to `bool_covariates_as_library = TRUE`, the default of `compute_posterior.eSVD()` and `compute_test_per_gene()` and what `eSVD()` passes. No existing test changed with it; the stage-by-stage rates change for direct callers, and `NEWS.md` says so.
- Decision (Claude, within Kevin's "refuse nonsensical inputs"): one internal `.check_posterior_args()` in `R/posterior.R`, called first by `compute_posterior.default` and by `compute_test_per_gene`, so both refuse the same settings with identical messages. It replaces the matrix path's bare `!a | !b` `stopifnot`.
- Resolved: `nuisance_lower_quantile = NULL` now means no floor on the matrix path too. It had computed `quantile(x, probs = NULL)`, which is `numeric(0)`, and died later with "missing value where TRUE/FALSE needed".
- Finding (`/code-review`): after the per-gene refusal, `eSVD(bool_adjust_covariates = TRUE)` with the default library flag would fail only after the whole fit. `eSVD()` now runs the same check on entry (T-PGENE-05 asserts that the intermediate file is never written).
- Open, accepted: an object built before this change whose `param` records both flags TRUE on the diet path is now refused by `recompute_pvalue()`; that combination was never meaningful.
- Tests written first and seen failing: T-PGENE-02 to -05 and T-NUIS-09. T-NUIS-08 was taken (`bool_use_log`) and T-POST-14 is reserved in `CRAN_READINESS.md` §0.4; the NULL-quantile case on the matrix path is covered inside T-PGENE-04 rather than a separate T-POST test.
- Verification: suite under `NOT_CRAN=true` 2110 pass / 0 fail / 0 skip across 298 blocks; `R CMD check --as-cran` `Status: 2 NOTEs`, tests `[ FAIL 0 | SKIP 20 | PASS 1608 ]`; `RoxygenNote` rewrite reverted again. Nothing is committed.

### [2026-09-30] (Session 21 — the `_claude` suffix dropped from `R/` and `tests/testthat/`)
- Kevin asked for the suffix to be removed from the package folders; 5 files in `R/` and 22 in `tests/testthat/` had it, none elsewhere in the package.
- Decision (approved in plan mode): six drafted test files had the name of an older test file (`compute_test_per_gene`, `compute_test_statistic`, `gamma_rate`, `initialization`, `nuisance`, `posterior`), so each was appended to the older file under a divider rather than given a new name. The only line removed was `context("Test compute_test_per_gene (claude)")`.
- Scope: the 26 `*_claude.*` files under `additional_context/` keep the suffix, since Kevin named the package folders only. Open: whether they lose it too, and whether the `_claude` convention in the master `CLAUDE.md` still stands for new files.
- `TEST_RUN_REPORT.md` is a dated record, so its filenames were left as they were and a note was added to its banner instead.
- Found while renaming: a comment in the posterior tests placed `.library_column_oracle()` in the cap test file; it lives in `helper-fixtures.R`. Corrected.
- Found: the tree was clean at the start of the session, so the changes of sessions 18 to 20 are committed (`b0dd8f2`); `CLAUDE_kevin.md` had said otherwise.
- Verification: suite under `NOT_CRAN=true` 2110 pass / 0 fail / 0 warnings / 0 skip, 298 blocks in 28 files; `R CMD check --as-cran` `Status: 3 NOTEs`, tests `[ FAIL 0 | WARN 0 | SKIP 20 | PASS 1608 ]`. The third NOTE is "unable to verify current time" (no time server reached), not a package matter. `RoxygenNote` rewrite reverted again. Nothing is committed.
