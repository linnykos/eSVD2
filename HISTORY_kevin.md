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
