# eSVD2 → CRAN: proposed unit-test suite

**Status:** proposal for review, started 2026-08-27. This is the companion to
`additional_context/CRAN_READINESS.md`: that document says *what is wrong*, this
one says *what we would have to assert to know it is right*. Every test in §3 of
`CRAN_READINESS.md` has been absorbed here and expanded; that section now points
back at this file.

**Nothing here has been written yet.** The intent is that Kevin reads the list
first, strikes the tests that are not worth their cost, flags the ones whose
expected answer he disagrees with, and only then do we write code and regenerate
fixtures. Several tests below cannot be written until an open question is
answered — those carry a **⚠ BLOCKED** tag naming the decision.

**Review pass 1 — 2026-08-29.** Kevin answered 15 of the 16 open questions
in-line as `[[KZL: ...]]`. Those answers are now folded into the affected tests
and marked **✅ RESOLVED (Q-…)**; the `[[KZL: ...]]` text is kept verbatim
alongside so the decision and its author stay visible. That pass also added two
pieces of scope that were not in this document before:

- **§2.15** — a test group that settles whether `Rmpfr` is needed, rather than
  asserting the equivalence and moving on (the old T-PVAL-07).
- **§2.16** — a new `gene_status` output: all-zero genes are removed before
  fitting and reinserted at the end. This is a *feature*, so §2.16 pins its
  specification first and the tests second.

**Review pass 2 — 2026-08-29.** Kevin answered the remaining eight, and pointed
at an existing wrapper — `esvd_helper.R` in the `WAS2CODE_REPO` — to be brought
into the package. That wrapper answers Q-TSTAT-2 and adds three cohort-level
filters that were not in this plan at all, so it gets its own section, **§2.17**.
**Every question from passes 1 and 2 is now closed**; §9.2 holds only the five
new ones that reading `esvd_helper.R` raised.

**Review pass 3 — 2026-08-29.** All five Q-COH questions answered, including a
five-step ordered specification for `eSVD_helper()` (quoted in §2.17.2). The
drop-vs-stop conflict is settled — **drop, with a warning**.

**Review pass 4 — 2026-08-29. Every question in this document is now closed.**
The last two fixed the call graph: **drop donors first, then check the cohort**
(Q-COH-6), and **all filtering, gene-status labelling, temporary removal and
reinsertion live in `eSVD_helper()`, with `eSVD()` erroring on anything that gets
through** (Q-COH-7). The second supersedes Q-STATUS-1's placement and moves
`gene_status` out of `eSVD()` — §2.16 and §2.17 are rewritten to match. It also
turns Q-STATUS-2's "BH on analyzed genes only" from a convention into a
structural guarantee, which is a real improvement on what this plan originally
proposed. §9.5 lists three things to confirm while implementing; nothing blocks.

## How to read this

Every test has an ID (`T-<AREA>-nn`) so it can be referred to in review comments,
and three fields:

- **Assert** — the executable claim.
- **Oracle** — where the expected answer comes from. This is the field worth
  scrutinising: a test whose oracle is "whatever the code currently returns" is a
  *snapshot*, not a correctness test, and is marked **[snapshot]**. Tests whose
  oracle is an independent computation are marked **[oracle]**; tests whose
  oracle is a mathematical invariant are **[invariant]**.
- **Why** — the failure mode it catches. If this field cannot be filled in, the
  test should be cut.

Status tags:

- **[new]** — no coverage today.
- **[strengthen]** — a test exists but does not actually constrain the answer.
- **[regression]** — pins a specific defect from `CRAN_READINESS.md`; must fail
  before the fix and pass after.
- **⚠ BLOCKED** — needs a decision from Kevin first.
- **✅ RESOLVED (Q-…)** — decided in the 2026-08-29 review pass; the decision and
  its consequences are stated at the test, and summarised in §9.1.

Evidence tags on the *claims* (not the tests):

- **[verified]** — run in R in the session that wrote this file. Note these were
  run against the **installed** `eSVD2 1.0.1.2`, not the working tree's
  `1.0.1.07`; anything version-sensitive is called out.
- **[inspection]** — read from source.

**Recurring parameterization warning.** The paper's over-dispersion `γ_j` is the
Gamma **scale**; the code's `nuisance_vec` is the Gamma **rate** `β = 1/γ_j`.
Larger `nuisance_vec` means *less* over-dispersed. Several proposed tests below
assert a *direction*, and their expected sign depends on this. They are marked
where it matters.

---

## 0. Conventions the suite should adopt before any test is written

These are decisions about the harness, not tests. They are listed first because
every test below assumes them.

| # | Convention | Rationale |
|---|---|---|
| C-01 | `Config/testthat/edition: 3` in `DESCRIPTION`; drop `context()` from all 12 files | 3e is required for `expect_no_error`/`expect_snapshot`; `context()` is deprecated |
| C-02 | Fixtures loaded via `testthat::test_path("fixtures", "...")`, never `load("../assets/...")` | the current form depends on the working directory and breaks outside `test_check()` |
| C-03 | One shared `helper-fixtures.R` exposing `small_esvd_obj()`, `tiny_counts()`, `feasible_point(family)` | ~40 of the tests below need the same 3 objects; building them inline triples check time |
| C-04 | Every stochastic test sets a seed **inside** the `test_that` block, and states the seed's role in a comment | a seed chosen to make a flaky test pass is a bug in disguise; say so out loud |
| C-05 | Numerical comparisons use `expect_equal(..., tolerance = ...)`, never `sum(abs(a-b)) <= tol` | `abs(sum(differences))` lets per-element errors of opposite sign cancel (this is exactly §1.5's bug); `expect_equal` is element-wise |
| C-06 | `MASS::mvrnorm` removed from the suite in favour of `matrix(rnorm(...), ...) %*% chol(Sigma)` | removes the undeclared-dependency WARNING (§2.9) rather than papering over it with a `Suggests` entry |
| C-07 | C++ tests that need `numDeriv` are guarded with `testthat::skip_if_not_installed("numDeriv")` | `numDeriv` is in `Suggests`; CRAN runs with `_R_CHECK_FORCE_SUGGESTS_=false` on some machines |
| C-08 | Total suite runtime budget: **≤ 90 s** on one core | CRAN's whole-check budget is 10 min and the vignette needs most of it |
| C-09 | `covr::package_coverage()` in CI with a floor, set to whatever the suite achieves on the day it lands, then ratcheted | prevents silent coverage loss on later edits |

**✅ RESOLVED (C-10 / Q-VAL-1).** [[KZL: Let's actually use regexp()'s for
checking all the error messages]] — **no snapshot tier at all.** Every error and
warning assertion in §5 and §6 uses
`expect_error(..., regexp = "<short stable fragment>")` /
`expect_warning(..., regexp = ...)`. Add this as convention **C-11**:

| # | Convention | Rationale |
|---|---|---|
| C-11 | No `tests/testthat/_snaps/`. Message assertions use `expect_error`/`expect_warning`/`expect_message` with a `regexp` fragment chosen to be the *stable* part of the message (the offending argument name, the failing value), never the full sentence | Kevin's call. Keeps message edits a one-step change; the cost is that a message can be reworded around the fragment without failing, which is acceptable because the fragment is the part carrying the information |

**One test is a casualty of this.** T-CPP-LOAD-06 proposed a `[snapshot]` of
`test_data_loader()`'s printed debug output — there is no non-snapshot way to
assert that usefully. Recommendation: **delete `test_data_loader()`** (it is an
exported, undocumented debug printer, already flagged in `CRAN_READINESS.md`
§2.3) rather than write a weaker test for it. T-CPP-LOAD-06 is struck.

---

## 1. Fixtures — what the tests need, and what has to be regenerated

This section exists because you asked about the synthetic data. **Most of the
work in this plan is fixture work, not assertion work.** The current
`tests/assets/synthetic_data.RData` cannot support the suite below for three
independent reasons: it is 10.1 MB against a 5 MB tarball limit (§2.8), it stores
a *fully fitted* `eSVD_obj` so tests that should exercise the pipeline instead
read its output back, and it has only one shape (2000 × 150, Poisson, 2 factors,
20 individuals) so no degenerate case is reachable.

### 1.1 Proposed fixture set

| Fixture | Shape | Built by | Used by |
|---|---|---|---|
| `F-TINY` | 120 cells × 20 genes, 6 individuals (3 case / 3 control), 1 factor + 1 numeric covariate, Poisson–Gamma per `data_generation.R`'s model | script, saved with `compress = "xz"` | the bulk of §2 — anything that needs a real but fast end-to-end object |
| `F-SMALL` | 400 cells × 40 genes, 8 individuals (4/4) | same script, different seed | §4 invariants that need enough individuals for a stable Welch df |
| `F-DEGEN` | a family of *hand-built* degenerate inputs, not simulated: an all-zero gene, a constant gene, a gene with one non-zero count, an individual with one cell, a group with one individual, a collinear covariate pair, a factor level with spaces and parentheses | inline in the test files, not stored | all of §5, plus T-CPP-GAM-03..06 |
| `F-NULL` | `generate_null(cell_per_person = 25, num_genes = 120, num_individuals = 8)` under a fixed seed | called at test time, not stored | T-PROP-06 (null calibration), T-ESVD-01 |
| `F-IDENT` | keep `identification1.rda` as-is | already small | existing `.identification` underflow test |

`F-TINY` and `F-SMALL` should store **only** `dat`, `covariates`, `metadata`,
`nuisance_vec`, `true_cc_status` — the raw inputs. The eleven-object
`synthetic_data.RData` (fitted `eSVD_obj`, `library_mat`, `gamma_mat`,
`nat_mat_nolib`, `session_info`) goes away. Tests that need a fitted object build
one in a helper from `F-TINY`, which is both smaller and a stronger test.

**✅ RESOLVED (Q-FIX-1)** [[KZL: Yes, this is a good test, and I don't think it
will be degenerate.]] — build the fitted object at test time, and `F-TINY` stays
at 120 × 20 with `k = 2`.

Agreed on the reasoning: after the covariate columns are orthogonalized out the
factorization fits 120·2 + 20·2 = 280 parameters against 2400 counts, and
`.reparameterize` needs rank 2 out of 20 gene rows, which is comfortable. The
one rider that stands: **the recovery assertions (T-NUIS-05 / T-PROP-07) still go
on `F-SMALL`**, not because 120 × 20 is degenerate but because ±30% recovery of a
Gamma rate needs cells, and 120 is thin for that. Everything structural runs on
`F-TINY`.

**Two new constraints on the fixtures**, from the 2026-08-29 decisions:

1. `F-TINY` must contain **at least two all-zero genes** (§2.16), which leaves 18
   analyzed genes out of 20. If Kevin raises the gene count, keep the all-zero
   genes as a fixed *count*, not a fraction.
2. Every individual in `F-TINY`, `F-SMALL` and `F-NULL` must have **≥ 3 cells**
   (§2.12, T-ESVD-09), or `eSVD()` will now refuse to run on them. At 6
   individuals × 120 cells this is automatic; it becomes a real constraint on
   any hand-built fixture.

### 1.2 Fixture provenance

`tests/assets/data_generation.R` stays and is updated in the same commit as the
fixtures — it is the only record of what the numbers mean. Drop
`devtools::session_info()` from it (that is the sole reason `devtools` is in
`Suggests`, §2.6) and record `R.version.string` and package versions with
`utils::sessionInfo()` in a comment instead.

### 1.3 Known-truth fixtures for the C++ layer

`F-DERIV`: for each of the 7 families, a **feasible** point `(XC, YZ, s, gamma)`
plus a 5 × 4 count matrix, stored as a plain list. The families have different
domains — `curved_gaussian` needs `θ > 0`, `exponential` and `neg_binom` need
`θ < 0`, the rest are unconstrained — so a single random point does *not* work
for all seven and the current suite's silence on six of them is partly explained
by that. This fixture is a prerequisite for all of §3.2 and should be built by a
helper `feasible_point(family)` rather than stored.

**Also required by F-DERIV: a non-`NA` `gamma`.** `opt_esvd.default`'s default is
`nuisance_vec = rep(NA, ncol(input_obj))`, and four of the seven families consume
`gamma` — so with the default the objective is `NA`. [verified] against the
installed build:

```
gaussian         objfn_all_r with gamma=NA -> NA
curved_gaussian                            -> NA
neg_binom                                  -> NA
neg_binom2                                 -> NA
poisson                                    -> 3.663431   (gamma unused)
bernoulli                                  -> -0.01228382 (gamma unused)
```

and `opt_esvd.default(family = "gaussian")` with the default `nuisance_vec` fails
with `Error: missing value where TRUE/FALSE needed` [verified]. See T-OPT-05.

---

## 2. R-level tests, in pipeline order

### 2.1 `format_covariates()` — `test_format_covariates.R` **[new file]**

Zero coverage today, and every downstream matrix depends on its column names and
ordering.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-FMT-01 | `colnames()[1] == "Intercept"`, `colnames()[2] == "Log_UMI"`, factor indicators last and in `levels()` order | [invariant] the contract every downstream `which(colnames(...) == ...)` relies on | column *position* is load-bearing in `compute_posterior.default` (`covariates[,-library_idx]`) |
| T-FMT-02 | a 3-level factor yields exactly 2 indicator columns, and the **first** level is the one dropped | [oracle] `stats::model.matrix(~ g)` gives the same 2 columns | **the roxygen says "dropping the last level"; the code drops the first** — [verified]: `factor(c("a","a","b","b","c","c"))` yields columns `g_b`, `g_c`. One of the two must change |
| T-FMT-03 | `variables_enumerate_all = "g"` keeps all 3 levels | [invariant] documented behaviour | the branch is untested and silently changes the rank of `C` |
| T-FMT-04 | `Log_UMI == log(Matrix::rowSums(dat))` exactly, for both `matrix` and `dgCMatrix` input | [oracle] direct recomputation | the paper's library-size definition; a `log1p` vs `log` slip here shifts every fitted `z` |
| T-FMT-05 | a variable named in `rescale_numeric_variables` has `sd == 1` and its mean is **unchanged** (not centered) with `bool_center = FALSE` | [invariant] the paper is explicit about not centering | centering here would break the interpretation of the intercept |
| T-FMT-06 | a numeric variable **not** named in `rescale_numeric_variables` passes through untouched | [oracle] identity | the roxygen says "rescales all the numerical variables"; the code rescales only the named ones. Second doc/code mismatch in this function |
| T-FMT-07 | `nrow(dat) != nrow(covariate_df)` gives an informative error | [invariant] | currently a bare `stopifnot` |
| T-FMT-08 | a single-level factor gives an informative error naming the variable | [invariant] | there is a `stopifnot(length(levels(vec)) > 1)` with no message; `eSVD()` works around it by pre-filtering, so the error is user-facing only for direct callers |
| T-FMT-09 | a factor level containing a space and parentheses (`"ASD (severe)"`) survives into `colnames()` **verbatim**, un-mangled | [invariant] | this is the input that breaks `reparameterization_esvd_covariates()` downstream (§1.3 of the readiness doc). Asserting it here localises the bug to the right function |

### 2.2 `initialize_esvd()` — `test_initialization.R` **[strengthen]**

Two tests exist; both assert `res[,"Intercept"] == 0`, i.e. only the
`bool_intercept = FALSE` branch.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-INIT-01 | `bool_intercept = TRUE`: `z_mat[,"Intercept"]` is **not** all zero and equals `glmnet`'s `a0` at the smallest lambda | [oracle] refit one gene with `glmnet::glmnet` directly | the branch is completely untested; a `c(0, ...)` vs `c(a0, ...)` slip is invisible |
| T-INIT-02 | `offset_variables = NULL`: `z_mat[,"Log_UMI"]` is *estimated*, not pinned to 1 | [oracle] contrast against the `offset_variables = "Log_UMI"` run | untested branch; `offset_vec <- NULL` changes the glmnet call shape |
| T-INIT-03 | `offset_variables = c("Log_UMI", "Age")`: the offset is `rowSums` of both columns and both coefficients come back as 1 | [oracle] direct recomputation of `Matrix::rowSums(covariates[,offset_variables])` | the multi-offset path is reachable and untested; note it forces *both* coefficients to 1, which is a modelling assumption worth stating |
| T-INIT-04 | the corner case `ncol(covariates_tmp) <= 1` takes the `stats::glm` branch and agrees with the `glmnet` branch to ~1e-3 when both are applicable | [oracle] the two branches are each other's oracle | the `glm` fallback has different penalisation (none) — the agreement tolerance is the honest measure of how different |
| T-INIT-05 | `dat` with `NA` entries: `is.matrix(dat)` path zeroes them; the `dgCMatrix` path does **not** | [inspection] `dat[is.na(dat)] <- 0` is guarded by `is.matrix(dat)` | **✅ RESOLVED (Q-INIT-1)** [[KZL: Yes, the NA should be 0'd as well]]. The test asserts *symmetry*: both the `matrix` and the `dgCMatrix` path zero `NA`s, so `initialize_esvd()` never hands an `NA` downstream. This is a **fix plus a test**, not a test alone — move `dat[is.na(dat)] <- 0` out from under the `is.matrix(dat)` guard and give the sparse case its own `dgCMatrix`-safe form (`dat@x[is.na(dat@x)] <- 0`, then `Matrix::drop0`). Note the consequence for §3.1: the C++ `Flag::na` machinery is then unreachable *through `initialize_esvd`*, but stays reachable for direct `opt_esvd`/`objfn_all_r` callers, so T-CPP-LOAD-02/04 remain worth writing |
| T-INIT-06 | `k > ncol(dat)`, `k = 0`, non-integer `k` each give an informative error | [invariant] | `stopifnot` today |
| T-INIT-07 | `covariates` without an `"Intercept"` column errors | [invariant] | there is a `stopifnot`; test it |
| T-INIT-08 | `metadata_individual` not a factor errors | [invariant] | there is a `stopifnot`; test it |
| T-INIT-09 | `length(metadata_individual) != nrow(dat)` errors **at `initialize_esvd`**, not three functions later in `compute_test_statistic` | [invariant] | currently unchecked at entry — the failure surfaces much later with an unrelated message |
| T-INIT-10 | the returned object has exactly the documented names and `class == "eSVD"`, and `latest_Fit == "fit_Init"` | [snapshot] | exists today; keep |

### 2.3 SVD helpers — `test_svd_safe.R` **[new file]**

`.svd_safe` / `.svd_in_sequence` / `.compute_matrix_mean` / `.compute_matrix_sd`
have no direct tests; only the first link of the `irlba → RSpectra → base::svd`
chain is ever exercised.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-SVD-01 | on a well-conditioned dense matrix, `.svd_safe(K = k)` reproduces `base::svd`'s top-`k` singular values to 1e-8 and `method == "irlba"` | [oracle] `base::svd` | establishes the happy path actually is the happy path |
| T-SVD-02 | the returned `u`/`v` carry `rownames(mat)` / `colnames(mat)` | [invariant] | downstream `rownames(x_mat) <- rownames(dat)` relies on this not already being wrong |
| T-SVD-03 | feed a matrix that makes `irlba` warn (e.g. rank < K, or K within 2 of `min(dim)`), assert `method == "RSpectra"` | [invariant] the fallback contract | the fallback chain has never run in a test; a broken second link would only show up on a user's data |
| T-SVD-04 | force both to fail (mock or a pathological matrix), assert `method == "base"` and the result is still a valid SVD | [invariant] | same |
| T-SVD-05 | `.compute_matrix_sd(mat, sd_vec = TRUE)` on a `dgCMatrix` equals `matrixStats::colSds(as.matrix(mat))` | [oracle] dense recomputation | **✅ RESOLVED (Q-SVD-1) — target CRAN, drop the dependency** [[KZL: Let's simply include the relevant sparseMatrixStats function directly into our package so we don't need to depend on it. Just remember to cite the function]]. The test is unchanged and is what makes either implementation safe. **⚠ NEW (Q-SVD-3): vendor, or reimplement?** — see the note below the table; I recommend reimplementing rather than vendoring, and the branch turns out to be dead code, which makes it a small call either way |
| T-SVD-06 | `.compute_matrix_mean(mat, mean_vec = 0.5)` — a length-1 *numeric* — is silently treated as `TRUE` and returns `colMeans` | [inspection] `if(mean_vec)` coerces `0.5` to `TRUE` | **✅ RESOLVED (Q-SVD-2) — reject scalars** [[KZL: Please make only `TRUE`/`FALSE`/`NULL`/full-vector accepted]]. The test flips from "assert the silent coercion" to "assert the error": `.compute_matrix_mean(mat, mean_vec = 0.5)` and `.compute_matrix_sd(mat, sd_vec = 0.5)` each error with a message naming the argument and the accepted forms. Note the guard has to test `is.logical()`, **not** `length(x) == 1` — `mean_vec = 0.5` and `mean_vec = TRUE` are both length 1, which is precisely why the current code cannot tell them apart. Also covered by T-VAL-29/30 |
| T-SVD-07 | `check_stability = TRUE, K = 3` does not run the stability check (guard is `K > 5`), and `K = 10` does | [inspection] | the guard uses `&` on scalars (§1.8); the test pins the intended semantics before `&&` is substituted |

**⚠ NEW (Q-SVD-3) — now resolved, see below: vendoring `sparseMatrixStats::colSds` is more work and more
CRAN paperwork than rewriting it, and the branch is dead.** Three facts worth
having before choosing:

1. **There is exactly one call site**, `R/reparameterization.R:334`, inside
   `.compute_matrix_sd()`. It fires only when the caller passes `sd_vec = TRUE`
   *and* `mat` is a `dgCMatrix`/`dgTMatrix`. The only in-package caller,
   `.initialize_residuals()`, passes `sd_vec = NULL` and a **dense**
   `residual_mat`. **The sparse branch never executes in the package as
   shipped** [inspection]. Whatever replaces it, no eSVD2 result changes.
2. **Vendoring is not free on CRAN.** `sparseMatrixStats` is MIT-licensed, so
   copying is permitted, but CRAN then requires its copyright holder added to
   `Authors@R` with role `cph` and the origin recorded in `inst/COPYRIGHTS`
   (`LICENSE` alone is not enough). And `colSds` there is *C++*, so vendoring
   means a new `src/` file, not a new `R/` function.
3. **The `Matrix`-only rewrite is four lines** and needs no attribution beyond a
   roxygen `@references`. The one thing it gives up is numerical stability:
   `sqrt((colSums(mat^2) − n·colMeans(mat)^2)/(n−1))` is the naive two-pass form
   and loses precision on columns with large means, which is exactly what the
   C++ version avoids.

**Recommended decision rule, so this does not need another round trip:** write
T-SVD-05 first, with `matrixStats::colSds(as.matrix(mat))` as the oracle and a
**relative tolerance of 1e-10**, on a fixture that deliberately includes a column
with a large mean and small variance. Then put the four-line `Matrix` rewrite
under it. If it passes, ship the rewrite and cite `sparseMatrixStats` in
`@references`. If it fails, vendor the C++ and do the `cph` paperwork.

**✅ RESOLVED (Q-SVD-3)** [[KZL: Yes, I think your plan is great. Let's do this
before vendor the C++]] — the decision rule stands: **try the four-line `Matrix`
rewrite first, vendor only if T-SVD-05 fails at 1e-10.** Given that the sparse
branch never executes in the package as shipped, my expectation is that it
passes and no C++ is vendored, so no third-party `cph` entry lands in
`Authors@R`. If it does fail, that is worth knowing before the paperwork.

A fourth option worth naming: `rescale`/`sd_vec` in `.svd_safe()` is never used
by any caller either. **Deleting the parameter** removes the dependency, the
dead branch, T-SVD-05 and T-SVD-06 all at once.

### 2.4 `opt_esvd()` — `test_optimization.R` **[strengthen]**

One test, Poisson only, `covariates` always present.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-OPT-01 | `opt_esvd.default(covariates = NULL)`: runs, returns `z_mat = NULL`, and `loss` is monotone | [invariant] | the no-covariate path works only because `.opt_esvd_format_matrices()` short-circuits on `is.null(covariates)`; nothing tests it |
| T-OPT-02 | `all(diff(loss) < 0)` for **every** family, not just `fit_First` under Poisson | [invariant] the optimizer's defining property | monotone descent is the one thing an alternating Newton scheme must guarantee; a sign error in any `family_*.cpp` breaks it |
| T-OPT-03 | run twice on identical input → `identical()` output | [invariant] the paper's determinism claim: *"different practitioners using our method would necessarily obtain the same resulting fit"* | guards against any future thread/atomic non-determinism in the C++ |
| T-OPT-04 | `offset_variables` columns of `z_mat` are **bit-identical** before and after | [invariant] | exists for `opt_yz` at the helper level; assert it at the `opt_esvd` level too, where `fixed_cols` is computed by name matching and could silently match nothing |
| T-OPT-05 | `opt_esvd.default(family = "gaussian")` with the **default** `nuisance_vec = rep(NA, p)` gives an informative error naming the family and the missing nuisance | [verified] currently `Error: missing value where TRUE/FALSE needed` | **[regression]** four of seven families are unusable at their documented defaults and the error names neither the family nor the parameter. This is arguably a §1-class defect that `CRAN_READINESS.md` does not list |
| T-OPT-06 | for each of the 6 non-Poisson families, with a *valid* `nuisance_vec`, `opt_esvd.default` runs to completion and reduces the loss | [invariant] | six documented, reachable, entirely untested code paths |
| T-OPT-07 | `max_iter = 1` returns after one iteration with a length-1 `loss` (no convergence test is attempted) | [invariant] | the `if(i >= 2)` guard; an off-by-one here would index `losses[0]` |
| T-OPT-08 | `.opt_esvd_setup_z_mat(covariates = NULL, ...)` returns `NULL` and `(covariates, z_init = NULL)` returns a `p × ncol(covariates)` zero matrix | [invariant] | the function returns the value of an `if`/`else` whose branches are assignments — correct by accident (§1.8) |
| T-OPT-09 | `verbose = 2` runs without error | [invariant] | **[regression]** §1.2: `print("Residual of loss: ", resid)` throws `invalid printing digits 0`. The most verbose mode of the core optimizer has never been executed |
| T-OPT-10 | `sum(!is.na(input_obj)) == 0` errors with a message about an all-`NA` matrix | [invariant] | there is a `stopifnot`; also exercises the C++ `"all elements in Ai are NA"` path from R |

### 2.5 `reparameterization_esvd_covariates()` — `test_reparameterization.R` **[strengthen]**

`.identification()` and `.reparameterize()` are well covered. The exported
covariate-orthogonalization step — the paper's Step 1 — is not covered at all.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-REP-01 | after the call, `Ŷ X̂ᵀ + Ẑ Cᵀ` is unchanged to 1e-8 | [invariant] the paper asserts this twice: *"the predictive power of our factorization did not change"* | **the single most valuable assertion about this function.** Nothing checks it today |
| T-REP-02 | after the call, `crossprod(x_mat, covariates) ≈ 0` | [invariant] the paper's Step 1 identifiability condition | the function's entire purpose |
| T-REP-03 | after the call, `crossprod(x_mat)/n` and `crossprod(y_mat)/p` are diagonal **and equal** | [invariant] Step 2 | covered for `.reparameterize` in isolation; not for the composed function |
| T-REP-04 | `omitted_variables` columns of `z_mat` are unchanged, and `x_mat` is **not** orthogonalized against them | [invariant] | `eSVD()` passes `c("Log_UMI", batch_vars)` here; if the argument silently did nothing, every published fit would differ |
| T-REP-05 | a covariate column whose name is non-syntactic (`"Diagnosis_ASD (severe)"`) either works or errors informatively — never `NA`-contaminates `z_mat` | [invariant] | **[regression]** §1.3: `as.data.frame()` applies `make.names()`, so `names(coef_vec)` carries the mangled name and `z_mat[j, names(coef_vec)]` is a subscript error. Reachable from ordinary factor levels; T-FMT-09 is its upstream half |
| T-REP-06 | a deliberately collinear covariate pair gives an informative error naming the collinear columns | [invariant] | **[regression]** §1.3: `lm` drops aliased terms, `coef()` returns `NA`, the `NA` propagates into `z_mat` → `μ̂` → every posterior, with no warning. The paper *explicitly* warns that one-hot individual vectors make `C` collinear, so this is a documented user footgun |
| T-REP-07 | `.identification()` with a `sym_prod` having a negative eigenvalue errors or clamps, rather than returning `NaN` via `sqrt()` | [invariant] | §1.8: the function already `warning()`s about rank deficiency and then proceeds anyway |
| T-REP-08 | `.reparameterize` on rank-deficient `x_mat` (duplicated column) **warns and proceeds**: `expect_warning(regexp = "rank")`, and the returned object is still a valid (unreparameterized) fit | **✅ RESOLVED (Q-REP-1)** [[KZL: Let's warn (and print via an appropriate verbose level) and proceed]] | `opt_esvd.default` wraps it in a `tryCatch` that falls back to the unreparameterized matrices *silently* — so a rank-deficient fit is currently indistinguishable from a good one |
| T-REP-09 | the warning **survives `opt_esvd`**: a rank-deficient fit driven through `opt_esvd(verbose = 0)` still emits it | [invariant] | **the half that makes T-REP-08 worth anything.** A `warning()` raised inside `.reparameterize` is swallowed today by `opt_esvd.default`'s `tryCatch(..., warning = ...)`. Warning at the inner level and catching it at the outer level leaves the user exactly as uninformed as before |
| T-REP-10 | at `verbose >= 1` the diagnostic **names the offending columns**; at `verbose = 0` only the bare `warning()` fires | [invariant] Kevin's "print via an appropriate verbose level" | two channels, deliberately: the `warning()` is unconditional so it shows up in `R CMD check` and in non-interactive runs, and the verbose message carries the detail that would be noise otherwise |

### 2.6 `estimate_nuisance()` — `test_nuisance.R` **[strengthen]**

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-NUIS-01 | `bool_use_log = TRUE` and `FALSE` agree to ~1e-3 on a well-conditioned fixture | [oracle] the two are each other's oracle — `gamma_rate` and `exp(log_gamma_rate)` estimate the same quantity by different routes | a free equivalence test, same idea as §1.5. Currently the `bool_use_log = TRUE` path has zero coverage |
| T-NUIS-02 | `.nuisance_in_sequence()` failures are **countable**: a helper returns how many genes fell through to `0` | [invariant] | today the function returns `0` and warns only when `verbose > 0`; the `0` is then clamped to `min_val`, so a silent total failure across all genes looks identical to a successful fit |
| T-NUIS-03 | the returned vector carries `colnames(input_obj)` as names | [invariant] | `compute_posterior` sweeps by position and `report_results` names by gene; a names/position mismatch here mislabels every result |
| T-NUIS-04 | all values `>= min_val` and finite, including on the `mean_mat` containing `Inf` fixture | [invariant] | exists; keep |
| T-NUIS-05 | on `generate_data()` output with a known `nuisance_param_vec`, the estimate recovers the truth within a stated tolerance | [oracle] the simulation truth | **✅ RESOLVED (Q-NUIS-1) — ±30% relative on the rate** [[KZL: Yes, let's use 30%]]. Written as `expect_equal(est, truth, tolerance = 0.3)` **gene by gene**, not on the mean of the vector — the mean would let a few wildly wrong genes hide behind the rest. Two riders: (i) run it on `F-SMALL`, not `F-TINY`, per Q-FIX-1; (ii) 30% at `n = 400` is a *loose* bar, so if the suite ever passes it only barely, that is a finding about the estimator, not about the tolerance — record the observed spread in a comment when the test lands |
| T-NUIS-06 | `estimate_nuisance.eSVD(bool_covariates_as_library = TRUE)` and `FALSE` give different answers, and the `TRUE` branch's `library_idx` matches `compute_posterior`'s | [invariant] | the two functions build `library_size_variables` with *nearly* identical but not shared code — `nuisance.R` uses `c(...)` where `posterior.R` uses `unique(c(...))`. If a variable were ever listed twice the index sets would diverge |
| T-NUIS-07 | dimension mismatch between `input_obj`, `mean_mat`, `library_mat` errors informatively | [invariant] | `stopifnot` today |

### 2.7 `compute_posterior()` — `test_posterior.R` **[strengthen]**

Currently three assertions, all about `dim()`.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-POST-01 | `posterior_var_mat == posterior_mean_mat / SplusBeta` exactly | [invariant] Eq. 14 | the cheapest possible check that the two matrices are consistent |
| T-POST-02 | all entries of both matrices are finite and strictly positive | [invariant] | a single `Inf` here propagates to `teststat_vec` and then to `locfdr` (§1.1's mechanism) |
| T-POST-03 | with `pseudocount = 0` and one all-zero gene, the posterior mean equals `Alpha/SplusBeta` and is **not** zero | [oracle] direct recomputation | the prior is what keeps an unobserved gene finite; if `Alpha` were dropped the gene would silently read as exactly 0 |
| T-POST-04 | dense `matrix` and `dgCMatrix` input give identical results | [invariant] | `compute_posterior.default` densifies internally; the test pins that the densification is faithful |
| T-POST-05 | `bool_return_components = TRUE` returns `numerator_mat`/`denominator_mat` and `mean == num/denom`, `var == num/denom^2` | [invariant] | untested branch, and the components are the documented debugging hook |
| T-POST-06 | `bool_adjust_covariates = TRUE` runs, and `bool_adjust_covariates = TRUE, bool_covariates_as_library = TRUE` errors | [invariant] the `stopifnot(!a \| !b)` | the experimental branch is documented, reachable and entirely untested; the mutual exclusion has no test either |
| T-POST-07 | `nuisance_vec` keeps `is.numeric() && !is.matrix()` and keeps its `names()` on **both** branches of the `bool_stabilize_underdispersion` guard | [invariant] | **[regression]** §1.4: `scale()` returns an `n × 1` matrix and moves `names()` to `rownames()` [verified]. Latent today, a bug on the next edit |
| T-POST-08 | the stabilization fires when `mean(log10(nuisance_vec)) > 0` — i.e. when the cohort is on average **under**-dispersed, `nuisance_vec` being the rate | **✅ RESOLVED (Q-POST-1) — the code is correct** [[KZL: The code is right, prose is wrong]] | §1.4. The change is therefore a **documentation fix, not a code fix**: rewrite the roxygen to say "when the global mean over-dispersion `γ = 1/nuisance_vec` is less than 1", and land it in the same commit as this test. Add the rate-vs-scale sentence to `?compute_posterior` and `?estimate_nuisance` while there — its absence is why this looked like a bug |
| T-POST-09 | `library_min` actually binds: with `library_min = 1e6` every entry of `SplusBeta` is `>= 1e6` | [invariant] | pins the clamp, which is the parameter whose default *disagrees* between the two pipelines (§1.5) |
| T-POST-10 | `alpha_max` binds: with `alpha_max = 1`, `Alpha <= nuisance_vec` element-wise | [invariant] | same reasoning |
| T-POST-11 | `nuisance_lower_quantile = 0.5` floors half the genes at the median | [oracle] direct `quantile()` recomputation | untested parameter |

### 2.8 `compute_test_statistic()` — `test_compute_test_statistic.R` **[strengthen]**

The helpers are well tested. The top-level statistic is checked only by a
weak power-ish comparison against `true_cc_status`.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-TSTAT-01 | on a fixture with **one cell per individual** and zero posterior variance, `teststat_vec` equals `stats::t.test(x, y, var.equal = FALSE)$statistic` gene by gene | [oracle] `stats::t.test` | in that degenerate case the mixture-Gaussian machinery reduces exactly to Welch's t. This is the only available *external* oracle for the headline statistic, and it is strong |
| T-TSTAT-02 | an individual appearing in **both** case and control errors | [invariant] | there is a `stopifnot`; test it |
| T-TSTAT-03 | `compute_test_statistic()` called **directly** with an individual having one cell still returns a finite statistic (no new guard at this level) | **✅ RESOLVED (Q-TSTAT-1) — guard at the wrapper, not here** [[KZL: There should be preprocessing check that stops the eSVD function early if there's individuals with 2 or less cells]] | the check goes in `eSVD()` as T-ESVD-09/10. **It must not go in `compute_test_statistic()`** — see the note below the table |
| T-TSTAT-04 | a group with exactly **one** individual gives a real error, not a `NaN` df | [invariant] | `n1 - 1 = 0` makes the Welch denominator `0/0`. `CRAN_READINESS.md` §3.5 already flags this: *"deserves a real error message, not a `NaN`"* |
| T-TSTAT-05 | `case_mean`/`control_mean`/`teststat_vec` all carry `colnames(posterior_mean_mat)` | [invariant] | `report_results` builds its `genes` column from `names(teststat_vec)` |
| T-TSTAT-06 | permuting the row order of `posterior_mean_mat` and `individual_vec` together leaves `teststat_vec` unchanged | [invariant] exchangeability | catches any accidental positional (rather than name-based) indexing in the averaging matrix |
| T-TSTAT-07 | `.construct_averaging_matrix` with an `idx_list` entry of length 0 errors, rather than producing an all-zero row that silently averages to 0 | [invariant] | reachable via an unused factor level; the resulting all-zero row biases the group mean with no signal |

**Where the ≤ 2-cell check must *not* live.** Kevin's answer to Q-TSTAT-1 is a
preprocessing check that stops early when any individual has 2 or fewer cells.
Putting that `stopifnot` inside `compute_test_statistic()` would **delete
T-TSTAT-01 and T-DF-01** — both are built on a one-cell-per-individual fixture,
because that is the degenerate case where the mixture-Gaussian machinery
collapses exactly to Welch's t and `stats::t.test` becomes an external oracle.
That oracle is the only independent check on the package's headline statistic,
and it is worth more than the redundancy of a second guard.

**Recommendation: the check lives in `eSVD()` only** (T-ESVD-09/10).
`compute_test_statistic()` keeps accepting whatever it is given, which is also
consistent with how the rest of the pipeline is layered — the `.default` methods
compute, the wrapper validates.

### 2.9 `.compute_df()` and `compute_pvalue()` — `test_compute_pvalue.R` **[new file]**

`compute_pvalue` produces the package's headline output and has **zero** tests.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-DF-01 | on the one-cell-per-individual fixture, `.compute_df()` equals `stats::t.test(...)$parameter` | [oracle] `stats::t.test` | same reduction as T-TSTAT-01; an external oracle for Welch–Satterthwaite |
| T-DF-02 | `.compute_df()` and `compute_test_statistic()` are computed from the **same** group variances | [invariant] | §1.6: `.compute_df()` re-derives `case_individuals`, the averaging matrix, and both group variances from scratch. Two copies of the same code drift; this test is what makes the refactor safe |
| T-DF-03 | `.compute_df()` errors informatively when `input_obj$dat` is absent | [invariant] | §1.6: the recomputation is why `compute_pvalue` needs `dat`, which conflicts with `bool_diet = TRUE` in `eSVD()`. The error should say so |
| T-PVAL-01 | with `teststat_vec` containing `t = 40, df = 18`, `gaussian_teststat` is **finite** | [verified] `qnorm(pt(40, 18))` is `Inf`; `qnorm(pt(-40, 18, log.p = TRUE), log.p = TRUE)` is `-8.915293` | **[regression]** §1.1. This is the highest-value single assertion in the whole plan: it is the defect that changes published-style results |
| T-PVAL-02 | with that same input, `pvalue_list$method == "locfdr"` | [invariant] | §1.1's real damage: one `Inf` makes `locfdr` error, the `tryCatch` swallows it, and the empirical null silently degrades to `.multtest_simple()` — the estimator whose own comment disclaims it. **Requires `method` to be stored in `pvalue_list`, which it is not today** |
| T-PVAL-03 | `log10pvalue` is monotone decreasing in `\|gaussian_teststat - null_mean\|` | [invariant] | a p-value that is not monotone in the statistic is definitionally broken; cheap and catches sign errors in the mirroring branch |
| T-PVAL-04 | `10^(-log10pvalue)` lies in `[0, 1]` for every gene | [invariant] | the mirroring construction `null_mean - (x - null_mean)` guarantees `2*pnorm(...) <= 1`; assert it |
| T-PVAL-05 | `fdr_vec == stats::p.adjust(pvalue_vec, "BH")` | [oracle] `p.adjust` | exists inside `test_report_results.R`; move it here where it belongs |
| T-PVAL-06 | `compute_pvalue`'s internally recomputed `log10pvalue_vec` equals `multtest()`'s `logpvalue_vec` | [invariant] | the same 6-line computation is written twice, once in `multtest.R` and once in `compute_pvalue.R` (and a third time in `compute_test_per_gene.R`). Either they agree or one is wrong |
| ~~T-PVAL-07~~ | **superseded by §2.15** (T-MPFR-01..08) | — | the original assertion — `Rmpfr::pnorm` and `stats::pnorm` are `identical()` on double input — shows the current calls are pointless but not that arbitrary precision is unnecessary *in principle*. Kevin asked for the stronger question to be settled, so §2.15 does it against a genuine 200-bit MPFR oracle |

### 2.10 `multtest()` — `test_multtest.R` **[new file]**

Zero coverage. This is the paper's Type-1-error control mechanism.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-MT-01 | `.multtest_locfdr()` on a clean N(0,1) sample of 1000 recovers `null_mean ≈ 0`, `null_sd ≈ 1` within a stated tolerance | [oracle] the simulation truth | establishes the primary estimator works before testing what happens when it doesn't |
| T-MT-02 | `.multtest_truncatedGauss()` on the same sample recovers the same, and `.multtest_simple()` too | [oracle] same | three estimators of one quantity — they are each other's oracle |
| T-MT-03 | the fallback chain is entered **in order**, and `method` reports which one ran: force `locfdr` to fail → `method == "truncated_mle"`; force both → `"simple"` | [invariant] | §1.1's silent-degradation mechanism. Note `.multtest_locfdr()` catches **warnings** as well as errors, so any `locfdr` warning — which is not rare — downgrades the estimator |
| T-MT-04 | a non-finite entry in `teststat_vec` errors at entry rather than propagating | [invariant] | §1.1 fix step 3: `stopifnot(all(is.finite(teststat_vec)))`. Today an `Inf` reaches `mean`/`sd` in `.multtest_truncatedGauss` and yields `NA` |
| T-MT-05 | `.multtest_truncatedGauss()` never evaluates `dnorm`/`pnorm` at `sigma0 <= 0` | [invariant] | **[regression]** §1.7: only `theta` is bounded; Nelder–Mead is free to step `sigma0` negative and `pnorm(0, 0, -1)` is `NaN` [verified]. Not yet reproduced in a real run — the test should therefore assert the *guard*, by calling the objective directly with `sigma0 = -1` and expecting `Inf` |
| T-MT-06 | `optim_res$convergence` is checked | [invariant] | currently ignored entirely |
| T-MT-07 | `.multtest_simple()` on input whose middle 90% is constant returns `null_sd == 0`, and the caller handles it | [invariant] | `pnorm(x, sd = 0)` is a step function; every p-value becomes 0 or 1 |
| T-MT-08 | `observed_quantile` is respected: `c(0.25, 0.75)` uses a strictly smaller subset than `c(0.05, 0.95)` | [invariant] | untested parameter |

### 2.11 `compute_test_per_gene()` — `test_compute_test_per_gene.R` **[strengthen]**

The existing equivalence test is `abs(sum(res1$teststat_vec - res2$teststat_vec)) <= 1e-3`,
which per-gene errors of opposite sign cancel out of. `CRAN_READINESS.md` calls
the strengthened version *"the single highest-value test in the package"* and I
agree: it pins one full implementation of the posterior/test/p-value stack
against another.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-PG-01 | `expect_equal(res1$teststat_vec, res2$teststat_vec, tolerance = 1e-8)` — element-wise | [oracle] the matrix pipeline | **[strengthen]** replaces the cancelling assertion |
| T-PG-02 | the same for `case_mean`, `control_mean`, `pvalue_list$df_vec`, `pvalue_list$gaussian_teststat`, `pvalue_list$log10pvalue`, `pvalue_list$fdr_vec`, `null_mean`, `null_sd` | [oracle] same | the current test compares one of eight outputs |
| T-PG-03 | the two paths agree when both are called with **default** arguments | [inspection] they currently do not | **[regression]** §1.5: `compute_posterior.eSVD` defaults `library_min = 0.1`, `compute_test_per_gene` defaults `1e-2`. `eSVD()` passes it explicitly, which is why this has gone unnoticed. **✅ RESOLVED (Q-PG-1) — `library_min = 0.1` everywhere** [[KZL: Yes, use 0.1]]. So `compute_test_per_gene`'s default changes from `1e-2` to `0.1`, `eSVD()`'s explicit pass-through stays, and T-PG-03 becomes a permanent guard against the two defaults drifting apart again. Document in both roxygen blocks that the value is shared |
| T-PG-04 | equivalence holds across the parameter grid: `bool_adjust_covariates ∈ {T,F}` × `bool_covariates_as_library ∈ {T,F}` × `bool_stabilize_underdispersion ∈ {T,F}` × `pseudocount ∈ {0, 1}` | [oracle] same | each boolean is implemented twice, in two different shapes (matrix sweep vs scalar multiply). 16 cells, cheap on `F-TINY` |
| T-PG-05 | `compute_test_per_gene` does **not** write `posterior_mean_mat`/`posterior_var_mat` onto the object | [invariant] its documented memory contract | the whole reason the function exists |
| T-PG-06 | on a fixture where `p = 1`, both paths still work | [invariant] | `Matrix::colMeans` on a 1-column matrix vs `mean()` on a vector — the `drop = FALSE` discipline differs between the two implementations |

### 2.12 `eSVD()` — `test_eSVD.R` **[new file]**

The main user-facing wrapper. 200 lines orchestrating the whole pipeline,
exported by accident via `exportPattern`, undocumented, and untested.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-ESVD-01 | one end-to-end run on a small `SeuratObject` returns an object with `teststat_vec`, `case_mean`, `control_mean`, `pvalue_list` | [invariant] | zero coverage on the one function a new user is most likely to call |
| T-ESVD-02 | `bool_diet = TRUE` and `bool_diet = FALSE` give the same `teststat_vec` and `log10pvalue` to 1e-8 | [oracle] each other | the two branches take *different code paths* (`compute_test_per_gene` vs the three-function chain), so this is T-PG-01 lifted to the wrapper — and it is the assertion that would catch a `library_min` drift reappearing at the wrapper level |
| T-ESVD-03 | `bool_diet = TRUE` drops `dat`, `covariates`, `fit_*` from the returned object | [invariant] documented behaviour | it is the memory contract; also confirms `compute_pvalue`'s dependence on `dat` (T-DF-03) does not silently break the diet path |
| T-ESVD-04 | `batch_var_prefix` non-`NULL`: the matching covariate columns end up in `omitted_variables` and their `z_mat` columns survive reparameterization unchanged | [invariant] | T-REP-04 lifted to the wrapper; `grep(batch_var_prefix, ...)` is a *regex* on user input, so a prefix containing `.` or `(` matches more than intended |
| T-ESVD-05 | `alpha_max = NULL` derives `2*max(dat@x)` and errors if that is `<= 0` | [invariant] | `eSVD_obj$dat@x` assumes a `dgCMatrix`; a dense `dat` would error obscurely |
| T-ESVD-06 | `intermediate_save` writes a loadable file at each of the 4 checkpoints | [invariant] | untested I/O; also the mechanism §4.1 worries about for pointer serialization. Use `withr::local_tempfile()` |
| T-ESVD-07 | `case_control_levels` in the wrong order flips the sign of every test statistic | [invariant] | the argument is documented as "Control and then Case"; getting it backwards is a silent, total inversion of the result. Worth an explicit test so the behaviour is at least *defined* |
| T-ESVD-08 | a categorical variable with one level after `droplevels` is dropped rather than erroring | [invariant] | `categorical_vars_subset` filters these out — this is the workaround for T-FMT-08 and it deserves a test |
| T-ESVD-09 | an individual with ≤ 2 cells is handled before any fitting, so one bad donor out of forty does not cost a full fit first | **✅ Q-COH-3 / Q-COH-7 — the helper drops with a warning; `eSVD()` errors** | pass 1 said *stop early*, `esvd_helper.R` *drops*; both survive, at different levels. The drop's own behaviour is now §2.17's (T-COH-01/02); this row pins the *middle* layer, and is the same assertion as T-COH-11's first half |
| T-ESVD-10 | an individual with exactly 3 cells is accepted | [invariant] the boundary | pins which side of `< 3` the threshold sits on — off-by-one here silently excludes valid donors. Unaffected by Q-COH-3 |
| T-ESVD-11 | `eSVD()` is **explicitly exported** and reachable as `eSVD2::eSVD()` | [invariant] | **[regression] and a live landmine.** `NAMESPACE` has no `export(eSVD)` — the package's main user-facing function is exported *only* by the blanket `exportPattern("^[[:alpha:]]+")` [verified]. Q-CPP-1 removes that line, which would silently un-export `eSVD` itself. One `@export` tag, but it has to land in the same commit |

### 2.13 `generate_data()` / `generate_null()`

`generate_data` has good moment-matching tests for poisson and neg_binom.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-GEN-01 | moment matching for the remaining families: `gaussian`, `curved_gaussian`, `exponential`, `neg_binom2`, `bernoulli` | [oracle] the closed-form mean/variance of each | 5 of 7 families untested in the generator that every other test's fixture depends on |
| T-GEN-02 | `generate_data` errors when `nuisance_param_vec = NA` and the family requires it | [invariant] | there is a `stopifnot(!nuisance_param_na)`; test it per family |
| T-GEN-03 | `generate_data` errors on an infeasible `nat_mat` (e.g. positive `θ` for `exponential`) | [invariant] `family$feasibility` | exercises the R-side feasibility wrapper, which is otherwise untested |
| T-GEN-04 | `generate_null()` returns the documented names, dimensions, and column names: `covariates` is `n × 5` with `c("Intercept","Log_UMI","CC","Sex","Age")`, `obs_mat` is `n × p` `dgCMatrix`, `metadata_individual` is a factor of length `n` | [invariant] the documented contract | exported and documented; the paper's Type-1-error claims rest on it and nothing checks its shape |
| T-GEN-05 | `generate_null()` marks exactly 10 genes as truly DE (columns 1–10) and the rest null | [invariant] the function's internal design | this is the fact that T-PROP-06 (null calibration) depends on; if the DE block ever moves, the calibration test silently becomes meaningless |
| T-GEN-06 | `num_individuals` odd, or `< 2`, errors — `rep(c(0,1), each = s/2)` silently misbehaves otherwise | [inspection] | `s/2` is not integer-checked |

### 2.14 Utilities — `test_utils.R` **[new file]**

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-UTIL-01 | `.mult_mat_vec(m, v) == m %*% diag(v)` and `.mult_vec_mat(v, m) == diag(v) %*% m` | [verified] holds today | the whole point of these functions is to be a fast identity; assert the identity |
| T-UTIL-02 | `.nonzero_col()` returns `numeric(0)` for an all-zero column, and matches `which(mat[,j] != 0)` otherwise, for both `bool_value` settings | [oracle] dense recomputation | hand-rolled `@p`/`@i`/`@x` pointer arithmetic; an off-by-one is invisible until it isn't |
| T-UTIL-03 | `.combine_two_named_lists()` preserves a `NULL` *value* under its name (the `TEMP_NAME` trick), and does not overwrite existing names | [invariant] | subtle and untested; it is how `param` accumulates across the pipeline, and a dropped `NULL` there means `.get_object(which_fit = "param")` fails later with `"what_obj is not found" is not TRUE` |
| T-UTIL-04 | `.get_object()` returns the right thing for **every** `what_obj` branch | [invariant] | the accessor every pipeline function routes through; 14 branches, zero tests |
| T-UTIL-05 | an unrecognized `what_obj` errors with a message naming the key | [verified] currently `stopifnot("what_obj is not found")` errors with `"what_obj is not found" is not TRUE` — a character is not a logical | §1.8 |
| T-UTIL-06 | `.get_object(what_obj = "z_mat1", which_fit = "fit_First")` errors (the `stopifnot(which_fit == "initial_Reg")` guard) | [invariant] | untested guards |
| T-UTIL-07 | `report_results()` on an incomplete object emits a `message()` and returns `invisible(NULL)` | [invariant] | the happy path is tested; the guard branch is not |
| T-UTIL-08 | `report_results()`'s `logFC` equals `log2(case_mean/control_mean)` and its `pvalue` equals `10^(-log10pvalue)` | [oracle] direct recomputation | one line each, and they are the numbers the user actually reads |
| T-UTIL-09 | `fisher_test()` matches `stats::fisher.test(..., alternative = "greater")$p.value` on several 2×2 configurations | [oracle] `stats::fisher.test` | the existing test asserts one hard-coded number to 1e-2. An external oracle is strictly better and free |
| T-UTIL-10 | `fisher_test(verbose = 1)` actually emits its diagnostic | [inspection] it does not — the `paste0` result is computed and discarded | trivial but it is a whole branch that does nothing |
| T-UTIL-11 | `fisher_test()` errors when `set1_genes ⊄ all_genes` | [invariant] | `stopifnot` today |
| T-UTIL-12 | `print.esvd_data_loader()` prints and returns invisibly | [invariant] | an exported S3 method with no test |

### 2.15 Is `Rmpfr` needed? — `test_tail_precision.R` **[new file]**

Kevin asked for tests that settle this rather than assert around it. The old
T-PVAL-07 only claimed `Rmpfr::pnorm` and `stats::pnorm` are `identical()` on
double input, which shows the *current calls* are pointless but not that
arbitrary precision is unnecessary in principle. The group below answers the
stronger question: **is there any regime this package reaches where a double
loses information that MPFR would keep?**

Everything tagged **[verified]** here was run in this session against the
installed `Rmpfr` on R 4.5.1 / macOS. The grid is
`z ∈ {−5, −10, −20, −37, −38.5, −40, −100, −1000, −1e4}`, chosen to straddle the
double underflow edge at `z ≈ −38.5`.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-MPFR-01 | `Rmpfr::pnorm(z, 0, 1, log.p = TRUE)` is `identical()` to `stats::pnorm(z, log.p = TRUE)` over the grid, and `class()` of the result is `"numeric"` — never `"mpfr"` | [verified] `identical()` is `TRUE` at every grid point; class is `numeric` | the four call sites (`multtest.R` ×4, `compute_pvalue.R` ×2, `compute_test_per_gene.R` ×2) pass plain doubles, so `Rmpfr`'s S4 method dispatches straight back to `stats`. **This is the "buys nothing" claim, and it is only half the argument** |
| T-MPFR-02 | double `stats::pnorm(z, log.p = TRUE)` agrees with a *genuine* 200-bit computation, `Rmpfr::pnorm(Rmpfr::mpfr(z, precBits = 200), log.p = TRUE)`, to **relative 1e-15** over the whole grid | [oracle] MPFR at 200 bits — the strongest oracle available | **the decisive test.** [verified] max relative error `7.3e-17` (at `z = −37`), and `4.1e-17` at `z = −1e4`. Doubles are correct to the last bit in log space across every regime this package can reach, so extra precision has nothing to buy *provided the computation stays in log space* |
| T-MPFR-03 | `2 * Rmpfr::pnorm(-40, 0, 1, log.p = FALSE)` is **exactly `0`**, identically to `2 * stats::pnorm(-40)` | [verified] both are `0` | **the failure `Rmpfr` was presumably added to prevent, and does not.** `multtest()`'s `pvalue_vec` is computed with `log.p = FALSE`, so it underflows to zero for `\|z\| > 38.5` — with or without `Rmpfr`. The zeros that reach `p.adjust` are caused by the *parameterization*, not by the precision |
| T-MPFR-04 | `stats::p.adjust(<an mpfr vector>, "BH")` silently coerces to double, so an mpfr p-value below `1e-308` comes back as exactly `0` | [verified] `p.adjust(2*Rmpfr::pnorm(Rmpfr::mpfr(c(-40,-30,-2), 200)), "BH")` returns `c(0, 1.472e-197, 0.0455)` — the first element has lost its value | **closes the argument.** Even a correctly-written MPFR pipeline cannot deliver its precision to the user, because the last step is `p.adjust`, which is double-only. Using `Rmpfr` *properly* would mean reimplementing BH in mpfr — for a difference that changes no rejection at any FDR threshold anyone uses |
| T-MPFR-05 | `stats::pnorm(z, log.p = TRUE)` matches the Mills-ratio asymptotic expansion `−z²/2 − log\|z\| − log(2π)/2 + log1p(−z⁻² + 3z⁻⁴ − 15z⁻⁶ + 105z⁻⁸)` to relative 1e-12 for `\|z\| ≥ 20` | [oracle] the analytic series — **no package dependency at all** | [verified] relative error `4.4e-13` at `z = −20` and `≤ 1.2e-16` for `\|z\| ≥ 40`. **This is the one to keep permanently**: it pins double-precision adequacy after `Rmpfr` is gone, needs nothing in `Suggests`, and costs microseconds |
| T-MPFR-06 | on end-to-end `F-NULL` and `F-SMALL` fits, `max(\|gaussian_teststat − null_mean\|/null_sd)` is recorded and asserted to lie inside the range T-MPFR-02/05 cover (say `< 200`) | [invariant] | the precision argument above is *conditional on the reachable range*. This test is what makes it stay true — if a future change starts producing `\|z\| = 10⁶`, this fails and the argument gets re-examined instead of silently expiring |
| T-MPFR-07 | with §1.1's `log.p` fix in place, `log10pvalue` for `t = 40, df = 18` is finite and equals the 200-bit MPFR computation of the same quantity to 1e-12 | [oracle] MPFR | ties T-PVAL-01 to this group: **`log.p = TRUE` is the fix, `Rmpfr` is not.** Asserting both in one place stops the two from being confused again |
| T-MPFR-08 | two genes with distinct `pvalue_list$log10pvalue` above 308 get distinct values in whatever `report_results()` returns | [invariant] | **the one place information is genuinely lost today, and MPFR is not the cure.** `report_results()` returns `pvalue = 10^(-log10pvalue)`, which underflows to `0` for `log10pvalue > 308`, so every strongly-DE gene reports `p = 0` and cannot be ranked — while `pvalue_list$log10pvalue` holds the distinction at full double precision the whole time. **The fix is to add a `log10pvalue` column to `report_results()`**, which is one line |

**Verdict, stated so it can be argued with.** Doubles in log space are exact to
the last bit over every regime this package reaches (T-MPFR-02, T-MPFR-05);
`Rmpfr` as currently called is a no-op (T-MPFR-01); it does not prevent the one
real underflow (T-MPFR-03); and its precision could not reach the user even if
it were used correctly, because BH is double-only (T-MPFR-04). **Recommendation:
drop `Rmpfr` from `Imports`, replace all eight calls with `stats::pnorm`, and
drop the install caveat from the README.** The information that *is* lost —
unrankable extreme genes in `report_results()` — is fixed by exposing
`log10pvalue`, not by adding precision (T-MPFR-08).

**⚠ NEW (Q-MPFR-1) — now resolved, see below: do T-MPFR-01..04 and T-MPFR-07 ship?** They need `Rmpfr`,
which is the dependency they exist to remove — keeping them means moving `Rmpfr`
from `Imports` to `Suggests` and guarding each with
`skip_if_not_installed("Rmpfr")`, which is a real cost for tests that answer a
question once. **Recommendation:** run T-MPFR-01..04 as a one-off decision
experiment, paste the output into `CRAN_READINESS.md` §2.2 as the evidence, and
ship only T-MPFR-05, 06 and 08 — none of which need `Rmpfr`. T-MPFR-07 keeps its
value with `stats::pnorm`-based expected values instead of MPFR ones.

**✅ RESOLVED (Q-MPFR-1)** [[KZL: Yes, let's do this outside of the package.]] —
so **T-MPFR-01, 02, 03, 04 never enter `tests/`.** They run once as a script,
their output is pasted into `CRAN_READINESS.md` §2.2 as the evidence for removing
the dependency, and `Rmpfr` leaves `DESCRIPTION` entirely — not `Imports` to
`Suggests`, but gone. Ships: **T-MPFR-05, T-MPFR-06, T-MPFR-08**, all three
`Rmpfr`-free. T-MPFR-07 ships with `stats::pnorm` expected values.

Keep the one-off script under `additional_context/` (which is `.Rbuildignore`d),
not `tests/`, so it is reproducible without being a dependency.

### 2.16 Gene status — `test_gene_status.R` **[new file, new feature]**

This is not a test of existing behaviour; it is a **new output**, so the
specification has to be pinned before the tests mean anything. Kevin's
instruction, verbatim:

> [[KZL: In general, in eSVD, let's have add a new entry to the eSVD called the
> "status" condition for each gene. Genes with status = 1 are analyzed as per
> usual. Genes with status = 2 are genes with all zeros. They are temporarily
> REMOVED prior to the initialization fitting of eSVD, and then finally put back
> into the matrix factorization at the VERY end (after p-values are calculated).
> The matrix factorization puts these genes back into the original position with
> NAs in the factorization, NAs for the test statistics, and pvalues of 1.]]

**Why this is worth doing — the current behaviour, [verified].** An all-zero gene
today goes all the way through the pipeline and misbehaves at three separate
places, none of which is visible to the user:

- `initialize_esvd()` → `glmnet::glmnet(family = "poisson")` on an all-zero
  response **warns**: `"from glmnet C++ code (error code -1); Convergence for 1th
  lambda value not reached after maxit=100000 iterations; solutions for larger
  lambdas returned"`, and returns the previous lambda's coefficients. On a real
  dataset with hundreds of all-zero genes this is hundreds of warnings, which is
  how it has stayed unnoticed — they are indistinguishable from noise.
- `estimate_nuisance()` → `gamma_rate` returns `9.706e-4` while `log_gamma_rate`
  returns exactly `-10`, its clamp boundary — a 21× disagreement (T-CPP-GAM-03).
- `compute_test_statistic()` → both group means are 0, the mixture variance is 0,
  and the Welch statistic is `0/0`.

#### 2.16.1 Specification

| Item | Decision |
|---|---|
| **Name and type** | `input_obj[["gene_status"]]`, a named **`factor`** over the *original* genes with `levels = c("analyzed", "all_zero")` — **exactly two levels, no more** (✅ Q-STATUS-4 [[KZL: Let's not add more status for now. Let's keep it at 1 and 2]]). `as.integer()` gives Kevin's 1 and 2 exactly; the factor is kept over a bare integer only because a printed object then says `"all_zero"` rather than an opaque `2`. With the enum frozen at two levels this is a readability preference, not a design necessity — say the word and it becomes an integer |
| **Length and names** | `length(gene_status) == ncol(dat)` *as passed in*, `names(gene_status) == colnames(dat)`, order preserved. It is the only record of the original gene set once the fitting object has fewer columns |
| **Definition of `all_zero`** | `Matrix::colSums(dat, na.rm = TRUE) == 0`. `na.rm` is consistent with Q-INIT-1 (NAs are zeroed on entry anyway), so a gene that is all `NA`-and-zero is `all_zero` |
| **Where the filter runs** | inside **`eSVD_helper()`**, at step 3 of the sequence in §2.17.2 — after the donor drop and the cohort checks, before `eSVD()` is called (✅ Q-COH-7 [[KZL: let's do all the filtering and gene status labeling/temporarily-removal, in eSVD_helper]]). **This supersedes Q-STATUS-1's placement**, which had `eSVD()` filtering |
| **What `eSVD()` does** | **errors**, naming the genes. It never sees an all-zero gene in normal use; the error is for a direct caller, and for catching a helper bug loudly rather than silently (✅ Q-COH-7) |
| **What `initialize_esvd()` does** | **errors** too (✅ Q-STATUS-1 [[KZL: Inside initialize_esvd(), let's error and leave filtering to `eSVD()`.]]). Three levels now test the same condition, so they must share one predicate — T-GS-17 |
| **Where the reinsertion runs** | in **`eSVD_helper()`**, after `eSVD()` returns (✅ Q-COH-7). Not in `eSVD()`'s `"Finalizing"` block, and therefore with no `bool_diet` interaction to reason about |
| **What gets `NA`** | every per-gene *estimate*: `y_mat` and `z_mat` rows, `nuisance_vec`, `posterior_mean_mat` / `posterior_var_mat` columns, `teststat_vec`, `case_mean`, `control_mean`, `pvalue_list$df_vec`, `pvalue_list$gaussian_teststat` |
| **What does not get `NA`** | `pvalue_list$log10pvalue` → `0` (i.e. p = 1) and `pvalue_list$fdr_vec` → `1`, both set explicitly (✅ Q-STATUS-3 [[KZL: The fdr_vec should be manually set to 1 as well]]). Rationale: a downstream `which(fdr < 0.05)` must not have to think about `NA` |
| **What is not padded at all** | `x_mat` — it is cells × k and has no gene dimension. Guards against a blanket "insert NA rows everywhere" implementation (T-GS-11) |
| **Multiple testing** | the empirical null (`multtest()`) and the BH adjustment are computed on **analyzed genes only** (✅ Q-STATUS-2 [[KZL: Empirical null and BH should be only done to analyzed genes.]]). Under Q-COH-7's placement this is now **structural rather than conventional**: `eSVD()` never receives an all-zero gene, so it cannot include one — see below |

**The multiple-testing point deserves to be stated explicitly, because getting it
wrong is silent and consequential.** If the status-2 genes were included in
`p.adjust(..., "BH")`, they would enlarge `n` in the `p · n / rank` formula and
**raise the adjusted p-value of every genuinely DE gene**. A cell type with 20000
measured genes of which 8000 are all-zero would lose roughly 40% of its effective
power, with no warning and no visible cause. The same applies to `locfdr`'s
empirical null, which estimates `null_mean`/`null_sd` from the *distribution* of
test statistics — a pile of degenerate genes distorts it.

Q-COH-7 makes this impossible by construction rather than by discipline, which is
the strongest form the guarantee can take. T-GS-09 and T-GS-10 stay as
architecture checks — they would fail loudly if a later refactor moved the filter
back inside `eSVD()` and got the exclusion wrong.

**Ordering: the donor drop must run before `gene_status` is computed**, and both
now live in `eSVD_helper()` (steps 1 and 3 of §2.17.2's sequence). A gene
expressed *only* in a dropped low-cell donor becomes all-zero **after** that
drop. If `gene_status` were computed first it would mark that gene `analyzed`,
and it would then go through the whole pipeline with an all-zero count vector —
exactly the failure §2.16 exists to prevent, reintroduced by ordering alone.
T-COH-07 pins the order. With both steps inside one function this is now an
internal invariant rather than a contract between two exported functions, which
is a much cheaper thing to keep true.

**Two things that removal deliberately does and does not change**, both worth a
test because both are easy to get wrong:

- **`Log_UMI` and `alpha_max` are unaffected.** `Log_UMI = log(rowSums(dat))` and
  dropping all-zero columns changes no row sum; `alpha_max = 2*max(dat@x)` reads
  stored non-zeros, and an all-zero gene stores none. Both are exact identities
  (T-GS-03, T-GS-04). If either changed, the filter would be silently altering
  the model for the retained genes.
- **The factorization for retained genes *does* change, and that is expected.**
  `.initialize_residuals()` takes the SVD of `log1p(A) − ZᵀC`. An all-zero column
  of `A` is **not** a zero column of that matrix — it is `−ZᵀC`, which is
  generally large and negative. So all-zero genes have been actively pulling on
  the factorization of every other gene. Removing them changes `x_mat`,
  and therefore changes results for analyzed genes. **This is a behaviour change,
  not a no-op, and the release note has to say so.** T-GS-05 pins the new
  definition; no test should assert the old numbers.

#### 2.16.2 Tests

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-GS-01 | `gene_status` is a factor of length `ncol(dat)` with `names() == colnames(dat)`, `levels == c("analyzed","all_zero")`, and `as.integer()` giving 1/2 | [invariant] the contract above | the contract every other test in this table reads from |
| T-GS-02 | the boundary of the definition: an all-zero gene → 2; a gene with a **single** non-zero count → 1; a gene that is all-zero **within one arm** but not overall → 1 | [invariant] | the two near-misses. The second is common in real data and must *not* be filtered — it is often the most interesting gene in the dataset |
| T-GS-03 | `Log_UMI` from `format_covariates(dat)` and from `format_covariates(dat[, status == 1])` are `identical()` | [oracle] the `rowSums` identity | if this fails the filter is changing the model for retained genes |
| T-GS-04 | `2*max(dat@x)` is identical before and after removal | [oracle] same reasoning | `eSVD()` derives `alpha_max` this way when it is `NULL` |
| T-GS-05 | every analyzed-gene output of `eSVD_helper(obj)` equals the corresponding output of `eSVD(obj[status == "analyzed", ])` element-wise to 1e-8 | [oracle] the direct filtered run **is** the definition of correctness here | **still the load-bearing test of the feature, but it now proves a narrower and more useful thing.** Under Q-COH-7 `eSVD()` never sees an all-zero gene, so this no longer asks "does the pipeline ignore them" — it asks "does the reinsertion disturb the retained entries", which is the only remaining way the feature can be wrong |
| T-GS-06 | `y_mat`, `z_mat`, `nuisance_vec`, `teststat_vec`, `case_mean`, `control_mean`, `pvalue_list$df_vec`, `pvalue_list$gaussian_teststat` are `NA` at exactly the status-2 positions and finite everywhere else | [invariant] the spec | one assertion per output; a loop over names keeps it to three lines |
| T-GS-07 | `pvalue_list$log10pvalue == 0` and `pvalue_list$fdr_vec == 1` at status-2 positions | [invariant] Kevin's "pvalues of 1" | deliberately *not* `NA` — see Q-STATUS-3 |
| T-GS-08 | `names(teststat_vec)` equals the **original** `colnames(dat)` in the original order, all-zero genes included | [invariant] | "put back into the original position" is the whole requirement; a reinsertion that appends at the end would pass T-GS-06 and fail here |
| T-GS-09 | appending 50 all-zero genes leaves `fdr_vec[status == "analyzed"]` **exactly** unchanged | [oracle] the run without them | the inflated-`n` failure. Q-COH-7 makes this structural — `eSVD()` cannot see the extra genes — so this is demoted from "catches a silent power loss" to **an architecture check**: it fails loudly if a later refactor moves the filter back inside `eSVD()` and gets the exclusion wrong |
| T-GS-10 | the same appending leaves `pvalue_list$null_mean` and `null_sd` unchanged | [oracle] same | the empirical-null half of the same point, same demotion |
| T-GS-11 | `x_mat` is `nrow(dat) × k` and contains no `NA` | [invariant] | guards against padding the cell dimension by symmetry with the gene dimension |
| T-GS-12 | `report_results()` returns one row per **original** gene in the original order, with `logFC = NA`, `pvalue = 1`, `pvalue_adj = 1` for status-2 genes, and drops none of them | [invariant] | `logFC = log2(NA/NA)` gives `NA` for free; the row must still be present, so a user's gene list round-trips |
| T-GS-13 | `eSVD_helper(bool_diet = TRUE)` and `bool_diet = FALSE` produce the same `gene_status` and the same reinserted outputs | [oracle] each other | T-ESVD-02 lifted to the helper. **Note Q-COH-7 makes this cheaper than it was**: because reinsertion happens after `eSVD()` returns, `compute_test_per_gene` does *not* need to learn about status — it simply never receives an all-zero gene |
| T-GS-14 | every gene all-zero → informative error naming the analyzed count, raised before `initialize_esvd()` | [invariant] | otherwise a zero-column matrix reaches `glmnet` |
| T-GS-15 | `k > sum(status == "analyzed")` errors naming both numbers | [invariant] | **a `k` that was legal before filtering can become illegal after it.** `k = 30` is `eSVD()`'s default and a cell type with 40 genes measured and 15 empty now fails where it used to run |
| T-GS-16 | `gene_status` survives an `intermediate_save` write/read round trip | [invariant] | **the `bool_diet` half of this is moot under Q-COH-7** — `gene_status` is attached after `eSVD()` returns, so the teardown cannot collect it. Kept for the save/load path, which still runs inside `eSVD()` and still must not corrupt what the helper attaches afterwards |
| T-GS-17 | the helper's `all_zero` label, `eSVD()`'s error and `initialize_esvd()`'s error all agree on a boundary case — a gene that is all `NA`, and a gene that is all zero except one `NA` | [invariant] one shared predicate | **the cost of Q-COH-7's three-level defence.** "All zero" is now tested in three places; if they use three expressions rather than one `.which_all_zero()`, the helper starts passing genes that `eSVD()` rejects, and the user gets an internal error from a function they never called |

#### 2.16.3 Knock-on effects on tests elsewhere in this plan

- **T-CPP-GAM-03** keeps its recording test at the unit level — `gamma_rate` is
  still callable directly — but is no longer BLOCKED.
- **T-CPP-GAM-01** (the `gamma_rate` vs `log_gamma_rate` equivalence, currently
  commented out) restricts its grid to non-degenerate inputs, which §2.16
  guarantees the pipeline now supplies.
- **T-CPP-LOAD-03** (all-zero `dgCMatrix` column in the C++ loader) stays as-is:
  the filter is a wrapper-level policy, and a user calling `opt_esvd()` directly
  can still present one. Note under Q-COH-7 there are now *three* levels above it
  that would have refused (`eSVD_helper`, `eSVD`, `initialize_esvd`), so this test
  covers a genuinely last-resort path.
- **T-PROP-10** (gene permutation equivariance) must now hold *including* the
  status-2 genes — permuting the input genes must permute `gene_status`
  identically. Under Q-COH-7 it runs against `eSVD_helper()`, since that is where
  `gene_status` now lives.
- **T-PG-01..06** (`compute_test_per_gene` vs the matrix path) are **unaffected**,
  which they would not have been under Q-STATUS-1's placement: because the
  reinsertion happens after `eSVD()` returns, neither path ever handles a
  status-2 gene and neither needs to learn about `gene_status`.
- **§1.1 fixtures**: `F-TINY` gains two all-zero genes.

### 2.17 Cohort filtering — `test_cohort_filter.R` **[new file, imported code]**

[[KZL: Q-TSTAT-2: Let's actually put `.../Was2CODE/R/esvd_helper.R` into this
package. This will help address some of the questions]] — the file is
`R/esvd_helper.R` in the `WAS2CODE_REPO` (see `CLAUDE.md` → *External
Locations*). 60 lines. It is a thin wrapper that applies four cohort-level
filters and then calls `eSVD()` verbatim through `...`.

This does answer Q-TSTAT-2, and more besides — it brings **three filters this
plan had not considered at all**. Read as a description of *the file as
imported* — §2.17.2 below is the packaged target, which differs — it says:

| Threshold | Default | Rule | Action on failure |
|---|---|---|---|
| `min_cells_per_id` | 3 | donors with `< 3` cells | **drop those donors**, keep going |
| `min_cells` | 20 | `≤ 20` cells total, after the drop | `warning()` + `return(NA)` |
| `min_ids` | 4 | `≤ 4` donors, after the drop | `warning()` + `return(NA)` |
| `min_cells_casecontrol` | 20 | `≤ 20` cells in *either* arm | `warning()` + `return(NA)` |

all behind a single `bool_check_donors = TRUE` switch. **That switch is removed
in the packaged version** (✅ Q-COH-4), and a fifth threshold `min_ids_per_arm` is
added (✅ Q-COH-2).

#### 2.17.1 What it settles

- **Q-TSTAT-2's threshold**: `min_cells_per_id = 3`, tested as
  `indiv_count < 3`, i.e. donors with 1 or 2 cells go. That matches the pass-1
  instruction ("2 or less cells") exactly, and it is an **argument**, not
  hard-coded — as recommended.
- **The `NaN` Welch df of T-TSTAT-04** is partly guarded by `min_ids`, and the
  `±Inf` statistic of §1.1 is guarded by `min_cells_per_id`. Both were open
  hazards in this plan with no proposed owner; now they have one.
- It is **already in use** across the Was2CODE analyses, so these defaults are
  empirical rather than invented. That is worth more than my recommendation was.

#### 2.17.2 The packaged specification

All five questions from pass 2 are answered. Kevin's ordered spec, verbatim:

> [[KZL: Q-COH-1: Let's not actually return NA as in `.../Was2CODE/R/esvd_helper.R`.
> Let's see if this answers your question. When using eSVD_helper (the main
> function I would like most people to use my function): 1) check the total number
> of cells (`Seurat::Cells(seurat_obj)`), and if it fails, warning and
> `return(NA)`, 2) check the number of individuals
> (`unique(seurat_obj@meta.data[,id_var])`), and if it fails, warning and
> `return(NA)`, 3) check the number of cells in each class
> (`table(seurat_obj@meta.data[,case_control_var])`), and if it fails, warning and
> `return(NA)`, 4) remove any individuals with only 2 or less cells (only a
> warning), 5) check of genes with all zero's — if there are, apply the gene
> status conditions and proceed after temporarily removing those genes.]]

| # | Answer | Effect on this section |
|---|---|---|
| **Q-COH-1** | Keep `warning()` + `return(NA)`; no classed sentinel | my `eSVD_skipped` proposal is declined. The three rejections stay distinguishable only by their warning *text*, so each must name its own threshold — T-COH-05 asserts three distinct `regexp` fragments, which is what makes an analysis loop able to tell them apart at all |
| **Q-COH-2** | ✅ add `min_ids_per_arm`, default 2 | T-COH-04 becomes a passing test of new behaviour rather than a failing test of imported behaviour |
| **Q-COH-3** | ✅ **drop, with a warning** (step 4, "only a warning") | settles the pass-1/pass-2 conflict. T-ESVD-09 and T-VAL-31 can now be written. The warning must name the donors — see T-COH-02 |
| **Q-COH-4** | ✅ **`bool_check_donors` is removed entirely.** Every filter always runs; users change thresholds instead | better than the split switch I proposed — see the note on disabling below |
| **Q-COH-5** | Keep `eSVD_helper()` as the wrapper (it takes a `Seurat` object, which is how most users will arrive) **and** export `filter_cohort()` for the donor drop | two exported functions, not one. Raises Q-COH-7 below |

**Disabling a filter, now that the switch is gone.** This works and is worth
documenting explicitly, because it is the only remaining escape hatch: every
threshold compares with `<=` or `<`, so setting it to `0` disables it —
`length(Cells(obj)) <= 0` is `FALSE` for any non-empty object, and
`indiv_count < 0` is never `TRUE`. T-COH-12 pins it. **One caveat worth a
`stopifnot`:** `min_cells_per_id = 1` or `2` does *not* disable the drop, it
weakens it, and a surviving 1-cell donor gives zero within-individual variance, a
`±Inf` statistic, and the silent `locfdr` degradation of §1.1. Recommend
`stopifnot(min_cells_per_id == 0 || min_cells_per_id >= 3)` — off, or safe, but
not in between.

##### ✅ Q-COH-6 — drop first, then check

[[KZL: Q-COH-6: Drop first, then check.]] — as the imported file already does.
Every threshold therefore describes the cohort that is actually analyzed, and the
`min_ids_per_arm = 2` added by Q-COH-2 is evaluated against the donor count the
model will really see. T-COH-03 asserts this order.

##### ✅ Q-COH-7 — all preprocessing *and* postprocessing lives in `eSVD_helper()`; `eSVD()` refuses

> [[KZL: Q-COH-7: Yes, let's do all the filtering and gene status
> labeling/temporarily-removal, in eSVD_helper. Then, eSVD() errors on any of
> these corner cases that accidentally get remained. Then, in eSVD_helper, after
> eSVD() is run, put those temporarily removed genes back in the correct
> locations]]

This **supersedes Q-STATUS-1's placement** — `gene_status` moves out of `eSVD()`
and into `eSVD_helper()`. The resulting call graph is one clean rule, *the
wrapper filters, the pipeline refuses*, applied at every level:

| Level | Low-cell donors | All-zero genes |
|---|---|---|
| `eSVD_helper()` — Seurat in, the recommended entry point | **drops**, with a warning naming them (Q-COH-3) | **labels, removes, and reinserts after `eSVD()` returns** |
| `eSVD()` | **errors**, naming them | **errors**, naming them |
| `initialize_esvd()` | n/a | **errors** (✅ Q-STATUS-1) |

**The definitive order inside `eSVD_helper()`**, combining Q-COH-6 with the
five-step spec. Every §2.17 and §2.16 test is written against this sequence:

1. **Drop** donors with `< min_cells_per_id` cells; `warning()` naming them.
2. **Check** `min_cells`, `min_ids`, `min_ids_per_arm`, `min_cells_casecontrol`
   on the *filtered* object; on failure `warning()` + `return(NA)`.
3. **Label** `gene_status` from the donor-filtered counts, and remove the
   `all_zero` genes from the object.
4. **Call `eSVD()`** on what remains.
5. **Reinsert** the removed genes at their original positions, per §2.16.

**This is a better architecture than the one I proposed, for a reason worth
recording.** Under Q-STATUS-1's placement, "the empirical null and BH are
computed on analyzed genes only" (Q-STATUS-2) was a *convention* that
`compute_pvalue()` had to be careful to honour, and T-GS-09/T-GS-10 existed to
catch it being broken. Under this placement `eSVD()` **never sees** an all-zero
gene, so `multtest()` and `p.adjust()` cannot include one. The statistical
requirement becomes structural instead of conventional. T-GS-09/T-GS-10 stay, but
demoted from "catches a silent power loss" to "confirms the architecture holds" —
which is the right place for them to end up.

**Three consequences that change what gets written:**

- **`eSVD()`'s returned object has no `gene_status` field.** It is created by the
  helper and attached after `eSVD()` returns. Every §2.16 test re-points at
  `eSVD_helper()`; §2.16.1's placement rows are rewritten accordingly.
- **The removal is a `Seurat` *feature* subset, not a matrix subset.** Both
  `eSVD()` and `eSVD_helper()` take `seurat_obj`, so step 3 is
  `seurat_obj[keep_genes, ]` and step 1 is `seurat_obj[, keep_cells]`. Feature
  subsetting a `Seurat` object is not side-effect-free — variable features,
  `scale.data` and any reductions are affected — so the tests must assert that
  what `eSVD()` actually consumes (the `counts` layer and `meta.data`) survives
  intact. T-COH-13.
- **"All zero" is now defined in three places** — the helper labels it,
  `eSVD()` errors on it, `initialize_esvd()` errors on it. They must share one
  internal predicate (`.which_all_zero(dat)`), or the definitions drift and the
  helper starts passing genes that `eSVD()` rejects. T-GS-17.

**One point to confirm in passing, not a blocking question.** "Errors on any of
these corner cases" reads naturally as the *correctness* conditions — all-zero
genes, `≤ 2`-cell donors, and one donor in an arm (which is T-TSTAT-04's `NaN`
df). It should **not** include the four cohort *power* minima: a legitimate if
underpowered 15-cell run should still be possible through `eSVD()` directly, and
those thresholds are the helper's policy, not the model's. §2.17.3 is written that
way.

**Two mechanical points, no decision needed:**

- **`Seurat::Cells()` → `SeuratObject::Cells()`.** `Seurat` is `Suggests`-only
  and this is the one call that needs it; `SeuratObject` supplies `Cells` and is
  what `eSVD()` already uses. (Kevin's spec repeats `Seurat::Cells`, which is
  fine as pseudocode.)
- **`utils::getFromNamespace("eSVD", "eSVD2")` becomes a direct call**, and the
  `requireNamespace("eSVD2")` guard at the top goes away, once the function is
  inside the package.
- **Name**: the file defines `esvd_helper`; Kevin writes `eSVD_helper`. Pick one —
  `eSVD_helper` reads better beside `eSVD` and `eSVD2`.

#### 2.17.3 Tests

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-COH-01 | a donor with 1 cell and a donor with 2 cells are dropped; a donor with 3 is kept | [invariant] the `< min_cells_per_id` boundary | the threshold Q-TSTAT-2 asked for, pinned on both sides |
| T-COH-02 | the drop **warns and names the donors dropped** | [invariant] ✅ Q-COH-3 confirms warn-and-proceed | **fails against `esvd_helper.R` as written**, which drops silently. "Only a warning" settles that it warns; that it *names* them is the part still worth asserting — a cohort quietly losing donors between input and output is discovered at revision time |
| T-COH-03 | the donor drop runs **first**, and every count check is then evaluated on the *filtered* object | [inspection] ✅ Q-COH-6 [[KZL: Drop first, then check.]] | checking before dropping admits exactly the cohorts the checks exist to reject — 6 donors of which 2 have 1 cell passes a 6-donor check, then becomes 4, and `min_ids_per_arm` ends up checked against a count the model never sees |
| T-COH-04 | a 4-case / 1-control cohort is **rejected** by `min_ids_per_arm = 2`, not run | [invariant] ✅ Q-COH-2 | the file as imported would run it: 5 donors passes the pooled `min_ids > 4`, then `.compute_df()` yields `n2 - 1 = 0` and a `NaN` df. This is T-TSTAT-04 reached through the front door, and the new argument is what closes it |
| T-COH-05 | each of the five filters fires on its own boundary and not on its neighbour's, and each of the three rejections emits a **distinct** warning naming its own threshold | [invariant] ✅ Q-COH-1 keeps the bare `NA` | with the return value fixed at `NA`, the warning text is the *only* channel telling a caller which filter fired. Three `expect_warning(regexp = ...)` fragments, per C-11. This is what makes an analysis loop over cell types diagnosable |
| T-COH-06 | `min_cells_per_id` of `1` or `2` is **rejected at entry** — off (`0`) or safe (`>= 3`), never in between | [invariant] ✅ Q-COH-4 removed the switch, so thresholds are the only escape hatch | a surviving 1-cell donor gives zero within-individual variance, `±Inf`, and the silent `locfdr` degradation of §1.1. With `bool_check_donors` gone this is the only way back into that state, and a `stopifnot` closes it |
| T-COH-07 | **ordering**: on a fixture where gene *g* is non-zero **only** in a 2-cell donor, `g` ends up `all_zero` in `gene_status` — i.e. the donor drop ran first | [invariant] | the §2.16 ordering constraint. If the order inverts, `g` is marked `analyzed` and goes through the pipeline with an all-zero count vector, which is the whole failure mode §2.16 removes |
| T-COH-08 | the filters change nothing when no threshold binds: a clean cohort gives results `identical()` to calling `eSVD()` directly | [oracle] the direct call | the wrapper must be a no-op on good data, or every existing result silently moves |
| T-COH-09 | `...` reaches `eSVD()`: a non-default `k` or `library_min` passed to the wrapper arrives intact | [invariant] | the passthrough is the whole rest of the function, and a typo'd argument name would be swallowed by `...` in silence |
| T-COH-10 | all **five** thresholds — `min_cells_per_id`, `min_cells`, `min_ids`, `min_ids_per_arm`, `min_cells_casecontrol` — are reachable as arguments, documented, and each demonstrably changes an outcome | [invariant] | five user-facing knobs; `R CMD check` catches an undocumented one, nothing catches a knob that has no effect. `bool_check_donors` must be **gone**, not deprecated |
| T-COH-11 | `eSVD()` called directly **errors** on a 1-cell donor and on an all-zero gene, naming them; and does **not** error merely for being underpowered (15 cells, 3 donors) | [invariant] ✅ Q-COH-7 | the correctness/power split. Without the errors, `eSVD()` stays reachable on the §1.1 `±Inf` path through the documented API the vignette uses. With them applied too broadly, a legitimate small run becomes impossible outside the helper |
| T-COH-13 | after the helper's cell and feature subsetting, the `counts` layer and `meta.data` that `eSVD()` actually consumes are intact, and `Log_UMI` recomputed from the subset matches the identity in T-GS-03 | [oracle] direct recomputation | **Q-COH-7 makes both filters `Seurat` object subsets**, `obj[, keep_cells]` and `obj[keep_genes, ]`. Feature subsetting a `Seurat` object touches variable features, `scale.data` and reductions; none of those matter to `eSVD()`, but a subset that quietly dropped or reordered the counts layer would |
| T-COH-14 | the five-step order of §2.17.2 holds end to end: a fixture where a gene is non-zero only in a 2-cell donor **and** the cohort is borderline on `min_ids` exercises steps 1, 2 and 3 in sequence | [invariant] | T-COH-03 and T-COH-07 each pin one adjacent pair; this pins the whole chain, which is where an ordering regression would actually land |
| T-COH-12 | setting a threshold to `0` disables that filter | [invariant] the `<=` / `<` comparisons | ✅ Q-COH-4 makes this the only escape hatch, so it needs to be a tested contract rather than an emergent property of the comparison operators |

**Fixture note — mostly defused by Q-COH-5.** `F-TINY` (6 donors × 120 cells,
3 per arm) passes `min_ids > 4` by one and passes `min_cells_casecontrol` at 60
cells per arm, but `min_ids_per_arm = 2` against 3 donors per arm leaves no
margin at all. Because Q-COH-5 keeps `eSVD_helper()` a *sibling* of `eSVD()`
rather than a layer inside it, **every §2 test that calls `eSVD()` is unaffected**
— which is the main reason that answer is convenient. The §2.17 tests build their
own purpose-built cohorts, and only T-COH-08 (the no-op check) needs a fixture
that clears all five thresholds comfortably. Worth sizing that one deliberately
rather than reusing `F-TINY`.

---

### 2.18 Log2 fold change and its standard error — `test_compute_log_fold_change_claude.R` **[new file, new feature]**

Added 2026-09-28 with eSVD2 1.1.0. Unlike the rest of this plan, this section
was written *with* the tests, not before them: Kevin asked for the feature and
the tests together and approved the list below as a plan-mode plan. All 20
(T-LFC-01 to -19, with -08b) are implemented and green (306 expectations). The statistic is

```
log2fc    = log2(case_mean / control_mean)
log2fc_se = (1 / ln 2) * sqrt(case_var / (n1 * case_mean^2) + control_var / (n0 * control_mean^2))
```

with `n1`, `n0` the number of *individuals* and `case_var`, `control_var` the
mixture variances the Welch statistic already divides by.

**Fixtures**, all built in the test file: a hand-built pair of posterior
matrices (3 case / 5 control individuals, unequal cells, shuffled rows); a
donor-level pair (identical cells, zero posterior variance); F-SMALL fitted
through both `opt_esvd` rounds; `generate_null()` cohorts (8 individuals, 20
cells each, 40 genes) fitted the same way.

**A. The formula**

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-LFC-01 | the four new vectors equal an explicit per-individual loop; `log2fc_se_vec * log(2)` is the natural-log SE | [oracle] independent loop, centred variance, no averaging matrix | the headline numbers, and the scale |
| T-LFC-02 | SE equals `sqrt(g' Σ g)` with `g` from `numDeriv::grad` of `log2(a/b)` | [oracle] numerical delta method | an algebra slip in the closed form cannot also be in a numerical gradient |
| T-LFC-03 | on donor-level data `case_var = (n-1)/n * stats::var(.)`, the log2 SE is the textbook form, and with equal arms the linear SE is `stats::t.test()$stderr * sqrt((n-1)/n)` | [oracle] `stats::var`, `stats::t.test` | the only external oracle; pins the no-Bessel decision for the SE as T-TSTAT-01a does for the statistic |
| T-LFC-04 | `(case_mean - control_mean) / sqrt(case_var/n1 + control_var/n0)` equals `teststat_vec` | [invariant] | the stored variances are the ones the statistic used |
| T-LFC-05 | swapping arms negates `log2fc_vec`, keeps the SE, swaps the variances | [invariant] | |
| T-LFC-06 | means × c and variances × c² leave both unchanged, c from 0.01 to 100 | [invariant] | a log ratio has no units; catches an unsquared mean |
| T-LFC-07 | k× the cells leaves the SE unchanged; k× the individuals divides it by `sqrt(k)` | [invariant] | the unit of replication is the individual; a cell-level SE fails the first half |
| T-LFC-08 | a non-positive arm mean gives `NA` in both, one warning, other genes untouched | decision D4 | `log2()` would return `NaN` / `-Inf` silently |
| T-LFC-08b | at the helper: a `NaN` variance, a negative mean beside an `NA` variance, a negative variance, one `NA` among valid inputs, a zero mean and an infinite mean all give `NA_real_` in both outputs (never `NaN` or `Inf`) under one warning counting 6; a gene whose four inputs are all `NA` is silent | decision D4 | added after code review, which found an SE of `NaN` and a guard skipped by a single `NA` |
| T-LFC-09 | every new vector is a plain named numeric vector carrying the gene names | [invariant] | |

**B. The plumbing**

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-LFC-10 | `compute_test_per_gene` matches the matrix path on the four vectors (1e-8) and records the individuals in `param` | [oracle] the other implementation | two implementations, one oracle for free |
| T-LFC-11 | `compute_log_fold_change()` matches a recomputation from the posterior matrices via `rowsum()`; idempotent | [oracle] | |
| T-LFC-12 | `report_results()` has `logFC_se` directly after `logFC`, equal to the stored vectors, finite and positive; a gene blanked in the stored vectors is `NA` in both columns | [oracle] | the numbers the user reads; `logFC` is read from the object and not recomputed |
| T-LFC-13 | both `bool_diet` paths of `eSVD()` agree; `compute_log_fold_change()` works on a diet object; `eSVD_helper()` gives `NA` at exactly the all-zero genes (positions 3 and 12) | [invariant] | reinsertion by position, and the reason the variances are stored |
| T-LFC-14 | an object without `case_var` is refused by `compute_log_fold_change()` by name; `report_results()` returns `logFC_se = NA` with a warning | decision D5 | objects saved before 1.1.0 |
| T-LFC-19 | after a stale `param` is planted, `compute_test_statistic()` and `compute_test_per_gene()` both refresh the recorded individuals, and `compute_log_fold_change()` then reproduces the SE | [invariant] | added after code review: `.combine_two_named_lists` never overwrites, so a rerun on a changed cohort left the old individuals and the SE was divided by the wrong number |

**C. Does it behave as a standard error**

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-LFC-15 | the between-individual part of the SE matches a bootstrap over individuals (B = 2000, fit held fixed) within 10%, on the genes with between-individual CV ≤ 0.2. `skip_on_cran()` | [oracle] bootstrap | the linearization and the 1/n scaling, at 4 individuals per arm |
| T-LFC-16 | planted genes have the right sign; on genes whose nuisance estimate has not diverged, the generator's truth is within 3 SE. `skip_on_cran()` | [oracle] generator truth from `nat_mat` | recovery |
| T-LFC-17 | over 30 refitted `generate_null()` cohorts, on the exactly-null genes: `sqrt(mean(log2fc²) / mean(se²)) ≤ 1.2` and ±2 SE covers 0 at ≥ 0.831, for all genes and for the not-diverged genes. `skip_on_cran()` | [oracle] repeated sampling | the only test that does not hold the fit fixed |
| T-LFC-18 | setting every nuisance rate to 1e6 shrinks the within part of the SE by more than 100× and leaves the SE equal to its between part | [invariant] posterior variance is `mu / r` for large `r` | the mechanism behind the restriction in T-LFC-16 and -17 |

**Teeth.** Fourteen deliberate breakages were applied one at a time in a
scratch copy of the final code; every one turned at least one test red.

| Breakage | Tests that go red |
|---|---|
| drop `1 / ln 2` from the SE | 01 02 03 08b 11 15 17 18 |
| divide by cells, not individuals | 01 02 03 07 10 11 13 15 16 17 18 19 |
| means not squared | 01 02 03 06 08b 11 15 18 |
| arm means swapped in the SE | 01 02 03 11 15 18 |
| natural-log fold change | 01 08b 11 |
| per-gene path stores the wrong variance | 10 13 19 |
| reinsertion does not pad the new vectors | 13 |
| no guard on invalid genes | 08 08b |
| Bessel-corrected variance | 01 03 07 10 11 13 15 19 |
| `report_results()` reports the natural-log SE | 12 |
| `report_results()` recomputes `logFC` from the arm means | 12 |
| a single `NA` input counts as padded | 08b |
| recorded individuals not refreshed, matrix path | 19 |
| recorded individuals not refreshed, per-gene path | 19 |

T-LFC-04, -05, -09 and -14 were not turned red by any of the fourteen; they
guard properties none of these breakages touches (the stored variances being
the statistic's, the arm swap, the names, the pre-1.1.0 object).

**What section C found.** These are properties of the statistic, measured on
the toy cohorts, and they are why two of the tests are narrower than first
planned.

1. *The within term is most of the SE*: a median of 98% of `log2fc_se²` on
   F-SMALL and 95% on `generate_null()`. The SE is a median of 7.8 and 4.3
   times the fit-fixed bootstrap SD.
2. *Against refitting it is about right.* Over 30 refitted cohorts, on the 88%
   of genes with an ordinary nuisance estimate: calibration ratio 0.774,
   coverage of ±2 SE 0.959.
3. *Where the nuisance estimate diverges the SE is anti-conservative.* On the
   other 12% (estimated rate near 1e7, true rates 0.1 to 10): calibration
   ratio 1.538, coverage 0.357. On F-SMALL, 5 of the 7 such genes miss their
   truth by more than 3 SE, the planted gene 2 by 12.5. This is the
   downstream face of Q10 in `CRAN_READINESS.md` (the nuisance blow-up).
4. *The delta method understates the between part at high CV*: by up to 20%
   at a between-individual CV near 1 with 4 individuals per arm. Inside
   CV ≤ 0.2 the bootstrap / delta ratio is 0.98 to 1.03.
5. *The depth adjustment shifts every fold change.* On F-SMALL the 35 null
   genes have a median estimate of −0.21: the five planted genes raise a case
   cell's total count by 2^0.23. Standard errors of 0.3 to 0.6 hide it.

**Corrections made to the tests after their first run**, recorded because a
test changed after it failed deserves a second look:

- T-LFC-15 first applied its 10% tolerance to every gene. The tolerance had
  been derived for CV ≤ 0.2 and is now asserted there. Its second assertion,
  "the reported SE is never below the bootstrap SD", was removed: it is not a
  property of the statistic (finding 4 with a vanishing within term), and
  passed only through its 5% allowance.
- T-LFC-16 first failed on one gene, the planted gene 2 of F-SMALL, whose
  nuisance estimate had diverged (finding 3). The 3-SE assertion is now made
  on the genes whose estimate has not diverged. The sign assertion is made on
  every planted gene.
- T-LFC-18 was added, to state the mechanism as a property of the model.
- After `/code-review` of the diff: T-LFC-08b and T-LFC-19 were added and
  T-LFC-12 extended, each for a defect the review found in the first draft
  of the code (see their *Why*); T-LFC-15 and T-LFC-16 were marked
  `skip_on_cran()`, because they put thresholds on the output of an iterative
  fit that is known to amplify last-digit differences. **Three of the tests
  that say most about the SE (15, 16, 17) therefore run only with
  `NOT_CRAN=true`**, which `devtools::test()` sets and a bare `R CMD check`
  does not.

**Questions for Kevin.**

- **Q-LFC-1.** Finding 3: leave the SE as it is and document it (done in
  `?compute_log_fold_change`), or act on it? Options: repair the nuisance
  divergence upstream (Q10), which removes the cause; or flag such genes in
  `report_results()`.
- **Q-LFC-2.** Decision D4 makes `compute_test_statistic.default()` warn on a
  non-positive arm mean. One existing test, T-TSTAT-06, feeds it mean-zero
  Gaussian matrices and now asserts that warning. Keep the warning, or return
  `NA` silently from the matrix method?
- **Q-LFC-3.** Names: `log2fc_vec` / `log2fc_se_vec` on the object, `logFC` /
  `logFC_se` in `report_results()`.

The comparison against DESeq2, dreamlet and NEBULA is
`lfc-se-comparison_2026-09-28_claude.R` in this folder. It is not a test.

## 3. C++ backend tests

Restricting to correctness. Everything here is driven from R through the
`RcppExports` bindings; no C++ test framework is proposed.

**⚠ BLOCKED (Q-CPP-1) — now resolved, see below.** §2.3 of the readiness doc proposes *un-exporting*
`objfn_Xi_r`, `grad_Xi_r`, `hessian_Xi_r`, `feas_*_r`, `data_loader`,
`esvd_family`, `gamma_rate`, `log_gamma_rate` — but the entire C++ test strategy
below calls them. The resolution is that they become internal (`:::`) rather than
exported, and the tests use `:::`. This is normal for a package's own test suite
and does not weaken the tests. Confirm before we start.

**✅ RESOLVED (Q-CPP-1)** [[KZL: Yes, let's do this]] — the bindings become
internal and every test in §3 calls them as `eSVD2:::objfn_Xi_r(...)` etc. Two
mechanical consequences: `exportPattern` in `NAMESPACE` has to go (it is what
exports them by accident today), and each function still needs a roxygen block
with `@noRd` so `R CMD check` does not complain about undocumented objects.

### 3.1 Data loader — `test_data_loader.R` **[new file]**

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-CPP-LOAD-01 | dense `double`, dense `integer`, and `dgCMatrix` versions of the *same* matrix produce identical `objfn_all_r` values | [oracle] each other | three loader implementations of one interface; this is the cheapest way to test all three at once |
| T-CPP-LOAD-02 | a matrix with `NA` entries: `objfn_all_r` skips them (the `Flag::na` path) and the value equals the same objective computed in R over the non-`NA` entries only | [oracle] an R reimplementation | the `num_non_na` normalisation divides by a count the R side never sees; if the count were wrong, every loss would be scaled wrongly *and monotonicity would still hold*, so T-OPT-02 cannot catch it |
| T-CPP-LOAD-03 | a `dgCMatrix` column with **zero** stored non-zeros iterates correctly (the `m_nnz < 1` early return in `SparseVecIterator::value`) | [oracle] the dense equivalent | an all-zero gene is common in real single-cell data and hits a code path nothing else does |
| T-CPP-LOAD-04 | a `dgCMatrix` with an explicitly-stored `NA` in the value vector is flagged `na`, not `regular` | [oracle] the dense equivalent | `SparseVecIterator::value` checks `NumericVector::is_na(val)` after the index match; the interaction of "stored NA" with `m_innerpos` advancement is intricate |
| T-CPP-LOAD-05 | `data_loader()` on an unsupported type errors **at construction** with a message naming the type | [verified] it does not: on a `dgeMatrix` and on a `character` matrix it returns without error, and the failure surfaces at first use as `Error: external pointer is not valid` | `data_loader.cpp` only reaches its `Rcpp::stop("unsupported matrix type")` for a non-S4 non-matrix. An S4 that is not `dgCMatrix`, or a dense matrix that is neither integer nor numeric, leaves `loader = nullptr` and returns a null `XPtr`. A `dgeMatrix` is entirely plausible user input |
| T-CPP-LOAD-06 | `test_data_loader()` output for a small mixed matrix is stable | [snapshot] | a debug printer; a snapshot test is proportionate. **Or delete the function** — it is currently exported and undocumented (§2.3) |

### 3.2 Analytic derivatives — `test_family_derivatives.R` **[new file]**

`numDeriv` is already in `Suggests` and has never been used. **This is the
highest-value C++ test available**: an error in any `family_*.cpp` derivative
produces a silently mis-converged fit with no symptom at all — the loss still
decreases, the optimizer still terminates, and the answer is wrong.

For each of the 7 families, at a feasible point from `feasible_point(family)`:

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-CPP-FAM-01 | `grad_Xi_r(x, ...) ≈ numDeriv::grad(function(v) objfn_Xi_r(v, ...), x)`, tolerance 1e-6 | [oracle] numerical differentiation | 7 families × 1 assertion |
| T-CPP-FAM-02 | `hessian_Xi_r(x, ...) ≈ numDeriv::hessian(...)`, tolerance 1e-5 | [oracle] same | 7 families |
| T-CPP-FAM-03 | `grad_YZj_r` vs `numDeriv::grad` over the **free** coordinates only (`YZind`) | [oracle] same | 7 families; the `YZind` subsetting is itself a source of index bugs |
| T-CPP-FAM-04 | `hessian_YZj_r` vs `numDeriv::hessian` over `YZind` | [oracle] same | 7 families; also confirms `subset_matrix()`'s hand-rolled `ind[j]*n` pointer arithmetic |
| T-CPP-FAM-05 | derivatives are correct at a point where the data column contains a **zero** and an **`NA`** | [oracle] same | every family implements `log_prob_single` **twice** — a general overload and an `Aij == 0` special case. The zero overload is only reached through the sparse loader, so a discrepancy between the two overloads shows up only on sparse input. Poisson, gaussian, neg_binom, neg_binom2, bernoulli, exponential, curved_gaussian all have this duplication |
| T-CPP-FAM-06 | `objfn_all_r` equals the sum of per-column `objfn_YZj_r` values, up to the documented normalisation | [invariant] | `objfn_all_r` divides by `total_non_na` while `objfn_YZj_r` divides by its own column's count. If those disagree, `opt_esvd`'s convergence threshold is being compared against a differently-scaled loss than the one the inner optimizer minimises |
| T-CPP-FAM-07 | the `l2pen` terms: with `l2penx = 0` the gradient equals the unpenalized one; with `l2penx = c` the difference is exactly `2*c*XCi[1:k]` | [oracle] closed form | penalties are added in four places (`objfn`, `grad`, `hessian`, `direction`) and must agree |
| T-CPP-FAM-08 | `direction_Xi`'s returned `grad` equals `grad_Xi_r`, and its `direction` equals `-solve(hessian_Xi_r, grad_Xi_r)` when the Hessian is PD | [invariant] | `direction_*` duplicates `grad_*` and `hessian_*` for speed. Two implementations, one answer — a free equivalence test. **Note the scalings differ between the LLT branch and the fallback branch** (`direc = -g/non_na` vs `-H⁻¹g`); assert whichever is intended |
| T-CPP-FAM-09 | `feas_Xi_r`/`feas_YZj_r` return `TRUE` for the 4 unconstrained families always, and correctly reject out-of-domain `θ` for `curved_gaussian` (`θ > 0`), `exponential` and `neg_binom` (`θ < 0`) | [invariant] the declared `domain()` | feasibility is what keeps the line search inside the parameter space; a wrong sign silently permits `log(-θ)` of a positive number |
| T-CPP-FAM-10 | `esvd_family("nonsense")` errors | [invariant] | there is an `Rcpp::stop`; test it |
| T-CPP-FAM-11 | `esvd_family(f)$domain` and `$feas_always` match the family's documented domain, for all 7 | [invariant] | these fields are read by `generate_data` and by the R-side `feasibility()` closure |
| T-CPP-FAM-12 | the R-side `.dat_to_nat.*` / `.nat_to_canon()` round-trip: `nat_to_canon(dat_to_nat(x)) ≈ x` where the family admits it (`gaussian`, `poisson` up to `tol`, `curved_gaussian`, `exponential`, `neg_binom2`) | [invariant] | 7 families × 2 conversions, entirely untested. Note `bernoulli` and `neg_binom` do **not** round-trip and the test must say so explicitly rather than skip them silently |

### 3.3 Constrained Newton — `test_constrained_newton.R` **[new file]**

Not directly exported; tested through `opt_x`/`opt_yz`.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-CPP-NEWT-01 | on a strictly convex quadratic-like problem (gaussian family, `l2pen > 0`), one `opt_x` call lands within 1e-6 of the closed-form ridge solution | [oracle] closed form | the only place a Newton solver's *answer* can be checked exactly |
| T-CPP-NEWT-02 | when the line search fails, the failure is **reported** — not merely warned about — and `opt_esvd` does not report convergence | [invariant] | **[regression]** §4.2: `line_search` returns `step = 0.0` and an unchanged iterate. `opt_esvd.default`'s convergence test then sees no change and `break`s, *reporting convergence*. **A failed optimization is currently indistinguishable from a converged one.** Requires the status flag from §4.2's fix |
| T-CPP-NEWT-03 | under `options(warn = 2)`, a line-search failure does not corrupt state or crash | [invariant] | §4.2: `Rcpp::warning` becomes an error and `longjmp`s past the destructors of live `MatrixXd`/`NumericVector` objects. The fix (return a flag, warn from R) makes this test trivially pass; without the fix the test is the demonstration |
| T-CPP-NEWT-04 | a singular Hessian takes the gradient-descent fallback and still decreases the objective | [invariant] | `compute_direction` falls back to `-g` when `Eigen::LLT` reports failure; untested |
| T-CPP-NEWT-05 | the early-exit path (`‖grad‖ <= eps_rel·max(1,‖x‖)` at iteration 0) returns the input unchanged | [invariant] | starting from the optimum should be a no-op |

### 3.4 `opt_x` / `opt_yz` — `test_optimization_helper.R` **[strengthen]**

The existing test is good: it checks that `C` is unchanged, that the fixed column
of `Z` is unchanged, and that the loss decreases twice. Extensions:

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-CPP-OPT-01 | the same, for all 7 families with valid `gamma` | [invariant] | Poisson only today |
| T-CPP-OPT-02 | `inplace = TRUE` modifies `XC_init` in place; `inplace = FALSE` leaves it untouched and returns a new matrix | [invariant] | the default differs between the C++ signature (`inplace = true`) and the R wrapper (`inplace = FALSE`) — a caller relying on either has a 50% chance of surprise |
| T-CPP-OPT-03 | `fixed_cols = integer(0)` (nothing fixed) and `fixed_cols = 1:ncol(YZ)` (everything fixed, so `YZind` is empty) both behave sanely | [invariant] | the empty-`YZind` case makes `subset_vector`/`subset_matrix` produce 0-length objects and `Eigen::LLT` factor a 0×0 matrix |
| T-CPP-OPT-04 | `s` (library multiplier) containing a zero or a negative value errors, rather than producing `log(0) = -Inf` inside `compute_fn_si` | [invariant] | Poisson's `fn_si = log(si)`; `si = 0` gives `-Inf` and every `exp(fn_si + theta)` becomes 0, so the gene contributes nothing and the loss looks *better*. Nothing validates `s` |
| T-CPP-OPT-05 | `opt_x` called twice in a row is idempotent to 1e-6 (already at the optimum) | [invariant] | a cheap convergence sanity check |

### 3.5 `gamma_rate` / `log_gamma_rate` — `test_gamma_rate.R` **[strengthen]**

The existing test checks that `gamma_rate` returns one finite non-negative number
per gene on one well-behaved configuration. The `log_gamma_rate` half is
commented out.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-CPP-GAM-01 | `gamma_rate(x, mu, s)` and `exp(log_gamma_rate(x, mu, s))` agree to ~1e-4 across a grid of `(mu, s)` regimes | [oracle] each other — two routes to the same MLE | **the commented-out half of the existing test.** This is a second free equivalence test, exactly like §1.5's. It is commented out, which suggests it was failing; finding out why is a prerequisite. **We now have a strong candidate for why** — see T-CPP-GAM-03: on a degenerate (all-zero) gene the two routines disagree by a factor of 20, because `log_gamma_rate` saturates at its clamp while `gamma_rate` does not. Restrict this test's grid to non-degenerate `(x, mu, s)` and it should pass; §2.16's filter is what guarantees the pipeline never feeds it a degenerate one |
| T-CPP-GAM-02 | the estimate maximizes the likelihood: `objfn(b̂) >= objfn(b)` for `b` on a grid around `b̂` | [oracle] the R implementation of `objfn` that already exists **in the comment block of `src/gamma_rate.cpp` lines 58–75** | the comment ships a complete R reference implementation of the objective, gradient and Hessian. Turning it into a test is nearly free and gives an independent oracle for both C++ routines |
| T-CPP-GAM-03 | record what the two routines *actually* return for an all-zero gene, and assert they are finite — but **do not** try to make them agree | **✅ RESOLVED (Q-GAM-1) — by removing the gene, not by fixing the estimator.** [[KZL: In general, in eSVD, let's have add a new entry to the eSVD called the "status" condition for each gene. Genes with status = 1 are analyzed as per usual. Genes with status = 2 are genes with all zeros. They are temporarily REMOVED prior to the initialization fitting of eSVD, and then finally put back into the matrix factorization at the VERY end (after p-values are calculated). The matrix factorization puts these genes back into the original position with NAs in the factorization, NAs for the test statistics, and pvalues of 1.]] Full specification and test list in **§2.16** | The degenerate answer is now unreachable through `eSVD()`, but `gamma_rate` stays exported-as-internal, so the unit-level behaviour still deserves one recording test. **[verified] on `x = rep(0, 100)`, `mu = rep(1e-8, 100)`, `s = rep(1, 100)`: `gamma_rate` returns `9.706e-4`; `log_gamma_rate` returns exactly `-10`, which is its clamp boundary.** `exp(-10) = 4.5e-5` — the two disagree by a factor of 21, and the second is saturated. Under the old behaviour that value fed straight into `.nuisance_in_sequence()`, which accepts any finite positive number without checking, and then into every posterior for that gene |
| T-CPP-GAM-04 | a gene with a single non-zero count behaves | [invariant] | same regime, less degenerate |
| T-CPP-GAM-05 | very large `mu` with tiny `s` does not silently return a bracket endpoint | [invariant] | **[regression]** §4.3: if the bracket loop exhausts `max_try` without satisfying `[l(ub)]'' <= 0`, `ub` is used anyway and Boost returns a *bound* rather than a root. The function returns it as if it were an MLE |
| T-CPP-GAM-06 | `log_gamma_rate`'s clamping is visible: when the true `log(β)` lies outside `[lower, upper] = [-10, 10]`, the function returns the boundary — assert that it does, and that the caller can tell | [inspection] the two early returns `return -lb` / `return -ub` | a saturated estimate is silently indistinguishable from a converged one, and `exp(10) ≈ 22026` then propagates into every posterior for that gene |
| ~~T-CPP-GAM-07~~ | **struck** — no convergence status will be added | **✅ RESOLVED (Q-GAM-2) — no** [[KZL: We don't need a convergence status for gamma_rate]] | Two consequences to accept knowingly, neither fatal now that §2.16 removes the input that provoked them. **(i)** T-CPP-GAM-05 and T-CPP-GAM-06 survive unchanged — a bracket endpoint and a clamp boundary are *observable in the return value* (they are exactly `lower`/`upper`), so the tests can assert them even though the caller cannot act on them. **(ii)** `.nuisance_in_sequence()`'s `log_gamma_rate` fallback stays permanently unreachable: it is triggered by a non-finite `gamma_rate` result, and a saturated result is finite. That is now a deliberate choice rather than an oversight — say so in a code comment so the next reader does not "fix" it |
| T-CPP-GAM-08 | `length(x) != length(mu)` or `!= length(s)` errors | [invariant] | the C++ reads `m_n = x.length()` and indexes `mu`/`s` with it — **a shorter `mu` is an out-of-bounds read**, not an error |

### 3.6 Pointer hygiene — `test_external_pointers.R` **[new file]**

`CRAN_READINESS.md` §4.1 predicts a **segfault** when a serialized `XPtr` is
dereferenced. **That prediction is wrong for the current Rcpp** — [verified]:
`saveRDS`/`readRDS` round-tripping either `esvd_family()` or `data_loader()`
output and then calling `feas_Xi_r()` or `objfn_all_r()` gives a clean R error,
`Error: external pointer is not valid`, from Rcpp's own checked `XPtr(SEXP)`
constructor. So the risk is smaller than the readiness doc states. The tests are
still worth writing, because the guarantee is Rcpp's and not ours, and because
the *message* is unhelpful.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-CPP-PTR-01 | `readRDS(saveRDS(esvd_family("poisson")))` then `feas_Xi_r()` gives an **informative** error naming the cause ("was this object restored from disk?"), not `external pointer is not valid` | [verified] the current message | locks in a guarantee we currently get by accident from a dependency |
| T-CPP-PTR-02 | the same for a restored `data_loader()` result through `objfn_all_r()` | [verified] same | |
| T-CPP-PTR-03 | `eSVD()`'s `intermediate_save` output can be `load()`ed and the pipeline resumed | [invariant] | this is the *practical* reason anyone would hit T-CPP-PTR-01/02. If the saved object contains no pointers (as inspection suggests), the test documents that fact so a future edit that adds one is caught |
| T-CPP-PTR-04 | `data_loader()` output survives a `gc()` and remains usable | [invariant] | `XPtr<DataLoader>(loader, true)` registers a finalizer; the `DenseDataLoader` holds an `Eigen::Ref` to the **R matrix's memory**, so if the R matrix is garbage-collected while the loader lives, the loader dangles. Construct a loader from a temporary, drop the R reference, `gc()`, then use the loader |

T-CPP-PTR-04 is the one I would prioritise: it is the only genuinely dangerous
pointer issue I can see, and unlike §4.1's it is not covered by Rcpp's checks.

---

## 4. Cross-cutting invariants (property-based)

These encode the paper's own claims. Each would catch a whole class of
regression that no single-function test can.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-PROP-01 | **Reparameterization preserves predictions.** `Ŷ X̂ᵀ + Ẑ Cᵀ` unchanged to 1e-8 across `reparameterization_esvd_covariates()` | [invariant] the paper states it twice | = T-REP-01, listed again because it should run on every fixture and after every fit stage, not just once |
| T-PROP-02 | **Orthogonality.** `X̂ᵀC ≈ 0` after Step 1; `X̂ᵀX̂/n` and `ŶᵀŶ/p` diagonal and equal after Step 2 | [invariant] | = T-REP-02/03 |
| T-PROP-03 | **Monotone loss** for every `opt_esvd` call and every family | [invariant] | = T-OPT-02 |
| T-PROP-04 | **Determinism** — `identical()` across two runs | [invariant] | = T-OPT-03 |
| T-PROP-05 | **Posterior sanity** — strictly positive, finite, and `var = mean/SplusBeta` exactly | [invariant] Eq. 14 | = T-POST-01/02 |
| T-PROP-06 | **Null calibration.** On `generate_null()`'s genuinely-null genes (11..p), p-values are approximately uniform: `ks.test(p_null, "punif")$p.value > 0.01` under a fixed seed | [oracle] the simulation truth | **this is the paper's Type-1-error claim as an executable test, and it is what would have caught §1.1.** It is also the single most expensive test in the plan and the most likely to be flaky |
| T-PROP-07 | **Recovery.** On `generate_data()` output with a known `nat_mat` and `nuisance_param_vec`, `estimate_nuisance()` recovers the truth within a stated tolerance | [oracle] the simulation truth | = T-NUIS-05 |
| T-PROP-08 | **Scale equivariance.** Multiplying every count by a constant `c` and adding `log(c)` to the `Log_UMI` covariate leaves `teststat_vec` unchanged | [invariant] the library-size model | **✅ RESOLVED (Q-PROP-1) — my judgement, per [[KZL: I don't have a strong opinion on this corner case of what happens with the clamps. Use your judgement and make sure it's well documented]].** Split into two tests, so the property and its boundary are both pinned rather than one being fudged into the other. **T-PROP-08a:** for `c ∈ {2, 10}` with `library_min` and `alpha_max` held at values that provably do not bind, `teststat_vec` is unchanged to `tolerance = 1e-6`. **T-PROP-08b:** for `c` large enough that `alpha_max` binds (`c = 1e4` on `F-TINY`), `teststat_vec` **does** change, and the test asserts that it changes — with a comment naming the clamp responsible. Rationale for asserting the failure rather than widening the tolerance until it passes: scale equivariance is a property of the *model*, and the clamps are a deliberate numerical safeguard that trades it away; a test that hid the trade would be asserting something the code does not do. **Documentation obligation:** `?compute_posterior` and `?eSVD` gain one sentence — "results are invariant to a global rescaling of the counts only while `library_min` and `alpha_max` do not bind; both are absolute thresholds, not relative ones" |
| T-PROP-09 | **Case/control label symmetry.** Swapping the case and control labels negates `teststat_vec` exactly and leaves `log10pvalue` unchanged | [invariant] | the test is two-sided; if it were not symmetric, direction of effect would bias significance |
| T-PROP-10 | **Gene permutation equivariance.** Permuting the columns of `dat` permutes every per-gene output identically | [invariant] | catches positional-vs-named indexing bugs anywhere in the stack, of which there are several candidate sites (`library_idx`, `case_control_idx`, `nuisance_vec` sweeps) |

**⚠ BLOCKED (Q-PROP-2) on T-PROP-06 — now resolved, see below.** A KS test with a fixed seed is a
one-sample check of a claim that is statistical, not deterministic. The options
are (a) fixed seed + loose threshold, accepting that it tests one draw; (b) a
handful of seeds with a "at least k of m pass" rule; (c) move it out of the CRAN
suite into a slow/CI-only tier. My recommendation is (c) — it is the most
valuable test in the plan and also the one most likely to make a CRAN check
flaky on a machine we do not control. Kevin should decide.

**✅ RESOLVED (Q-PROP-2) — option (a), fixed seed with a loose threshold**
[[KZL: Let's do (a) actually]]. So T-PROP-06 ships in the normal CRAN suite. Three
things it then has to do to earn that place, none optional:

1. **State the seed's role out loud** (convention C-04). The comment must say the
   seed was fixed *before* the threshold was chosen, not tuned until the test
   passed — otherwise the test asserts nothing.
2. **Threshold `> 0.01`, not `> 0.05`.** With ~110 null genes a KS p-value below
   0.01 under the true null happens once in a hundred runs; at 0.05 it happens
   once in twenty, and CRAN runs the suite on a dozen platforms.
3. **Budget it.** This is the most expensive test in the plan and C-08 caps the
   whole suite at 90 s. If `F-NULL` at 200 cells × 120 genes cannot fit, shrink
   `F-NULL` rather than the threshold.

Note the risk Kevin is accepting by choosing (a) over (c): a genuine
platform-specific numerical difference in `locfdr` will surface as a CRAN check
failure rather than a CI failure. That is a real possibility and the reason (c)
was recommended — but it is also the scenario in which we would *want* to know.

---

## 5. Input validation and error messages

The principle: **every exported function gets one test that malformed input
produces a useful error rather than a downstream `NA`.** `stopifnot()` is used
pervasively and its messages are unhelpful — `"what_obj is not found" is not TRUE`
is the worst example but not the only one.

These are cheap and mostly mechanical, so they are listed in one table.

| ID | Function | Malformed input | Expected |
|---|---|---|---|
| T-VAL-01 | `initialize_esvd` | `k > ncol(dat)` | error naming `k` and `ncol(dat)` |
| T-VAL-02 | `initialize_esvd` | unnamed `dat` | error or documented tolerance |
| T-VAL-03 | `initialize_esvd` | `covariates` lacking `"Intercept"` | error naming the column |
| T-VAL-04 | `initialize_esvd` | `metadata_individual` not a factor | error naming the argument |
| T-VAL-05 | `initialize_esvd` | `lambda` outside `[1e-4, 1e4]` | error naming the bound |
| T-VAL-06 | `initialize_esvd` | `"Intercept"` in `offset_variables` | error (there is a `stopifnot`) |
| T-VAL-07 | `opt_esvd.default` | `offset_variables` not in `colnames(covariates)` | error naming the missing ones |
| T-VAL-08 | `opt_esvd.default` | `nrow(x_init) != nrow(dat)` | error naming the mismatch |
| T-VAL-09 | `compute_posterior.default` | `colnames(covariates) != colnames(z_mat)` | error (there is a `stopifnot`) |
| T-VAL-10 | `compute_posterior.default` | negative `nuisance_vec` | error (there is a `stopifnot(all(nuisance_vec >= 0))`) |
| T-VAL-11 | `compute_posterior.eSVD` | `bool_adjust_covariates` and `bool_covariates_as_library` both `TRUE` | error |
| T-VAL-12 | `compute_test_statistic` | individual in both arms | error |
| T-VAL-13 | `compute_test_statistic` | one individual per arm | **error, not `NaN` df** |
| T-VAL-14 | `compute_test_statistic` | `length(individual_vec) != nrow(mat)` | error |
| T-VAL-15 | `compute_pvalue` | object without `teststat_vec` | error naming what is missing |
| T-VAL-16 | `compute_pvalue` | object without `dat` (the `bool_diet` case) | error explaining `.compute_df` needs it |
| T-VAL-17 | `format_covariates` | `nrow` mismatch | error |
| T-VAL-18 | `format_covariates` | single-level factor | error naming the variable |
| T-VAL-19 | `fisher_test` | `set1_genes ⊄ all_genes` | error |
| T-VAL-20 | `fisher_test` | non-character input | error |
| T-VAL-21 | `esvd_family` | unknown family string | error listing the 7 valid names |
| T-VAL-22 | `data_loader` | `dgeMatrix` / `character` matrix | **error at construction** — see T-CPP-LOAD-05 |
| T-VAL-23 | `gamma_rate` | length mismatch among `x`, `mu`, `s` | error — see T-CPP-GAM-08 |
| T-VAL-24 | `generate_data` | infeasible `nat_mat` for the family | error |
| T-VAL-25 | `.get_object` | unrecognized `what_obj` | error naming the key |
| T-VAL-26 | `eSVD` | `case_control_levels` not length 2 | error |
| T-VAL-27 | `eSVD` | duplicated entries in `categorical_vars` | error (there is a `stopifnot`) |
| T-VAL-28 | `reparameterization_esvd_covariates` | `fit_name` not in `names(input_obj)` | error naming the available fits |
| T-VAL-29 | `.compute_matrix_mean` | `mean_vec = 0.5` (a length-1 numeric) | error naming `mean_vec` and listing `TRUE`/`FALSE`/`NULL`/full-vector — see Q-SVD-2 |
| T-VAL-30 | `.compute_matrix_sd` | `sd_vec = 0.5` | same, for `sd_vec` |
| T-VAL-31 | `eSVD_helper` / `eSVD` | an individual with ≤ 2 cells | a `warning()` naming the dropped donors in `eSVD_helper` (✅ Q-COH-3); an **error** naming them in `eSVD` (✅ Q-COH-7). Both messages must name the donors — `esvd_helper.R` currently drops silently |
| T-VAL-32 | `eSVD` | every gene all-zero | error naming the analyzed-gene count (`0`), not a zero-column matrix into `initialize_esvd` — see T-GS-14 |
| T-VAL-33 | `eSVD` / `initialize_esvd` | `k` greater than the number of **analyzed** genes | error naming both `k` and the analyzed count, since a legal `k` can become illegal after §2.16's filter — see T-GS-15 |
| T-VAL-34 | `eSVD` | `SeuratObject` not installed | a `requireNamespace()` guard naming the package and how to install it, not `there is no package called 'SeuratObject'` — see §7. `esvd_helper.R` already models the right pattern for `eSVD2` itself |

**⚠ BLOCKED (Q-VAL-1) — now resolved, see below.** Writing these means committing to specific error
*messages*, which means either snapshot tests (C-10) or `expect_error(regexp=)`
with a fragment. Recommendation: `expect_error(..., regexp = "<key phrase>")` with
a short, stable fragment — snapshots for 28 messages is too brittle.

**✅ RESOLVED (Q-VAL-1) — `regexp`** [[KZL: Yes, let's use regexp]]. Formalized as
convention **C-11** in §0. Every row in this table is written as
`expect_error(<call>, regexp = "<fragment>")`, where the fragment is the
*identifying* part of the message — the argument name, the offending value, the
missing column — and never the surrounding sentence. Practical consequence:
writing these 33 tests means writing 33 error messages, since most of the
`stopifnot()` calls today produce no message worth matching on. That work is the
bulk of §5, not the assertions.

---

## 6. Verbose branches

Every `verbose` level of every function, smoke-tested under
`expect_no_error()` / `expect_no_warning()`. This is the cheapest section in the
plan and it catches §1.2, a hard crash sitting in the core optimizer's most
verbose mode that no current test would ever reach.

| ID | Target | Levels |
|---|---|---|
| T-VERB-01 | `opt_esvd.default` | 0, 1, 2, 3 — **level 2 currently throws** (§1.2) |
| T-VERB-02 | `initialize_esvd` | 0, 1, 2 |
| T-VERB-03 | `estimate_nuisance.default` | 0, 1, 2 |
| T-VERB-04 | `compute_test_statistic.default` | 0, 1 |
| T-VERB-05 | `compute_test_per_gene` | 0, 1, 2, 3, 4 |
| T-VERB-06 | `reparameterization_esvd_covariates` | 0, 1 |
| T-VERB-07 | `eSVD` | 0, 1, 2 |
| T-VERB-08 | `opt_x` / `opt_yz` | 0, 1, 2, 3 (level 3 enables the C++ Newton trace) |
| T-VERB-09 | `fisher_test` | 0, 1 — see T-UTIL-10, the branch currently does nothing |

Once §5.3 replaces `print()`/`cat()` with `message()` package-wide, these become
`expect_message()` tests with a fragment, which is strictly better: they then
assert that the diagnostic *says something*, not merely that it does not crash.

---

## 7. Regression tests, indexed against `CRAN_READINESS.md`

Each of these must **fail before the fix and pass after**. If one passes before
the fix, either the fix or the test is wrong. This table is the acceptance
criterion for §6 of the readiness doc.

| Defect | Test(s) | Must fail before fix? |
|---|---|---|
| §1.1 `qnorm(pt())` saturates to `+Inf` | T-PVAL-01, T-PVAL-02, T-MT-03, T-MT-04 | yes — [verified] `qnorm(pt(40,18)) == Inf` |
| §1.2 `opt_esvd(verbose = 2)` throws | T-OPT-09, T-VERB-01 | yes — [verified] |
| §1.3 non-syntactic covariate names / rank deficiency | T-REP-05, T-REP-06, T-FMT-09 | yes for name mangling; the collinearity half needs a fixture that actually goes rank-deficient |
| §1.4 `scale()` type change | T-POST-07 | yes — [verified] `scale()` returns a matrix and drops `names()` |
| §1.4 doc/code direction mismatch | T-POST-08 | **✅ Q-POST-1: the code is correct** — so this is a *documentation* defect, not a numerical one. T-POST-08 pins the code's direction and the roxygen is rewritten to match |
| §1.5 two-pipeline drift | T-PG-01..04, T-ESVD-02 | T-PG-03 yes (defaults differ); T-PG-01 no (they agree today when called correctly) — it is a *guard*, not a bug demonstration |
| §1.6 `.compute_df` recomputation | T-DF-02, T-DF-03 | no — it is correct today, just duplicated. These are guards for the refactor |
| §1.7 unconstrained `sigma0` | T-MT-05, T-MT-06 | no — a latent hazard, not a demonstrated failure. Test the guard, not a reproduction |
| §1.8 `stopifnot("what_obj is not found")` | T-UTIL-05 | yes — [verified] |
| §1.8 `.identification` negative eigenvalues | T-REP-07 | needs a fixture that produces one; may not be constructible |
| §1.8 `&` vs `&&` scalar guards | T-SVD-07 | no — behaviour is identical for scalars today; the test pins semantics before the change |
| §2.1 `sparseMatrixStats` | T-SVD-05 | no — write the test, then swap the implementation under it |
| §2.2 `Rmpfr` buys nothing | T-PVAL-07 | no — it documents the equivalence that licenses removal |
| §2.9 `MASS` undeclared | C-06 (remove it, don't declare it) | n/a |
| §4.1 external pointers | T-CPP-PTR-01..04 | **partially** — [verified] Rcpp already errors cleanly, contra the readiness doc's segfault prediction. T-CPP-PTR-04 (dangling `Eigen::Ref` after `gc()`) is the live risk |
| §4.2 `Rcpp::warning` under `warn = 2`; silent line-search failure | T-CPP-NEWT-02, T-CPP-NEWT-03 | yes for T-CPP-NEWT-02: a failed optimization currently reports convergence |
| §4.3 `gamma_rate` bracket search | T-CPP-GAM-05, T-CPP-GAM-07 | needs a fixture that exhausts `max_try`; may require constructing `(x, mu, s)` adversarially |
| **new** four families unusable at default `nuisance_vec` | T-OPT-05 | yes — [verified]. Not in the readiness doc; should be added to §1 |
| **new** `format_covariates` drops first level, roxygen says last | T-FMT-02 | yes — [verified]. Documentation defect |
| **new** `format_covariates` rescales only named variables, roxygen says all | T-FMT-06 | yes. Documentation defect |
| **new** `data_loader` returns a null pointer for `dgeMatrix` | T-CPP-LOAD-05, T-VAL-22 | yes — [verified] |
| **new** `gamma_rate` out-of-bounds read on short `mu`/`s` | T-CPP-GAM-08 | yes by inspection; verify under ASan |
| **new** all-zero gene makes `glmnet` warn and return the wrong lambda's fit | T-GS-01..16 | yes — [verified] `"Convergence for 1th lambda value not reached after maxit=100000"`. Not in the readiness doc; fixed by §2.16's filter rather than by a `glmnet` change |
| **new** `gamma_rate` and `log_gamma_rate` disagree 21× on an all-zero gene, the latter saturated at its clamp | T-CPP-GAM-03 | yes — [verified] `9.706e-4` vs `exp(-10) = 4.5e-5`. The most likely reason T-CPP-GAM-01 was commented out |
| **new** `Rmpfr` does not prevent the `pvalue_vec` underflow it appears to guard, and `p.adjust` is double-only | T-MPFR-01..08 | yes — [verified]. Strengthens `CRAN_READINESS.md` §2.2 from "buys nothing" to "cannot buy anything" |
| **new** `report_results()` reports `p = 0` for every gene with `log10pvalue > 308`, making extreme genes unrankable | T-MPFR-08 | yes — [verified] by construction. One-line fix: expose `log10pvalue` |
| **new** removing `exportPattern` would un-export `eSVD()` itself | T-ESVD-11 | yes — [verified] `NAMESPACE` has no `export(eSVD)`. A one-tag fix, but it has to land with the Q-CPP-1 change or the package's headline function disappears |
| **new** `SeuratObject` is used unguarded from `Suggests` | T-VAL-34 | yes — [verified] `eSVD()` calls `SeuratObject::LayerData()` at `R/eSVD.R:38`, `SeuratObject` is `Suggests`-only, and there is **not one `requireNamespace()` call anywhere in `R/`**. CRAN requires conditional use of a suggested package; this is a check failure waiting, and it is in neither working document |
| **new** `esvd_helper.R`'s `min_ids` does not prevent the `NaN` Welch df | T-COH-04 | yes by inspection — 4 case + 1 control passes a pooled `min_ids > 4` on 5 donors and still gives `n2 - 1 = 0`. See Q-COH-2 |

---

## 8. What this plan deliberately does *not* propose

Listed so the omissions are choices rather than oversights.

- **Performance/benchmark tests.** Explicitly out of scope per the current goal.
- **A C++ unit-test framework** (`testthat`'s `cpp11test`, Catch). Everything
  above is reachable from R through the existing bindings; adding a C++ test
  harness would add a `LinkingTo` dependency for no coverage gain.
- **Snapshot tests on numbers.** Locking in a fitted `x_mat` would break on every
  BLAS/LAPACK difference across CRAN's platforms. Numbers are tested against
  oracles and invariants only.
- **Vignette output tests.** The vignettes need external downloads (§2.7) and
  their fate is undecided.
- **Tests for `oldcode/`.** Not shipped.
- **`.rotate()` and `.l2norm()`** in `R/utils.R`. `.rotate()` appears to be dead
  code — grep finds no caller. **Recommend deleting it** rather than testing it;
  if it is kept, it needs one test.
- **Round-trip tests for `bernoulli` and `neg_binom`** `dat_to_nat`/`nat_to_canon`.
  These conversions are deliberately lossy (`bernoulli` maps to ±1 regardless of
  input). T-CPP-FAM-12 asserts this explicitly instead of pretending otherwise.

---

## 9. Open questions for Kevin, collected

### 9.1 Resolved in the 2026-08-29 review pass 1

Kept for the record, since each decision is now load-bearing for a test's
expected value. "Consequence" is what the answer costs beyond the test itself.

| # | Answer | Consequence |
|---|---|---|
| Q-INIT-1 | Zero `NA`s on the sparse path too | a fix, not just a test: move `dat[is.na(dat)] <- 0` out from under the `is.matrix()` guard |
| Q-SVD-1 | Target CRAN; absorb the one `sparseMatrixStats` function | opens **Q-SVD-3** (vendor vs. reimplement) — see §2.3 |
| Q-SVD-2 | Reject a length-1 numeric | guard must be `is.logical()`, not `length() == 1`; adds T-VAL-29/30 |
| Q-REP-1 | Warn and proceed, with a verbose-gated detail message | needs T-REP-09: `opt_esvd`'s `tryCatch` must stop swallowing the warning, or the decision has no effect |
| Q-NUIS-1 | ±30% relative on the rate | asserted gene-by-gene, on `F-SMALL` not `F-TINY` |
| Q-POST-1 | The code is right; the prose is wrong | a roxygen rewrite, plus the rate-vs-scale sentence added to `?compute_posterior` and `?estimate_nuisance` |
| Q-TSTAT-1 | Error early when an individual has ≤ 2 cells | the check goes in `eSVD()` **only** — putting it in `compute_test_statistic()` would delete T-TSTAT-01 and T-DF-01, the only external `t.test` oracle in the plan. Adds T-ESVD-09/10, T-VAL-31 |
| Q-PG-1 | `library_min = 0.1` everywhere | changes `compute_test_per_gene`'s default from `1e-2` |
| Q-CPP-1 | Bindings become internal; tests use `:::` | `exportPattern` in `NAMESPACE` has to go; each binding needs `@noRd` |
| Q-GAM-1 | Neither — filter the gene out instead; new `gene_status` output | a **new feature**, specified in §2.16, 16 new tests |
| Q-GAM-2 | No convergence status | T-CPP-GAM-05/06 survive (a boundary is observable in the return value); `.nuisance_in_sequence()`'s `log_gamma_rate` fallback stays unreachable *by decision* — comment it as such |
| Q-PROP-1 | My judgement, documented | split into T-PROP-08a (holds, moderate `c`) and T-PROP-08b (fails when `alpha_max` binds, and the test asserts the failure) |
| Q-PROP-2 | Option (a): fixed seed, loose threshold, in the CRAN suite | threshold `> 0.01` not `> 0.05`; the seed's role stated in a comment; counts against the 90 s budget |
| Q-VAL-1 / C-10 | `regexp` fragments; **no snapshot tier** | formalized as convention C-11. Strikes T-CPP-LOAD-06 — recommend deleting `test_data_loader()` instead |

### 9.2 Resolved in the 2026-08-29 review pass 2

| # | Answer | Consequence |
|---|---|---|
| Q-FIX-1 | `F-TINY` stays at 120 × 20, `k = 2`; not degenerate | recovery assertions (T-NUIS-05, T-PROP-07) still go on `F-SMALL` — a cell-count issue, not a degeneracy one |
| Q-SVD-3 | Try the four-line `Matrix` rewrite first; vendor the C++ only if T-SVD-05 fails | on the expected outcome, no third-party `cph` entry lands in `Authors@R` |
| Q-MPFR-1 | Run the `Rmpfr` tests outside the package | T-MPFR-01..04 never enter `tests/`; `Rmpfr` leaves `DESCRIPTION` entirely rather than moving to `Suggests`. The one-off script lives in `additional_context/` |
| Q-STATUS-1 | `initialize_esvd()` errors; `eSVD()` filters | one definition of "all zero", in one place |
| Q-STATUS-2 | Empirical null and BH on analyzed genes only | T-GS-09/10 enforce it |
| Q-STATUS-3 | `fdr_vec` set to 1 at status-2 positions | matches `log10pvalue = 0`; no `NA` reaches a threshold comparison |
| Q-STATUS-4 | Two levels only, `analyzed` / `all_zero` | with the enum frozen, factor-vs-integer is now only a readability preference |
| Q-TSTAT-2 | Import `esvd_helper.R` from `WAS2CODE_REPO` | brings three cohort filters this plan had not considered; specified as **§2.17**, and raises Q-COH-1..5 below |

### 9.3 Resolved in the 2026-08-29 review pass 3

| # | Answer | Consequence |
|---|---|---|
| Q-COH-1 | Keep `warning()` + `return(NA)`; no classed sentinel | the warning *text* becomes the only channel identifying which filter fired, so each rejection must name its own threshold (T-COH-05) |
| Q-COH-2 | Add `min_ids_per_arm`, default 2 | closes the `NaN` Welch df reachable through the front door (T-COH-04) |
| Q-COH-3 | **Drop low-cell donors, with a warning** | settles the pass-1/pass-2 conflict. T-ESVD-09 and T-VAL-31 can be written; the warning must *name* the donors (T-COH-02) |
| Q-COH-4 | Remove `bool_check_donors` entirely; all filters always run, thresholds are the knob | better than the split switch I proposed. Makes "set the threshold to 0" the only escape hatch, so that becomes a tested contract (T-COH-12) plus a `stopifnot` against the unsafe middle (T-COH-06) |
| Q-COH-5 | Keep `eSVD_helper()` as the wrapper **and** export `filter_cohort()` | two exported functions. Conveniently leaves every §2 test that calls `eSVD()` unaffected by the thresholds — but raises Q-COH-7 |

### 9.4 Resolved in the 2026-08-29 review pass 4

| # | Answer | Consequence |
|---|---|---|
| Q-COH-6 | **Drop first, then check** | every threshold describes the cohort actually analyzed; `min_ids_per_arm` is evaluated against a real count (T-COH-03) |
| Q-COH-7 | **All filtering, gene-status labelling, temporary removal and reinsertion live in `eSVD_helper()`. `eSVD()` errors on anything that gets through** | **supersedes Q-STATUS-1's placement.** `gene_status` moves out of `eSVD()`. Rewrites §2.16.1's four placement rows, re-points T-GS-05/09/10/13/16 at the helper, and adds T-GS-17, T-COH-13, T-COH-14 |

**No open questions remain.** Every decision needed to start writing is made.

### 9.5 Three things to confirm while implementing, not before

These do not block anything and each has a stated default in the text. Flagged
only so they are not discovered as surprises.

| # | Item | Default taken |
|---|---|---|
| 1 | What exactly `eSVD()` errors on — "any of these corner cases" | the three **correctness** conditions (all-zero gene, `≤ 2`-cell donor, one donor in an arm), **not** the four cohort *power* minima. A legitimately underpowered 15-cell run stays possible through `eSVD()` directly; those thresholds are the helper's policy, not the model's. T-COH-11 |
| 2 | `min_cells_per_id` of `1` or `2` | rejected at entry by a `stopifnot` — off (`0`) or safe (`>= 3`), never the middle, which weakens the drop rather than disabling it and puts a 1-cell donor back on the §1.1 path. T-COH-06 |
| 3 | Function name | `eSVD_helper` (Kevin's spelling), not the file's `esvd_helper` — it sits better beside `eSVD` and `eSVD2` |

**Two mechanical findings that need no decision, only a commit:**

- `NAMESPACE` has no `export(eSVD)` — the main user-facing function is exported
  only by `exportPattern("^[[:alpha:]]+")`, which Q-CPP-1 removes. Add the
  `@export` tag in the same commit (T-ESVD-11).
- `eSVD()` calls `SeuratObject::LayerData()` with `SeuratObject` in `Suggests`
  and **no `requireNamespace()` guard anywhere in `R/`**. CRAN requires
  conditional use of a suggested package (T-VAL-34).

---

## 10. Suggested order of work, once the list is agreed

Deliberately different from `CRAN_READINESS.md` §6, because tests come *before*
the fixes they guard.

0. **Harness** (§0): edition 3, `test_path()`, helpers, drop `context()` and
   `MASS`. No new assertions yet; the existing suite must still pass.
1. **Fixtures** (§1): regenerate `F-TINY`/`F-SMALL`, update
   `data_generation.R`, verify the tarball drops under 5 MB. This unblocks
   everything and also clears `CRAN_READINESS.md` §2.8.
2. **The two free equivalence tests** — T-PG-01..04 (matrix vs per-gene) and
   T-CPP-GAM-01 (`gamma_rate` vs `log_gamma_rate`). Both compare an existing
   implementation against an existing implementation, so they need no new
   oracle, and together they cover the posterior/test/p-value stack and the
   nuisance estimator. **These become the safety net for every later change.**
3. **Numerical-gradient tests** (§3.2). These make the C++ safe to touch, which
   is a precondition for §4.2 and §4.3's status-flag work.
4. **The §1.1 regression tests** (T-PVAL-01/02, T-MT-03/04) and the `log.p` fix.
   Two lines of implementation, and it is the defect that changes results.
   **Fold §2.15 in here** — T-MPFR-01..04 run once as the decision experiment,
   `Rmpfr` comes out of `Imports` in the same commit as the `log.p` fix, and
   T-MPFR-05/06/08 ship. The two are the same piece of work: both are about what
   happens in the extreme tail, and doing them together means touching the eight
   `pnorm` call sites once instead of twice.
5. **Cohort filtering** (§2.17) **then the `gene_status` feature** (§2.16), in
   that order and in one step, because two ordering constraints run between them:
   the donor drop must happen before `gene_status` is computed (T-COH-07), and
   the count checks must happen after the donor drop (T-COH-03, Q-COH-6). The
   cheapest way to guarantee both is to build them together — and under Q-COH-7
   they are literally the same function, steps 1–3 and 5 of §2.17.2. Build
   `eSVD_helper()` top to bottom in that order, with `eSVD()`'s three errors
   (T-COH-11) written *first*, so a mistake in the helper fails loudly rather
   than silently.
   The `gene_status` half comes deliberately after step 2, because
   T-GS-05 — "the run with all-zero genes equals the run without them" — is only
   meaningful once T-PG-01..04 have established that the two pipelines agree with
   each other in the first place. Order within the step: spec, then T-GS-01/02
   (the classifier), then the filter, then T-GS-03/04 (removal changes nothing it
   shouldn't), then the reinsertion and the rest. **T-GS-09/10 before shipping** —
   they are the two that are silent when wrong.
6. **Validation and verbose tiers** (§5, §6). Mechanical, high count, cheap, and
   §6 catches §1.2 immediately. Now 33 rows, not 28.
7. **The invariants** (§4), hardest last. T-PROP-06 is no longer gated (Q-PROP-2
   resolved) but is still the one to watch against the 90 s budget.
8. **Coverage floor in CI** (C-09).

---

## Appendix: coverage scorecard

Where the suite stands today, and what the plan above would add. "Functions" here
means functions defined under `R/` plus the `RcppExports` bindings.

| Area | Functions | Tested today | After this plan |
|---|---|---|---|
| `format_covariates.R` | 1 | 0 | 1 |
| `initialization.R` | 4 | 2 (partial) | 4 |
| `optimization.R` / `_helper.R` | 7 | 3 (Poisson only) | 7 (all 7 families) |
| `reparameterization.R` | 8 | 4 | 8 |
| `nuisance.R` | 4 | 2 | 4 |
| `posterior.R` | 3 | 1 (dims only) | 3 |
| `compute_test_statistic.R` | 6 | 5 | 6 |
| `compute_pvalue.R` | 2 | 0 | 2 |
| `multtest.R` | 4 | 0 | 4 |
| `compute_test_per_gene.R` | 1 | 1 (weak) | 1 (strong) |
| `eSVD.R` | 1 | 0 | 1 |
| `esvd_family.R` | 9 | 0 | 9 |
| `generate_data.R` / `generate_null.R` | 3 | 2 (2 of 7 families) | 3 (7 of 7) |
| `utils.R` / `data_management.R` / `data_loader.R` | 9 | 0 | 9 |
| `fisher_test.R` / `report_results.R` | 2 | 2 (weak oracles) | 2 (external oracles) |
| `src/` bindings | 16 | 4 | 16 |
| `gene_status` (new, §2.16) | — | — | 17 tests |
| tail precision (new, §2.15) | — | — | 8 tests (4 shipped, 4 run once) |
| cohort filtering (new, §2.17) | 2 (1 imported, 1 new) | 0 | 14 tests |

Roughly: **~40 `test_that` blocks today, ~242 proposed** (~200 in the original
draft, plus 17 for `gene_status`, 14 for cohort filtering, 4 shipped for tail
precision, 6 new validation rows and a handful of splits). The single largest
block is still §3.2's 12 assertions × 7 families, which is also the cheapest to
write per unit of risk removed.

Three of the additions are worth more than their line count suggests:
**T-GS-05** makes the reinsertion provably non-disturbing, **T-COH-04** catches a
`NaN` df that the imported filter looks like it prevents and does not, and
**T-GS-17** guards the one cost of Q-COH-7's three-level defence — "all zero" is
now decided in three places and must be decided by one predicate.

Worth noting what Q-COH-7 *removed*: T-GS-09/T-GS-10 dropped from load-bearing to
architecture checks, and `compute_test_per_gene` no longer needs to know about
`gene_status` at all. A better call graph paid for itself in tests.
