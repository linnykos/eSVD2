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

**Open question for Kevin (C-10):** should the suite gain a `tests/testthat/_snaps/`
tier at all? Snapshot tests are the cheapest way to lock in error *messages*
(§5 below has ~20 of them) but they make every intentional message edit a
two-step change. Recommendation: yes for error messages, no for numbers.

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

**⚠ BLOCKED (Q-FIX-1).** Building the fitted object at test time costs runtime
(`initialize_esvd` runs a glmnet Poisson ridge per gene). At 20 genes that is
cheap; the question is whether `F-TINY` is *large enough to be a meaningful test*
of `opt_esvd` convergence. My recommendation is yes for structural assertions and
no for the recovery assertions (T-PROP-07), which should use `F-SMALL`. Kevin
should confirm 120 × 20 is not so small that the eSVD factorization is degenerate
at `k = 2`.

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
| T-INIT-05 | `dat` with `NA` entries: `is.matrix(dat)` path zeroes them; the `dgCMatrix` path does **not** | [inspection] `dat[is.na(dat)] <- 0` is guarded by `is.matrix(dat)` | an asymmetry between dense and sparse input that no test covers. **⚠ BLOCKED (Q-INIT-1): is the sparse path supposed to zero NAs too?** The C++ loader has a whole `Flag::na` machinery, so the intent may be to keep them |
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
| T-SVD-05 | `.compute_matrix_sd(mat, sd_vec = TRUE)` on a `dgCMatrix` equals `matrixStats::colSds(as.matrix(mat))` | [oracle] dense recomputation | this is the *only* `sparseMatrixStats` call site (§2.1); the test is what makes it safe to replace with `Matrix` primitives. **⚠ BLOCKED (Q-SVD-1) on the CRAN-vs-Bioconductor decision** — but write the test first, then swap the implementation under it |
| T-SVD-06 | `.compute_matrix_mean(mat, mean_vec = 0.5)` — a length-1 *numeric* — is silently treated as `TRUE` and returns `colMeans` | [inspection] `if(mean_vec)` coerces `0.5` to `TRUE` | a caller passing a scalar centering constant gets column centering instead. Either document that only `TRUE`/`FALSE`/`NULL`/full-vector are accepted, or reject scalars. **⚠ BLOCKED (Q-SVD-2): which?** |
| T-SVD-07 | `check_stability = TRUE, K = 3` does not run the stability check (guard is `K > 5`), and `K = 10` does | [inspection] | the guard uses `&` on scalars (§1.8); the test pins the intended semantics before `&&` is substituted |

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
| T-REP-08 | `.reparameterize` on rank-deficient `x_mat` (duplicated column) — assert the documented behaviour, whatever we decide it is | **⚠ BLOCKED (Q-REP-1): should this warn, error, or silently proceed?** | `opt_esvd.default` wraps it in a `tryCatch` that falls back to the unreparameterized matrices *silently* — so a rank-deficient fit is currently indistinguishable from a good one |

### 2.6 `estimate_nuisance()` — `test_nuisance.R` **[strengthen]**

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-NUIS-01 | `bool_use_log = TRUE` and `FALSE` agree to ~1e-3 on a well-conditioned fixture | [oracle] the two are each other's oracle — `gamma_rate` and `exp(log_gamma_rate)` estimate the same quantity by different routes | a free equivalence test, same idea as §1.5. Currently the `bool_use_log = TRUE` path has zero coverage |
| T-NUIS-02 | `.nuisance_in_sequence()` failures are **countable**: a helper returns how many genes fell through to `0` | [invariant] | today the function returns `0` and warns only when `verbose > 0`; the `0` is then clamped to `min_val`, so a silent total failure across all genes looks identical to a successful fit |
| T-NUIS-03 | the returned vector carries `colnames(input_obj)` as names | [invariant] | `compute_posterior` sweeps by position and `report_results` names by gene; a names/position mismatch here mislabels every result |
| T-NUIS-04 | all values `>= min_val` and finite, including on the `mean_mat` containing `Inf` fixture | [invariant] | exists; keep |
| T-NUIS-05 | on `generate_data()` output with a known `nuisance_param_vec`, the estimate recovers the truth within a stated tolerance | [oracle] the simulation truth | the only test that says the estimator is *right* rather than merely finite. **⚠ BLOCKED (Q-NUIS-1): what tolerance is acceptable?** At `n = 120` cells per gene I would expect roughly ±30% on the rate — Kevin should set the number he is willing to defend |
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
| T-POST-08 | the stabilization fires in the documented direction | **⚠ BLOCKED (Q-POST-1)** | §1.4: the roxygen says "when the global mean over-dispersion is less than 1"; the code fires when `mean(log10(nuisance_vec)) > 0`. Because `nuisance_vec` is the *rate*, the code may be right and the prose wrong. **Only Kevin can say which was intended, and the test's expected value depends entirely on the answer** |
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
| T-TSTAT-03 | an individual with exactly one cell: the within-individual variance is 0 and the result is still finite (or errors informatively) | **⚠ BLOCKED (Q-TSTAT-1): finite-but-huge, or error?** | with one cell the mixture variance can go to 0 and `t = ±Inf`, which is precisely the input that breaks `locfdr` (§1.1) |
| T-TSTAT-04 | a group with exactly **one** individual gives a real error, not a `NaN` df | [invariant] | `n1 - 1 = 0` makes the Welch denominator `0/0`. `CRAN_READINESS.md` §3.5 already flags this: *"deserves a real error message, not a `NaN`"* |
| T-TSTAT-05 | `case_mean`/`control_mean`/`teststat_vec` all carry `colnames(posterior_mean_mat)` | [invariant] | `report_results` builds its `genes` column from `names(teststat_vec)` |
| T-TSTAT-06 | permuting the row order of `posterior_mean_mat` and `individual_vec` together leaves `teststat_vec` unchanged | [invariant] exchangeability | catches any accidental positional (rather than name-based) indexing in the averaging matrix |
| T-TSTAT-07 | `.construct_averaging_matrix` with an `idx_list` entry of length 0 errors, rather than producing an all-zero row that silently averages to 0 | [invariant] | reachable via an unused factor level; the resulting all-zero row biases the group mean with no signal |

### 2.9 `.compute_df()` and `compute_pvalue()` — `test_compute_pvalue.R` **[new file]**

`compute_pvalue` produces the package's headline output and has **zero** tests.

| ID | Assert | Oracle | Why |
|---|---|---|---|
| T-DF-01 | on the one-cell-per-individual fixture, `.compute_df()` equals `stats::t.test(...)$parameter` | [oracle] `stats::t.test` | same reduction as T-TSTAT-01; an external oracle for Welch–Satterthwaite |
| T-DF-02 | `.compute_df()` and `compute_test_statistic()` are computed from the **same** group variances | [invariant] | §1.6: `.compute_df()` re-derives `case_individuals`, the averaging matrix, and both group variances from scratch. Two copies of the same code drift; this test is what makes the refactor safe |
| T-DF-03 | `.compute_df()` errors informatively when `input_obj$dat` is absent | [invariant] | §1.6: the recomputation is why `compute_pvalue` needs `dat`, which conflicts with `bool_diet = TRUE` in `eSVD()`. The error should say so |
| T-PVAL-01 | with `teststat_vec` containing `t = 40, df = 18`, `gaussian_teststat` is **finite** | [verified] `qnorm(pt(40, 18))` is `Inf`; `qnorm(pt(-40, 18, log.p = TRUE), log.p = TRUE)` is `-8.915293` | **[regression]** §1.1. This is the highest-value single assertion in the whole plan: it is the defect that changes published-style results |
| T-PVAL-02 | with that same input, `pvalue_list$method == "locfdr"` | [invariant] | §1.1's real damage: one `Inf` makes `locfdr` error, the `tryCatch` swallows it, and the empirical null silently degrades to `.multtest_simple()` — the estimator whose own comment disclaims it. **Requires `method` to be stored in `pvalue_list`, which it is not today** |
| T-PVAL-03 | `log10pvalue` is monotone decreasing in `|gaussian_teststat - null_mean|` | [invariant] | a p-value that is not monotone in the statistic is definitionally broken; cheap and catches sign errors in the mirroring branch |
| T-PVAL-04 | `10^(-log10pvalue)` lies in `[0, 1]` for every gene | [invariant] | the mirroring construction `null_mean - (x - null_mean)` guarantees `2*pnorm(...) <= 1`; assert it |
| T-PVAL-05 | `fdr_vec == stats::p.adjust(pvalue_vec, "BH")` | [oracle] `p.adjust` | exists inside `test_report_results.R`; move it here where it belongs |
| T-PVAL-06 | `compute_pvalue`'s internally recomputed `log10pvalue_vec` equals `multtest()`'s `logpvalue_vec` | [invariant] | the same 6-line computation is written twice, once in `multtest.R` and once in `compute_pvalue.R` (and a third time in `compute_test_per_gene.R`). Either they agree or one is wrong |
| T-PVAL-07 | `Rmpfr::pnorm` and `stats::pnorm` are `identical()` on the inputs actually used | [verified] `Rmpfr::pnorm(-40, log.p=TRUE)` returns a base `numeric` equal to `stats::pnorm(-40, log.p=TRUE)` | §2.2: this test is what licenses dropping a GMP/MPFR system dependency. Write it, then delete it with the dependency |

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
| T-PG-03 | the two paths agree when both are called with **default** arguments | [inspection] they currently do not | **[regression]** §1.5: `compute_posterior.eSVD` defaults `library_min = 0.1`, `compute_test_per_gene` defaults `1e-2`. `eSVD()` passes it explicitly, which is why this has gone unnoticed. **⚠ BLOCKED (Q-PG-1): align on 0.1 (matching `compute_posterior`) — confirm** |
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

---

## 3. C++ backend tests

Restricting to correctness. Everything here is driven from R through the
`RcppExports` bindings; no C++ test framework is proposed.

**⚠ BLOCKED (Q-CPP-1).** §2.3 of the readiness doc proposes *un-exporting*
`objfn_Xi_r`, `grad_Xi_r`, `hessian_Xi_r`, `feas_*_r`, `data_loader`,
`esvd_family`, `gamma_rate`, `log_gamma_rate` — but the entire C++ test strategy
below calls them. The resolution is that they become internal (`:::`) rather than
exported, and the tests use `:::`. This is normal for a package's own test suite
and does not weaken the tests. Confirm before we start.

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
| T-CPP-GAM-01 | `gamma_rate(x, mu, s)` and `exp(log_gamma_rate(x, mu, s))` agree to ~1e-4 across a grid of `(mu, s)` regimes | [oracle] each other — two routes to the same MLE | **the commented-out half of the existing test.** This is a second free equivalence test, exactly like §1.5's. It is commented out, which suggests it was failing; finding out why is a prerequisite |
| T-CPP-GAM-02 | the estimate maximizes the likelihood: `objfn(b̂) >= objfn(b)` for `b` on a grid around `b̂` | [oracle] the R implementation of `objfn` that already exists **in the comment block of `src/gamma_rate.cpp` lines 58–75** | the comment ships a complete R reference implementation of the objective, gradient and Hessian. Turning it into a test is nearly free and gives an independent oracle for both C++ routines |
| T-CPP-GAM-03 | an all-zero count vector for a gene returns something finite and documented | **⚠ BLOCKED (Q-GAM-1): what is the right answer?** | the MLE is degenerate; the code will return *something*. `.nuisance_in_sequence()` accepts any finite positive value without checking |
| T-CPP-GAM-04 | a gene with a single non-zero count behaves | [invariant] | same regime, less degenerate |
| T-CPP-GAM-05 | very large `mu` with tiny `s` does not silently return a bracket endpoint | [invariant] | **[regression]** §4.3: if the bracket loop exhausts `max_try` without satisfying `[l(ub)]'' <= 0`, `ub` is used anyway and Boost returns a *bound* rather than a root. The function returns it as if it were an MLE |
| T-CPP-GAM-06 | `log_gamma_rate`'s clamping is visible: when the true `log(β)` lies outside `[lower, upper] = [-10, 10]`, the function returns the boundary — assert that it does, and that the caller can tell | [inspection] the two early returns `return -lb` / `return -ub` | a saturated estimate is silently indistinguishable from a converged one, and `exp(10) ≈ 22026` then propagates into every posterior for that gene |
| T-CPP-GAM-07 | `gamma_rate` and `log_gamma_rate` return a convergence status | **⚠ BLOCKED (Q-GAM-2)** | §4.3's fix. `.nuisance_in_sequence()` already has a `log_gamma_rate` fallback — it just never learns it is needed. The test cannot be written until the status exists |
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
| T-PROP-08 | **Scale equivariance.** Multiplying every count by a constant `c` and adding `log(c)` to the `Log_UMI` covariate leaves `teststat_vec` unchanged | [invariant] the library-size model | a strong, cheap end-to-end property that no individual function test implies. **⚠ BLOCKED (Q-PROP-1): is this actually a property of the model as implemented?** The `library_min` clamp and `alpha_max` clamp both break it for extreme `c`. I believe it holds for moderate `c`; Kevin should confirm the model intends it |
| T-PROP-09 | **Case/control label symmetry.** Swapping the case and control labels negates `teststat_vec` exactly and leaves `log10pvalue` unchanged | [invariant] | the test is two-sided; if it were not symmetric, direction of effect would bias significance |
| T-PROP-10 | **Gene permutation equivariance.** Permuting the columns of `dat` permutes every per-gene output identically | [invariant] | catches positional-vs-named indexing bugs anywhere in the stack, of which there are several candidate sites (`library_idx`, `case_control_idx`, `nuisance_vec` sweeps) |

**⚠ BLOCKED (Q-PROP-2) on T-PROP-06.** A KS test with a fixed seed is a
one-sample check of a claim that is statistical, not deterministic. The options
are (a) fixed seed + loose threshold, accepting that it tests one draw; (b) a
handful of seeds with a "at least k of m pass" rule; (c) move it out of the CRAN
suite into a slow/CI-only tier. My recommendation is (c) — it is the most
valuable test in the plan and also the one most likely to make a CRAN check
flaky on a machine we do not control. Kevin should decide.

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

**⚠ BLOCKED (Q-VAL-1).** Writing these means committing to specific error
*messages*, which means either snapshot tests (C-10) or `expect_error(regexp=)`
with a fragment. Recommendation: `expect_error(..., regexp = "<key phrase>")` with
a short, stable fragment — snapshots for 28 messages is too brittle.

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
| §1.4 doc/code direction mismatch | T-POST-08 | **⚠ BLOCKED (Q-POST-1)** |
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

Every **⚠ BLOCKED** tag above, in one place. Each blocks at least one test.

| # | Question | Blocks |
|---|---|---|
| Q-FIX-1 | Is 120 cells × 20 genes at `k = 2` a non-degenerate eSVD problem? Should `F-TINY` be bigger? | all of §2's fixtures |
| Q-INIT-1 | Should the sparse `dat` path zero `NA`s the way the dense path does, or does the C++ `Flag::na` machinery mean `NA`s are meant to survive? | T-INIT-05 |
| Q-SVD-1 | CRAN or Bioconductor — i.e. drop `sparseMatrixStats` or keep it? (`CRAN_READINESS.md` §2.1) | T-SVD-05, and most of the packaging work |
| Q-SVD-2 | Should `.compute_matrix_mean`/`.compute_matrix_sd` reject a length-1 *numeric* (currently silently coerced to `TRUE`)? | T-SVD-06 |
| Q-REP-1 | Should `.reparameterize` on a rank-deficient input warn, error, or proceed? `opt_esvd` currently swallows the failure in a `tryCatch` and silently keeps the unreparameterized matrices | T-REP-08 |
| Q-NUIS-1 | What recovery tolerance for `estimate_nuisance` is defensible at `n ≈ 120` cells/gene? | T-NUIS-05, T-PROP-07 |
| Q-POST-1 | Which direction is `bool_stabilize_underdispersion` meant to fire? The roxygen and the code disagree, and because `nuisance_vec` is the *rate*, the code may well be right | T-POST-08 |
| Q-TSTAT-1 | An individual with exactly one cell: finite-but-huge statistic, or an error? | T-TSTAT-03 |
| Q-PG-1 | Align `library_min` on 0.1 (matching `compute_posterior`)? | T-PG-03 |
| Q-CPP-1 | Confirm the raw Rcpp bindings become internal (`:::`) rather than disappearing, so the C++ tests can call them | all of §3 |
| Q-GAM-1 | What should `gamma_rate` return for an all-zero gene? | T-CPP-GAM-03 |
| Q-GAM-2 | Add a convergence status to `gamma_rate`/`log_gamma_rate`? | T-CPP-GAM-07 |
| Q-PROP-1 | Is scale equivariance (T-PROP-08) an intended property of the model as implemented? | T-PROP-08 |
| Q-PROP-2 | Where does the null-calibration test live — CRAN suite, or a slow/CI-only tier? | T-PROP-06 |
| Q-VAL-1 | `expect_error(regexp = ...)` fragments vs snapshots for the 28 validation messages? | all of §5 |
| C-10 | Snapshot tier at all? | §5, T-CPP-LOAD-06 |

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
5. **Validation and verbose tiers** (§5, §6). Mechanical, high count, cheap, and
   §6 catches §1.2 immediately.
6. **The invariants** (§4), hardest last, with T-PROP-06 gated on Q-PROP-2.
7. **Coverage floor in CI** (C-09).

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

Roughly: **~40 `test_that` blocks today, ~200 proposed.** The single largest
block is §3.2's 12 assertions × 7 families, which is also the cheapest to write
per unit of risk removed.
