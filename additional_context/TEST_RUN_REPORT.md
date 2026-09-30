# eSVD2 test suite — first run against unmodified code

> **Dated record (2026-08-29 to 2026-09-01), not updated since.** Every
> finding here was resolved in versions 1.0.2 to 1.2.0. The current state of
> the suite is at the top of `UNIT_TEST_PLAN.md` (2014 expectations passing
> under 1.2.0), and the current state of the package is in §0 of
> `CRAN_READINESS.md`. Finding 1.1 (`gamma_rate` could not exceed the library
> size) was fixed in Part 6.1. Its consequence, a diverging rate, is what
> the cap of 1.2.0 addresses (`OVERDISPERSION_BRAINSTORM.md`).

**Run 2026-08-29**, R 4.5.1 / macOS (Darwin 23.2.0), `devtools::load_all()` on
the working tree (version 1.0.1.07). **No eSVD2 R or C++ code was changed**, as
instructed — every failure below is a statement about the package as it stands.

```
Rscript -e 'devtools::load_all("."); source("tests/testthat/helper-fixtures.R");
            testthat::test_dir("tests/testthat", filter = "_claude")'
```

## Headline numbers

| | Pass | Fail | Error | Skip |
|---|---|---|---|---|
| **New suite** (`test_*_claude.R`, 15 files) | **388** | **87** | 1 | 29 |
| **Existing suite** (12 files, untouched) | **79** | 0 | 0 | 0 |

The existing suite still passes, which was the point of writing the new files
alongside rather than editing them: the baseline is undisturbed and every
failure below is new information.

The 87 failures are **31 distinct assertions** — several are loops over genes or
families (T-NUIS-05 alone contributes 40, one per gene). The 29 skips are the
`gene_status` and cohort-filtering tests, which are written against functions
that do not exist yet.

## Where the files are

New files carry the `_claude` suffix per the repo convention, which also keeps
them from colliding with the existing `test_initialization.R`, `test_nuisance.R`
and so on.

| File | Plan section | Pass | Fail | Skip |
|---|---|---|---|---|
| `helper-fixtures.R` | §0 C-03, §1 | — | — | — |
| `test_family_derivatives_claude.R` | §3.2 | 53 | 0 | 0 |
| `test_utils_claude.R` | §2.14 | 99 | 2 | 0 |
| `test_validation_claude.R` | §2.2, §2.11, §3.1, §5, §6 | 29 | 10 | 1 |
| `test_gamma_rate_claude.R` | §3.5 | 28 | 5 | 0 |
| `test_svd_safe_claude.R` | §2.3 | 26 | 2 | 0 |
| `test_posterior_claude.R` | §2.7 | 24 | 0 | 0 |
| `test_invariants_claude.R` | §4, §2.4, §2.5 | 23 | 6 | 1 |
| `test_compute_pvalue_claude.R` | §2.9 | 22 | 9 | 0 |
| `test_multtest_claude.R` | §2.10 | 22 | 6 | 0 |
| `test_format_covariates_claude.R` | §2.1 | 22 | 1 | 0 |
| `test_tail_precision_claude.R` | §2.15 | 16 | 1 | 0 |
| `test_compute_test_statistic_claude.R` | §2.8 | 12 | 1 | 0 |
| `test_nuisance_claude.R` | §2.6 | 9 | 42 | 0 |
| `test_gene_status_claude.R` | §2.16 | 2 | 0 | 15 |
| `test_cohort_filter_claude.R` | §2.17 | 1 | 2+1err | 12 |

**A fixture note.** `helper-fixtures.R` builds F-TINY and F-SMALL in code with a
fixed seed rather than shipping `.rda` files. That is a departure from an
earlier draft of §1 and it buys three things: the tarball carries no fixture
bytes (CRAN_READINESS.md §2.8 is about a 10.1 MB `.RData`), the provenance of
every number is readable code, and a test needing a *fitted* object exercises
the pipeline instead of reading its own answer back. F-TINY does **not** carry
the two all-zero genes §1.1 asks for — with `gene_status` unimplemented they
would make every unrelated test fail for the same reason and drown the signal;
`.tiny_counts_with_zero_genes()` supplies them to §2.16 alone.

---

# Part 1 — NEW findings, in neither working document

Five of these are substantive. Two bear directly on the paper's Type-1-error
claim.

## 1.1 `gamma_rate()` cannot return a rate above the library size ⚠ serious

**`src/gamma_rate.cpp` line 181:**

```cpp
double ub = Rcpp::max(s);          // upper bound of the search bracket
```

`s` is the **library size**. The bracket-refinement loop that follows only ever
multiplies `ub` by 0.5 — it never grows it. So the estimate is capped by the
largest library size, however strong the evidence.

Measured on 4000 draws from the eSVD hierarchical model, `mu = 5`, `s = 1`:

| true β | `gamma_rate` | `exp(log_gamma_rate)` |
|---|---|---|
| 0.10 | 0.0985 | 0.0985 |
| 0.60 | 0.5975 | 0.5975 |
| 0.90 | 0.8825 | 0.8825 |
| 1.00 | 0.9432 | 0.9432 |
| **1.20** | **0.9996** | 1.1555 |
| **1.50** | **0.9995** | 1.7248 |
| **3.00** | **0.9995** | 2.9250 |
| **10.00** | **0.9996** | 10.2225 |

`gamma_rate` tracks the truth up to β = 1 and then pins at ≈0.9995 forever.
`log_gamma_rate` is correct throughout. Rescaling `s` from 1 to 10 moves the
answer, which is the proof that the ceiling is `max(s)` rather than a property
of the likelihood (T-CPP-GAM-05b). The R objective from the comment block of
`gamma_rate.cpp` confirms independently that the likelihood at β = 3 is strictly
higher than at the returned value (T-CPP-GAM-02b) — the cap is losing
likelihood, not finding a different optimum.

**Why it matters.** `estimate_nuisance()` defaults to `bool_use_log = FALSE`, so
it calls `gamma_rate`. `nuisance_vec` is the Gamma **rate** β = 1/γ, so any gene
whose over-dispersion γ is below 1 — any gene *better* behaved than γ = 1 — is
assigned β ≈ max(s) instead of its true, larger value. `nuisance_vec` enters the
posterior denominator directly, so posterior variances for those genes are too
large, their test statistics too small, and power is lost on exactly the
well-behaved genes.

This is a sharper form of what the plan anticipated as T-CPP-GAM-05 ("returns a
bound rather than a root"). It is also almost certainly **why T-CPP-GAM-01 was
commented out** in the existing `test_gamma_rate.R`: above the cap the two
routines cannot agree, because one of them is pinned.

- Fails: **T-NUIS-05** (all 40 genes), **T-CPP-GAM-05**, **T-CPP-GAM-05b**,
  **T-NUIS-08**, **T-CPP-GAM-01b** (passes, documenting the disagreement).
- Suggested fix: search on the log scale, or set `ub` from the data
  (e.g. a method-of-moments estimate) and let the bracket grow as well as shrink.

## 1.2 The test statistic omits the Bessel correction ⚠ bears on Type-1 error

`.compute_mixture_gaussian_variance()` computes

```r
mean(var) + mean(m^2) - mean(m)^2
```

which is the **population** variance of the mixture — it divides by *n*, with no
Bessel correction. That is then divided by `n1`/`n2` to form the Welch
denominator, and the **same** variances feed `.compute_df()`'s
Welch–Satterthwaite degrees of freedom.

Welch's *t* uses the **sample** variance, which divides by *n − 1*. On the
one-cell-per-individual fixture — where the posterior variances are zero and the
two should coincide exactly — they differ by exactly √(n/(n−1)) per arm:

```
teststat_eSVD = teststat_Welch × sqrt(n/(n-1))        (n1 = n2 = n)
```

verified to 1e-8 across every gene (T-TSTAT-01a passes, pinning the relation).

The statistic is therefore **inflated** — 9.5% at 6 individuals per arm, 5.4% at
10, 2.6% at 20 — and the resulting p-values are **anti-conservative**, while
being referred to a *t* distribution that assumes the sample variance.

Whether this is intended is a modelling question I cannot settle: the mixture
genuinely *is* a population object, so computing its population variance is
defensible; using that in a statistic compared against a *t* distribution is the
step that does not follow. **This needs your judgement.**

- Fails: **T-DF-01** (6 genes). Passes: **T-TSTAT-01a**, which documents the
  exact relationship.

## 1.3 `.multtest_simple()` underestimates the null sd by 21% ⚠

On 1000 draws from exactly N(0, 1):

| estimator | null mean | null sd | |
|---|---|---|---|
| `locfdr` | +0.0065 | 0.9866 | accurate |
| `truncated_mle` | +0.0671 | 1.1079 | sd 12% too **high** (conservative) |
| `simple` | +0.0104 | 0.7901 | sd 21% too **low** (anti-conservative) |

`.multtest_simple()` is the **last** fallback, and it is the one whose own source
comment disclaims it. Dividing by a null sd 21% too small inflates every z-score
by ~27%.

Combine with CRAN_READINESS.md §1.1 and the consequence is concrete: one strongly
DE gene produces an `Inf`, `locfdr` errors, the `tryCatch` swallows it,
`.multtest_truncatedGauss` also fails on the `Inf`, and the run silently lands on
`.multtest_simple` **for every gene in the dataset**. This quantifies the cost of
§1.1's silent degradation, which was previously described only qualitatively.

Separately: **at 200 genes the fallback fires on ordinary data**, because
`.multtest_locfdr()` catches *warnings* as well as errors and `locfdr` warns at
modest gene counts. It then returns null mean 0.70 and sd 2.57 on N(0,1) data.

- Fails: **T-MT-02**, **T-MT-03**. Passes: **T-MT-02a**, which pins the biases.

## 1.4 `format_covariates()` "rescaling" is not standardization

`scale(x, center = FALSE, scale = TRUE)` divides by the **root-mean-square**,
not the standard deviation — those coincide only for a centered column. For
`Age` with mean 40 and sd 8 the RMS is ≈40.8, so the rescaled column has sd
≈0.2, not 1.

The roxygen says the function "rescales all the numerical variables (but does
not center them)". A reader will take that to mean unit variance. It does not.

The property it *does* deliver is real and worth keeping: the result is
invariant to the unit the covariate was recorded in (years vs. months give
identical columns — T-FMT-05b passes). But the doc should say what it does.

- **T-FMT-05** pins the actual semantics and passes. No failure; recorded here
  because it is a documentation defect nobody had noticed.

## 1.5 Q-SVD-3 is answered by the test, and the answer is "neither option"

The decision rule from §2.3 was: try the four-line `Matrix` rewrite, vendor the
C++ only if T-SVD-05 fails at 1e-10.

**The naive rewrite fails, and badly.**
`sqrt((colSums(x²) − n·colMeans(x)²)/(n−1))` suffers catastrophic cancellation on
a column with a large mean and small variance, goes negative under the square
root, and returns **NaN** (T-SVD-05b, passes — it asserts the NaN).

But vendoring C++ is still unnecessary. **T-SVD-05c ships a stable `Matrix`-only
form that matches `matrixStats::colSds` to 1e-10** on the same pathological
fixture, by summing squared deviations over the stored non-zeros and correcting
for the implied zeros:

```r
sum_i (x_i - m)^2 = sum_{stored} (x_i - m)^2 + (n - nnz) * m^2
```

**Recommendation: ship T-SVD-05c's implementation.** No `cph` entry in
`Authors@R`, no new `src/` file, and the branch is dead code anyway.

## 1.6 Not root-caused: scale equivariance fails larger than expected

**T-PROP-08a** asserts that doubling every count and adding log(2) to `Log_UMI`
leaves `teststat_vec` unchanged, with `library_min` and `alpha_max` set so they
cannot bind. It fails by up to **38% relative**, all in the same direction:

```
original: 2.3878 1.4830 1.8988 2.7037
scaled:   3.1489 2.0528 2.5480 3.6261
```

**I have not root-caused this and am not claiming a defect.** The most likely
explanation is my test rather than the code: it reuses the `nuisance_vec`
estimated from the original data instead of re-fitting, so the rescaling is not
carried through the whole model. A conclusive version needs a full re-fit on the
scaled data, which is slower and which I have not written. Flagged so it is not
mistaken for a confirmed finding.

---

# Part 2 — Failures the plan predicted

These are the tests written specifically to pin a known defect. All fail, as
designed.

| Test | Defect | Source |
|---|---|---|
| **T-PVAL-01** | `qnorm(pt(40, 18))` is `+Inf`; needs `log.p = TRUE` | §1.1 |
| **T-PVAL-02** | `pvalue_list` has no `method` field, so the silent fallback is invisible | §1.1 |
| **T-MT-04** | a non-finite `teststat_vec` is absorbed rather than rejected | §1.1 |
| **T-MT-06** | `optim`'s `convergence` code is ignored entirely | §1.7 |
| **T-VERB-01** | `opt_esvd(verbose = 2)` throws | §1.2 |
| **T-OPT-05** | `gaussian`, `curved_gaussian`, `neg_binom`, `neg_binom2` all die at their documented defaults (`nuisance_vec = rep(NA, p)`) — 4 of 7 families | new in session 2 |
| **T-PG-03** | `library_min` defaults disagree: `compute_test_per_gene` 1e-2 vs `compute_posterior.eSVD` 0.1 | §1.5 |
| **T-INIT-05** | the sparse path does not zero `NA`s; `initialize_esvd` on a `dgCMatrix` with an `NA` **errors inside glmnet** | Q-INIT-1 |
| **T-SVD-06** | `mean_vec = 0.5` is silently coerced to `TRUE` | Q-SVD-2 |
| **T-REP-08/09** | `.reparameterize` on a rank-deficient input **errors inside `eigen()`** rather than warning and proceeding | Q-REP-1 |
| **T-NUIS-02** | nuisance failures are not countable | §2.6 |
| **T-TSTAT-07** | `.construct_averaging_matrix` accepts an empty index set | §2.8 |
| **T-CPP-LOAD-05** | `data_loader()` on a `dgeMatrix` or character matrix returns a null pointer without error | §4.1 |
| **T-CPP-GAM-08** | a short `mu`/`s` is an out-of-bounds read, not an error | §3.5 |
| **T-MPFR-08** | `report_results()` reports `p = 0` for every gene with `log10pvalue > 308` | §2.15 |
| **T-ESVD-11** | `NAMESPACE` has no `export(eSVD)` — only the blanket `exportPattern` | §7 |
| **T-VAL-34** | `SeuratObject` used from `Suggests` with no `requireNamespace()` guard anywhere in `R/` | §7 |
| **T-UTIL-05, T-VAL-15/16, T-VAL-21, T-VAL-28, T-DF-03, T-FMT-08** | `stopifnot()` messages name nothing the caller can act on | §5, §1.8 |
| **T-UTIL-10** | `fisher_test(verbose = 1)` computes a `paste0` and discards it | §2.14 |

## What passed, and is worth knowing

- **All 28 analytic gradient and Hessian checks pass, across all seven
  families** (`grad_Xi_r`, `hessian_Xi_r`, `grad_YZj_r`, `hessian_YZj_r` against
  `numDeriv`). This was the plan's largest single block and its highest-value
  C++ test. **The derivatives are correct.** `numDeriv` has now earned its place
  in `Suggests` after years of sitting unused.
- The l2 penalty enters the gradient as `2·l2pen·x / p`, and `objfn_all_r` is
  exactly the mean over cells of `objfn_Xi_r` — the normalizations are
  consistent (T-CPP-FAM-07, T-CPP-FAM-08).
- **All 24 posterior tests pass**, including `var = mean/denominator` exactly,
  the `bool_return_components` hook, both stabilization branches keeping
  `nuisance_vec`'s type and names, and every clamp binding as documented.
- **Q-POST-1 confirmed**: the stabilization fires when
  `mean(log10(nuisance_vec)) > 0`. The code is right and the roxygen is wrong,
  exactly as you said (T-POST-08 passes).
- The three loader implementations agree, including on an all-zero sparse column
  (T-CPP-LOAD-01, T-CPP-LOAD-03).
- Reparameterization preserves predictions and orthogonalizes `X` against `C`
  to 1e-6 (T-PROP-01, T-PROP-02, T-REP-03).

---

# Part 3 — Not yet implemented (29 skips)

`test_gene_status_claude.R` (15 skips) and `test_cohort_filter_claude.R`
(12 skips) are written **first**, against the §2.16 and §2.17 specifications, so
the features can be built to them. They skip cleanly on
`exists("eSVD_helper")` / `exists("filter_cohort")` rather than producing a wall
of "could not find function".

Four of their tests need no implementation and already run:

- **T-GS-03** ✓ removing all-zero genes leaves `Log_UMI` identical.
- **T-GS-04** ✓ removing them leaves `alpha_max = 2·max(dat)` identical.
- **T-COH-11** ✗ `eSVD()` does **not** error on a 1-cell donor or an all-zero
  gene — the Q-COH-7 guard is not there yet, as expected.
- **T-COH-13** ✗ my own test bug: `seurat_obj[genes, ]` errors with "the default
  assay will be removed"; `subset(seurat_obj, features = )` is the supported
  route. Worth knowing before the helper is written, since Q-COH-7 makes both
  filters Seurat *object* subsets.

The one **error** (not failure) is in `test_cohort_filter_claude.R`'s
`.cohort_seurat()` helper at very small cell counts, where Seurat refuses to
build the object. It only affects tests that are skipped anyway; I will fix it
when the feature lands.

---

# Part 4 — What I would do next

Ordered by value, not by effort.

1. **`gamma_rate`'s `max(s)` ceiling** (1.1). It silently distorts
   `nuisance_vec` for a whole class of genes and nothing else in the package
   would ever reveal it. `log_gamma_rate` already works, so a one-line change of
   `estimate_nuisance`'s default is a stopgap while the bracket is fixed.
2. **Decide on the Bessel correction** (1.2). This is yours to call, and it
   changes every reported p-value.
3. **The `log.p = TRUE` fix** (§1.1) together with the `Rmpfr` removal — see
   `RMPFR_REPORT.md`. Same call sites, same commit. Add `method` to
   `pvalue_list` at the same time so T-PVAL-02 can pass.
4. **`.multtest_simple`'s 21% bias** (1.3). Even with §1.1 fixed, the fallback
   fires on ordinary 200-gene datasets.
5. **`@export` on `eSVD()` and a `requireNamespace()` guard** — two lines, and
   the first is a landmine under the Q-CPP-1 change.
6. **The error messages** (§5). Mechanical, high count, and it is what makes the
   remaining 8 validation failures pass.
7. Then `gene_status` and the cohort filter, against the 27 tests already
   waiting for them.

**Nothing here has been fixed.** All 15 new test files, `helper-fixtures.R`, the
`Rmpfr` experiment script and `RMPFR_REPORT.md` are the only additions; `R/`,
`src/`, `DESCRIPTION` and `NAMESPACE` are untouched.

---

# Part 5 — Self-audit corrections (2026-08-29, same day)

I re-checked every test after the first write-up. Ten are defective, and two
statements in Parts 1–4 above are **wrong**. Corrections first.

## 5.1 Corrections to what I reported

**(a) T-PROP-06 never ran.** The null-calibration test — the paper's
Type-1-error claim, and the most valuable single test in the plan — carries
`skip_on_cran()`, which skips under a plain `Rscript` because `NOT_CRAN` is
unset. It was counted as a skip, not a pass. **Run with `NOT_CRAN=true` it
PASSES.** The headline 388/87/29 therefore understates coverage by the one test
that matters most. Either drop `skip_on_cran()` (Q-PROP-2 chose option (a), the
CRAN suite, so it should not be there at all) or always run with `NOT_CRAN` set.

**(b) The Q-REP-1 claim in Part 2 is wrong.** I wrote that `.reparameterize`
"errors inside `eigen()` rather than warning and proceeding". It does not.
`.reparameterize` on a duplicated column **already warns and proceeds**, which is
exactly the behaviour Q-REP-1 asked for:

```
WARNING: Detecting rank defficiency in reparameterization step
errored: FALSE
```

The real defect is one level up and is *sharper* than what I described:
`reparameterization_esvd_covariates()` on a genuinely finite rank-deficient fit
warns from `.identification`, **proceeds anyway, and then dies in `eigen()`**
with "infinite or missing values in 'x'". So `.identification` is producing
non-finite values and the warning is not preventing anything — which is precisely
the concern T-REP-07 was written for, confirmed at a real call site rather than a
synthetic one.

## 5.2 Finding 1 re-verified, and the characterization tightened

My first sweep passed mismatched `s` values (data generated at `s = 1`, then
estimated with `s = 10`), which is a misspecified model and not evidence of
anything. Redone with data generated **consistently** for each `s`, true rate
β = 3, μ = 5, n = 4000:

| max(s) | `gamma_rate` | R MLE | `exp(log_gamma_rate)` |
|---|---|---|---|
| 0.25 | **0.24989** | 3.29650 | 3.29649 |
| 0.50 | **0.49977** | 3.53434 | 3.53434 |
| 1.00 | **0.99952** | 3.04782 | 3.04781 |
| 2.00 | **1.99992** | 3.07958 | 3.07958 |
| 4.00 | 2.98141 | 2.98142 | 2.98141 |
| 10.00 | 2.96497 | 2.96496 | 2.96497 |

The finding is **confirmed and now exact**: `gamma_rate` returns
**`min(MLE, max(s))`** — it saturates at `max(s)` to four significant figures
whenever the MLE exceeds it, and is correct otherwise. `log_gamma_rate` equals
the independent R MLE in every row. Finding 1 stands as written.

## 5.3 Ten defective tests of my own

| Test | Defect |
|---|---|
| **T-CPP-GAM-05b** | Invalid: reuses data generated at `s = 1` while passing `s = 10`. Must regenerate per `s`. (The *finding* is fine — 5.2 proves it properly.) |
| **T-REP-08** | Name says "warns and still returns a valid fit"; I replaced `expect_warning` with a try/error check, so it passes **without checking the warning at all**. |
| **T-TSTAT-04** | Tests the wrong function. A one-individual arm does **not** give NaN in `compute_test_statistic` — it returns finite values. The NaN concern is real but lives in `.compute_df`: `compute_pvalue` on such a cohort dies with "missing values and NaN's not allowed if 'na.rm' is FALSE". Test should move. |
| **T-REP-07** | Passes through a vacuous `expect_true(TRUE)` branch when the call errors. Asserts nothing about the package. |
| **T-NUIS-06** | Second assertion is a tautology: `length(unique(x)) == length(unique(unique(x)))`, where `x` was built with `unique()`. Always true. |
| **T-MT-05** | Name promises a `sigma0 <= 0` guard test; the body asserts a fact about `stats::pnorm` and that `null_sd > 0` on clean data. Never exercises the guard. |
| **T-CPP-FAM-12** | Asserts `is.function(family_obj$feasibility)` — unrelated to the bernoulli lossiness its name claims. |
| **T-VAL-34, T-ESVD-11** | Read package source through `test_path("..", "..")`, which will not resolve under `R CMD check` (tests run from an installed copy). They will silently **skip** there rather than fail. |
| **T-COH-03, T-COH-04** | `expect_true(all(is.na(res)) || isTRUE(res$rejected))` — `all(is.na(NA))` is `TRUE`, so once the function returns `NA` for *any* reason these pass trivially, not for the intended one. |

None of these produce a false *finding* — the defects make tests weak or
mis-aimed, not wrong about the package. But six of them currently **pass**, and a
test that passes without constraining anything is worse than no test.

## 5.4 Confirmed sound

- **Finding 2 (Bessel correction) verified exactly at four sample sizes**:
  observed ratio to `stats::t.test` is `sqrt(n/(n-1))` to 1e-10 at n = 3, 5, 10
  and 25 individuals per arm. Not an artefact of one fixture.
- **The gradient tests are not vacuous.** Gradient magnitudes run 0.04–1.37 and
  Hessians 0.08–5.31 across the seven families, against a 1e-5 absolute
  tolerance. Worth noting that for `neg_binom` (|grad| = 0.043) and `bernoulli`
  (0.054) an absolute 1e-5 is a *relative* ~2e-4, which is looser than it looks;
  a relative tolerance would be tighter, though all seven pass comfortably today.

---

# Part 6 — After Kevin's decisions (2026-08-29)

**442 pass / 36 fail / 29 skip**, up from 388/87/29. The existing 79-test suite
is still green. Three code changes, all authorized.

| | Before | After |
|---|---|---|
| New suite | 388 pass / 87 fail | **442 pass / 36 fail** |
| Existing suite | 79 / 0 | 79 / 0 |

## 6.1 `gamma_rate`'s bracket — FIXED (`src/gamma_rate.cpp`)

Kevin: *"Can you let the bound grow if needed?"*

Four changes to the root search:

1. **Grow the upper bound** while `[l(ub)]' > 0`, so the root is bracketed from
   the right instead of being clamped at `max(s)`.
2. **Shrink the lower bound** until `[l(lb)]' > 0`, so it is bracketed from the
   left too.
3. **Start Newton from the geometric mean** of the bracket rather than a fixed
   `1.0`, which was only ever adequate because the bracket was capped.
4. **Raise the iteration cap** from 10 to 200. Boost bisects whenever a Newton
   step leaves the bracket, so this costs nothing on easy inputs.

The old retry loop that shrank `lb` on `[l(b*)]'' > 0` is gone; it was a
workaround for a bad bracket and, with a correct one, it pulled the search off
the root.

**Verified against an independent `stats::optimize` of the R objective** at 20
combinations — `max(s)` in {0.25, 1, 2, 4, 10} × true β in {0.3, 1, 3, 10}:

- Before: **18 of 20 missed**, the worst by 84%.
- After: **0 of 20 miss**; all agree to ≤ 6e-5 relative.

`gamma_rate` and `log_gamma_rate` now agree to **1e-12 or better** across
β = 0.1 … 10, where they previously differed by a factor of 20 above β = 1. That
is almost certainly why the equivalence half of the existing `test_gamma_rate.R`
was commented out — it is now writable, as T-CPP-GAM-01b.

## 6.2 The min-cells and per-arm guards — ADDED

Kevin: *"the correct behavior in all these functions regarding the p-value and
test statistics is that there should be an error."*

New internal `.check_cohort_is_testable()` in `R/compute_test_statistic.R`,
called from `compute_test_statistic.default()`, `.compute_df()` and
`compute_test_per_gene()` — all three, so the two pipelines cannot diverge on
exactly the inputs the guard exists for. Two conditions:

- **An arm with fewer than 2 individuals** → error. This is the genuine
  numerical failure: the Welch denominator contains `(v/n)²/(n−1)`, so `n = 1`
  gives **df = 0** and every `stats::pt()` returns NaN. It used to surface as
  *"missing values and NaN's not allowed"* from inside `multtest`, naming nothing.
- **An individual with fewer than `min_cells_per_individual` cells** (default 3)
  → error naming the individuals and their counts.

`min_cells_per_individual = 0` disables the second check, and is the documented
escape hatch. Worth recording that the second is a **policy** rather than a
numerical necessity: I verified that a 1-cell donor produces perfectly finite
statistics, even with zero posterior variances. The refusal states what the
method is willing to be asked, which is Kevin's call and now explicit.

**The external oracle survived.** T-TSTAT-01a needed one cell per individual,
which the new guard forbids. Rebuilt with **3 identical cells per donor and zero
posterior variance** — the reduction to Welch's t depends on the zero variance
and the identical cells, not the cell count. Still exact to 1e-10.

## 6.3 The Bessel correction — RESOLVED, not a defect

Kevin: *"Do not use a Bessel correction."* So the population variance is
deliberate and `.compute_mixture_gaussian_variance()` is correct. **T-DF-01 is
deleted** — asserting equality with `stats::t.test` was simply the wrong oracle.
T-TSTAT-01a keeps the contract in its true form:

```
teststat_eSVD = teststat_Welch * sqrt(n/(n-1))     per arm
```

verified exact at n = 3, 5, 10 and 25. Six failures gone.

## 6.4 T-PROP-08 (scale equivariance) — CUT

Kevin: *"Let's cut this."* The reason is recorded in the test file so it is not
proposed again: multiplying counts by *c* is not the same as sequencing *c*
times as deep. `Poisson(c·ℓ·λ)` has mean **and** variance `c·ℓ·μ`; multiplying
observed counts by *c* gives mean `c·ℓ·μ` but variance `c²·ℓ·μ`. The rescaled
data no longer satisfies the model's mean–variance relationship, so there is no
reason for the statistic to be preserved. A genuine Poisson invariance would be
*thinning*, which is stochastic and not worth its seed discipline here.

## 6.5 Two things the fix exposed

**(a) `estimate_nuisance` now computes the true MLE — and Q-NUIS-1's ±30% is a
statement about the estimator, not about our code.** Measured:

| cells/gene | median rel. err | 90th pct | fraction over 30% |
|---|---|---|---|
| 400 | 0.257 | 0.858 | 0.45 |
| 2000 | 0.107 | 0.256 | 0.07 |
| 10000 | 0.053 | 0.157 | 0.00 |

The Gamma rate is **weakly identified when over-dispersion is small** — the
likelihood is flat in β once counts look near-Poisson, so a single gene's MLE can
sit far from the truth while still being the MLE. Gene-by-gene ±30% needs
~10 000 cells per gene, which no fixture here can afford inside the 90 s budget.

T-NUIS-05 is therefore split: **T-NUIS-05** asserts our code reproduces an
independent `stats::optimize` of the R objective (exact, deterministic — this is
the assertion that would have caught the bracket cap), and **T-NUIS-05b** checks
the *median* relative error across genes, which is stable. Both pass.

**(b) `log_gamma_rate` is now the routine with the binding constraint.** Its
bounds are `log(β) ∈ [−10, 10]`, i.e. β ≤ 22026. On the fixture, one gene of
twenty is essentially Poisson (over-dispersion ≈ 2.6e−7, MLE ≈ 3.8e6) and the
log route returns exactly `exp(10) = 22026.47` for it while `gamma_rate` now
returns the MLE. That is a deliberate clamp, not a bug — but the two are **not
interchangeable at the extremes**, and the asymmetry has flipped. T-NUIS-08 now
asserts agreement off the clamp and that the log route is the one that saturates.

## 6.6 The ten defective tests — all repaired

T-CPP-GAM-05b (regenerates data per `s`), T-REP-08 (asserts the warning again —
`.reparameterize` **already** warns and proceeds correctly, so my earlier
"errors inside `eigen()`" claim was wrong at that level), T-TSTAT-04 (moved to
the real failure and split into 04/04b), T-REP-07 (no vacuous branch),
T-NUIS-06 (compares the two index sets instead of a tautology), T-MT-05 (renamed
to what it asserts, with the untestable gap stated), T-CPP-FAM-12 (asserts the
bernoulli lossiness), T-COH-03/04 (match the identifying warning rather than
`all(is.na(NA))`).

Still outstanding and unfixed: **T-VAL-34 / T-ESVD-11** read package source via
`test_path("..", "..")` and will silently skip under `R CMD check`. They need a
different mechanism.

## 6.7 The 36 remaining failures

All are the plan's predicted findings, unchanged in character from Part 2 —
`log.p`/§1.1 (3), the `multtest` estimators (4), error messages that name
nothing (8), `data_loader` on a `dgeMatrix`, `exportPattern`/`SeuratObject`,
`verbose = 2`, the un-run `gene_status` and cohort work (2 + 1 error), and
T-OPT-05's four families that die at their documented defaults. Nothing in
this list is new, and nothing is a test defect.

---

# Part 7 — After the code fixes (2026-09-01)

Kevin reviewed the suite and asked for the code to be fixed to it. Result:

| | Before (Part 6) | After |
|---|---|---|
| Whole suite | 522 pass / 37 fail / 28 skip | **640 pass / 0 fail / 0 skip** |

The 28 skips were the `gene_status` and cohort-filter features, now
implemented in `R/eSVD_helper_claude.R` (`filter_cohort()`, `eSVD_helper()`,
`.reinsert_genes()`) with `.which_all_zero()` in `R/utils.R` as the one shared
predicate. Everything else is in the file it was always in.

## 7.1 Code changes, by defect

| Defect (test) | Fix |
|---|---|
| §1.1 `qnorm(pt())` saturates (T-PVAL-01) | `.t_to_gaussian()` in `compute_pvalue.R`: log-scale composition mirrored through zero, used by both pipelines |
| `method` discarded (T-PVAL-02) | stored in `pvalue_list`; `multtest()` also **warns** whenever it falls below `locfdr` |
| non-finite input absorbed (T-MT-04) | `multtest()` errors at entry |
| `.multtest_simple` sd 21% low (T-MT-02) | moment-match the truncated sample to a normal truncated at the same null quantiles |
| `.multtest_truncatedGauss` sd 2.57 at 200 genes (T-MT-03) | **the defect was in the model, not the optimizer**: Efron's `theta` was a free parameter, dropping the constraint `p0 <= 1` that ties the window count to the null mass. Re-parameterized as `(delta0, log sigma0, p0)` under `L-BFGS-B` with `p0 in [1e-4, 1]`; `convergence` returned (T-MT-06) |
| `Rmpfr` | gone from `DESCRIPTION` and every call site; `sparseMatrixStats` replaced by `.sparse_col_sds()` (T-SVD-05c's form) |
| `report_results` ties at `p = 0` (T-MPFR-08) | a `log10pvalue` column; `pvalue` documented as underflowing |
| `verbose = 2` throws (T-VERB-01) | `print(paste0())` |
| four families die at defaults (T-OPT-05) | `nuisance_vec = NULL` means `rep(1, p)`, documented as a placeholder |
| sparse NAs not zeroed (T-INIT-05) | `dat@x[is.na(dat@x)] <- 0` |
| `library_min` defaults differ (T-PG-03) | `compute_test_per_gene` now `0.1` |
| `data_loader` null pointer (T-CPP-LOAD-05) | `Rcpp::stop` when no branch built a loader |
| `gamma_rate` out-of-bounds read (T-CPP-GAM-08) | length check in both C++ routines |
| rank-deficient fit dies in `eigen()` (T-REP-09) | `.identification` floors eigenvalues at `tol` and proceeds after its warning |
| nuisance failures uncountable (T-NUIS-02) | `.estimate_nuisance_matrix()` returns the count; stored as `param$nuisance_num_failed`; warns when positive |
| `length(x) == 1` coerces `0.5` to `TRUE` (T-SVD-06) | `is.logical()` guards naming the argument |
| empty index set (T-TSTAT-07) | error naming the position |
| error messages naming nothing (T-DF-03, T-VAL-15/16/21/28, T-UTIL-05, T-FMT-08) | `stop()` with the offending value |
| `eSVD` not exported, `SeuratObject` unguarded (T-ESVD-11, T-VAL-34) | full roxygen block with `@export`; `requireNamespace()`; `@exportPattern` removed from `zzz.R` |
| `eSVD()` refuses nothing (T-COH-11) | errors on all-zero genes, `k > ncol`, an individual in both arms, and `.check_cohort_is_testable()` before any fitting |
| §2.4, §2.5 mechanical blockers | `LICENSE` is the two-line stub; `override` on the five virtuals |

**`bool_diet = TRUE` now keeps the final fit.** `eSVD()` used to `NULL` all
three fits, including `fit_Second`, leaving `latest_Fit` dangling and the
object usable only by `report_results()`. T-GS-11 reads `x_mat` under the
default `bool_diet`, and §2.16.1's reinsertion spec pads `y_mat`/`z_mat`
rows — both presume the final fit survives. It is `n x k` plus `p x (k + r)`,
small next to what the diet removes. **Kevin's call to keep or revert**; the
revert is one line plus `bool_diet = FALSE` in T-GS-11.

## 7.2 Two tests changed because they could not pass, and one fixture

- **T-PVAL-01** recomputed the naive `qnorm(pt())` *inline* and asserted it
  was finite — it never called the package. Now calls `.t_to_gaussian()` and
  additionally pins `8.915293` and exactness where the naive form is finite.
- **T-PVAL-01b** asserted `fdr_vec` finite after passing an `Inf` to
  `multtest()`, while T-MT-04 (and the readiness doc's fix step 3) demand an
  error on that input. Its own comment said "rejected at entry"; it now
  asserts the error.
- **T-MT-02a** pinned the *buggy* estimator values ("so a change is caught")
  and directly contradicted T-MT-02. Rewritten to pin the corrected values.
  **T-MT-08** asserted narrower window → strictly smaller sd, which is the
  bias itself; now asserts both windows give ≈1 and differ.
- **T-MPFR-08** asserted `res$pvalue[1] != res$pvalue[2]` for `log10pvalue`
  400 vs 800 — impossible for any double column. Now asserts on the new
  `log10pvalue` column and that `pvalue` still ties.
- **T-COH-04** built "4 case / 1 control" by indexing `levels()` positionally;
  with ten donors the levels sort `indiv_1, indiv_10, indiv_2, ...`, so it got
  2 controls and 3 cases and the filter correctly accepted it. Donors are now
  named explicitly.
- **`helper-fixtures.R`**: gene names are `gene1` not `gene_1`, because
  `SeuratObject::CreateSeuratObject()` rewrites underscores to dashes, which
  would have failed every name-based `gene_status` assertion and was the
  actual cause of T-COH-13's `subset()` error.

## 7.3 New finding: the pipeline was not deterministic, and is sensitive

Two runs of `eSVD()` on identical input differed. The initialization differed
by 2e-9 (`irlba` starts from a random vector drawn from the user's RNG), and
by the end of the pipeline the test statistics differed by **up to 8** on
the strongest genes (6.70 vs 7.41; 13.99 vs 14.29) while agreeing to 1e-5 on
the rest. Fixed for reproducibility: `.svd_start_vector()` gives `irlba` and
`RSpectra` a deterministic start (a golden-ratio sequence through `qnorm`, no
RNG state touched); two runs now agree exactly, which is what made T-COH-08
and T-GS-05/09/10 pass.

**The amplification is not fixed by that, only hidden.** A 1e-9 perturbation
of the start moving a Welch statistic by 8 means the alternating optimization
plus reparameterization does not land on a well-defined point for the
strongly DE genes — most likely a flat direction between the latent
factors and the free covariate coefficients in the second fit. Worth a look
before submission; T-OPT-03's "deterministic" only ever tested `opt_esvd`
from a fixed start.

## 7.4 Code review of the diff (`/code-review high`, same day)

Ten candidates, verified by independent agents before the review stopped.

| Verdict | Finding | Status |
|---|---|---|
| CONFIRMED | `compute_test_per_gene` accepted `min_cells_per_individual` and never used it | fixed: calls `.check_cohort_is_testable()` |
| CONFIRMED | `eSVD_helper(min_cells_per_id = 0)` did not forward `min_cells_per_individual = 0`, so `eSVD()`'s default of 3 still errored | fixed: forwarded unless the caller passes it |
| CONFIRMED | `multtest()` on a tiny panel: `.multtest_simple` gives `null_sd = NA`, and the result was NA p-values behind a warning | fixed: errors naming the gene count |
| CONFIRMED | `eSVD()`'s unchanged `grep(id_var, colnames(covariates))` drops any covariate whose name contains `id_var` (`donor_age` for `id_var = "donor"`) | fixed: exact `<id_var>_<level>` names |
| PLAUSIBLE | `compute_test_statistic.eSVD` on a diet object fails in a `stopifnot` naming nothing | fixed: guard naming `dat` and `bool_diet` |
| CONFIRMED | T-SVD-05c defined a local `.sparse_col_sds` that shadowed the package's | fixed: test now calls the package function |
| CONFIRMED | **`gamma_rate` without the cap returns 1e4–6e7 for half the fixture's genes (true rates 2.3–7.6); `mean(log10) = 4.55`, so `bool_stabilize_underdispersion` rescales every gene by `10^-4.55`** | **open — Kevin** |
| CONFIRMED | `.multtest_locfdr` catches `locfdr`'s routine "f(z) misfit" warning, which fires on large heavy-tailed gene sets while `mlest` is still valid; real datasets may land on `truncated_mle` (pre-existing behaviour, now with a warning) | **open — Kevin**; T-MT-03 pins it |
| CONFIRMED | the case/control-individual derivation exists in four places with different error messages | open — refactor to one helper |
| CONFIRMED | `eSVD_helper` keeps a transposed count copy alive across `eSVD()`; `.reinsert_genes` re-allocates even when nothing was removed | open — efficiency, deferred by policy |

Suite after the six fixes: **644 pass / 0 fail / 0 skip**.
