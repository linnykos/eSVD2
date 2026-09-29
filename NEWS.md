# eSVD2 1.1.0

First CRAN submission.

## New features

* Every analysis now reports the log2 fold change between the case and the
  control individuals and a standard error for it. `report_results()` gains
  the column `logFC_se` beside `logFC`, and the `eSVD` object gains
  `log2fc_vec`, `log2fc_se_vec`, `case_var` and `control_var`, from both
  `compute_test_statistic()` and `compute_test_per_gene()` and therefore from
  `eSVD()` and `eSVD_helper()` under either setting of `bool_diet`.
* The standard error is the delta-method standard error of
  `log2(case_mean / control_mean)`, built from the two per-arm variances the
  Welch statistic already divides by. It is divided by the number of
  individuals and not the number of cells, which puts it on the scale and the
  unit of replication of the standard errors of 'DESeq2', 'dreamlet' and
  'NEBULA'. `?compute_log_fold_change` states what it accounts for and what
  it does not; in particular it is unreliable for a gene whose estimated
  nuisance rate is very large.
* New exported function `compute_log_fold_change()`, which recomputes the two
  vectors on an `eSVD` object. It needs neither the counts nor the posterior
  matrices, so it works on an object built with `bool_diet = TRUE`.

## Changes that can affect existing code

* `report_results()` returns six columns where it returned five. Code that
  reads its columns by name is unaffected; code that reads them by position
  is not.
* `compute_test_statistic()` on matrices returns a list of seven elements
  where it returned three; the first three are unchanged and in the same
  order. It warns, and returns `NA` for the fold change and its standard
  error, for a gene whose case or control mean is not positive. A posterior
  mean from `compute_posterior()` is always positive, so this affects only
  matrices supplied by the caller.
* An `eSVD` object saved by an earlier version stores the arm means but not
  the arm variances. `report_results()` still works on it, with `logFC_se`
  set to `NA` and a warning; `compute_log_fold_change()` refuses it. Rerun
  the test step, or `eSVD()` for an object built with `bool_diet = TRUE`.
* `compute_test_per_gene()` now records the individuals of each arm in
  `param`, as `compute_test_statistic()` does, and both replace that record
  when they are rerun. Before, a rerun on a changed cohort kept the
  individuals of the first run.

# eSVD2 1.0.2

Development version, not released. Relative to the GitHub release 1.0.0:

## Correctness

* `compute_pvalue()` and `compute_test_per_gene()` map Welch statistics to
  Gaussian statistics on the log scale, so a strongly up-regulated gene no
  longer produces an infinite statistic that silently pushed the empirical
  null onto its crudest fallback estimator.
* `multtest()` now reports which empirical-null estimator ran (`method`) and
  warns whenever `locfdr` could not be used. The truncated-Gaussian fallback
  is Efron's estimator with the `p0 <= 1` constraint restored; the simple
  moment estimator is corrected for truncation.
* The per-gene overdispersion MLE (`gamma_rate`) no longer caps its search
  at the largest library size; both `estimate_nuisance()` methods warn when
  a gene's estimate fails and falls to `min_val`.
* `initialize_esvd()` refuses all-zero genes, zeroes `NA` counts in sparse
  matrices too, and keeps the library-size offset in the single-covariate
  GLM fallback.
* `eSVD()` refuses individuals present in both arms, arms with fewer than
  two individuals, and individuals with fewer than
  `min_cells_per_individual` cells; the new `eSVD_helper()` and
  `filter_cohort()` apply the cohort filters, remove all-zero genes, and
  reinsert them with a `gene_status` record.
* The iterative SVD solvers start from a fixed vector, so identical inputs
  give identical fits without consuming the user's RNG stream.
* `opt_esvd()` warns when any Newton line search failed instead of raising
  the warning from inside C++.
* `bool_stabilize_underdispersion` keeps `nuisance_vec` a named numeric
  vector on both branches.

## Packaging

* `Rmpfr` and `sparseMatrixStats` are no longer dependencies.
* Only the intended API is exported; the raw Rcpp bindings are internal.
* `eSVD()` is documented.
* The two ASD tutorials are pkgdown articles rather than vignettes, since
  they need multi-gigabyte downloads.

# eSVD2 1.0.0

GitHub release accompanying Lin, Qiu and Roeder (2024),
<doi:10.1186/s12859-024-05724-7>.
