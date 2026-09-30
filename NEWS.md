# eSVD2 1.2.0

## The nuisance rate is capped

* `estimate_nuisance()` bounds each gene's nuisance rate (the Gamma rate, the
  reciprocal of the over-dispersion) at `cap_multiplier` times the median over
  cells of the gene's library size, `min(MLE, cap_multiplier * median library
  size)`, with `cap_multiplier = 10` by default. The bound is applied in R
  after the maximum-likelihood estimate; `gamma_rate()` is unchanged.
* A gene at the boundary (no finite maximum-likelihood rate) is set to the
  bound itself, and not to the smaller of the bound and the value at which
  the iterations of the estimate stopped. That value is arbitrary (about 1e7
  on one route, `exp(10)` on the other) and is below the bound when the
  library size is in the thousands, as it is when the library excludes the
  intercept.
* The floor `min_val` is a multiple of the median library size too, so that
  no gene's rate is above its bound: the rate is `max(min(MLE, cap_multiplier
  * m), min_val * m)`, with `m` the gene's median library size, and `min_val`
  must be below `cap_multiplier`. In 1.1.0 the floor was an absolute number;
  it held a gene whose fit was degenerate (median library size below
  `min_val / cap_multiplier`) above its bound. A gene whose estimation fails
  gets the floor.
* Why: a gene whose counts are no more variable around the fit than a Poisson
  has no finite maximum-likelihood rate. In 1.1.0 such a gene received a rate
  near 1e7 (or exactly `exp(10)`, when the first estimation route did not
  return), its posterior followed the fit, and its test statistic was
  inflated. The bound is a calibration device, and the documentation says so.
* **Results change.** The statistics, p-values, fold changes and standard
  errors of every analysis differ from 1.1.0, for all genes and not only the
  capped ones: the posterior divides every rate by the geometric mean of the
  rates when that mean is above 1 (`bool_stabilize_underdispersion`), and
  the empirical null is fitted to all genes. `cap_multiplier = Inf` gives
  the rates of 1.1.0.
  `cap_multiplier = 1` is close to, but not the same as, the version that
  accompanied Lin, Qiu and Roeder (2024), whose bound was the largest library
  size of the gene.
* `eSVD()` has the new argument `cap_multiplier`, and `eSVD_helper()` passes
  it on.

## New features

* Every gene has a status, stored as the factor `nuisance_status` beside
  `nuisance_vec`: `estimated`, `capped` (a finite maximum-likelihood rate
  above the bound), `boundary` (no finite maximum-likelihood rate exists) or
  `failed`. The rate before the cap is kept as `nuisance_mle_vec`, and
  `param` records `nuisance_cap_multiplier`, `nuisance_num_capped` and
  `nuisance_num_boundary`.
* New exported function `recompute_pvalue()`: redo the posterior, the test
  statistic and the p-values at another `cap_multiplier` on a fitted object,
  with every other setting as it was and without fitting or estimating
  anything again. An object built with `bool_diet = TRUE` has no counts; pass
  the Seurat object the analysis was run on, and the counts and covariates
  are rebuilt from it and checked against the fit. The check compares
  summaries recorded at the fit (each gene's mean count, each covariate's
  column sum, and sums weighted by the position of the cell, which change
  when values are exchanged between cells); it notices counts or metadata
  that were edited since, and it is not a proof of equality.
* New exported function `plot_nuisance()`: each gene's rate, in units of its
  library size, against its mean count, its -log10 p-value or its sparsity,
  with the cap drawn and chosen genes labeled.
* New exported function `plot_fitted_vs_observed()`: the fitted count against
  the observed count for each cell and gene, the points more than `num_sd`
  standard deviations from the fit in red, and optionally each point's
  interval of `num_sd` standard deviations. Both plots need 'ggplot2'
  ('ggrepel' for the first), which are suggested and not required.

## Changes that can affect existing code

* `report_results()` returns seven columns where it returned six; the new
  one is `nuisance_status`.
* The fit on an `eSVD` object has five more per-gene vectors
  (`nuisance_mle_vec`, `nuisance_library_median_vec`, `nuisance_status`,
  `gene_mean_count_vec`, `gene_sparsity_vec`), and `param` has more entries:
  the ones above, the variable names `eSVD()` was called with (`esvd_*`),
  `test_min_cells_per_individual`, and, from `compute_test_per_gene()`, the
  settings of the posterior (`posterior_*`).
* `estimate_nuisance()` and `compute_posterior()` overwrite their entries of
  `param` when they are called again. They used to keep the entries of the
  first call, which `recompute_pvalue()` would then have repeated the
  posterior with.
* An object built with 1.1.0 or earlier works with `report_results()`
  (`nuisance_status` is `NA`) and is refused by `recompute_pvalue()` and the
  two plots, which need what 1.2.0 stores.
* `estimate_nuisance()` on an `eSVD` object now defaults to
  `bool_covariates_as_library = TRUE`, as `compute_posterior()` and
  `compute_test_per_gene()` do and as `eSVD()` has always passed it. Before,
  a stage-by-stage analysis run with defaults estimated the nuisance rate on
  one library size and computed the posterior on another. `eSVD()` and
  `eSVD_helper()` are unaffected; a direct call to `estimate_nuisance()`
  without the argument gives different rates. Pass
  `bool_covariates_as_library = FALSE` for the old default.
* `compute_test_per_gene()` refuses the settings `compute_posterior()`
  refuses, with the same messages: `bool_adjust_covariates = TRUE` together
  with `bool_covariates_as_library = TRUE` (it used to run), a boolean that
  is not one `TRUE` or `FALSE`, and a non-positive or missing `alpha_max`,
  `library_min` or `pseudocount`, or a `nuisance_lower_quantile` outside
  `[0, 1]`. Both functions refuse these by name.
* `nuisance_lower_quantile = NULL` means no floor in `compute_posterior()`,
  as it already did in `compute_test_per_gene()`. It used to empty the vector
  of rates and fail with an unrelated error.

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
