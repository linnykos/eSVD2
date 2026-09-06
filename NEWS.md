# eSVD2 1.0.2

First CRAN submission. Relative to the GitHub release 1.0.0:

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
