# eSVD2 1.2.0

The first version intended for CRAN. Everything below is relative to
1.0.1.07, the version on the `master` branch (April 2026) and the last one
installed from GitHub. The development versions 1.0.2 and 1.1.0 in between
were never released; their changes are folded in here.

## Results differ from 1.0.1.07

Every analysis gives different statistics, p-values and fold changes, for
all genes. The reasons, from largest to smallest effect:

* **The nuisance rate is capped.** `estimate_nuisance()` bounds each gene's
  Gamma rate (the reciprocal of the over-dispersion) at `cap_multiplier`
  times the median of the gene's library size, with `cap_multiplier = 10`
  by default, and floors it at `min_val` times the same median. A gene
  whose counts are no more variable than a Poisson has no finite
  maximum-likelihood rate; 1.0.1.07 gave such a gene the largest library
  size (see the next item), its posterior followed the fit, and its
  statistic was inflated. The cap is a calibration device, and the
  documentation says so. `cap_multiplier = Inf` removes it. `eSVD()` and
  `eSVD_helper()` take the new argument `cap_multiplier`.
* **The rate is now the maximum-likelihood estimate.** `gamma_rate()`
  truncated its search at the largest library size of the gene and returned
  that bound as if it were the estimate whenever the true maximizer lay
  above it. The bracket now grows until the maximizer is inside it.
* **A degraded empirical null is no longer silent.** Welch statistics are
  mapped to Gaussian statistics on the log scale, so one strongly
  up-regulated gene no longer produces an infinite statistic that pushed
  `locfdr` onto its crudest fallback for the whole data set. `multtest()`
  reports which estimator ran (`pvalue_list$method`) and warns whenever it
  is not `locfdr`. The truncated-Gaussian fallback is Efron's estimator
  with the constraint `p0 <= 1` restored (without it, 200 draws from
  `N(0, 1)` gave a null standard deviation of 2.6), and the moment
  estimator is corrected for truncation.
* **Two defaults changed.** `compute_test_per_gene()` floors the library
  size at `library_min = 0.1`, as `eSVD()` and `compute_posterior()` always
  did, instead of `0.01`. `estimate_nuisance()` on an `eSVD` object defaults
  to `bool_covariates_as_library = TRUE`, as `compute_posterior()` does and
  as `eSVD()` has always passed it; a stage-by-stage analysis run with
  defaults used to estimate the rate on one library size and compute the
  posterior on another.
* **The fit is deterministic.** The iterative SVD solvers start from a fixed
  vector, so identical inputs give identical fits and the user's random
  seed is not consumed. Two runs of 1.0.1.07 on the same data differed by
  up to 0.006 in the Welch statistic.

## New functions

* `eSVD_helper()` wraps `eSVD()` with the cohort filters of the paper, and
  `filter_cohort()` applies them alone: minimum cells per individual, per
  arm and overall, and minimum individuals per arm. All-zero genes are
  removed before the fit and reinserted afterwards with `NA` results; the
  outcome of every filter is recorded in `gene_status` and `param`.
* `compute_log_fold_change()` gives each gene a log2 fold change between
  the case and control individuals and a delta-method standard error for
  it, on the unit of replication of the individual (as 'DESeq2', 'dreamlet'
  and 'NEBULA' report it). Both are computed by `compute_test_statistic()`
  and `compute_test_per_gene()` and stored as `log2fc_vec` and
  `log2fc_se_vec`. `?compute_log_fold_change` says what the standard error
  does and does not account for.
* `recompute_pvalue()` redoes the posterior, the statistic and the p-values
  at another `cap_multiplier` on a fitted object, without fitting again. An
  object built with `bool_diet = TRUE` has no counts; pass the Seurat object
  the analysis was run on, and the counts and covariates are rebuilt and
  checked against summaries recorded at the fit.
* `plot_nuisance()` draws each gene's rate, in units of its library size,
  against its mean count, its -log10 p-value or its sparsity, with the cap
  drawn. `plot_fitted_vs_observed()` draws the fitted against the observed
  count per cell and gene, marking the points more than `num_sd` model
  standard deviations from the fit. Both need 'ggplot2' (and 'ggrepel'),
  which are suggested, not required.
* `multtest()` and `reparameterization_esvd_covariates()` are exported by
  name, with documentation and examples.

## Changes that can affect existing code

* `report_results()` returns seven columns where it returned four: `genes`,
  `logFC`, `logFC_se`, `log10pvalue`, `pvalue`, `pvalue_adj` and
  `nuisance_status`. Code that reads columns by position breaks.
* Every gene has a `nuisance_status` (`estimated`, `capped`, `boundary` or
  `failed`), stored as a factor beside `nuisance_vec` in the fit, with the
  uncapped estimate in `nuisance_mle_vec`. The fit and `param` carry more
  entries (gene means and sparsity, the library medians, the counts of
  capped and boundary genes, the variable names `eSVD()` was called with,
  the posterior settings, and the individuals of each arm).
* `compute_test_statistic()`, `compute_pvalue()` and
  `compute_test_per_gene()` take `min_cells_per_individual` (default 3),
  which `eSVD()` passes through.
* `opt_esvd()` on matrices takes `nuisance_vec = NULL` (meaning 1 for every
  gene) instead of a vector of `NA`; four of the seven families consumed
  the `NA` inside the objective and failed in the line search.
* `estimate_nuisance()`, `compute_posterior()`, `compute_test_statistic()`
  and `compute_test_per_gene()` overwrite their entries of `param` when
  rerun, so a rerun on a changed cohort or setting is recorded.
* An object saved by 1.0.1.07 still works with `report_results()`
  (`logFC_se` and `nuisance_status` are `NA`, with a warning) and is
  refused by `compute_log_fold_change()`, `recompute_pvalue()` and the two
  plots, which need what 1.2.0 stores.
* `compute_pvalue()` on an object built with `eSVD(bool_diet = TRUE)` stops
  with a message saying the p-values were already computed by
  `compute_test_per_gene()`; it used to fail inside `.compute_df()`.

## Correctness fixes

Each of these was a silent wrong answer in 1.0.1.07.

* `initialize_esvd()` zeroes `NA` counts in sparse matrices too (a single
  `NA` in a `dgCMatrix` errored inside 'glmnet'), keeps the library-size
  offset in the single-covariate `stats::glm` fallback (a design with only
  the case-control indicator was initialized without depth adjustment),
  and refuses all-zero genes instead of mis-fitting them.
* `reparameterization_esvd_covariates()` regresses on the design matrix
  directly instead of through `stats::lm(x ~ ., data = ...)`, which
  rewrote non-syntactic covariate names such as `"Diagnosis_ASD (severe)"`
  (a subscript error) and returned `NA` coefficients for aliased columns
  that propagated into `z_mat` and every posterior. The identification step
  floors an eigenvalue at or below `tol` instead of inverting it (an exact
  zero gave `Inf`, a rounding-error negative gave `NaN`).
* `opt_esvd()` counts failed Newton line searches in C++ and raises one R
  warning at the end, instead of `Rcpp::warning()` from inside a C++ frame
  (which can skip destructors under `options(warn = 2)`), and stops with a
  message when the objective becomes non-finite. Its verbose output printed
  nothing useful.
* `gamma_rate()` and `log_gamma_rate()` check that `x`, `mu` and `s` have
  the same length (a shorter input was an out-of-bounds read).
* The C++ data loader rejects an unsupported matrix type on entry instead
  of returning a null pointer that failed later as "external pointer is
  not valid".
* `bool_stabilize_underdispersion` keeps `nuisance_vec` a named numeric
  vector on both branches.
* Column standard deviations of a sparse matrix are computed from the
  stored non-zeros without cancellation; the two-pass formula returned
  `NaN` for a column with a large mean and a small variance.
* `compute_pvalue()` stops when an arm has one individual (the Welch
  degrees of freedom were zero and every p-value `NaN`), and
  `multtest()` stops when no estimator can fit a null and warns when the
  null has no spread.
* `fisher_test(verbose = TRUE)` prints its summary; the string was built
  and discarded.

## Input validation

New errors that name the cause, in place of failures far downstream:

* A rank-deficient design is refused by `initialize_esvd()` and by the
  reparameterization, naming the aliased columns. Including an indicator
  for every individual, or for every level of a factor, causes this. The
  2024 ASD design of the articles is rank deficient (`Seqbatch` and
  `region` are determined by `Capbatch`); 1.0.1.07 fit it anyway.
* `eSVD()` refuses a `batch_var_prefix` matching no variable, variables
  missing from the Seurat object, a `case_control_var` whose values are not
  the two `case_control_levels`, `k` above the number of genes, all-zero
  genes, individuals in both arms, an arm with fewer than two individuals,
  and individuals with fewer than `min_cells_per_individual` cells.
* The posterior's settings are validated once, for `compute_posterior()`,
  `compute_test_per_gene()` and `eSVD()` alike: a boolean that is not one
  `TRUE` or `FALSE`, `bool_adjust_covariates = TRUE` together with
  `bool_covariates_as_library = TRUE`, and a non-positive `alpha_max`,
  `library_min` or `pseudocount`. `nuisance_lower_quantile = NULL` means no
  floor on both paths; it used to empty the vector of rates.
* `initialize_esvd()` refuses `NA` or non-finite covariates and a `k` that
  is not a positive integer; `opt_esvd()` refuses `NA` in any input;
  `format_covariates()` names a factor with one level; the family
  constructor names an unknown family; `.svd_safe()` refuses a `mean_vec`
  or `sd_vec` of the wrong type or length.

## Packaging

* `DESCRIPTION`: the title and description are rewritten for CRAN and cite
  the paper; 'Rmpfr' and 'sparseMatrixStats' are no longer needed (the
  log-scale p-values and the sparse standard deviations replace them);
  'Rcpp' moved from `Depends` to `Imports`; 'stats' and 'utils' are
  imported; 'devtools', 'EnhancedVolcano' and 'npregfast' leave `Suggests`
  and 'ggplot2', 'ggrepel', 'rmarkdown' and 'withr' join it; `biocViews`
  is dropped and `BugReports` added.
* `LICENSE` is the two-line CRAN form for MIT, and `NEWS.md` and
  `cran-comments.md` exist.
* `NAMESPACE` exports only the functions marked `@export`; the pattern
  export that exposed every non-dot object, including the raw Rcpp bindings
  that take external pointers, is gone.
* Every exported function has documentation and a runnable example, and
  about 35 misleading comments and roxygen descriptions in `R/` and `src/`
  are corrected.
* The test suite grows from 12 files and 31 blocks to 29 files, 298 blocks
  and about 2,100 expectations, with fixtures regenerated by
  `tests/assets/data_generation.R` (the shipped fixture shrinks from 10 MB
  to 0.4 MB).
* The vignette `eSVD2.Rmd` runs without 'Seurat' or 'ggplot2' (those
  chunks are skipped) and gains sections on `report_results()`, the
  over-dispersion diagnostics and `recompute_pvalue()`. The two ASD
  tutorials, which need multi-gigabyte downloads, are pkgdown articles under
  `vignettes/articles/` and were rerun under 1.2.0 through `eSVD_helper()`,
  with the diagnostic plots and a comparison of caps.
* `.Rbuildignore` keeps the collaboration material (`additional_context/`,
  the `CLAUDE*.md` and `HISTORY_*.md` notes, `.githooks/`, `oldcode/`) out
  of the tarball; it stays on GitHub.

# eSVD2 1.0.1.07

The `master` branch on GitHub (April 2026): small fixes on top of 1.0.0.

# eSVD2 1.0.0

GitHub release accompanying Lin, Qiu and Roeder (2024),
<doi:10.1186/s12859-024-05724-7>.
