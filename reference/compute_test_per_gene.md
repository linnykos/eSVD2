# Gene-wise posterior + test statistic + p-value without storing posterior matrices

This function reproduces the behavior of
[`compute_posterior()`](https://linnykos.github.io/eSVD2/reference/compute_posterior.md),
[`compute_test_statistic()`](https://linnykos.github.io/eSVD2/reference/compute_test_statistic.md),
and
[`compute_pvalue()`](https://linnykos.github.io/eSVD2/reference/compute_pvalue.md),
but it computes everything one gene at a time without ever allocating an
n x p posterior mean/variance matrix.

## Usage

``` r
compute_test_per_gene(
  input_obj,
  alpha_max = 1000,
  bool_adjust_covariates = FALSE,
  bool_covariates_as_library = TRUE,
  bool_stabilize_underdispersion = TRUE,
  library_min = 0.01,
  nuisance_lower_quantile = 0.01,
  pseudocount = 0,
  verbose = 0
)
```

## Arguments

- input_obj:

  eSVD object after nuisance estimation.

- alpha_max:

  Maximum value of the prior mean (same role as in
  `compute_posterior.eSVD`); default 1e3.

- bool_adjust_covariates:

  Boolean; if TRUE, adjust the prior by confounding covariates (same as
  in `compute_posterior.default`).

- bool_covariates_as_library:

  Boolean; if TRUE, include non–case-control covariates in the
  covariate-adjusted library size.

- bool_stabilize_underdispersion:

  Boolean; if TRUE, mean-center log10(nuisance_vec) when it suggests
  under-dispersion.

- library_min:

  Minimum value for the covariate-adjusted library size.

- nuisance_lower_quantile:

  Lower quantile at which to floor nuisance_vec.

- pseudocount:

  Numeric; additional count added to each entry in the count matrix when
  forming the posterior.

- verbose:

  Integer; controls printed messages.

## Value

The input `eSVD` object with:

- `teststat_vec` — Welch t-statistics (length p).

- `case_mean` — case Gaussian means per gene.

- `control_mean` — control Gaussian means per gene.

- `pvalue_list` — list with `df_vec`, `fdr_vec`, `gaussian_teststat`,
  `log10pvalue`, `null_mean`, `null_sd`.
