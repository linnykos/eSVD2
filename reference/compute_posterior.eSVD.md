# Compute posterior according to Gamma-Poisson model for eSVD object

The posterior is computed based on whatever `input_obj$latest_Fit` is
set to.

## Usage

``` r
# S3 method for class 'eSVD'
compute_posterior(
  input_obj,
  alpha_max = 1000,
  bool_adjust_covariates = F,
  bool_covariates_as_library = T,
  bool_return_components = F,
  bool_stabilize_underdispersion = T,
  library_min = 0.1,
  nuisance_lower_quantile = 0.01,
  pseudocount = 0,
  ...
)
```

## Arguments

- input_obj:

  `eSVD` object outputed from `opt_esvd.eSVD`.

- alpha_max:

  Maximum value of numerator when computing posterior, default is `1e3`.

- bool_adjust_covariates:

  Boolean to adjust the numerator in the posterior by the donor
  covariates, default is `FALSE`. This parameter is experimental, and we
  have not yet encountered a scenario where it is useful to be set to be
  `TRUE`.

- bool_covariates_as_library:

  Boolean to include the donor covariates effects in the adjusted
  library size, default is `TRUE`

- bool_return_components:

  Boolean to return the numerator and denominator of the posterior terms
  as well (which will themselves by matrices that are cell-by-gene
  matrices), default is `FALSE`

- bool_stabilize_underdispersion:

  Boolean to stabilize the over-dispersion parameter, specifically to
  rescale all the over-dispersions the global mean over-disperion is
  less than 1, default is `TRUE`

- library_min:

  All covariate-adjusted library size smaller than this value are set to
  this value, default is 0.1.

- nuisance_lower_quantile:

  All the nuisance values that are smaller than this quantile are set to
  this quantile, default is 0.01

- pseudocount:

  The additional count that is added to the count matrix, default is 0.

- ...:

  Additional parameters.

## Value

`eSVD` object with `posterior_mean_mat` and `posterior_var_mat` appended
to the list in `input_obj[[input_obj[["latest_Fit"]]]]`.
