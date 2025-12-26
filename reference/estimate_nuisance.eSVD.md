# Estimate nuisance values for eSVD objects (i.e., over-dispersion)

Assumes a Gamma-Poisson model where the mean and variance are
proportionally related.

## Usage

``` r
# S3 method for class 'eSVD'
estimate_nuisance(
  input_obj,
  bool_covariates_as_library = F,
  bool_library_includes_interept = T,
  bool_use_log = F,
  min_val = 1e-04,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  `eSVD` object outputed from `opt_esvd.eSVD`. Specifically, the
  nuisance parameters will be estimated based on the fit in
  `input_obj[[input_obj[["latest_Fit"]]]]`.

- bool_covariates_as_library:

  Boolean to adjust the numerator in the posterior by the donor
  covariates, default is `FALSE`. This parameter is experimental, and we
  have not yet encountered a scenario where it is useful to be set to be
  `TRUE`.

- bool_library_includes_interept:

  Boolean if the intercept term from the eSVD matrix factorization
  should be included in the calculation for the covariate-adjusted
  library size, default is `TRUE`.

- bool_use_log:

  Boolean if the nuisance (i.e., over-dispersion) parameter should be
  estimated on the log scale, default is `FALSE`.

- min_val:

  Minimum value of the nuisance parameter.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

`eSVD` object with `nuisance_vec` appended to the list in
`input_obj[[input_obj[["latest_Fit"]]]]`.
