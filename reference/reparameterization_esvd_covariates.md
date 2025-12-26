# Reparameterize eSVD object

Reparameterize eSVD object

## Usage

``` r
reparameterization_esvd_covariates(
  input_obj,
  fit_name,
  omitted_variables = NULL,
  verbose = 0
)
```

## Arguments

- input_obj:

  `eSVD` object, after either
  [`initialize_esvd()`](https://linnykos.github.io/eSVD2/reference/initialize_esvd.md)
  or
  [`opt_esvd()`](https://linnykos.github.io/eSVD2/reference/opt_esvd.md)

- fit_name:

  The name of the fit in `input_obj` that you wish to reparameterize.
  This should be in `names(input_obj)`

- omitted_variables:

  Either `NULL` (the default) or variables in `input_obj$covariates`
  that should not be reparameterized

- verbose:

  Integer.

## Value

`eSVD` object after adjusting the fit in `fit_name`.
