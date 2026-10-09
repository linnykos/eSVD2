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

## Examples

``` r
set.seed(10)
sim <- generate_null(cell_per_person = 15, num_genes = 40,
                     num_individuals = 8)
esvd_obj <- initialize_esvd(dat = sim$obs_mat,
                            covariates = sim$covariates,
                            metadata_individual = sim$metadata_individual,
                            case_control_variable = "CC",
                            bool_intercept = TRUE,
                            k = 2,
                            lambda = 0.1)
esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                               fit_name = "fit_Init",
                                               omitted_variables = "Log_UMI")
# x_mat is now orthogonal to the retained covariates
max(abs(crossprod(esvd_obj$fit_Init$x_mat,
                  esvd_obj$covariates[, "Sex"])))
#> [1] 1.361133e-13
```
