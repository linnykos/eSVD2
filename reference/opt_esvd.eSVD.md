# Optimize eSVD for eSVD objects

Optimize eSVD for eSVD objects

## Usage

``` r
# S3 method for class 'eSVD'
opt_esvd(
  input_obj,
  fit_name = "fit_First",
  fit_previous = "fit_Init",
  l2pen = 0.1,
  max_iter = 100,
  offset_variables = NULL,
  tol = 1e-06,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  `eSVD` object output from `initialize_esvd` (or an earlier
  `opt_esvd`).

- fit_name:

  String for the name of that will become the current fit when storing
  the results in `input_obj`.

- fit_previous:

  String for the name of the previous fit that this function will grab
  the initialization values from.

- l2pen:

  Small positive number for the amount of penalization for both the
  cells' and the genes' latent vectors as well as the coefficients.

- max_iter:

  Positive integer for number of iterations.

- offset_variables:

  A vector of strings depicting which column names in
  `input_obj$covariates` be treated as an offset during the optimization
  (i.e., their coefficients will not change throughout the
  optimization).

- tol:

  Relative tolerance of the stopping rule: iteration stops once the
  objective changes by at most `tol` times
  `max(1, |previous objective|)`.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

`eSVD` object with added elements with name to whatever `fit_name` was
set to.

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
esvd_obj <- opt_esvd(input_obj = esvd_obj,
                     max_iter = 5,
                     offset_variables = setdiff(colnames(esvd_obj$covariates), "CC"),
                     fit_name = "fit_First",
                     fit_previous = "fit_Init")
esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                               fit_name = "fit_First",
                                               omitted_variables = "Log_UMI")
esvd_obj$latest_Fit
#> [1] "fit_First"
esvd_obj$fit_First$loss
#> [1] -0.7940044 -0.8235692 -0.8328533 -0.8387476 -0.8419992
```
