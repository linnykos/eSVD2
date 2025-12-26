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

  `eSVD` object outputed from `apply_initial_threshold`.

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
  `input_obj$covariate` be treated as an offset during the optimization
  (i.e., their coefficients will not change throughout the
  optimization).

- tol:

  Small positive number to differentiate between zero and non-zero.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

`eSVD` object with added elements with name to whatever `fit_name` was
set to.
