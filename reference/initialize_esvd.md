# Initialize eSVD

For each gene, this function estimates two ridge-regression penalized
GLMs (using the Poisson model) – one using the `case_control_variable`
and one without, and both sets of coefficients as well as the p-value
(according to a deviance test) is returned. This p-value is on the
log10-scale.

## Usage

``` r
initialize_esvd(
  dat,
  covariates,
  metadata_individual,
  bool_intercept = F,
  case_control_variable = NULL,
  k = 30,
  lambda = 0.01,
  library_size_variable = "Log_UMI",
  offset_variables = "Log_UMI",
  metadata_case_control = NULL,
  verbose = 0
)
```

## Arguments

- dat:

  Dataset (either `matrix` or `dgCMatrix`) where the \\n\\ rows
  represent cells and \\p\\ columns represent genes. The rows and
  columns of the matrix should be named.

- covariates:

  `matrix` object with \\n\\ rows with the same rownames as `dat` where
  the columns represent the different covariates. Notably, this should
  contain only numerical columns (i.e., all categorical variables should
  have already been split into numerous indicator variables), and all
  the columns in `covariates` will (strictly speaking) be included in
  the eSVD matrix factorization model.

- metadata_individual:

  `factor` vector of length \\n\\ that denotes which cell originates
  from which individual.

- bool_intercept:

  Boolean on whether or not an intercept will be included as a
  covariate.

- case_control_variable:

  A string of the column name of `covariates` which depicts the
  case-control status of each cell. Notably, this should be a binary
  variable where a `1` is hard-coded to describe case, and a `0` to
  describe control.

- k:

  Number of latent dimensions.

- lambda:

  Penalty of the `mixed_effect_variables` when using
  [`glmnet::glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)
  to initialize the coefficients.

- library_size_variable:

  A string of the variable name (which must be in `covariates`) of which
  variable denotes the sequenced (i.e., observed) library size.

- offset_variables:

  A vector of strings depicting which column names in `covariate` will
  be set to have a coefficient of `1` automatically (i.e., there will be
  no estimation of their coefficient).

- metadata_case_control:

  (Optional) vector of length \\n\\ with values strictly 0 or 1 that
  denotes if a cell is from cases or controls. By default, this is set
  to `NULL` since the code will extract this information from
  `covariates`.

- verbose:

  Integer

## Value

`eSVD` object with elements `dat`, `covariates`, `initial_Reg` and
`param`
