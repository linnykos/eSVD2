# Format covariates

Mainly, this method splits the categorical variables (which should be
\`factor\` variables) into indicator variables (i.e., one-hot encoding),
dropping the last level, and then rescales all the numerical variables
(but does not center them), and computes the `"Log_UMI"` (i.e., log
total counts) for each cell. `"Log_UMI"` is added as its own column.

## Usage

``` r
format_covariates(
  dat,
  covariate_df,
  bool_center = FALSE,
  rescale_numeric_variables = NULL,
  variables_enumerate_all = NULL
)
```

## Arguments

- dat:

  Dataset (either `matrix` or `dgCMatrix`) where the \\n\\ rows
  represent cells and \\p\\ columns represent genes. The rows and
  columns of the matrix should be named.

- covariate_df:

  `data.frame` where each row represents a cell, and the columns are the
  different categorical or numerical variables that you wish to adjust
  for

- bool_center:

  Boolean if the numerical variables should be centered around zero,
  default is `FALSE`

- rescale_numeric_variables:

  A vector of strings denoting the column names in `covariate_df` that
  are numerical and you wish to rescale

- variables_enumerate_all:

  If not `NULL`, this allows you to control specifically which `factor`
  variables in `covariate_df` you would like to split into indicators.
  By default, this is `NULL`, meaning all the `factor` variables are
  split into indicators

## Value

a `matrix` with the same number of rows as `dat`
