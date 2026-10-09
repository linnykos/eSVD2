# Format covariates

Mainly, this method splits the categorical variables (which should be
\`factor\` variables) into indicator variables (i.e., one-hot encoding),
dropping the first level of each factor as the reference, then rescales
the numerical variables named in `rescale_numeric_variables` (only
those; the others are left as they are), and computes the `"Log_UMI"`
(i.e., log total counts) for each cell. `"Log_UMI"` is added as its own
column.

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

  Names of `factor` variables in `covariate_df` that get an indicator
  for every level, with no reference level dropped. Every `factor` is
  split into indicators either way; by default (`NULL`) each drops its
  first level. A full set of indicators is collinear with the intercept,
  which `initialize_esvd` refuses.

## Value

a `matrix` with the same number of rows as `dat`, whose first two
columns are `"Intercept"` and `"Log_UMI"`, followed by the numerical
variables and then one indicator column `"<variable>_<level>"` per
retained factor level.

## Details

The rescaling divides each named column by its root-mean-square
(`scale(x, center = FALSE, scale = TRUE)`), not by its standard
deviation, so the result does not have unit variance unless the column
is also centered (`bool_center = TRUE`). What it does deliver is
invariance to the unit the covariate was recorded in: age in years and
age in months give identical columns.

## Examples

``` r
set.seed(10)
dat <- matrix(stats::rpois(60 * 5, lambda = 3), nrow = 60, ncol = 5)
dimnames(dat) <- list(paste0("cell", 1:60), paste0("gene", 1:5))
covariate_df <- data.frame(
  CC = factor(rep(c("control", "case"), each = 30),
              levels = c("control", "case")),
  Sex = factor(rep(c("F", "M"), times = 30)),
  Age = stats::rnorm(60, mean = 40, sd = 10)
)
covariates <- format_covariates(dat = dat,
                                covariate_df = covariate_df,
                                rescale_numeric_variables = "Age")
colnames(covariates)
#> [1] "Intercept" "Log_UMI"   "Age"       "CC_case"   "Sex_M"    
utils::head(covariates, 3)
#>       Intercept  Log_UMI       Age CC_case Sex_M
#> cell1         1 2.079442 0.9020923       0     0
#> cell2         1 2.197225 0.9373963       0     1
#> cell3         1 2.639057 1.2885881       0     0
```
