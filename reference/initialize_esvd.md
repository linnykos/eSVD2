# Initialize eSVD

For each gene, this function fits one ridge-penalized Poisson GLM of the
counts on `covariates` (with `offset_variables` as offsets;
[`stats::glm`](https://rdrr.io/r/stats/glm.html) when at most one
covariate is left to estimate) to obtain `z_mat`, and then takes a
rank-`k` SVD of `log1p(dat)` minus the fitted covariate effects to
obtain `x_mat` and `y_mat`. The result is stored as `fit_Init`.

## Usage

``` r
initialize_esvd(
  dat,
  covariates,
  metadata_individual,
  bool_intercept = FALSE,
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

  Ridge penalty used by
  [`glmnet::glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)
  on every covariate that is neither the intercept nor an offset; must
  lie in \\\[10^{-4}, 10^4\]\\.

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

`eSVD` object with elements `dat`, `covariates`, `param`, `fit_Init` (an
`eSVD_Fit` with `x_mat`, `y_mat`, `z_mat`), `latest_Fit` (`"fit_Init"`),
`case_control` and `individual`

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
names(esvd_obj)
#> [1] "dat"          "covariates"   "param"        "fit_Init"     "latest_Fit"  
#> [6] "case_control" "individual"  
dim(esvd_obj$fit_Init$x_mat)
#> [1] 120   2
```
