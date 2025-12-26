# Compute posterior according to Gamma-Poisson model for matrices and sparse matrices.

Compute posterior according to Gamma-Poisson model for matrices and
sparse matrices.

## Usage

``` r
# Default S3 method
compute_posterior(
  input_obj,
  case_control_variable,
  covariates,
  esvd_res,
  library_size_variable,
  nuisance_vec,
  alpha_max = 1000,
  bool_adjust_covariates = F,
  bool_covariates_as_library = T,
  bool_library_includes_interept = T,
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

  Dataset (either `matrix` or `dgCMatrix`) where the \\n\\ rows
  represent cells and \\p\\ columns represent genes. The rows and
  columns of the matrix should be named.

- case_control_variable:

  A string of the column name of `covariates` which depicts the
  case-control status of each cell. Notably, this should be a binary
  variable where a `1` is hard-coded to describe case, and a `0` to
  describe control.

- covariates:

  `matrix` object with \\n\\ rows with the same rownames as `dat` where
  the columns represent the different covariates. Notably, this should
  contain only numerical columns (i.e., all categorical variables should
  have already been split into numerous indicator variables).

- esvd_res:

  Output of `opt_esvd.default`, notably with elements `x_mat`, `y_mat`
  and `z_mat`

- library_size_variable:

  A string of the variable name (which must be in `covariates`) of which
  variable denotes the sequenced (i.e., observed) library size.

- nuisance_vec:

  Vector of non-negative numerics of length `ncol(input_obj)`, such as
  the output of `estimate_nuisance.default`.

- alpha_max:

  Maximum value of numerator when computing posterior, default is `1e3`.

- bool_adjust_covariates:

  Boolean to adjust the numerator in the posterior by the donor
  covariates, default is `FALSE`. This parameter is experimental, and we
  have not yet encountered a scenario where it is useful to be set to be
  `TRUE`.

- bool_covariates_as_library:

  Boolean to include the donor covariates effects in the adjusted
  library size, default is `TRUE`.

- bool_library_includes_interept:

  Boolean if the intercept term from the eSVD matrix factorization
  should be included in the calculation for the covariate-adjusted
  library size, default is `TRUE`.

- bool_return_components:

  Boolean to return the numerator and denominator of the posterior terms
  as well (which will themselves by matrices that are cell-by-gene
  matrices), default is `FALSE`.

- bool_stabilize_underdispersion:

  Boolean to stabilize the over-dispersion parameter, specifically to
  rescale all the over-dispersions the global mean over-disperion is
  less than 1, default is `TRUE`.

- library_min:

  All covariate-adjusted library size smaller than this value are set to
  this value, default is 0.1.

- nuisance_lower_quantile:

  All the nuisance values that are smaller than this quantile are set to
  this quantile.

- pseudocount:

  The additional count that is added to the count matrix, default is 0.

- ...:

  Additional parameters.

## Value

A `list` of elements `posterior_mean_mat` and `posterior_var_mat`
