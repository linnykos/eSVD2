# Compute posterior according to Gamma-Poisson model for eSVD object

The posterior is computed based on whatever `input_obj$latest_Fit` is
set to.

## Usage

``` r
# S3 method for class 'eSVD'
compute_posterior(
  input_obj,
  alpha_max = 1000,
  bool_adjust_covariates = FALSE,
  bool_covariates_as_library = TRUE,
  bool_return_components = FALSE,
  bool_stabilize_underdispersion = TRUE,
  library_min = 0.1,
  nuisance_lower_quantile = 0.01,
  pseudocount = 0,
  ...
)
```

## Arguments

- input_obj:

  `eSVD` object output from `opt_esvd.eSVD`.

- alpha_max:

  Upper bound on the fitted mean without the library, \\\mu\_{ji}\\,
  before it forms the prior \\\mu\_{ji}\beta_j\\, default is `1e3`.

- bool_adjust_covariates:

  Boolean to adjust the numerator in the posterior by the donor
  covariates, default is `FALSE`. This parameter is experimental, and we
  have not yet encountered a scenario where it is useful to be set to be
  `TRUE`.

- bool_covariates_as_library:

  Boolean to include the donor covariates effects in the adjusted
  library size, default is `TRUE`

- bool_return_components:

  Boolean to return the numerator and denominator of the posterior terms
  as well (themselves cell-by-gene matrices), default is `FALSE`

- bool_stabilize_underdispersion:

  Boolean, default `TRUE`. `nuisance_vec` holds Gamma *rates*, which
  share the units of the library size, so a geometric mean above 1 means
  less over-dispersion than a unit rate in those units (not
  under-dispersion: a finite rate is always over-dispersed). When that
  happens (`mean(log10(nuisance_vec)) > 0`), every gene's rate, capped
  ones included, is divided by the one geometric mean so that
  `mean(log10(nuisance_vec))` becomes 0; otherwise the rates are left
  alone

- library_min:

  All covariate-adjusted library size smaller than this value are set to
  this value, default is 0.1.

- nuisance_lower_quantile:

  All the nuisance values that are smaller than this quantile are set to
  this quantile, default is 0.01; one number in \\\[0, 1\]\\, or `NULL`
  for no floor.

- pseudocount:

  The additional count that is added to the count matrix, default is 0.

- ...:

  Additional parameters.

## Value

`eSVD` object with `posterior_mean_mat` and `posterior_var_mat` appended
to the list in `input_obj[[input_obj[["latest_Fit"]]]]`.

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
esvd_obj <- estimate_nuisance(input_obj = esvd_obj)
esvd_obj <- compute_posterior(input_obj = esvd_obj,
                              alpha_max = 2 * max(sim$obs_mat))
dim(esvd_obj$fit_First$posterior_mean_mat)
#> [1] 120  40
```
