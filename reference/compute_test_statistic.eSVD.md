# Compute test statistics for eSVD object

Compute test statistics for eSVD object

## Usage

``` r
# S3 method for class 'eSVD'
compute_test_statistic(
  input_obj,
  min_cells_per_individual = 3,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  `eSVD` object output from `compute_posterior.eSVD`.

- min_cells_per_individual:

  Minimum number of cells an individual must contribute; see
  `compute_test_statistic.default`.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

`eSVD` object with added elements, each a named numeric vector with one
entry per gene: `teststat_vec` (the Welch statistic), `case_mean` and
`control_mean` (the mean over the individuals of each arm of their mean
posterior expression), `case_var` and `control_var` (the variance of
each arm, the two quantities the statistic divides by), and `log2fc_vec`
and `log2fc_se_vec` (the log2 fold change and its standard error; see
`compute_log_fold_change`). The individuals of each arm and
`min_cells_per_individual` are recorded in `param`.

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
esvd_obj <- compute_test_statistic(input_obj = esvd_obj)
utils::head(esvd_obj$teststat_vec)
#>    gene_1    gene_2    gene_3    gene_4    gene_5    gene_6 
#> -4.710563 -4.069603  1.434486  4.243108 -2.545816  5.770229 
```
