# Compute p-values

Converts the Welch statistics in `input_obj$teststat_vec` to Gaussian
statistics via their Welch-Satterthwaite degrees of freedom, fits an
empirical null with `multtest`, and returns two-sided p-values on the
\\-\log\_{10}\\ scale with Benjamini-Hochberg adjustment.

## Usage

``` r
compute_pvalue(input_obj, min_cells_per_individual = 3, verbose = 0, ...)
```

## Arguments

- input_obj:

  `eSVD` object output from `compute_test_statistic`.

- min_cells_per_individual:

  Minimum number of cells an individual must contribute; see
  `compute_test_statistic`.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

`eSVD` object with added element `"pvalue_list"`, a list with `df_vec`,
`fdr_vec`, `gaussian_teststat`, `log10pvalue` (\\-\log\_{10}\\ of the
two-sided p-value), `method` (which empirical-null estimator ran; see
`multtest`), `null_mean` and `null_sd`.

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
esvd_obj <- compute_pvalue(input_obj = esvd_obj)
esvd_obj$pvalue_list$method
#> [1] "locfdr"
utils::head(esvd_obj$pvalue_list$fdr_vec)
#>    gene_1    gene_2    gene_3    gene_4    gene_5    gene_6 
#> 0.4332071 0.4531028 0.9764068 0.4332071 0.8221153 0.4332071 
```
