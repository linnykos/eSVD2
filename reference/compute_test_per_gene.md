# Gene-wise posterior + test statistic + p-value without storing posterior matrices

This function reproduces the behavior of
[`compute_posterior()`](https://linnykos.github.io/eSVD2/reference/compute_posterior.md),
[`compute_test_statistic()`](https://linnykos.github.io/eSVD2/reference/compute_test_statistic.md),
and
[`compute_pvalue()`](https://linnykos.github.io/eSVD2/reference/compute_pvalue.md),
but it computes everything one gene at a time without ever allocating an
n x p posterior mean/variance matrix. It refuses the same settings as
`compute_posterior`, among them `bool_adjust_covariates = TRUE` together
with `bool_covariates_as_library = TRUE`.

## Usage

``` r
compute_test_per_gene(
  input_obj,
  alpha_max = 1000,
  bool_adjust_covariates = FALSE,
  bool_covariates_as_library = TRUE,
  bool_stabilize_underdispersion = TRUE,
  library_min = 0.1,
  min_cells_per_individual = 3,
  nuisance_lower_quantile = 0.01,
  pseudocount = 0,
  verbose = 0
)
```

## Arguments

- input_obj:

  eSVD object after nuisance estimation.

- alpha_max:

  Upper bound on the fitted mean without the library, before it forms
  the prior (same role as in `compute_posterior.eSVD`); default 1e3.

- bool_adjust_covariates:

  Boolean; if TRUE, divide the posterior numerator (count plus prior) by
  the exponentiated effect of the confounding covariates (same as in
  `compute_posterior.default`).

- bool_covariates_as_library:

  Boolean; if TRUE, include non–case-control covariates in the
  covariate-adjusted library size.

- bool_stabilize_underdispersion:

  Boolean; if TRUE, mean-center log10(nuisance_vec) when its mean is
  above 0, i.e. when the geometric mean of the rates is above 1 (less
  over-dispersion than a unit rate; see `compute_posterior.eSVD`).

- library_min:

  Minimum value for the covariate-adjusted library size; default 0.1,
  the same as `compute_posterior.eSVD`.

- min_cells_per_individual:

  Minimum number of cells an individual must contribute; see
  `compute_test_statistic`.

- nuisance_lower_quantile:

  Lower quantile at which to floor nuisance_vec, one number in \[0, 1\];
  `NULL` means no floor.

- pseudocount:

  Numeric; additional count added to each entry in the count matrix when
  forming the posterior.

- verbose:

  Integer; controls printed messages.

## Value

The input `eSVD` object with:

- `teststat_vec` — Welch t-statistics (length p).

- `case_mean` — case Gaussian means per gene.

- `control_mean` — control Gaussian means per gene.

- `case_var`, `control_var` — the variance of each arm, the two
  quantities the statistic divides by.

- `log2fc_vec`, `log2fc_se_vec` — the log2 fold change and its standard
  error; see `compute_log_fold_change`.

- `pvalue_list` — list with `df_vec`, `fdr_vec`, `gaussian_teststat`,
  `log10pvalue`, `method`, `null_mean`, `null_sd`; see `compute_pvalue`.

The individuals of each arm are recorded in `param`, as
`compute_test_statistic` does, and so are the settings of the posterior,
under the names `compute_posterior` uses (`posterior_alpha_max` and so
on). Calling the function again overwrites them.

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
esvd_obj <- compute_test_per_gene(input_obj = esvd_obj,
                                  alpha_max = 2 * max(sim$obs_mat))
utils::head(esvd_obj$teststat_vec)
#>    gene_1    gene_2    gene_3    gene_4    gene_5    gene_6 
#> -4.710563 -4.069603  1.434486  4.243108 -2.545816  5.770229 
# no cell-by-gene posterior matrices were stored
names(esvd_obj$fit_First)
#>  [1] "x_mat"                       "y_mat"                      
#>  [3] "z_mat"                       "loss"                       
#>  [5] "nuisance_vec"                "nuisance_mle_vec"           
#>  [7] "nuisance_library_median_vec" "nuisance_status"            
#>  [9] "gene_mean_count_vec"         "gene_sparsity_vec"          
```
