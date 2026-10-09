# Report the DE results from eSVD

Report the DE results from eSVD

## Usage

``` r
report_results(input_obj)
```

## Arguments

- input_obj:

  `eSVD` object output from `compute_pvalue`.

## Value

a data frame with one row per gene (in the order of
`input_obj$teststat_vec`) and columns `genes`, `logFC` (\\\log_2\\ of
the case-to-control ratio of posterior means), `logFC_se` (the standard
error of `logFC`, also on the \\\log_2\\ scale; see
`compute_log_fold_change` for what it is a standard error of and how it
compares with those of 'DESeq2', 'dreamlet' and 'NEBULA'), `log10pvalue`
(\\-\log\_{10}\\ of the two-sided p-value), `pvalue`, `pvalue_adj`
(Benjamini-Hochberg) and `nuisance_status` (a factor: what
`estimate_nuisance` did with the gene's over-dispersion, one of
`estimated`, `capped`, `boundary`, `failed`; `NA` for a gene that was
not analyzed, and for every gene of an object built before version
1.2.0). The fold change and the statistic of a `capped` or `boundary`
gene reflect the cap on its rate, and those of a `failed` gene the
floor. `pvalue` underflows to `0` for any gene with `log10pvalue` above
about 308, so genes should be ranked by `log10pvalue`, which keeps the
distinction. `logFC / logFC_se` is not the statistic behind `pvalue`,
which is computed on the linear scale and calibrated with an empirical
null. `logFC_se` is `NA`, with a warning, for an object built before
version 1.1.0, which does not store it.

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
result_df <- report_results(esvd_obj)
utils::head(result_df[order(result_df$log10pvalue, decreasing = TRUE), ])
#>           genes     logFC  logFC_se log10pvalue     pvalue pvalue_adj
#> gene_6   gene_6  1.597469 0.2859912    1.701819 0.01986921  0.4332071
#> gene_1   gene_1 -1.315935 0.2917377    1.393840 0.04037943  0.4332071
#> gene_8   gene_8 -1.238795 0.2844472    1.372885 0.04237547  0.4332071
#> gene_4   gene_4  1.636802 0.3983898    1.363304 0.04332071  0.4332071
#> gene_2   gene_2 -1.752031 0.4640375    1.191151 0.06439447  0.4531028
#> gene_10 gene_10  1.280898 0.3991789    1.167712 0.06796542  0.4531028
#>         nuisance_status
#> gene_6        estimated
#> gene_1        estimated
#> gene_8        estimated
#> gene_4        estimated
#> gene_2        estimated
#> gene_10       estimated
```
