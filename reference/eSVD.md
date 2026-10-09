# Run the full eSVD-DE pipeline on a Seurat object

The end-to-end wrapper: `format_covariates`, `initialize_esvd` and two
rounds of `opt_esvd` (the case-control coefficient held fixed, then
freed), each of the three followed by
`reparameterization_esvd_covariates`, `estimate_nuisance`, and then
either the fused `compute_test_per_gene` (`bool_diet = TRUE`) or
`compute_posterior`, `compute_test_statistic` and `compute_pvalue` in
sequence.

## Usage

``` r
eSVD(
  batch_var_prefix,
  case_control_levels,
  case_control_var,
  categorical_vars,
  id_var,
  numerical_vars,
  seurat_obj,
  alpha_max = NULL,
  bool_adjust_covariates = FALSE,
  bool_covariates_as_library = TRUE,
  bool_diet = TRUE,
  bool_intercept = TRUE,
  bool_library_includes_interept = TRUE,
  bool_stabilize_underdispersion = TRUE,
  bool_use_log = FALSE,
  cap_multiplier = 10,
  intermediate_save = NULL,
  k = 30,
  l2pen = 0.1,
  lambda = 0.1,
  library_min = 0.1,
  max_iter = 100,
  min_cells_per_individual = 3,
  pseudocount = 0,
  tol = 1e-06,
  verbose = 0
)
```

## Arguments

- batch_var_prefix:

  `NULL`, or a regular expression matched against the names of
  `categorical_vars`; the matching covariates are treated as batch
  variables and their coefficients are not reparameterized.

- case_control_levels:

  Character vector of length 2 giving the two levels of
  `case_control_var`, **control first, then case**. Reversing them flips
  the sign of every test statistic.

- case_control_var:

  Name of the column in `seurat_obj@meta.data` holding the case-control
  status.

- categorical_vars:

  Character vector of column names in `seurat_obj@meta.data` to adjust
  for as factors (one-hot encoded, first level dropped), or `NULL`. A
  variable with a single level is dropped.

- id_var:

  Name of the column in `seurat_obj@meta.data` naming each cell's
  individual (donor).

- numerical_vars:

  Character vector of column names in `seurat_obj@meta.data` to adjust
  for as numerical covariates, or `NULL`. They are rescaled by
  `format_covariates`.

- seurat_obj:

  A `Seurat` object whose `"RNA"` assay has a `"counts"` layer (genes by
  cells) and whose `meta.data` holds the variables above.

- alpha_max:

  Maximum value of the prior mean in the posterior; `NULL` (the default)
  uses twice the largest count.

- bool_adjust_covariates:

  See `compute_posterior`.

- bool_covariates_as_library:

  See `estimate_nuisance` and `compute_posterior`.

- bool_diet:

  If `TRUE` (the default), compute the posterior, test statistics and
  p-values one gene at a time without allocating cell-by-gene posterior
  matrices, and drop `dat`, `covariates` and the intermediate fits
  `fit_Init` and `fit_First` from the returned object. The final fit
  (`fit_Second`: `x_mat`, `y_mat`, `z_mat`, `nuisance_vec`) is kept
  either way.

- bool_intercept:

  Whether the per-gene GLMs include an intercept.

- bool_library_includes_interept:

  See `estimate_nuisance`.

- bool_stabilize_underdispersion:

  See `compute_posterior`.

- bool_use_log:

  See `estimate_nuisance`.

- cap_multiplier:

  One positive number, default `10`: no gene's nuisance rate may exceed
  this many times the median of its library size, and `Inf` means no
  cap. See `estimate_nuisance`. The factorization does not depend on it,
  so `recompute_pvalue` can change it afterwards without fitting again.

- intermediate_save:

  `NULL`, or a file path at which the object is saved after each major
  stage.

- k:

  Number of latent dimensions; must not exceed the number of genes.

- l2pen:

  L2 penalty for `opt_esvd`.

- lambda:

  Ridge penalty for `initialize_esvd`.

- library_min:

  See `compute_posterior`.

- max_iter:

  Maximum iterations for each `opt_esvd` call.

- min_cells_per_individual:

  Minimum number of cells every individual must contribute; an
  individual with fewer is an error. `0` disables the check. See
  `compute_test_statistic`.

- pseudocount:

  See `compute_posterior`.

- tol:

  Convergence tolerance for `opt_esvd`.

- verbose:

  Integer; `0` is silent.

## Value

An `eSVD` object with elements `teststat_vec`, `case_mean`,
`control_mean`, `case_var`, `control_var`, `log2fc_vec` and
`log2fc_se_vec` (the log2 fold change and its standard error; see
`compute_log_fold_change`), `pvalue_list` (see `compute_pvalue`),
`param`, `case_control`, `individual`, `latest_Fit` and the fit it names
(an `eSVD_Fit` with `x_mat`, `y_mat`, `z_mat`, `nuisance_vec` and the
other per-gene vectors of `estimate_nuisance`), plus `dat`,
`covariates`, the earlier fits and the posterior matrices when
`bool_diet = FALSE`. `param` records the settings of every stage (the
`opt_*` entries are those of the first `opt_esvd` call, because a later
call does not overwrite an entry that is already there), the number of
genes whose nuisance rate was capped (`nuisance_num_capped`) and the
arguments of this call that name the variables (as `esvd_*`), which is
what lets `recompute_pvalue` and `plot_fitted_vs_observed` rebuild the
counts and covariates from `seurat_obj` alone. `report_results` turns
the object into a data frame, and `plot_nuisance` and
`plot_fitted_vs_observed` draw its diagnostics.

## Details

**This function refuses inputs the model cannot handle**, rather than
filtering them: a gene whose counts are all zero, an individual with
fewer than `min_cells_per_individual` cells, or an arm with fewer than
two individuals each produce an error naming the offending genes or
individuals. The recommended entry point, `eSVD_helper`, applies the
cohort filters and removes all-zero genes before calling this function,
and reinserts the removed genes afterwards with a `gene_status` record.
It does *not* refuse an underpowered cohort (few cells or few
individuals overall); those minima are `eSVD_helper`'s policy.

The model is a Gamma-Poisson hierarchy: \\A\_{ji} \mid \lambda\_{ji}
\sim \mathrm{Poisson}(\ell\_{ji} \lambda\_{ji})\\ with \\\lambda\_{ji}
\sim \mathrm{Gamma}(\mathrm{mean} = \mu\_{ji}, \mathrm{var} = \gamma_j
\mu\_{ji})\\, where \\\mu\_{ji}\\ is a low-rank function of the cell and
gene latent vectors plus the case-control covariate, and \\\ell\_{ji}\\
is the covariate-adjusted sequencing depth. Note that the package stores
the Gamma *rate* \\\beta_j = 1/\gamma_j\\ as `nuisance_vec`, so larger
values mean *less* over-dispersion. Differential expression is a
difference in mean expression between the case and control individuals,
tested with a Welch statistic on the per-individual posterior means. See
Lin, Qiu and Roeder (2024)
[doi:10.1186/s12859-024-05724-7](https://doi.org/10.1186/s12859-024-05724-7)
.

## Examples

``` r
# \donttest{
if(requireNamespace("SeuratObject", quietly = TRUE)){
  set.seed(10)
  sim <- generate_null(cell_per_person = 15, num_genes = 40,
                       num_individuals = 8)
  meta_df <- data.frame(sim$covariates[, c("CC", "Sex", "Age")])
  meta_df$Sex <- factor(meta_df$Sex)
  meta_df$Individual <- sim$metadata_individual
  rownames(meta_df) <- rownames(sim$obs_mat)
  # Seurat rewrites "gene_1" as "gene-1" (with a warning); rename first
  count_mat <- Matrix::t(sim$obs_mat)
  rownames(count_mat) <- gsub("_", "", rownames(count_mat))
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = count_mat,
                                                 meta.data = meta_df)

  esvd_obj <- eSVD(batch_var_prefix = NULL,
                   case_control_levels = c("0", "1"),
                   case_control_var = "CC",
                   categorical_vars = "Sex",
                   id_var = "Individual",
                   numerical_vars = "Age",
                   seurat_obj = seurat_obj,
                   k = 2,
                   max_iter = 5)
  utils::head(report_results(esvd_obj))
}
#>       genes     logFC  logFC_se log10pvalue     pvalue pvalue_adj
#> gene1 gene1 -1.468000 0.2676634   1.4108375 0.03882956  0.5406965
#> gene2 gene2 -1.816972 0.4706980   1.0651751 0.08606467  0.5737645
#> gene3 gene3  1.050755 0.9300542   0.4481459 0.35633137  0.9944779
#> gene4 gene4  1.609953 0.4094298   1.3638726 0.04326408  0.5406965
#> gene5 gene5 -1.686198 0.5471413   0.8323053 0.14712780  0.7356390
#> gene6 gene6  1.764437 0.2837168   1.8455667 0.01427031  0.5406965
#>       nuisance_status
#> gene1       estimated
#> gene2       estimated
#> gene3       estimated
#> gene4       estimated
#> gene5       estimated
#> gene6       estimated
# }
```
