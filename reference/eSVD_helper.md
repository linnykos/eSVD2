# Run eSVD-DE on a Seurat object, with cohort filtering and gene status

The recommended entry point. It wraps `eSVD` in the five steps a real
dataset needs and `eSVD` itself refuses to perform:

1.  Drop individuals with fewer than `min_cells_per_id` cells, with a
    warning naming them (`filter_cohort`).

2.  Check `min_cells`, `min_ids`, `min_ids_per_arm` and
    `min_cells_casecontrol` on the filtered cohort; on failure, warn and
    return `NA` (`filter_cohort`).

3.  Label every gene's `gene_status` from the filtered counts and remove
    the genes whose counts are all zero. This runs *after* the drop,
    because a gene expressed only in a dropped individual becomes
    all-zero.

4.  Call `eSVD` on what remains; every further argument goes through
    `...`.

5.  Reinsert the removed genes at their original positions.

## Usage

``` r
eSVD_helper(
  batch_var_prefix,
  case_control_levels,
  case_control_var,
  categorical_vars,
  id_var,
  numerical_vars,
  seurat_obj,
  min_cells = 20,
  min_cells_casecontrol = 20,
  min_cells_per_id = 3,
  min_ids = 4,
  min_ids_per_arm = 2,
  verbose = 0,
  ...
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

- min_cells:

  Reject a cohort with this many cells or fewer, after the drop.

- min_cells_casecontrol:

  Reject a cohort in which either arm has this many cells or fewer.

- min_cells_per_id:

  Drop every individual with fewer than this many cells; `0` or at least
  `3`.

- min_ids:

  Reject a cohort with this many individuals or fewer, pooled across
  arms.

- min_ids_per_arm:

  Reject a cohort in which either arm has fewer than this many
  individuals. The default `2` is the smallest number for which the
  Welch degrees of freedom are defined.

- verbose:

  Integer; `0` is silent.

- ...:

  Further arguments to `eSVD`, such as `k`, `bool_diet`,
  `cap_multiplier` or `intermediate_save`.

## Value

Either `NA` (the cohort was rejected; see the warning), or the `eSVD`
object returned by `eSVD` with the removed genes reinserted as described
above and an added element `gene_status`.

## Details

**`gene_status`** is a named factor over the original genes, in their
original order, with levels `"analyzed"` (1) and `"all_zero"` (2). An
`all_zero` gene comes back with `NA` in every per-gene estimate
(`teststat_vec`, `case_mean`, `control_mean`, `case_var`, `control_var`,
`log2fc_vec`, `log2fc_se_vec`, `pvalue_list$df_vec`,
`pvalue_list$gaussian_teststat`, and the rows of `y_mat`, `z_mat`,
`nuisance_vec` and the other per-gene vectors of `estimate_nuisance` in
the final fit), a `pvalue_list$log10pvalue` of `0` (that is, a p-value
of 1) and a `pvalue_list$fdr_vec` of `1`. The empirical null and the
Benjamini-Hochberg adjustment are computed on the analyzed genes only,
which is guaranteed by construction: `eSVD` never sees an all-zero gene.

Removing all-zero genes changes the factorization of the retained genes,
because an all-zero column is not a zero column of the residual matrix
the initialization decomposes. Results for the retained genes therefore
differ from a run that kept the all-zero genes; they equal a run of
`eSVD` on the filtered object exactly.

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

  esvd_obj <- eSVD_helper(batch_var_prefix = NULL,
                          case_control_levels = c("0", "1"),
                          case_control_var = "CC",
                          categorical_vars = "Sex",
                          id_var = "Individual",
                          numerical_vars = "Age",
                          seurat_obj = seurat_obj,
                          k = 2,
                          max_iter = 5)
  table(esvd_obj$gene_status)
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
