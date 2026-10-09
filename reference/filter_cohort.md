# Filter a cohort before running eSVD-DE

Applies the cohort-level filters that `eSVD_helper` runs before `eSVD`,
in this order:

1.  **Drop** every individual with fewer than `min_cells_per_id` cells,
    with a warning naming them.

2.  **Check** the filtered cohort against `min_cells`, `min_ids`,
    `min_ids_per_arm` and `min_cells_casecontrol`; on the first failure,
    warn and return `NA`.

Dropping first means every threshold describes the cohort that will
actually be analyzed.

## Usage

``` r
filter_cohort(
  seurat_obj,
  id_var,
  case_control_var,
  min_cells = 20,
  min_cells_casecontrol = 20,
  min_cells_per_id = 3,
  min_ids = 4,
  min_ids_per_arm = 2,
  verbose = 0
)
```

## Arguments

- seurat_obj:

  A `Seurat` object.

- id_var:

  Column of `seurat_obj@meta.data` naming each cell's individual.

- case_control_var:

  Column of `seurat_obj@meta.data` holding the case-control status (two
  levels).

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

## Value

Either `NA` (the cohort was rejected; see the warning), or a list with
`dropped_individuals` (character vector, possibly empty) and
`seurat_obj` (the object with those individuals' cells removed).

## Details

`min_cells`, `min_ids` and `min_cells_casecontrol` reject at `<=`, and
`min_cells_per_id` and `min_ids_per_arm` at `<`, so setting one to `0`
disables it. `min_cells_per_id` is the one exception: a value of `1` or
`2` would not disable the drop but weaken it, and a surviving one-cell
individual has no within-individual variance, so it must be `0` (off) or
at least `3`.

The four rejections return the same bare `NA`, so the warning text is
what tells a caller which filter fired; each names its own threshold.

## Examples

``` r
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

  # drop two cells from one individual: with min_cells_per_id = 15 that
  # individual is removed, with a warning naming it
  seurat_small <- seurat_obj[, -(1:2)]
  res <- filter_cohort(seurat_obj = seurat_small,
                       id_var = "Individual",
                       case_control_var = "CC",
                       min_cells_per_id = 15)
  res$dropped_individuals
  ncol(res$seurat_obj)
}
#> Warning: dropping 1 individual(s) with fewer than `min_cells_per_id` = 15 cells: indiv_1 (13 cells)
#> [1] 105
```
