# Redo the test at another cap on the nuisance rate, without fitting again

The factorization (`x_mat`, `y_mat`, `z_mat`) is the expensive part of
an analysis and does not depend on the cap on the nuisance rate; see
`estimate_nuisance`. This function takes a fitted and tested `eSVD`
object, applies another `cap_multiplier` to the maximum-likelihood rates
stored on it, and repeats what comes after the rates: the posterior, the
test statistic, the log2 fold change and the p-values. Nothing is
estimated again, and every other setting is the one recorded in
`input_obj$param` when the object was built.

## Usage

``` r
recompute_pvalue(input_obj, cap_multiplier, seurat_obj = NULL, verbose = 0)
```

## Arguments

- input_obj:

  `eSVD` object that has been tested, by `eSVD`, `eSVD_helper`,
  `compute_test_per_gene` or `compute_pvalue`, with version 1.2.0 or
  later.

- cap_multiplier:

  One positive number: no gene's nuisance rate may exceed this many
  times the median of its library size. `Inf` means no cap. See
  `estimate_nuisance`.

- seurat_obj:

  `NULL`, or the `Seurat` object the analysis was run on. Required when
  `input_obj` has no counts (`input_obj$dat`) or covariates, and ignored
  otherwise.

- verbose:

  Integer; `0` is silent.

## Value

The `eSVD` object with `teststat_vec`, `case_mean`, `control_mean`,
`case_var`, `control_var`, `log2fc_vec`, `log2fc_se_vec` and
`pvalue_list` replaced, and, in the fit named by `latest_Fit`,
`nuisance_vec`, `nuisance_status` and the posterior matrices (if the
object had them) replaced. `param$nuisance_cap_multiplier`,
`param$nuisance_num_capped` and `param$nuisance_num_boundary` describe
the new cap. The factorization and `nuisance_mle_vec` are unchanged, so
the function can be called again at any other cap. An object that had no
counts is returned without them.

## Details

The posterior needs the counts and the covariates. An object built with
`eSVD(bool_diet = FALSE)`, or stage by stage, carries them. An object
built with `bool_diet = TRUE` (the default of `eSVD` and `eSVD_helper`)
does not, so pass `seurat_obj`, the Seurat object the analysis was run
on. The cells and genes of the fit are looked up in it **by name**, so
it may hold more cells and genes than were analyzed (the individuals and
the all-zero genes `eSVD_helper` removed), and the covariates are
rebuilt from its `meta.data` with the variable names `eSVD` recorded.
Summaries of the rebuilt counts and covariates are compared with those
recorded at the fit, and the function stops if they differ, for example
because the `meta.data` was edited since.

The result is what `eSVD` would have returned had it been called with
this `cap_multiplier`.
