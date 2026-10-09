# Estimate nuisance values for eSVD objects (i.e., over-dispersion)

Assumes a Gamma-Poisson model where the mean and variance are
proportionally related: \\\lambda\_{ji} \sim
\mathrm{Gamma}(\mathrm{mean} = \mu\_{ji}, \mathrm{var} = \gamma_j
\mu\_{ji})\\. The estimated `nuisance_vec` is the Gamma *rate* \\\beta_j
= 1/\gamma_j\\, the reciprocal of the over-dispersion \\\gamma_j\\ in
Lin, Qiu and Roeder (2024), so a **larger** value means **less**
over-dispersion.

## Usage

``` r
# S3 method for class 'eSVD'
estimate_nuisance(
  input_obj,
  bool_covariates_as_library = TRUE,
  bool_library_includes_interept = TRUE,
  bool_use_log = FALSE,
  cap_multiplier = 10,
  min_val = 1e-04,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  `eSVD` object output from `opt_esvd.eSVD`. Specifically, the nuisance
  parameters will be estimated based on the fit in
  `input_obj[[input_obj[["latest_Fit"]]]]`.

- bool_covariates_as_library:

  Boolean to include the donor covariates (every covariate other than
  the case-control variable) in the covariate-adjusted library size,
  default is `TRUE`, as in `compute_posterior.eSVD` and
  `compute_test_per_gene` and as `eSVD` passes it. Pass the same value
  to all of them, or the rate and the posterior use different libraries.

- bool_library_includes_interept:

  Boolean if the intercept term from the eSVD matrix factorization
  should be included in the calculation for the covariate-adjusted
  library size, default is `TRUE`.

- bool_use_log:

  Boolean if the nuisance (i.e., over-dispersion) parameter should be
  estimated on the log scale, default is `FALSE`.

- cap_multiplier:

  One positive number, default `10`: no gene's rate may exceed this many
  times the median over cells of its library size. `Inf` means no cap.

- min_val:

  Minimum value of the rate, as a multiple of the median over cells of
  the gene's library size (the units of `cap_multiplier`); must be below
  `cap_multiplier`. Default `1e-4`.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

`eSVD` object with the following named vectors, one entry per gene,
added to the list in `input_obj[[input_obj[["latest_Fit"]]]]`:
`nuisance_vec` (the rate after the floor and the cap, which is what
later stages read; `compute_posterior` still floors it at a quantile and
may rescale it), `nuisance_mle_vec` (the rate before the cap but after
the floor, so for a `failed` gene it is the floor and not an estimate),
`nuisance_library_median_vec` (the median over cells of the gene's
library size, so the cap is `cap_multiplier` times this),
`nuisance_status` (a factor with levels `estimated`, `capped`,
`boundary`, `failed`; see above), `gene_mean_count_vec` (the gene's mean
count over cells) and `gene_sparsity_vec` (the fraction of cells in
which its count is zero). `input_obj$param` records
`nuisance_cap_multiplier`, `nuisance_num_capped` (the number of genes
whose rate the cap replaced, the `boundary` genes among them included),
`nuisance_num_boundary` and `nuisance_num_failed`; calling the function
again overwrites them. `plot_nuisance` draws these.

## Details

**The cap.** Each gene's rate is the maximum-likelihood estimate,
bounded above: \\\hat\beta_j = \min\\\hat\beta_j^{\mathrm{MLE}},\\ c
\cdot \mathrm{median}\_i\\ \ell\_{ji}\\\\, with \\c\\ = `cap_multiplier`
and \\\ell\_{ji}\\ the covariate-adjusted library size of gene \\j\\ in
cell \\i\\. The rate has the units of the library size, which is why the
bound is a multiple of it. The bound is a calibration device and not an
estimate of the smallest over-dispersion in the data: a gene whose
counts are no more variable around the fit than a Poisson has no finite
maximum-likelihood rate, its posterior would collapse onto the fit, and
its test statistic would be inflated. The bound limits how much more
weight the posterior may give the fit for one gene than for the typical
gene. `cap_multiplier = Inf` removes it and gives the rates of version
1.1.0, except that the floor `min_val` is now a multiple of the library
size rather than an absolute number. `cap_multiplier = 1` is close to,
but not the same as, the bound of the version that accompanied Lin, Qiu
and Roeder (2024), which was the largest library size of the gene and
not a multiple of its median. `recompute_pvalue` changes the multiplier
on a fitted object without estimating anything again.

**The floor** `min_val` is in the same units: no gene's rate is below
`min_val` times the median of its library size. On the unit-free scale
of `plot_nuisance` every rate therefore lies between `min_val` and
`cap_multiplier`, and no gene is above its cap. `min_val` must be below
`cap_multiplier`.

**The status** of each gene records what happened to it:

- `estimated`:

  the gene keeps its maximum-likelihood rate.

- `capped`:

  a finite maximum-likelihood rate exists and is above the bound.

- `boundary`:

  no finite maximum-likelihood rate exists: the score statistic for
  over-dispersion, \\D_j = \sum_i \[(A\_{ji} - \ell\_{ji}\mu\_{ji})^2 -
  A\_{ji}\] / \mu\_{ji}\\, is not positive. Such a gene is set to the
  bound, whatever value the iterations of the estimate stopped at (which
  is what `nuisance_mle_vec` holds for it). With `cap_multiplier = Inf`
  it keeps that value.

- `failed`:

  both estimation routes failed and the rate is the floor, `min_val`
  times the median of its library size.

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
summary(esvd_obj$fit_First$nuisance_vec)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 3.446e-04 3.449e-01 1.358e+00 2.718e+00 3.305e+00 1.329e+01 
table(esvd_obj$fit_First$nuisance_status)
#> 
#> estimated    capped  boundary    failed 
#>        37         1         2         0 
esvd_obj$param$nuisance_num_capped
#> [1] 3
```
