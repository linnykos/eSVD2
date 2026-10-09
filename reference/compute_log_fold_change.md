# Compute the log2 fold change and its standard error

Reports, for every gene, the log2 fold change between the case and the
control individuals and a standard error for it, on the scale and with
the unit of replication that 'DESeq2', 'dreamlet' and 'NEBULA' report
theirs. `compute_test_statistic` and `compute_test_per_gene` call this
for you, so the two vectors are part of the ordinary output of `eSVD`;
it is exported for objects on which they need to be recomputed.

## Usage

``` r
compute_log_fold_change(input_obj, verbose = 0)
```

## Arguments

- input_obj:

  `eSVD` object output from `compute_test_statistic` or
  `compute_test_per_gene`. It must carry `case_mean`, `control_mean`,
  `case_var` and `control_var` (named numeric vectors, one entry per
  gene) and, in `param`, the individuals of each arm. The count matrix
  and the posterior matrices are not needed, so an object built with
  `eSVD(bool_diet = TRUE)` is accepted.

- verbose:

  Integer; `0` is silent.

## Value

The `eSVD` object with two named numeric vectors added or replaced, one
entry per gene in the order of `case_mean`: `log2fc_vec`, the log2 fold
change (positive when expression is higher in the cases), and
`log2fc_se_vec`, its standard error. A gene whose `case_mean` or
`control_mean` is `NA`, as for the all-zero genes reinserted by
`eSVD_helper`, is `NA` in both. A gene whose mean is not positive is
also `NA` in both, with a warning.

## Details

**The estimate.** Write \\\bar m\_{1j}\\ and \\\bar m\_{0j}\\ for
`case_mean` and `control_mean`: gene \\j\\'s posterior mean expression,
averaged over the cells of each individual and then over the individuals
of the arm. The fold change is \$\$\widehat{\beta}\_j = \log_2(\bar
m\_{1j} / \bar m\_{0j}).\$\$

**The standard error.** Write \\v\_{1j}\\ and \\v\_{0j}\\ for `case_var`
and `control_var`, the variances the Welch statistic divides by, and
\\n_1\\, \\n_0\\ for the number of case and control *individuals*.
Treating the two arms as independent, the delta method gives
\$\$\widehat{\mathrm{SE}}(\widehat{\beta}\_j) = \frac{1}{\ln 2}
\sqrt{\frac{v\_{1j}}{n_1 \bar m\_{1j}^2} + \frac{v\_{0j}}{n_0 \bar
m\_{0j}^2}}.\$\$ Multiply by `log(2)` for the standard error of the
natural-log fold change, which is the scale 'NEBULA' reports on. The
standard error of the difference \\\bar m\_{1j} - \bar m\_{0j}\\ itself,
the denominator of `teststat_vec`, is
`sqrt(case_var / n1 + control_var / n0)`.

**What it is a standard error of.** The target is the log2 ratio of the
arm-level mean expression, each arm mean being an average over
individuals of their covariate- and depth-adjusted expression. The
variance is divided by the number of individuals and never by the number
of cells, so sequencing more cells from the same individuals does not
shrink it. That is the property it shares with a pseudobulk analysis in
'DESeq2' or 'dreamlet' and with the standard error 'NEBULA' reports for
a subject-level predictor. Things to know when setting the numbers side
by side:

- It is a ratio of arithmetic means, as in 'DESeq2' and 'NEBULA'.
  'dreamlet' reports a difference of mean log expression, which is a
  different quantity whenever the two arms differ in spread.

- Each arm's variance is the variance of a mixture over individuals,
  \\v_g = \mathrm{mean}\_i(\bar V_i) + \mathrm{var}\_i(m_i)\\, where
  \\m_i\\ and \\\bar V_i\\ are individual \\i\\'s average posterior mean
  and variance. The first term is an average per-cell posterior
  variance, is not divided by the number of cells, and is usually most
  of \\v_g\\. The standard error is therefore larger, often several
  times larger, than the standard deviation of the fold change when
  individuals are resampled and the fit is held fixed.

- It treats the fitted factorization and the nuisance parameters as
  known, as 'DESeq2' treats its dispersions. This matters most for a
  gene whose estimated nuisance rate (`nuisance_vec`) is very large. The
  posterior of such a gene concentrates on the fitted prediction, the
  first term of \\v_g\\ tends to zero, and the standard error reflects
  only how much the fitted predictions differ between individuals. It
  can then be much smaller than the error in the fold change. Inspect
  `nuisance_vec` before relying on an unusually small standard error.

- Expression is adjusted for each cell's observed sequencing depth, so
  the fold change is relative to total expression. If the genes that
  differ between the arms make up a large share of the counts, every
  gene's fold change is shifted by the same amount in the opposite
  direction.

- `log2fc_vec / log2fc_se_vec` does not reproduce the p-value in
  `pvalue_list`. That p-value comes from the statistic on the linear
  scale, its Welch-Satterthwaite degrees of freedom, and an empirical
  null fitted across genes.

The fold change is a ratio of posterior means, which are shrunk toward
the fitted low-rank prediction; it is therefore a shrunken estimate,
closer in kind to the output of `lfcShrink` in 'DESeq2' than to an
unshrunken maximum likelihood estimate.

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
# compute_test_statistic has already stored both vectors
utils::head(cbind(log2fc = esvd_obj$log2fc_vec,
                  se = esvd_obj$log2fc_se_vec))
#>           log2fc        se
#> gene_1 -1.315935 0.2917377
#> gene_2 -1.752031 0.4640375
#> gene_3  1.295652 0.9232678
#> gene_4  1.636802 0.3983898
#> gene_5 -1.413173 0.5682952
#> gene_6  1.597469 0.2859912
# the standard error on the natural-log scale
utils::head(esvd_obj$log2fc_se_vec * log(2))
#>    gene_1    gene_2    gene_3    gene_4    gene_5    gene_6 
#> 0.2022172 0.3216463 0.6399605 0.2761428 0.3939122 0.1982340 
esvd_obj <- compute_log_fold_change(input_obj = esvd_obj)
```
