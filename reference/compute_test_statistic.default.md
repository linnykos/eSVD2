# Compute test statistics for matrices

Compute test statistics for matrices

## Usage

``` r
# Default S3 method
compute_test_statistic(
  input_obj,
  posterior_var_mat,
  case_individuals,
  control_individuals,
  individual_vec,
  min_cells_per_individual = 3,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  Posterior mean matrix (a `matrix`) where the \\n\\ rows represent
  cells and \\p\\ columns represent genes. The rows and columns of the
  matrix should be named.

- posterior_var_mat:

  Posterior variance matrix (a `matrix`) where the \\n\\ rows represent
  cells and \\p\\ columns represent genes. The rows and columns of the
  matrix should be the same as those in `input_obj`.

- case_individuals:

  Vector of strings: the values of `individual_vec` that are the case
  individuals.

- control_individuals:

  Vector of strings: the values of `individual_vec` that are the control
  individuals.

- individual_vec:

  Vector of strings of length \\n\\ (i.e., the number of cells) that
  denote which cell originates from which individual.

- min_cells_per_individual:

  Minimum number of cells an individual must contribute for the test to
  be computed. Individuals with fewer are expected to have been dropped
  upstream by `eSVD_helper`; reaching here with one is an error rather
  than a silently noisy statistic. Set to `0` to disable the check.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

A list of named numeric vectors, each of length `ncol(input_obj)` and
carrying `colnames(input_obj)`: `teststat_vec`, `case_mean`,
`control_mean`, `case_var`, `control_var`, `log2fc_se_vec` and
`log2fc_vec`. The two variances are those of the mixture over
individuals that summarizes each arm, \\\mathrm{mean}\_i(\bar V_i) +
\mathrm{var}\_i(m_i)\\, where \\m_i\\ and \\\bar V_i\\ are individual
\\i\\'s average posterior mean and variance and the variance over
individuals divides by their number, with no Bessel correction. The
statistic is
`(case_mean - control_mean) / sqrt(case_var / n1 + control_var / n0)`
for `n1` case and `n0` control individuals. See
`compute_log_fold_change` for the last two.
