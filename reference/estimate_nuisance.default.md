# Estimate nuisance values for matrix or sparse matrices.

Each gene's value is the maximum-likelihood Gamma *rate* \\\beta_j =
1/\gamma_j\\ under \\A\_{ji} \mid \lambda\_{ji} \sim
\mathrm{Poisson}(\ell\_{ji} \lambda\_{ji})\\, \\\lambda\_{ji} \sim
\mathrm{Gamma}(\mu\_{ji} \beta_j, \beta_j)\\, with \\\mu\\ from
`mean_mat` and \\\ell\\ from `library_mat`. Larger values mean less
over-dispersion. The estimate is then bounded above by `cap_multiplier`
times the median of the gene's column of `library_mat`; see
`estimate_nuisance.eSVD` for why.

## Usage

``` r
# Default S3 method
estimate_nuisance(
  input_obj,
  mean_mat,
  library_mat,
  bool_use_log = FALSE,
  cap_multiplier = 10,
  min_val = 1e-04,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  Dataset (either `matrix` or `dgCMatrix`) where the \\n\\ rows
  represent cells and \\p\\ columns represent genes. The rows and
  columns of the matrix should be named.

- mean_mat:

  A `matrix` of \\n\\ rows and \\p\\ columns that represents the
  expected value of each entry.

- library_mat:

  A `matrix` of \\n\\ rows and \\p\\ columns that represents the library
  size of each entry.

- bool_use_log:

  Boolean if the nuisance (i.e., over-dispersion) parameter should be
  estimated on the log scale, default is `FALSE`.

- cap_multiplier:

  One positive number, default `10`: no gene's rate may exceed this many
  times the median of its column of `library_mat`. `Inf` means no cap.

- min_val:

  Minimum value of the rate, as a multiple of the median over cells of
  the gene's library size (the units of `cap_multiplier`); must be below
  `cap_multiplier`. Default `1e-4`.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

Numeric vector of length \\p\\ (named by `colnames(input_obj)` when
present) of Gamma rates, after the cap. A gene whose estimation fails on
both routes gets the floor, `min_val` times the median of its library
size, and a warning reports how many genes did. The rates before the cap
and each gene's status are returned by the `eSVD` method only.
