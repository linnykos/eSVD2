# Generate data

Generate data

## Usage

``` r
generate_data(
  nat_mat,
  family,
  nuisance_param_vec = NA,
  library_size_vec = 1,
  tol = 0.001
)
```

## Arguments

- nat_mat:

  An \\n\times p\\ matrix of natural parameters, where \\n\\ rows
  represent cells and \\p\\ columns represent genes.

- family:

  A character string, one of `"gaussian"`, `"exponential"`, `"poisson"`,
  `"neg_binom"`, `"curved_gaussian"`, and `"bernoulli"`.

- nuisance_param_vec:

  Either `NA` or a single numeric or a length-\\p\\ vector of numerics
  representing nuisance parameters (for `family = "neg_binom"` and
  `family = "curved_gausian"`). It is only required if
  `family %in% c("neg_binom", "curved_gaussian")`.

- library_size_vec:

  Either `NA` or a length-\\n\\ vector of numerics

- tol:

  Small positive value to determine the smallest possible value in the
  output matrix, useful for only `family = "curved_gaussian"`.

## Value

The generated data matrix
