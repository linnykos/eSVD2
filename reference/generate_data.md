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
  `"neg_binom"`, `"neg_binom2"`, `"curved_gaussian"`, and `"bernoulli"`.

- nuisance_param_vec:

  Either `NA` or a single numeric or a length-\\p\\ vector of numerics
  representing nuisance parameters. It is required if `family` is
  `"gaussian"`, `"curved_gaussian"`, `"neg_binom"` or `"neg_binom2"`,
  and ignored otherwise.

- library_size_vec:

  A single positive numeric (recycled to every cell; default `1`) or a
  length-\\n\\ vector of positive numerics.

- tol:

  Small positive value to determine the smallest possible value in the
  output matrix, applied only when `family` is `"gaussian"` or
  `"curved_gaussian"`.

## Value

The generated data matrix, of the same dimension as `nat_mat`

## Examples

``` r
set.seed(10)
x_mat <- matrix(stats::rnorm(20 * 2), nrow = 20, ncol = 2)
y_mat <- matrix(stats::rnorm(10 * 2), nrow = 10, ncol = 2)
nat_mat <- tcrossprod(x_mat, y_mat)
dat <- generate_data(nat_mat = nat_mat, family = "poisson")
dim(dat)
#> [1] 20 10
```
