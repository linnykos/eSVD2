# Function to reparameterize two matrices

test

## Usage

``` r
.reparameterize(x_mat, y_mat, equal_covariance)
```

## Arguments

- x_mat:

  matrix of dimension `n` by `k`

- y_mat:

  matrix of dimension `p` by `k`

- equal_covariance:

  boolean

## Value

list of two matrices

## Details

Designed to output matrices of the same dimension as `x_mat` and
`y_mat`, but linearly transformed so `x_mat %*% t(y_mat)` is preserved
but either `x_mat %*% t(x_mat)` is diagonal and equal to
`y_mat %*% t(y_mat)` (if `equal_covariance` is `FALSE`) or
`x_mat %*% t(x_mat)/nrow(x_mat)` is diagonal and equal to
`y_mat %*% t(y_mat)/nrow(y_mat)` (if `equal_covariance` is `TRUE`)
