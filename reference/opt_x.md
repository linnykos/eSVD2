# Optimize X given C, Y and Z

Optimize X given C, Y and Z

## Usage

``` r
opt_x(
  XC_init,
  YZ,
  k,
  loader,
  family,
  s,
  gamma,
  l2penx,
  verbose = 0,
  inplace = FALSE,
  ...
)
```

## Arguments

- XC_init:

  Initial value for the \`\[X C\]\` matrix, X \[n x k\], C \[n x r\]

- YZ:

  The \`\[Y Z\]\` matrix, Y \[p x k\], Z \[p x r\]

- k:

  Number of columns in X and Y

- loader:

  The data loader, typically returned by data_loader()

- family:

  A family object, typically returned by esvd_family()

- s:

  The library size vector, \[n x 1\]

- gamma:

  The nuisance parameter vector, \[p x 1\]

- l2penx:

  The l2 penalty parameter for X, a scalar

- verbose:

  Verbosity parameter

- inplace:

  Whether the input XC_init will be modified and returned

- ...:

  Additional parameters
