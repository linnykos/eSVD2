# Optimize Y and Z given X and C

Optimize Y and Z given X and C

## Usage

``` r
opt_yz(
  YZ_init,
  XC,
  k,
  fixed_cols,
  loader,
  family,
  s,
  gamma,
  l2peny,
  l2penz,
  verbose = 0,
  inplace = FALSE,
  ...
)
```

## Arguments

- YZ_init:

  Initial value for the \`\[Y Z\]\` matrix, Y \[p x k\], Z \[p x r\]

- XC:

  The \`\[X C\]\` matrix, X \[n x k\], C \[n x r\]

- k:

  Number of columns in X and Y

- fixed_cols:

  Which columns in YZ need to be fixed

- loader:

  The data loader, typically returned by data_loader()

- family:

  A family object, typically returned by esvd_family()

- s:

  The library size vector, \[n x 1\]

- gamma:

  The nuisance parameter vector, \[p x 1\]

- l2peny:

  The l2 penalty parameter for Y, a scalar

- l2penz:

  The l2 penalty parameter for Z, a scalar

- verbose:

  Verbosity parameter

- inplace:

  Whether the input XC_init will be modified and returned

- ...:

  Additional parameters
