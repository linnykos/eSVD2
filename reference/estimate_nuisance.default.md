# Estimate nuisance values for matrix or sparse matrices.

Estimate nuisance values for matrix or sparse matrices.

## Usage

``` r
# Default S3 method
estimate_nuisance(
  input_obj,
  mean_mat,
  library_mat,
  bool_use_log = F,
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

- min_val:

  Minimum value of the nuisance parameter.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

Vector of length \\p\\
