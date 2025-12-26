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

  Vector of strings representing the individuals in
  `metadata[,covariate_individual]` that are the case individuals.

- control_individuals:

  Vector of strings representing the individuals in
  `metadata[,covariate_individual]` that are the control individuals.

- individual_vec:

  Vector of strings of length \\n\\ (i.e., the number of cells) that
  denote which cell originates from which individual.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

A vector of test statistics of length `ncol(input_obj)`
