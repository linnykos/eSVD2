# Compute test statistics for eSVD object

Compute test statistics for eSVD object

## Usage

``` r
# S3 method for class 'eSVD'
compute_test_statistic(input_obj, verbose = 0, ...)
```

## Arguments

- input_obj:

  `eSVD` object outputed from `compute_posterior.eSVD`.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

`eSVD` object with added element `"teststat_vec"`
