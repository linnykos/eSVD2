# Compute p-values

Compute p-values

## Usage

``` r
compute_pvalue(input_obj, verbose = 0, ...)
```

## Arguments

- input_obj:

  `eSVD` object outputed from `compute_test_statistic`.

- verbose:

  Integer.

- ...:

  Additional parameters.

## Value

`eSVD` object with added element `"pvalue_list"`
