# Internal constructor for a distribution family

Internal constructor for a distribution family

## Usage

``` r
esvd_family(family)
```

## Arguments

- family:

  One of `"gaussian"`, `"curved_gaussian"`, `"exponential"`,
  `"poisson"`, `"neg_binom"`, `"neg_binom2"` or `"bernoulli"`.

## Value

An object of class `esvd_family`: a list with `name`, `domain` (the open
interval of admissible natural parameters), `feas_always`, the external
pointer `internal` used by the C++ objective, and the functions
`feasibility`, `dat_to_nat` and `nat_to_canon`.
