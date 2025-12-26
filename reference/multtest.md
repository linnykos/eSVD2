# Perform multiple-testing adjustment using Efron's empirical null

Perform multiple-testing adjustment using Efron's empirical null

## Usage

``` r
multtest(teststat_vec, observed_quantile = c(0.05, 0.95))
```

## Arguments

- teststat_vec:

  a vector of Gaussian-like test statistics

- observed_quantile:

  The quantiles of `teststat_vec` used to estimate the null distribution

## Value

a list with `fdr_vec`, `logpvalue_vec`, `method`, `null_mean`,
`null_sd`, and `pvalue_vec`
