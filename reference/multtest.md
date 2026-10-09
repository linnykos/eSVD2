# Perform multiple-testing adjustment using Efron's empirical null

Estimates the empirical null distribution \\N(\mu_0, \sigma_0^2)\\ of
the Gaussianized test statistics and converts every statistic to a
two-sided p-value under it, followed by a Benjamini-Hochberg adjustment.

## Usage

``` r
multtest(teststat_vec, observed_quantile = c(0.05, 0.95))
```

## Arguments

- teststat_vec:

  Named vector of finite Gaussian-like test statistics, one per gene.

- observed_quantile:

  Two numbers in `(0, 1)`, the lower and upper quantiles of
  `teststat_vec` that delimit the central subset used by the two
  fallback estimators.

## Value

A list with `fdr_vec` (BH-adjusted p-values), `logpvalue_vec`
(\\-\log\_{10}\\ of the two-sided p-values), `method` (which estimator
produced the null), `null_mean`, `null_sd`, and `pvalue_vec` (the
two-sided p-values). All vectors carry `names(teststat_vec)`.

## Details

Three estimators of the null are tried in order, each only when the
previous one fails:
[`locfdr::locfdr`](https://rdrr.io/pkg/locfdr/man/locfdr.html)
(`method = "locfdr"`), the truncated-Gaussian maximum-likelihood
estimator of Efron (2007, Equation 4.12) fitted to the statistics
between the two `observed_quantile` quantiles
(`method = "truncated_mle"`), and a moment-matching estimator on the
same central subset (`method = "simple"`). A warning is raised whenever
the chain falls below `locfdr`, so a degraded null is never silent; the
estimator that ran is returned as `method`.

## Examples

``` r
set.seed(10)
# 300 null statistics from N(0.2, 1.1^2) plus a few real signals
teststat_vec <- c(stats::rnorm(300, mean = 0.2, sd = 1.1), 6, -7, 8)
names(teststat_vec) <- paste0("gene", seq_along(teststat_vec))
res <- multtest(teststat_vec)
res$method
#> [1] "locfdr"
c(res$null_mean, res$null_sd)
#> [1] 0.1050933 1.0976984
which(res$fdr_vec < 0.05)
#> gene301 gene302 gene303 
#>     301     302     303 
```
