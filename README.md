---
output: github_document
---

# Purpose

This repository contains all the functions to perform eSVD-DE (for package `eSVD2`, version `1.0.0` as of June 16, 2024) and the downstream analysis, for the paper "eSVD-DE: cohort-wide differential expression in single-cell RNA-seq data using exponential-family embeddings". See the companion GitHub package https://github.com/linnykos/eSVD2_examples for all the analyses performed in the paper. (Note, the original analysis was performed on eSVD2 version 0.0.0.0071, and we are working to make the code in https://github.com/linnykos/eSVD2_examples to be compatibile with the latest version of eSVD2.)

This code was developed and tested primarily on R 4.3.2. on a 2023 Macbook Pro (macOS Sonoma 14.2.1) equipped with Apple M2 Max processor (32 Gb RAM).

<!-- badges: start -->
[![DOI:10.1186/s12859-024-05724-7](https://img.shields.io/badge/doi-10.1186/s12859--024--05724--7-firebrick.svg)](https://doi.org/10.1186/s12859-024-05724-7)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->
  
# Tutorials and vignettes

Please see https://linnykos.github.io/eSVD2/index.html for in-depth tutorials and vignettes on how to use eSVD-DE.

# Installation

This package can be installed through `devtools` in R.

```R
library("devtools")
devtools::install_github("linnykos/eSVD2")
```
The package itself depends on several packages. These include `irlba`, `glmnet`, `locfdr`, `Matrix`, `matrixStats`, `Rmpfr`, `RSpectra`, and `sparseMatrixStats`. See the last section of this README to see where (i.e., CRAN, Bioconductor, or GitHub) to download all such packages. We have noted that `Rmpfr` is sometimes tricky to install due to its required C++ libraries.

After installation of all the dependencies, the installation of the `eSVD2` package itself takes modest time (less than 10 minutes). The installation time mainly consists of time to compile the C++ code since the matrix factorization optimization was written using Rcpp for faster performance.

<details>
<summary>**Known installation issues and the solutions**</summary>

(Solution posted on December 24, 2023): If you come across the error,
```R
Error in irlba::irlba() : 
  function 'as_cholmod_sparse' not provided by package 'Matrix'
```
then it is likely you need to downgrade your version of `Matrix` to `1.6-1.1`. See https://github.com/bwlewis/irlba/issues/70#issuecomment-1826900769. Hence, in the R console,
```R
remove.packages("Seurat")
remove.packages("SeuratObject")
remove.packages("Matrix")
remotes::install_version("Matrix", version = "1.6-1.1")
remotes::install_version("SeuratObject", version = "5.0.0")
remotes::install_version("Seurat", version = "5.0.0")
```
</details>

# Small simulated dataset to demo the software

See https://linnykos.github.io/eSVD2/articles/eSVD2.html for the small demo on how to use eSVD2.

# Setup

The following shows the suggested package versions that the developer (GitHub username: linnykos) used when developing the eSVD2 package.

```R
> devtools::session_info()
─ Session info ──────────────────────────
 setting  value
 version  R version 4.3.2 (2023-10-31)
 os       macOS Sonoma 14.2.1
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/Chicago
 date     2024-06-16
 rstudio  2024.04.1+748 Chocolate Cosmos (desktop)
 pandoc   NA

─ Packages ──────────────────────────────
 package     * version date (UTC) lib source
 brio          1.1.5   2024-04-24 [1] CRAN (R 4.3.1)
 cachem        1.1.0   2024-05-16 [1] CRAN (R 4.3.3)
 cli           3.6.2   2023-12-11 [1] CRAN (R 4.3.1)
 devtools      2.4.5   2022-10-11 [1] CRAN (R 4.3.0)
 digest        0.6.35  2024-03-11 [1] CRAN (R 4.3.1)
 ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.3.0)
 eSVD2       * 1.0.0   2024-06-17 [1] local
 fastmap       1.2.0   2024-05-15 [1] CRAN (R 4.3.3)
 fs            1.6.4   2024-04-25 [1] CRAN (R 4.3.1)
 glue          1.7.0   2024-01-09 [1] CRAN (R 4.3.1)
 htmltools     0.5.8.1 2024-04-04 [1] CRAN (R 4.3.1)
 htmlwidgets   1.6.4   2023-12-06 [1] CRAN (R 4.3.1)
 httpuv        1.6.15  2024-03-26 [1] CRAN (R 4.3.1)
 later         1.3.2   2023-12-06 [1] CRAN (R 4.3.1)
 lifecycle     1.0.4   2023-11-07 [1] CRAN (R 4.3.1)
 magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.3.0)
 memoise       2.0.1   2021-11-26 [1] CRAN (R 4.3.0)
 mime          0.12    2021-09-28 [1] CRAN (R 4.3.0)
 miniUI        0.1.1.1 2018-05-18 [1] CRAN (R 4.3.0)
 pkgbuild      1.4.4   2024-03-17 [1] CRAN (R 4.3.1)
 pkgload       1.3.4   2024-01-16 [1] CRAN (R 4.3.1)
 profvis       0.3.8   2023-05-02 [1] CRAN (R 4.3.0)
 promises      1.3.0   2024-04-05 [1] CRAN (R 4.3.1)
 purrr         1.0.2   2023-08-10 [1] CRAN (R 4.3.0)
 R6            2.5.1   2021-08-19 [1] CRAN (R 4.3.0)
 Rcpp        * 1.0.12  2024-01-09 [1] CRAN (R 4.3.1)
 remotes       2.5.0   2024-03-17 [1] CRAN (R 4.3.1)
 rlang         1.1.4   2024-06-04 [1] CRAN (R 4.3.3)
 rstudioapi    0.16.0  2024-03-24 [1] CRAN (R 4.3.1)
 sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.3.0)
 shiny         1.8.1.1 2024-04-02 [1] CRAN (R 4.3.1)
 stringi       1.8.4   2024-05-06 [1] CRAN (R 4.3.1)
 stringr       1.5.1   2023-11-14 [1] CRAN (R 4.3.1)
 testthat    * 3.2.1.1 2024-04-14 [1] CRAN (R 4.3.1)
 urlchecker    1.0.1   2021-11-30 [1] CRAN (R 4.3.0)
 usethis       2.2.3   2024-02-19 [1] CRAN (R 4.3.1)
 vctrs         0.6.5   2023-12-01 [1] CRAN (R 4.3.1)
 xtable        1.8-4   2019-04-21 [1] CRAN (R 4.3.0)

 [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
```
