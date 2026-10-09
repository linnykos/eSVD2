---
output: github_document
---

# Purpose

This repository contains the `eSVD2` R package, which implements eSVD-DE, the method of the paper "eSVD-DE: cohort-wide differential expression in single-cell RNA-seq data using exponential-family embeddings" (Lin, Qiu and Roeder, *BMC Bioinformatics*, 2024). See the companion repository https://github.com/linnykos/eSVD2_examples for all the analyses performed in the paper. (Note, the original analysis was performed on eSVD2 version 0.0.0.0071, and we are working to make the code in https://github.com/linnykos/eSVD2_examples compatible with the latest version of eSVD2.)

The current version is `1.2.0`; see `NEWS.md` for what changed. This code was developed and tested primarily on R 4.3 and 4.5 on a 2023 Macbook Pro (macOS Sonoma 14.2.1) equipped with Apple M2 Max processor (32 Gb RAM).

<!-- badges: start -->
[![DOI:10.1186/s12859-024-05724-7](https://img.shields.io/badge/doi-10.1186/s12859--024--05724--7-firebrick.svg)](https://doi.org/10.1186/s12859-024-05724-7)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->
  
# Tutorials and vignettes

Please see https://linnykos.github.io/eSVD2/index.html for in-depth tutorials and vignettes on how to use eSVD-DE.

# Installation

The development version can be installed from GitHub with `remotes` (or `devtools`) in R.

```R
# install.packages("remotes")
remotes::install_github("linnykos/eSVD2")
```
(Once the package is accepted on CRAN, `install.packages("eSVD2")` will install the released version.)

The package itself depends on several packages. These include `irlba`, `glmnet`, `locfdr`, `Matrix`, `matrixStats`, `Rcpp`, and `RSpectra`, all on CRAN. `SeuratObject` is needed by `eSVD()` and `eSVD_helper()` to read the count matrix from a Seurat object; the lower-level functions (`initialize_esvd()`, `opt_esvd()`, ...) work on a plain `matrix` or `dgCMatrix` and do not need it.

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

The following shows the package versions that the developer (GitHub username: linnykos) used in October 2026 to run the vignette `eSVD2.Rmd` with version 1.2.0 of the eSVD2 package (version 1.0.0, in June 2024, was developed under R 4.3.2).

```R
> devtools::session_info()
─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.5.1 (2025-06-13)
 os       macOS Sonoma 14.2.1
 system   aarch64, darwin20
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/Los_Angeles
 date     2026-10-08
 pandoc   NA
 quarto   1.7.32 @ /usr/local/bin/quarto

─ Packages ───────────────────────────────────────────────────────────────────
 package          * version date (UTC) lib source
 abind              1.4-8   2024-09-12 [2] CRAN (R 4.5.0)
 cachem             1.1.0   2024-05-16 [2] CRAN (R 4.5.0)
 cli                3.6.6   2026-04-09 [2] CRAN (R 4.5.2)
 cluster            2.1.8.3 2026-07-30 [2] CRAN (R 4.5.2)
 codetools          0.2-20  2024-03-31 [2] CRAN (R 4.5.1)
 cowplot            1.2.0   2025-07-07 [2] CRAN (R 4.5.0)
 data.table         1.18.4  2026-05-06 [2] CRAN (R 4.5.2)
 deldir             2.0-4   2024-02-28 [2] CRAN (R 4.5.0)
 devtools         * 2.5.2   2026-04-30 [2] CRAN (R 4.5.2)
 dichromat          2.0-1   2026-07-22 [2] CRAN (R 4.5.2)
 digest             0.6.39  2025-11-19 [2] CRAN (R 4.5.2)
 dotCall64          1.2     2024-10-04 [2] CRAN (R 4.5.0)
 dplyr              1.2.1   2026-04-03 [2] CRAN (R 4.5.2)
 ellipsis           0.3.3   2026-04-04 [2] CRAN (R 4.5.2)
 eSVD2            * 1.2.0   2026-10-09 [1] local
 evaluate           1.0.5   2025-08-27 [2] CRAN (R 4.5.0)
 farver             2.1.2   2024-05-13 [2] CRAN (R 4.5.0)
 fastDummies        1.7.6   2026-04-22 [2] CRAN (R 4.5.2)
 fastmap            1.2.0   2024-05-15 [2] CRAN (R 4.5.0)
 fitdistrplus       1.2-6   2026-01-24 [2] CRAN (R 4.5.2)
 foreach            1.5.2   2022-02-02 [2] CRAN (R 4.5.0)
 fs                 2.1.0   2026-04-18 [2] CRAN (R 4.5.2)
 future             1.75.0  2026-07-20 [2] CRAN (R 4.5.2)
 future.apply       1.20.2  2026-02-20 [2] CRAN (R 4.5.2)
 generics           0.1.4   2025-05-09 [2] CRAN (R 4.5.0)
 ggplot2            4.0.3   2026-04-22 [2] CRAN (R 4.5.2)
 ggrepel            0.9.8   2026-03-17 [2] CRAN (R 4.5.2)
 ggridges           0.5.7   2025-08-27 [2] CRAN (R 4.5.0)
 glmnet             5.0     2026-05-04 [2] CRAN (R 4.5.2)
 globals            0.19.1  2026-03-13 [2] CRAN (R 4.5.2)
 glue               1.8.1   2026-04-17 [2] CRAN (R 4.5.2)
 goftest            1.2-3   2021-10-07 [2] CRAN (R 4.5.0)
 gridExtra          2.3.1   2026-06-25 [2] CRAN (R 4.5.2)
 gtable             0.3.6   2024-10-25 [2] CRAN (R 4.5.0)
 htmltools          0.5.9   2025-12-04 [2] CRAN (R 4.5.2)
 htmlwidgets        1.6.4   2023-12-06 [2] CRAN (R 4.5.0)
 httpuv             1.6.17  2026-03-18 [2] CRAN (R 4.5.2)
 httr               1.4.8   2026-02-13 [2] CRAN (R 4.5.2)
 ica                1.0-3   2022-07-08 [2] CRAN (R 4.5.0)
 igraph             2.3.3   2026-06-26 [2] CRAN (R 4.5.2)
 irlba              2.3.7   2026-01-30 [2] CRAN (R 4.5.2)
 iterators          1.0.14  2022-02-05 [2] CRAN (R 4.5.0)
 jsonlite           2.0.0   2025-03-27 [2] CRAN (R 4.5.0)
 KernSmooth         2.23-27 2026-08-12 [2] CRAN (R 4.5.2)
 knitr              1.51    2025-12-20 [2] CRAN (R 4.5.2)
 labeling           0.4.3   2023-08-29 [2] CRAN (R 4.5.0)
 later              1.4.8   2026-03-05 [2] CRAN (R 4.5.2)
 lattice            0.23-1  2026-08-12 [2] CRAN (R 4.5.2)
 lifecycle          1.0.5   2026-01-08 [2] CRAN (R 4.5.2)
 listenv            1.0.0   2026-06-22 [2] CRAN (R 4.5.2)
 lmtest             0.9-40  2022-03-21 [2] CRAN (R 4.5.0)
 locfdr             1.1-8   2015-07-15 [2] CRAN (R 4.5.0)
 magrittr           2.0.5   2026-04-04 [2] CRAN (R 4.5.2)
 MASS               7.3-66  2026-07-15 [2] CRAN (R 4.5.2)
 Matrix             1.7-6   2026-07-25 [2] CRAN (R 4.5.2)
 matrixStats        1.5.0   2025-01-07 [2] CRAN (R 4.5.0)
 memoise            2.0.1   2021-11-26 [2] CRAN (R 4.5.0)
 mime               0.13    2025-03-17 [2] CRAN (R 4.5.0)
 miniUI             0.1.2   2025-04-17 [2] CRAN (R 4.5.0)
 nlme               3.1-170 2026-07-15 [2] CRAN (R 4.5.2)
 otel               0.2.0   2025-08-29 [2] CRAN (R 4.5.0)
 parallelly         1.48.0  2026-06-29 [2] CRAN (R 4.5.2)
 patchwork          1.3.2   2025-08-25 [2] CRAN (R 4.5.0)
 pbapply            1.7-4   2025-07-20 [2] CRAN (R 4.5.0)
 pillar             1.11.1  2025-09-17 [2] CRAN (R 4.5.0)
 pkgbuild           1.4.8   2025-05-26 [2] CRAN (R 4.5.0)
 pkgconfig          2.0.3   2019-09-22 [2] CRAN (R 4.5.0)
 pkgload            1.5.3   2026-06-15 [2] CRAN (R 4.5.2)
 plotly             4.12.1  2026-07-22 [2] CRAN (R 4.5.2)
 plyr               1.8.9   2023-10-02 [2] CRAN (R 4.5.0)
 png                0.1-9   2026-03-15 [2] CRAN (R 4.5.2)
 polyclip           1.10-7  2024-07-23 [2] CRAN (R 4.5.0)
 progressr          1.0.0   2026-07-04 [2] CRAN (R 4.5.2)
 promises           1.5.0   2025-11-01 [2] CRAN (R 4.5.0)
 purrr              1.2.2   2026-04-10 [2] CRAN (R 4.5.2)
 R6                 2.6.1   2025-02-15 [2] CRAN (R 4.5.0)
 RANN               2.6.3   2026-08-24 [2] CRAN (R 4.5.2)
 RColorBrewer       1.1-3   2022-04-03 [2] CRAN (R 4.5.0)
 Rcpp               1.1.2   2026-07-05 [2] CRAN (R 4.5.2)
 RcppAnnoy          0.0.23  2026-01-12 [2] CRAN (R 4.5.2)
 RcppHNSW           0.7.0   2026-05-26 [2] CRAN (R 4.5.2)
 reshape2           1.4.5   2025-11-12 [2] CRAN (R 4.5.0)
 reticulate         1.46.0  2026-04-09 [2] CRAN (R 4.5.2)
 rlang              1.3.0   2026-07-05 [2] CRAN (R 4.5.2)
 ROCR               1.0-12  2026-01-23 [2] CRAN (R 4.5.2)
 RSpectra           0.16-2  2024-07-18 [2] CRAN (R 4.5.0)
 Rtsne              0.17    2023-12-07 [2] CRAN (R 4.5.0)
 S7                 0.2.2   2026-04-22 [2] CRAN (R 4.5.2)
 scales             1.4.0   2025-04-24 [2] CRAN (R 4.5.0)
 scattermore        1.2     2023-06-12 [2] CRAN (R 4.5.0)
 sctransform        0.4.3   2026-01-10 [2] CRAN (R 4.5.2)
 sessioninfo        1.2.4   2026-06-04 [2] CRAN (R 4.5.2)
 Seurat             5.5.1   2026-06-26 [2] CRAN (R 4.5.2)
 SeuratObject       5.4.0   2026-04-11 [2] CRAN (R 4.5.2)
 shape              1.4.6.1 2024-02-23 [2] CRAN (R 4.5.0)
 shiny              1.14.0  2026-06-21 [2] CRAN (R 4.5.2)
 sp                 2.2-3   2026-07-19 [2] CRAN (R 4.5.2)
 spam               2.11-4  2026-05-29 [2] CRAN (R 4.5.2)
 spatstat.data      3.1-9   2025-10-18 [2] CRAN (R 4.5.0)
 spatstat.explore   3.8-2   2026-07-27 [2] CRAN (R 4.5.2)
 spatstat.geom      3.8-2   2026-07-24 [2] CRAN (R 4.5.2)
 spatstat.random    3.5-1   2026-07-27 [2] CRAN (R 4.5.2)
 spatstat.sparse    3.2-0   2026-05-21 [2] CRAN (R 4.5.2)
 spatstat.univar    3.2-0   2026-05-18 [2] CRAN (R 4.5.2)
 spatstat.utils     3.2-4   2026-07-16 [2] CRAN (R 4.5.2)
 stringi            1.8.9   2026-08-04 [2] CRAN (R 4.5.2)
 stringr            1.6.0   2025-11-04 [2] CRAN (R 4.5.0)
 survival           3.8-11  2026-08-21 [2] CRAN (R 4.5.2)
 tensor             1.5.1   2025-06-17 [2] CRAN (R 4.5.0)
 tibble             3.3.1   2026-01-11 [2] CRAN (R 4.5.2)
 tidyr              1.3.2   2025-12-19 [2] CRAN (R 4.5.2)
 tidyselect         1.2.1   2024-03-11 [2] CRAN (R 4.5.0)
 usethis          * 3.2.1   2025-09-06 [2] CRAN (R 4.5.0)
 uwot               0.2.4   2025-11-10 [2] CRAN (R 4.5.0)
 vctrs              0.7.3   2026-04-11 [2] CRAN (R 4.5.2)
 viridisLite        0.4.3   2026-02-04 [2] CRAN (R 4.5.2)
 withr              3.0.3   2026-06-19 [2] CRAN (R 4.5.2)
 xfun               0.60    2026-07-09 [2] CRAN (R 4.5.2)
 xtable             1.8-8   2026-02-22 [2] CRAN (R 4.5.2)
 zoo                1.9-0   2026-07-31 [2] CRAN (R 4.5.2)

 [1] <local library holding eSVD2 1.2.0, installed from devel>
 [2] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
 * ── Packages attached to the search path.

──────────────────────────────────────────────────────────────────────────────
```
