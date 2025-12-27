# Neuron analysis for Layer 2/3 in ASD

For this vignette, we will be using a dataset that can be downloaded
from
<https://www.dropbox.com/scl/fi/csyz62obxx8xj5gl5sgp4/asd_seurat.RData?rlkey=d7dlhqjesunal52u8z2iomjsg&dl=0>.
This file is called `asd_seurat.RData`, and contains the following
variables:

- `asd`: The Seurat object we will be analyzing, containing 13,302 cells
  and 6,599 genes (all of which are deemed highly variable)
- `date_of_run`: The date when `asd` was constructed
- `session_info`: The session info when `asd` was constructed

Notably, this is **not** the version of the dataset used in the eSVD-DE
paper. Instead, the intent of this vignette to show you how to perform
an eSVD-DE analysis. That is, while the processed dataset in the paper
is different from this vignette’s processed dataset, the procedure we
showcase is here is intended as the typical usage of our method.

If you wish to see how this data was preprocessed starting from publicly
available data (or you cannot download the `.RData` file), please see
<https://linnykos.github.io/eSVD2/articles/asd-preprocess.html>.

Relevant citation:

    Velmeshev, D., Schirmer, L., Jung, D., Haeussler, M., Perez, Y., Mayer, S., ... & Kriegstein, A. R. (2019). Single-cell genomics identifies cell type–specific molecular changes in autism. Science, 364(6441), 685-689.

``` r
library(Seurat)
library(SeuratObject)
library(EnhancedVolcano)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
```

## Loading the data

We first load the data. We can use the
[`table()`](https://rdrr.io/r/base/table.html) function to roughly see
how many cells each donor contributes and how many donors were
classified to be `ASD` or `Control`.

``` r
load("asd_seurat.RData")
asd
table(asd$individual, asd$diagnosis)
```

    An object of class Seurat 
    6599 features across 13302 samples within 1 assay 
    Active assay: RNA (6599 features, 6599 variable features)
     2 layers present: counts, data

                ASD Control
    indiv_1823    0      72
    indiv_4341    0     733
    indiv_4849  428       0
    indiv_4899  349       0
    indiv_5144  202       0
    indiv_5163    0     145
    indiv_5242    0      81
    indiv_5278  944       0
    indiv_5294  449       0
    indiv_5387    0    1129
    indiv_5391    0     237
    indiv_5403   69       0
    indiv_5408    0     162
    indiv_5419  327       0
    indiv_5531  631       0
    indiv_5538    0     391
    indiv_5554    0     233
    indiv_5565  881       0
    indiv_5577    0     215
    indiv_5841  451       0
    indiv_5864  275       0
    indiv_5879    0     278
    indiv_5893    0     976
    indiv_5936    0     193
    indiv_5939 1016       0
    indiv_5945  422       0
    indiv_5958    0     826
    indiv_5976    0     106
    indiv_5978  264       0
    indiv_6032    0     289
    indiv_6033  528       0

## Using eSVD-DE

We will be testing for differential expression in the mean between the
15 donors with ASD and the 16 control donors, among the 6599 genes.

- **Step 0 (Formatting the covariates):** Before we use eSVD-DE, we want
  to reformat all the covariates. Typically, it is a good idea to set
  all the categorical covariates to be `factor` variables where the last
  level has the smallest size. This is because that last level will be
  the one that does not get its own indicator column during the one-hot
  encoding when using
  [`eSVD2::format_covariates()`](https://linnykos.github.io/eSVD2/reference/format_covariates.md).

``` r
gene_vec <- Seurat::VariableFeatures(asd[["RNA"]])
mat <- SeuratObject::LayerData(asd,
                               layer = "counts",
                               assay = "RNA",
                               features = gene_vec)
mat <- Matrix::t(mat)
covariate_df <- asd@meta.data[, c("individual",
                                  "region",
                                  "age",
                                  "sex",
                                  "diagnosis",
                                  "Seqbatch",
                                  "Capbatch")]

factor_variables <- c("region", "diagnosis", "sex", "Seqbatch", "Capbatch")
for (variable in factor_variables) {
  covariate_df[, variable] <- factor(covariate_df[, variable], levels = names(sort(table(covariate_df[, variable]), decreasing = T)))
}

covariate_df[, "diagnosis"] <- factor(covariate_df[, "diagnosis"], levels = c("Control", "ASD"))

covariates <- eSVD2::format_covariates(
  dat = mat,
  covariate_df = covariate_df,
  rescale_numeric_variables = "age"
)
```

- **Step 1 (Initializing the eSVD object):** First, we initialize the
  object of class `eSVD`, which we will call `eSVD_obj`. This is done
  using the
  [`eSVD2::initialize_esvd()`](https://linnykos.github.io/eSVD2/reference/initialize_esvd.md)
  function. All the following operations will use `eSVD_obj` for all
  calculations through the usage of this method. This is where we pass
  the main inputs `dat`, `covariates`, and `metadata_individual`. Notice
  that after
  [`eSVD2::initialize_esvd()`](https://linnykos.github.io/eSVD2/reference/initialize_esvd.md)
  (and later, after
  [`eSVD2::opt_esvd()`](https://linnykos.github.io/eSVD2/reference/opt_esvd.md))
  we always call
  [`eSVD2::reparameterization_esvd_covariates()`](https://linnykos.github.io/eSVD2/reference/reparameterization_esvd_covariates.md)
  to ensure that the fit remains identifiable.

Empirically, we have found it useful to set `omitted_variables` as the
`Log_UMI` and the batch variables. This means that when eSVD2
reparameterizes the variables, it does not fiddle with the coefficient
for `Log_UMI`. We found this beneficial since the coefficient for
specifically the covariate `Log_UMI` bears a special meaning (as it
impacts the normalization of the sequencing depth most directly), and we
do not wish to “mix” covariates correlated with `Log_UMI` too early. We
have also found it useful to keep the batch variables separate since
some datasets we have analyzed have donor covariates that are very
correlated with the batch variable, and keeping these two aspects (i.e.,
the biological effect of donor covariates, and the artifical effect of
the sequencing batches) seperate to be beneficial. (Overall, this is not
a strict recommendation, however, we have found this practice to be
beneficial in almost all instances we have used eSVD2, so we ourselves
have not tried deviating from this convention.)

On a single core of a 2023 Macbook Pro (Apple M2 Max, 32 Gb RAM), this
took roughly 7 minutes.

``` r
time_start1 <- Sys.time()
eSVD_obj <- eSVD2::initialize_esvd(
  dat = mat,
  covariates = covariates[, -grep("individual", colnames(covariates))],
  metadata_individual = covariate_df[, "individual"],
  case_control_variable = "diagnosis_ASD",
  bool_intercept = TRUE,
  k = 30,
  lambda = 0.1,
  metadata_case_control = covariates[, "diagnosis_ASD"],
  verbose = 1
)
time_end1 <- Sys.time()

omitted_variables <- colnames(eSVD_obj$covariates)[c(grep("Seqbatch", colnames(eSVD_obj$covariates)), grep("Capbatch", colnames(eSVD_obj$covariates)))]
eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Init",
  omitted_variables = c("Log_UMI", omitted_variables)
)
```

- **Step 2 (Fitting the eSVD object):** We now fit the eSVD. Typically,
  we recommend fitting twice (i.e., calling
  [`eSVD2::opt_esvd()`](https://linnykos.github.io/eSVD2/reference/opt_esvd.md)
  twice) due to the non-convex nature of the objective function. In the
  first call of
  [`eSVD2::opt_esvd()`](https://linnykos.github.io/eSVD2/reference/opt_esvd.md),
  hold all the coefficients to all the covariates aside from the
  case-control covariate (here, denoted as `"CC"` so we remove as much
  of the case-control covariate effect as possible. Then, in the second
  call of
  [`eSVD2::opt_esvd()`](https://linnykos.github.io/eSVD2/reference/opt_esvd.md),
  we allow all the coefficients to be optimized.

Note that in the last fit of
[`eSVD2::opt_esvd()`](https://linnykos.github.io/eSVD2/reference/opt_esvd.md),
we do not set any `omitted_variables`. This is our typical
recommendation – it is good to allow all coefficients to be
reparameterized right before the optimization is all done.

On a single core of a 2023 Macbook Pro (Apple M2 Max, 32 Gb RAM), for
this vignette’s dataset, the first fit (i.e., call of
[`eSVD2::opt_esvd()`](https://linnykos.github.io/eSVD2/reference/opt_esvd.md))
took roughly 7 minutes (for 11 iterations). The second fit took roughly
4 minutes (for 15 iterations).

``` r
time_start2 <- Sys.time()
eSVD_obj <- eSVD2::opt_esvd(
  input_obj = eSVD_obj,
  l2pen = 0.1,
  max_iter = 100,
  offset_variables = setdiff(colnames(eSVD_obj$covariates), "diagnosis_ASD"),
  tol = 1e-6,
  verbose = 1,
  fit_name = "fit_First",
  fit_previous = "fit_Init"
)
time_end2 <- Sys.time()
eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_First",
  omitted_variables = c("Log_UMI", omitted_variables)
)

time_start3 <- Sys.time()
eSVD_obj <- eSVD2::opt_esvd(
  input_obj = eSVD_obj,
  l2pen = 0.1,
  max_iter = 100,
  offset_variables = NULL,
  tol = 1e-6,
  verbose = 1,
  fit_name = "fit_Second",
  fit_previous = "fit_First"
)
time_end3 <- Sys.time()
eSVD_obj <- eSVD2::reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Second",
  omitted_variables = omitted_variables
)
```

- **Step 3 (Estimating the nuisance (i.e., over-dispersion)
  parameters):** After fitting the eSVD, we need to estimate the
  nuisance parameters. Since we are modeling the data via a negative
  binomial, this nuisance parameter is the over-dispersion parameter.
  This is done via the
  [`eSVD2::estimate_nuisance()`](https://linnykos.github.io/eSVD2/reference/estimate_nuisance.md)
  function.

On a single core of a 2023 Macbook Pro (Apple M2 Max, 32 Gb RAM), this
took roughly 9 minutes.

``` r
time_start4 <- Sys.time()
eSVD_obj <- eSVD2::estimate_nuisance(
  input_obj = eSVD_obj,
  bool_covariates_as_library = TRUE,
  bool_library_includes_interept = TRUE,
  bool_use_log = FALSE,
  verbose = 1
)
time_end4 <- Sys.time()
```

- **Step 4 (Computing the posterior):** With the over-disperison
  parameter estimated, we are now ready to compute the posterior mean
  and variance of each cell’s gene expression value. This posterior
  distribution is to balance out the fitted eSVD denoised value against
  the covariate-adjusted library size. This step ensures that even
  though we pooled all the genes to denoise the single-cell expression
  matrix, we do not later inflate the Type-1 error. This is done via the
  [`eSVD2::compute_posterior()`](https://linnykos.github.io/eSVD2/reference/compute_posterior.md)
  function.

``` r
eSVD_obj <- eSVD2::compute_posterior(
  input_obj = eSVD_obj,
  bool_adjust_covariates = FALSE,
  alpha_max = 2 * max(mat@x),
  bool_covariates_as_library = TRUE,
  bool_stabilize_underdispersion = TRUE,
  library_min = 0.1,
  pseudocount = 0
)
```

- **Step 5 (Computing the test statistics):** Now we can compute the
  test statistics using the
  [`eSVD2::compute_test_statistic()`](https://linnykos.github.io/eSVD2/reference/compute_test_statistic.md)
  function.

``` r
eSVD_obj <- eSVD2::compute_test_statistic(input_obj = eSVD_obj, verbose = 1)
```

- **Step 6 (Computing the p-values):** Lastly, we can compute the
  p-values via the
  [`eSVD2::compute_pvalue()`](https://linnykos.github.io/eSVD2/reference/compute_pvalue.md)
  function.

``` r
eSVD_obj <- eSVD2::compute_pvalue(input_obj = eSVD_obj)

save(eSVD_obj, session_info, date_of_run,
     file = "asd_esvd.RData")
```

## Plotting the p-values

To visualize the p-values, we can plot a volcano plot where the x-axis
is our estimated log fold change (based on the global case donor mean,
where each donor’s specific mean is first computed, minus the global
control donor mean) and the y-axis the negative log 10 p-value. The
dotted horizontal line denotes the FDR cutoff.

To give this volcano plot some biological interpretation, we show where
the genes found in Gandal et al. (see
<https://www.nature.com/articles/s41586-022-05377-7>) on the volcano
plot as well. That study was based on the bulk sequencing of brain
regions, where possibly many cell types are involved in the DEG
analysis. Nonetheless, we would expect a good number of DEGs from Gandal
et al. to be highly significant in our analysis, which is what our
volcano plot shows.

``` r
pvalue_vec <- 10 ^ (-eSVD_obj$pvalue_list$log10pvalue)
tmp_df <- data.frame(
  gene = names(eSVD_obj$pvalue_list$log10pvalue),
  lfc = log2(eSVD_obj$case_mean) - log2(eSVD_obj$control_mean),
  pvalue = pvalue_vec
)
pCutoff <- max(pvalue_vec[which(eSVD_obj$pvalue_list$fdr_vec <= 0.05)])

data("gandal_df", package = "eSVD2")
gandal_genes <- gandal_df[which(gandal_df[, "WholeCortex_ASD_FDR"] <= 0.05), "external_gene_name"]
gandal_genes <- intersect(gandal_genes, tmp_df$gene)

tmp_df2 <- tmp_df[c(which(!tmp_df$gene %in% gandal_genes),
                    which(tmp_df$gene %in% gandal_genes)), ]

colCustom <- sapply(tmp_df2$gene, function(gene) {
  ifelse(gene %in% gandal_genes, "coral", "navy")
})
names(colCustom)[colCustom == "coral"] <- "Gandal et al."
names(colCustom)[colCustom == "navy"] <- "others"

# png(
#   "volcano.png",
#   width = 1500,
#   height = 2000,
#   res = 300,
#   units = "px"
# )
EnhancedVolcano::EnhancedVolcano(
  tmp_df2,
  lab = tmp_df2$gene,
  selectLab = gandal_genes,
  colCustom = colCustom,
  x = "lfc",
  y = "pvalue",
  ylim = c(0, max(eSVD_obj$pvalue_list$log10pvalue)),
  labSize = 3,
  pCutoff = pCutoff,
  FCcutoff = 0,
  drawConnectors = TRUE
)
# graphics.off()
```

![Volcano plot the eSVD-DE analysis, highlighting the Gandal et al.
genes](https://github.com/linnykos/eSVD2/blob/master/vignettes/volcano.png?raw=true)

Volcano plot the eSVD-DE analysis, highlighting the Gandal et al. genes

Alternatively, you can make the volcano plot that specifically
highlights the housekeeping genes, which can act as a rough “negative
control”. That is, we would not expect many of the housekeeping genes to
be significantly different between case and control donors. This is what
we see in our volcano plot.

``` r
data("housekeeping_df", package = "eSVD2")
hk_genes <- housekeeping_df[, "Gene"]
hk_genes <- intersect(hk_genes, tmp_df$gene)

tmp_df2 <- tmp_df[c(which(!tmp_df$gene %in% hk_genes),
                    which(tmp_df$gene %in% hk_genes)), ]

colCustom <- sapply(tmp_df2$gene, function(gene) {
  ifelse(gene %in% hk_genes, "darkolivegreen1", "navy")
})
names(colCustom)[colCustom == "darkolivegreen1"] <- "Housekeeping genes"
names(colCustom)[colCustom == "navy"] <- "others"

# png(
#   "volcano_hk.png",
#   width = 1500,
#   height = 2000,
#   res = 300,
#   units = "px"
# )
EnhancedVolcano::EnhancedVolcano(
  tmp_df2,
  lab = tmp_df2$gene,
  selectLab = hk_genes,
  colCustom = colCustom,
  x = "lfc",
  y = "pvalue",
  ylim = c(0, max(eSVD_obj$pvalue_list$log10pvalue)),
  labSize = 3,
  pCutoff = pCutoff,
  FCcutoff = 0,
  drawConnectors = TRUE
)
# graphics.off()
```

![Volcano plot the eSVD-DE analysis, highlighting the housekeeping
genes](https://github.com/linnykos/eSVD2/blob/master/vignettes/volcano_hk.png?raw=true)

Volcano plot the eSVD-DE analysis, highlighting the housekeeping genes

## Setup

The following shows the suggested package versions that the developer
(GitHub username: linnykos) used when developing the eSVD2 package.

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.2 (2025-10-31)
    ##  os       Ubuntu 24.04.3 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language en
    ##  collate  C.UTF-8
    ##  ctype    C.UTF-8
    ##  tz       UTC
    ##  date     2025-12-27
    ##  pandoc   3.1.11 @ /opt/hostedtoolcache/pandoc/3.1.11/x64/ (via rmarkdown)
    ##  quarto   NA
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package          * version date (UTC) lib source
    ##  abind              1.4-8   2024-09-12 [1] RSPM
    ##  bslib              0.9.0   2025-01-30 [1] RSPM
    ##  cachem             1.1.0   2024-05-16 [1] RSPM
    ##  cli                3.6.5   2025-04-23 [1] RSPM
    ##  cluster            2.1.8.1 2025-03-12 [3] CRAN (R 4.5.2)
    ##  codetools          0.2-20  2024-03-31 [3] CRAN (R 4.5.2)
    ##  cowplot            1.2.0   2025-07-07 [1] RSPM
    ##  data.table         1.18.0  2025-12-24 [1] RSPM
    ##  deldir             2.0-4   2024-02-28 [1] RSPM
    ##  desc               1.4.3   2023-12-10 [1] RSPM
    ##  devtools           2.4.6   2025-10-03 [1] RSPM
    ##  digest             0.6.39  2025-11-19 [1] RSPM
    ##  dotCall64          1.2     2024-10-04 [1] RSPM
    ##  dplyr              1.1.4   2023-11-17 [1] RSPM
    ##  ellipsis           0.3.2   2021-04-29 [1] RSPM
    ##  EnhancedVolcano  * 1.28.2  2025-12-04 [1] Bioconduc~
    ##  evaluate           1.0.5   2025-08-27 [1] RSPM
    ##  farver             2.1.2   2024-05-13 [1] RSPM
    ##  fastDummies        1.7.5   2025-01-20 [1] RSPM
    ##  fastmap            1.2.0   2024-05-15 [1] RSPM
    ##  fitdistrplus       1.2-4   2025-07-03 [1] RSPM
    ##  fs                 1.6.6   2025-04-12 [1] RSPM
    ##  future             1.68.0  2025-11-17 [1] RSPM
    ##  future.apply       1.20.1  2025-12-09 [1] RSPM
    ##  generics           0.1.4   2025-05-09 [1] RSPM
    ##  ggplot2          * 4.0.1   2025-11-14 [1] RSPM
    ##  ggrepel          * 0.9.6   2024-09-07 [1] RSPM
    ##  ggridges           0.5.7   2025-08-27 [1] RSPM
    ##  globals            0.18.0  2025-05-08 [1] RSPM
    ##  glue               1.8.0   2024-09-30 [1] RSPM
    ##  goftest            1.2-3   2021-10-07 [1] RSPM
    ##  gridExtra          2.3     2017-09-09 [1] RSPM
    ##  gtable             0.3.6   2024-10-25 [1] RSPM
    ##  htmltools          0.5.9   2025-12-04 [1] RSPM
    ##  htmlwidgets        1.6.4   2023-12-06 [1] RSPM
    ##  httpuv             1.6.16  2025-04-16 [1] RSPM
    ##  httr               1.4.7   2023-08-15 [1] RSPM
    ##  ica                1.0-3   2022-07-08 [1] RSPM
    ##  igraph             2.2.1   2025-10-27 [1] RSPM
    ##  irlba              2.3.5.1 2022-10-03 [1] RSPM
    ##  jquerylib          0.1.4   2021-04-26 [1] RSPM
    ##  jsonlite           2.0.0   2025-03-27 [1] RSPM
    ##  KernSmooth         2.23-26 2025-01-01 [3] CRAN (R 4.5.2)
    ##  knitr              1.51    2025-12-20 [1] RSPM
    ##  later              1.4.4   2025-08-27 [1] RSPM
    ##  lattice            0.22-7  2025-04-02 [3] CRAN (R 4.5.2)
    ##  lazyeval           0.2.2   2019-03-15 [1] RSPM
    ##  lifecycle          1.0.4   2023-11-07 [1] RSPM
    ##  listenv            0.10.0  2025-11-02 [1] RSPM
    ##  lmtest             0.9-40  2022-03-21 [1] RSPM
    ##  magrittr           2.0.4   2025-09-12 [1] RSPM
    ##  MASS               7.3-65  2025-02-28 [3] CRAN (R 4.5.2)
    ##  Matrix             1.7-4   2025-08-28 [3] CRAN (R 4.5.2)
    ##  matrixStats        1.5.0   2025-01-07 [1] RSPM
    ##  memoise            2.0.1   2021-11-26 [1] RSPM
    ##  mime               0.13    2025-03-17 [1] RSPM
    ##  miniUI             0.1.2   2025-04-17 [1] RSPM
    ##  nlme               3.1-168 2025-03-31 [3] CRAN (R 4.5.2)
    ##  otel               0.2.0   2025-08-29 [1] RSPM
    ##  parallelly         1.46.0  2025-12-12 [1] RSPM
    ##  patchwork          1.3.2   2025-08-25 [1] RSPM
    ##  pbapply            1.7-4   2025-07-20 [1] RSPM
    ##  pillar             1.11.1  2025-09-17 [1] RSPM
    ##  pkgbuild           1.4.8   2025-05-26 [1] RSPM
    ##  pkgconfig          2.0.3   2019-09-22 [1] RSPM
    ##  pkgdown            2.2.0   2025-11-06 [1] RSPM
    ##  pkgload            1.4.1   2025-09-23 [1] RSPM
    ##  plotly             4.11.0  2025-06-19 [1] RSPM
    ##  plyr               1.8.9   2023-10-02 [1] RSPM
    ##  png                0.1-8   2022-11-29 [1] RSPM
    ##  polyclip           1.10-7  2024-07-23 [1] RSPM
    ##  progressr          0.18.0  2025-11-06 [1] RSPM
    ##  promises           1.5.0   2025-11-01 [1] RSPM
    ##  purrr              1.2.0   2025-11-04 [1] RSPM
    ##  R6                 2.6.1   2025-02-15 [1] RSPM
    ##  ragg               1.5.0   2025-09-02 [1] RSPM
    ##  RANN               2.6.2   2024-08-25 [1] RSPM
    ##  RColorBrewer       1.1-3   2022-04-03 [1] RSPM
    ##  Rcpp               1.1.0   2025-07-02 [1] RSPM
    ##  RcppAnnoy          0.0.22  2024-01-23 [1] RSPM
    ##  RcppHNSW           0.6.0   2024-02-04 [1] RSPM
    ##  remotes            2.5.0   2024-03-17 [1] RSPM
    ##  reshape2           1.4.5   2025-11-12 [1] RSPM
    ##  reticulate         1.44.1  2025-11-14 [1] RSPM
    ##  rlang              1.1.6   2025-04-11 [1] RSPM
    ##  rmarkdown          2.30    2025-09-28 [1] RSPM
    ##  ROCR               1.0-11  2020-05-02 [1] RSPM
    ##  RSpectra           0.16-2  2024-07-18 [1] RSPM
    ##  Rtsne              0.17    2023-12-07 [1] RSPM
    ##  S7                 0.2.1   2025-11-14 [1] RSPM
    ##  sass               0.4.10  2025-04-11 [1] RSPM
    ##  scales             1.4.0   2025-04-24 [1] RSPM
    ##  scattermore        1.2     2023-06-12 [1] RSPM
    ##  sctransform        0.4.2   2025-04-30 [1] RSPM
    ##  sessioninfo        1.2.3   2025-02-05 [1] RSPM
    ##  Seurat           * 5.4.0   2025-12-14 [1] RSPM
    ##  SeuratObject     * 5.3.0   2025-12-12 [1] RSPM
    ##  shiny              1.12.1  2025-12-09 [1] RSPM
    ##  sp               * 2.2-0   2025-02-01 [1] RSPM
    ##  spam               2.11-1  2025-01-20 [1] RSPM
    ##  spatstat.data      3.1-9   2025-10-18 [1] RSPM
    ##  spatstat.explore   3.6-0   2025-11-22 [1] RSPM
    ##  spatstat.geom      3.6-1   2025-11-20 [1] RSPM
    ##  spatstat.random    3.4-3   2025-11-21 [1] RSPM
    ##  spatstat.sparse    3.1-0   2024-06-21 [1] RSPM
    ##  spatstat.univar    3.1-5   2025-11-17 [1] RSPM
    ##  spatstat.utils     3.2-0   2025-09-20 [1] RSPM
    ##  stringi            1.8.7   2025-03-27 [1] RSPM
    ##  stringr            1.6.0   2025-11-04 [1] RSPM
    ##  survival           3.8-3   2024-12-17 [3] CRAN (R 4.5.2)
    ##  systemfonts        1.3.1   2025-10-01 [1] RSPM
    ##  tensor             1.5.1   2025-06-17 [1] RSPM
    ##  textshaping        1.0.4   2025-10-10 [1] RSPM
    ##  tibble             3.3.0   2025-06-08 [1] RSPM
    ##  tidyr              1.3.2   2025-12-19 [1] RSPM
    ##  tidyselect         1.2.1   2024-03-11 [1] RSPM
    ##  usethis            3.2.1   2025-09-06 [1] RSPM
    ##  uwot               0.2.4   2025-11-10 [1] RSPM
    ##  vctrs              0.6.5   2023-12-01 [1] RSPM
    ##  viridisLite        0.4.2   2023-05-02 [1] RSPM
    ##  withr              3.0.2   2024-10-28 [1] RSPM
    ##  xfun               0.55    2025-12-16 [1] RSPM
    ##  xtable             1.8-4   2019-04-21 [1] RSPM
    ##  yaml               2.3.12  2025-12-10 [1] RSPM
    ##  zoo                1.8-15  2025-12-15 [1] RSPM
    ## 
    ##  [1] /home/runner/work/_temp/Library
    ##  [2] /opt/R/4.5.2/lib/R/site-library
    ##  [3] /opt/R/4.5.2/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
