# Neuron preprocessing for eSVD-DE tutorial

## Downloading the data

The data we will be working with can be downloaded from
<https://cells.ucsc.edu/autism/rawMatrix.zip>. After downloading this,
you will get a file called `rawMatrix.zip`, which, after unzipping,
contains

- `barcodes.tsv`: A matrix with 104,559 rows (i.e., cells) and 1 column,
  denoting the cells’ barcodes
- `genes.tsv`: A matrix with 65,217 rows (i.e., genes) and 2 columns,
  denoting the ENSEMBL ID or Gene Symbols
- `matrix.mtx`: A sparse matrix with 65,217 rows and 104,559 columns,
  denoting the counts for each cell-gene pair
- `meta.txt`: A matrix with 104,559 rows (i.e., cells) and 16 columns,
  containing various covariates about each cell/donor

The original paper is from <https://pubmed.ncbi.nlm.nih.gov/31097668/>.

## Setting up in R

You will the following packages. `sparseMatrixStats` is from
Bioconductor.

``` r
library(Matrix)
library(Seurat)
library(sparseMatrixStats)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
```

## Reading in the data

First, we read in the data. For simplicity in this vignette’s
preprocessing, since there are some duplicated genes, we arbitrarily
remove any duplicates.

``` r
mat <- Matrix::readMM("/Users/kevinlin/Downloads/rawMatrix/matrix.mtx")
gene_mat <- read.csv("/Users/kevinlin/Downloads/rawMatrix/genes.tsv", sep = "\t", header = F)
cell_mat <- read.csv("/Users/kevinlin/Downloads/rawMatrix/barcodes.tsv", sep = "\t", header = F)
metadata <- read.csv("/Users/kevinlin/Downloads/rawMatrix/meta.txt", sep = "\t", header = T, stringsAsFactors = TRUE)
metadata$individual <- factor(paste0("indiv_", metadata$individual))

rownames(mat) <- gene_mat[,2]
colnames(mat) <- cell_mat[,1]

# remove duplicates
if(any(duplicated(rownames(mat)))){
  mat <- mat[-which(duplicated(rownames(mat))),]
}

rownames(metadata) <- metadata$cell
metadata <- metadata[,setdiff(colnames(metadata), "cell")]
```

## Creating Seurat object

Next, we create the Seurat object and filter only the cells that are in
Layer 2/3. (We will abbreviate this Seurat object as “asd” for “Autism
Spectrum Disorder”.)

``` r
asd <- Seurat::CreateSeuratObject(counts = mat,
                                  meta.data = metadata)

asd <- subset(asd, cluster == "L2/3")

asd[["percent.mt"]] <- Seurat::PercentageFeatureSet(asd, pattern = "^MT-")
```

## Preprocessing the Seurat object

We then do basic preprocessing. Mainly, this consists of first remove
genes that have no counts at all, then removing cells that have too many
counts, then removing genes with too many counts. Then, we aggregate
genes that are deemed highly variable, the original DEGs found in the
Velmeshev et al. paper, housekeeping genes, and SFARI genes. The union
of all these genes are the ones that we set as “highly variable,” as we
would like the eSVD-DE analysis later to use these genes.

Once we are done with this step, we are ready to perform the eSVD-DE
analysis (in the other vignette).

``` r
set.seed(10)
mat <- SeuratObject::LayerData(asd, layer = "counts", assay = "RNA")
mat <- Matrix::t(mat)

# remove genes with all 0's
gene_sum <- Matrix::colSums(mat)
if (any(gene_sum == 0)) {
  features <- SeuratObject::Features(asd)
  features <- setdiff(features, names(gene_sum)[which(gene_sum == 0)])
  asd <- subset(asd, features = features)
}

# remove cells with too high counts
keep_vec <- rep(TRUE, length(Seurat::Cells(asd)))
keep_vec[asd$nCount_RNA >= 60000] <- FALSE
asd[["keep"]] <- keep_vec
asd <- subset(asd, keep == TRUE)

# remove genes with too high counts
mat <- SeuratObject::LayerData(asd, layer = "counts", assay = "RNA")
mat <- Matrix::t(mat)
gene_total <- Matrix::colSums(log1p(mat))
gene_threshold <- log1p(10) * nrow(mat)
features <- colnames(mat)[which(gene_total <= gene_threshold)]
asd <- subset(asd, features = features)

# find the variable genes
asd <- Seurat::NormalizeData(asd)
asd <- Seurat::FindVariableFeatures(asd, selection.method = "vst", nfeatures = 5000)

data("velmeshev_gene_df", package = "eSVD2")
de_genes <- velmeshev_gene_df[velmeshev_gene_df[, "Cell type"] == "L2/3", "Gene name"]

data("housekeeping_df", package = "eSVD2")
hk_genes <- housekeeping_df[, "Gene"]

data("sfari_df", package = "eSVD2")
sfari_genes <- sfari_df[, "gene.symbol"]

gene_keep <- unique(c(
  de_genes,
  hk_genes,
  sfari_genes,
  Seurat::VariableFeatures(asd)
))
rm_idx <- grep("^MT", gene_keep) # remove mitochondrial genes
if(length(rm_idx) > 0) gene_keep <- gene_keep[-rm_idx]
gene_current <- SeuratObject::Features(asd)
gene_keep <- intersect(gene_current, gene_keep)
Seurat::VariableFeatures(asd) <- gene_keep
asd <- subset(asd, features = Seurat::VariableFeatures(asd))
```

With this done, we can save this minimal Seurat object for the eSVD-DE
vignette. This file will be roughly 97 Mb in size.

``` r
save(asd, date_of_run, session_info, file = "asd_seurat.RData")
```

## Visualizing the data

It is nice to visualize any dataset before we do any analysis, and
eSVD-DE is no different. Hence, we can do some simple preprocessing.

``` r
set.seed(10)
asd <- Seurat::ScaleData(asd, features = Seurat::VariableFeatures(asd))
asd <- Seurat::RunPCA(asd,
                      features = Seurat::VariableFeatures(object = asd),
                      verbose = FALSE)

set.seed(10)
asd <- Seurat::RunUMAP(asd, dims = 1:30)
```

``` r
p1 <- Seurat::DimPlot(
  asd,
  reduction = 'umap',
  group.by = 'diagnosis',
  label = TRUE,
  repel = TRUE,
  label.size = 2.5
) + Seurat::NoLegend()
p2 <- Seurat::DimPlot(
  asd,
  reduction = 'umap',
  group.by = 'individual',
  label = TRUE,
  repel = TRUE,
  label.size = 2.5
) + Seurat::NoLegend()
p1 + p2
```

![UMAP of the L2/3
neurons](https://github.com/linnykos/eSVD2/blob/master/vignettes/asd-umap.png?raw=true)

UMAP of the L2/3 neurons

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
    ##  date     2025-12-26
    ##  pandoc   3.1.11 @ /opt/hostedtoolcache/pandoc/3.1.11/x64/ (via rmarkdown)
    ##  quarto   NA
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package           * version date (UTC) lib source
    ##  abind               1.4-8   2024-09-12 [1] RSPM
    ##  bslib               0.9.0   2025-01-30 [1] RSPM
    ##  cachem              1.1.0   2024-05-16 [1] RSPM
    ##  cli                 3.6.5   2025-04-23 [1] RSPM
    ##  cluster             2.1.8.1 2025-03-12 [3] CRAN (R 4.5.2)
    ##  codetools           0.2-20  2024-03-31 [3] CRAN (R 4.5.2)
    ##  cowplot             1.2.0   2025-07-07 [1] RSPM
    ##  data.table          1.18.0  2025-12-24 [1] RSPM
    ##  deldir              2.0-4   2024-02-28 [1] RSPM
    ##  desc                1.4.3   2023-12-10 [1] RSPM
    ##  devtools            2.4.6   2025-10-03 [1] RSPM
    ##  digest              0.6.39  2025-11-19 [1] RSPM
    ##  dotCall64           1.2     2024-10-04 [1] RSPM
    ##  dplyr               1.1.4   2023-11-17 [1] RSPM
    ##  ellipsis            0.3.2   2021-04-29 [1] RSPM
    ##  evaluate            1.0.5   2025-08-27 [1] RSPM
    ##  farver              2.1.2   2024-05-13 [1] RSPM
    ##  fastDummies         1.7.5   2025-01-20 [1] RSPM
    ##  fastmap             1.2.0   2024-05-15 [1] RSPM
    ##  fitdistrplus        1.2-4   2025-07-03 [1] RSPM
    ##  fs                  1.6.6   2025-04-12 [1] RSPM
    ##  future              1.68.0  2025-11-17 [1] RSPM
    ##  future.apply        1.20.1  2025-12-09 [1] RSPM
    ##  generics            0.1.4   2025-05-09 [1] RSPM
    ##  ggplot2             4.0.1   2025-11-14 [1] RSPM
    ##  ggrepel             0.9.6   2024-09-07 [1] RSPM
    ##  ggridges            0.5.7   2025-08-27 [1] RSPM
    ##  globals             0.18.0  2025-05-08 [1] RSPM
    ##  glue                1.8.0   2024-09-30 [1] RSPM
    ##  goftest             1.2-3   2021-10-07 [1] RSPM
    ##  gridExtra           2.3     2017-09-09 [1] RSPM
    ##  gtable              0.3.6   2024-10-25 [1] RSPM
    ##  htmltools           0.5.9   2025-12-04 [1] RSPM
    ##  htmlwidgets         1.6.4   2023-12-06 [1] RSPM
    ##  httpuv              1.6.16  2025-04-16 [1] RSPM
    ##  httr                1.4.7   2023-08-15 [1] RSPM
    ##  ica                 1.0-3   2022-07-08 [1] RSPM
    ##  igraph              2.2.1   2025-10-27 [1] RSPM
    ##  irlba               2.3.5.1 2022-10-03 [1] RSPM
    ##  jquerylib           0.1.4   2021-04-26 [1] RSPM
    ##  jsonlite            2.0.0   2025-03-27 [1] RSPM
    ##  KernSmooth          2.23-26 2025-01-01 [3] CRAN (R 4.5.2)
    ##  knitr               1.51    2025-12-20 [1] RSPM
    ##  later               1.4.4   2025-08-27 [1] RSPM
    ##  lattice             0.22-7  2025-04-02 [3] CRAN (R 4.5.2)
    ##  lazyeval            0.2.2   2019-03-15 [1] RSPM
    ##  lifecycle           1.0.4   2023-11-07 [1] RSPM
    ##  listenv             0.10.0  2025-11-02 [1] RSPM
    ##  lmtest              0.9-40  2022-03-21 [1] RSPM
    ##  magrittr            2.0.4   2025-09-12 [1] RSPM
    ##  MASS                7.3-65  2025-02-28 [3] CRAN (R 4.5.2)
    ##  Matrix            * 1.7-4   2025-08-28 [3] CRAN (R 4.5.2)
    ##  MatrixGenerics    * 1.22.0  2025-10-29 [1] Bioconduc~
    ##  matrixStats       * 1.5.0   2025-01-07 [1] RSPM
    ##  memoise             2.0.1   2021-11-26 [1] RSPM
    ##  mime                0.13    2025-03-17 [1] RSPM
    ##  miniUI              0.1.2   2025-04-17 [1] RSPM
    ##  nlme                3.1-168 2025-03-31 [3] CRAN (R 4.5.2)
    ##  otel                0.2.0   2025-08-29 [1] RSPM
    ##  parallelly          1.46.0  2025-12-12 [1] RSPM
    ##  patchwork           1.3.2   2025-08-25 [1] RSPM
    ##  pbapply             1.7-4   2025-07-20 [1] RSPM
    ##  pillar              1.11.1  2025-09-17 [1] RSPM
    ##  pkgbuild            1.4.8   2025-05-26 [1] RSPM
    ##  pkgconfig           2.0.3   2019-09-22 [1] RSPM
    ##  pkgdown             2.2.0   2025-11-06 [1] RSPM
    ##  pkgload             1.4.1   2025-09-23 [1] RSPM
    ##  plotly              4.11.0  2025-06-19 [1] RSPM
    ##  plyr                1.8.9   2023-10-02 [1] RSPM
    ##  png                 0.1-8   2022-11-29 [1] RSPM
    ##  polyclip            1.10-7  2024-07-23 [1] RSPM
    ##  progressr           0.18.0  2025-11-06 [1] RSPM
    ##  promises            1.5.0   2025-11-01 [1] RSPM
    ##  purrr               1.2.0   2025-11-04 [1] RSPM
    ##  R6                  2.6.1   2025-02-15 [1] RSPM
    ##  ragg                1.5.0   2025-09-02 [1] RSPM
    ##  RANN                2.6.2   2024-08-25 [1] RSPM
    ##  RColorBrewer        1.1-3   2022-04-03 [1] RSPM
    ##  Rcpp                1.1.0   2025-07-02 [1] RSPM
    ##  RcppAnnoy           0.0.22  2024-01-23 [1] RSPM
    ##  RcppHNSW            0.6.0   2024-02-04 [1] RSPM
    ##  remotes             2.5.0   2024-03-17 [1] RSPM
    ##  reshape2            1.4.5   2025-11-12 [1] RSPM
    ##  reticulate          1.44.1  2025-11-14 [1] RSPM
    ##  rlang               1.1.6   2025-04-11 [1] RSPM
    ##  rmarkdown           2.30    2025-09-28 [1] RSPM
    ##  ROCR                1.0-11  2020-05-02 [1] RSPM
    ##  RSpectra            0.16-2  2024-07-18 [1] RSPM
    ##  Rtsne               0.17    2023-12-07 [1] RSPM
    ##  S7                  0.2.1   2025-11-14 [1] RSPM
    ##  sass                0.4.10  2025-04-11 [1] RSPM
    ##  scales              1.4.0   2025-04-24 [1] RSPM
    ##  scattermore         1.2     2023-06-12 [1] RSPM
    ##  sctransform         0.4.2   2025-04-30 [1] RSPM
    ##  sessioninfo         1.2.3   2025-02-05 [1] RSPM
    ##  Seurat            * 5.4.0   2025-12-14 [1] RSPM
    ##  SeuratObject      * 5.3.0   2025-12-12 [1] RSPM
    ##  shiny               1.12.1  2025-12-09 [1] RSPM
    ##  sp                * 2.2-0   2025-02-01 [1] RSPM
    ##  spam                2.11-1  2025-01-20 [1] RSPM
    ##  sparseMatrixStats * 1.22.0  2025-10-29 [1] Bioconduc~
    ##  spatstat.data       3.1-9   2025-10-18 [1] RSPM
    ##  spatstat.explore    3.6-0   2025-11-22 [1] RSPM
    ##  spatstat.geom       3.6-1   2025-11-20 [1] RSPM
    ##  spatstat.random     3.4-3   2025-11-21 [1] RSPM
    ##  spatstat.sparse     3.1-0   2024-06-21 [1] RSPM
    ##  spatstat.univar     3.1-5   2025-11-17 [1] RSPM
    ##  spatstat.utils      3.2-0   2025-09-20 [1] RSPM
    ##  stringi             1.8.7   2025-03-27 [1] RSPM
    ##  stringr             1.6.0   2025-11-04 [1] RSPM
    ##  survival            3.8-3   2024-12-17 [3] CRAN (R 4.5.2)
    ##  systemfonts         1.3.1   2025-10-01 [1] RSPM
    ##  tensor              1.5.1   2025-06-17 [1] RSPM
    ##  textshaping         1.0.4   2025-10-10 [1] RSPM
    ##  tibble              3.3.0   2025-06-08 [1] RSPM
    ##  tidyr               1.3.2   2025-12-19 [1] RSPM
    ##  tidyselect          1.2.1   2024-03-11 [1] RSPM
    ##  usethis             3.2.1   2025-09-06 [1] RSPM
    ##  uwot                0.2.4   2025-11-10 [1] RSPM
    ##  vctrs               0.6.5   2023-12-01 [1] RSPM
    ##  viridisLite         0.4.2   2023-05-02 [1] RSPM
    ##  xfun                0.55    2025-12-16 [1] RSPM
    ##  xtable              1.8-4   2019-04-21 [1] RSPM
    ##  yaml                2.3.12  2025-12-10 [1] RSPM
    ##  zoo                 1.8-15  2025-12-15 [1] RSPM
    ## 
    ##  [1] /home/runner/work/_temp/Library
    ##  [2] /opt/R/4.5.2/lib/R/site-library
    ##  [3] /opt/R/4.5.2/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
