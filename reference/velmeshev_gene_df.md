# Velmeshev et al. DEGs

From https://www.science.org/doi/10.1126/science.aav8130, data
supplement S4 (i.e., `aav8130_data-s4.xls`)

## Usage

``` r
velmeshev_gene_df
```

## Format

\## \`velmeshev_gene_df\` A data frame with 692 rows (genes) and 14
columns:

- Cell type:

  Cell type that the DEG was in

- gene ID:

  ENSG gene ID

- Gene name:

  Gene name

- Gene biotype:

  Biotype

- Fold change:

  See original paper for details

- Sample fold change:

  See original paper for details

- q value:

  See original paper for details

- correlation (bulk mRNA/bulkized nuclear RNA):

  See original paper for details

- Epilepsy DEG:

  boolean

- gene group:

  boolean

- SFARI gene:

  boolean

- Satterstrom:

  boolean

- Sanders:

  boolean

- cell type-specific expression:

  boolean

## Source

\<https://www.science.org/doi/10.1126/science.aav8130\>

## Examples

``` r
if (FALSE) { # \dontrun{
# How the data was loaded
tib <- readxl::read_xls("aav8130_data-s4.xls",
sheet = "ASD_DEGs")
velmeshev_gene_df <- as.data.frame(tib)
usethis::use_data(velmeshev_gene_df)
} # }
```
