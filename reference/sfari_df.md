# SFARI genes

From https://gene.sfari.org/database/human-gene/, after clicking on
"Download this dataset", accessed on June 15, 2024

## Usage

``` r
sfari_df
```

## Format

\## \`sfari_df\` A data frame with 1176 rows (genes) and 10 columns:

- status:

  See website for details

- gene.symbol:

  Gene symbol

- gene.name:

  Full gene name

- ensembl.id:

  Gene ID

- chromosome:

  Chromosome the gene is on

- genetic.category:

  See website for details

- gene.score:

  See website for details

- syndromic:

  See website for details

- eagle:

  See website for details

- number.of.reports:

  See website for details

## Source

\<https://www.science.org/doi/10.1126/science.aav8130\>

## Examples

``` r
if (FALSE) { # \dontrun{
# How the data was loaded
sfari_df <- read.csv("SFARI-Gene_genes_03-28-2024release_06-16-2024export.csv")
usethis::use_data(sfari_df)
} # }
```
