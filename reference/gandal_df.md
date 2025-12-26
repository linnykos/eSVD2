# Gandal et al. results

From https://www.nature.com/articles/s41586-022-05377-7, Supplementary
Data 3 (i.e., `41586_2022_5377_MOESM5_ESM.xlsx`)

## Usage

``` r
gandal_df
```

## Format

\## \`gandal_df\` A data frame with 24836 rows (genes) and 3 columns:

- external_gene_name:

  Gene symbol

- WholeCortex_ASD_logFC:

  Estimated log fold change

- WholeCortex_ASD_FDR:

  FDR q-value

## Source

\<https://www.nature.com/articles/s41586-022-05377-7\>

## Examples

``` r
if (FALSE) { # \dontrun{
# How the data was loaded
df <- openxlsx::read.xlsx("/Users/kevinlin/Downloads/41586_2022_5377_MOESM5_ESM.xlsx",
                          sheet = "DEGene_Statistics")
gandal_df <- df[,c("external_gene_name", "WholeCortex_ASD_logFC", "WholeCortex_ASD_FDR")]
usethis::use_data(gandal_df)
} # }
```
