# Hounkpe et al. housekeeping genes

From
https://academic.oup.com/nar/article/49/D1/D947/5871367#supplementary-data,
gkaa609_Supplemental_Files, specifically, `Supplementary_Table1.xlsx`)

## Usage

``` r
housekeeping_df
```

## Format

\## \`housekeeping_df\` A data frame with 1003 rows (genes) and 1
column:

- Gene:

  Gene name

## Source

\<https://academic.oup.com/nar/article/49/D1/D947/5871367\>

## Examples

``` r
if (FALSE) { # \dontrun{
# How the data was loaded
df <- openxlsx::read.xlsx("Supplementary_Table1.xlsx",
sheet = "Gene Model")
colnames(df) <- df[1,1]
housekeeping_df <- df[-1,,drop = FALSE]
usethis::use_data(housekeeping_df)
} # }
```
