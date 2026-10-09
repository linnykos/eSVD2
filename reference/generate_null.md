# Generate null data

This function does not allow much flexibility, and is meant primarily
for simple null simulations. Here, \`n\` (the number of cells) is equal
to \`cell_per_person\*num_individuals\`.

## Usage

``` r
generate_null(cell_per_person = 100, num_genes = 1000, num_individuals = 20)
```

## Arguments

- cell_per_person:

  The number of cells per individual

- num_genes:

  The number of genes in the simulation

- num_individuals:

  The number of individuals

## Value

a list with `covariates` (a \`matrix\` with \`n\` rows (cells) and 5
columns, named \`"Intercept"\`, \`"Log_UMI"\`, \`"CC"\`, \`"Sex"\` and
\`"Age"\`, ready for \`initialize_esvd\`; \`"CC"\` is the 0/1
case-control indicator), `df` (a \`data.frame\` with \`num_individuals\`
rows and the 6 columns \`"Intercept"\`, \`"Log_UMI"\`, \`"CC"\`,
\`"Sex"\`, \`"Age"\` and \`"Individual"\`), `metadata_individual` (a
\`factor\` of length \`n\` naming each cell's individual), `nat_mat` (a
\`n\` by \`p\` \`matrix\` of natural parameters), and `obs_mat` (a \`n\`
by \`p\` \`dgCMatrix\` of observed counts). The first 10 genes carry a
case-control effect; the rest are null.

## Examples

``` r
set.seed(10)
sim <- generate_null(cell_per_person = 15, num_genes = 40,
                     num_individuals = 8)
dim(sim$obs_mat)
#> [1] 120  40
table(sim$metadata_individual)
#> 
#> indiv_1 indiv_2 indiv_3 indiv_4 indiv_5 indiv_6 indiv_7 indiv_8 
#>      15      15      15      15      15      15      15      15 
```
