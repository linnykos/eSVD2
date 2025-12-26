# Generate null data

This function doesn't allow of much flexbility, and is meant primarily
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

a list with `covariates` (a \`matrix\` with \`n\` cells and 4 columns,
named \`"Intercept"\`, \`"Log_UMI"\`, \`"Sex"\`, and \`"Age"\`, which
will be used in \`eSVD2::initialize_esvd\`), `df` (a \`data.frame\` with
\`num_individuals\` rows and 5 columns, named named \`"Intercept"\`,
\`"Log_UMI"\`, \`"Sex"\`, \`"Age"\`, and \`"Individual"\`),
`metadata_individual` (a \`matrix\` with \`n\` rows and 1 column called
\`"Individual"\` for which cell originates from which individual)
`nat_mat` (a \`n\` by \`p\` \`matrix\` denoting the natural parameter
for each cell's gene expression), and `obs_mat` (a \`n\` by \`p\`
\`dgCMatrix\` denoting the observed count for each cell's gene
expression),
