# Fisher's exact test

Perform a Fisher's exact test. This formula is from
https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/

## Usage

``` r
fisher_test(set1_genes, set2_genes, all_genes, verbose = 0)
```

## Arguments

- set1_genes:

  A vector of characters of your resulting gene list

- set2_genes:

  A vector of characters of a gene list you are comparing against

- all_genes:

  A vector of characters of all the genes in your "universe". This
  should include all the genes in `set1_genes` and `set2_genes`.

- verbose:

  Numeric

## Value

A list of: `pvalue` for the p-value, `set_bg_len` for the number of
genes in `all_genes` but not in `set2_genes`, `set_overlap_len` for the
number of genes in both `set1_genes` and `set2_genes`, `set1_len` for
number of genes in `set1_genes`, and `set2_len` for number of genes in
`set2_genes`.
