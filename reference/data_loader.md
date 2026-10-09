# Internal data loader

Wraps a dense `matrix` or a `dgCMatrix` in a C++ iterator used by the
optimizer. The pointer is only valid in the R session that created it; a
restored copy (e.g. from `readRDS`) errors when used.

## Usage

``` r
data_loader(mat)
```

## Arguments

- mat:

  `matrix` (integer or double) or `dgCMatrix`, cells by genes.

## Value

An external pointer of class `esvd_data_loader`.
