# Optimize eSVD for matrices or sparse matrices.

Optimize eSVD for matrices or sparse matrices.

## Usage

``` r
# Default S3 method
opt_esvd(
  input_obj,
  x_init,
  y_init,
  z_init = NULL,
  covariates = NULL,
  family = "poisson",
  l2pen = 0.1,
  library_multipler = rep(1, nrow(input_obj)),
  max_iter = 100,
  nuisance_vec = rep(NA, ncol(input_obj)),
  offset_variables = NULL,
  tol = 1e-06,
  verbose = 0,
  ...
)
```

## Arguments

- input_obj:

  Dataset (either `matrix` or `dgCMatrix`) where the \\n\\ rows
  represent cells and \\p\\ columns represent genes. The rows and
  columns of the matrix should be named.

- x_init:

  Initial matrix of the cells' latent vectors that is \\n\\ rows and
  \\k\\ columns. The row names should be the same as `input_obj`.

- y_init:

  Initial matrix of the genes' latent vectors that is \\p\\ rows and
  \\k\\ columns. The row names should be the same as the column names of
  `input_obj`.

- z_init:

  Initial matrix of the genes' coefficient vectors that is \\p\\ rows
  and `ncol(covariates)` columns. The row names should be the same as
  the column names of `input_obj`, and the column names should be the
  same as `covariates`.

- covariates:

  `matrix` object with \\n\\ rows with the same rownames as `input_obj`
  where the columns represent the different covariates. Notably, this
  should contain only numerical columns (i.e., all categorical variables
  should have already been split into numerous indicator variables).

- family:

  String among `"gaussian"`, `"curved_gaussian"`, `"exponential"`,
  `"poisson"`, `"neg_binom"`, `"neg_binom2"`, or `"bernoulli"`. Notably,
  with exception of `"neg_binom2"`, all the other families are
  parameterized such that eSVD is fitting the dot product to be the
  canonical parameter of these expoential-family distributions. For
  `"neg_binom2"`, the dot product is the log-mean of the distribution
  (i.e., similar to the canonical parameterization of the Poisson
  family).

- l2pen:

  Small positive number for the amount of penalization for both the
  cells' and the genes' latent vectors as well as the coefficients.

- library_multipler:

  Vector of positive numerics of length \\n\\. It is the multiplier such
  that the variance of cell `i`'s entries is the mean of cell `i`'s
  entries times the square-root of cell `i`'s value in
  `library_multipler` (entry-wise). This is used as an alternative
  interpretation of how library-size affects a cell's gene expression
  (instead of using the library size as a covariate to be regressed
  out).

- max_iter:

  Positive integer for number of iterations.

- nuisance_vec:

  Vector of non-negative numerics (or `NA`'s) of length \\p\\,
  representing each gene's nuisance parameter when using an
  exponential-family distribution that requires one. It is used only
  when `family` is `"curved_gaussian"` or `"neg_binom"` or
  `"neg_binom2"`.

- offset_variables:

  A vector of strings depicting which column names in
  `input_obj$covariate` be treated as an offset during the optimization
  (i.e., their coefficients will not change throughout the
  optimization).

- tol:

  Small positive number to differentiate between zero and non-zero.

- verbose:

  Integer

- ...:

  Additional parameters

## Value

a `list` with elements `x_mat`, `y_mat`, `z_mat`, `library_multiplier`,
`loss`, `nuisance_vec` and `param`.
