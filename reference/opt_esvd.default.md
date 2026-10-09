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
  nuisance_vec = NULL,
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
  natural parameter of these exponential-family distributions (see the
  "Natural parameter" line at the top of each `src/family_*.cpp`). For
  `"neg_binom2"`, the dot product is the log-mean of the distribution
  (i.e., similar to the canonical parameterization of the Poisson
  family).

- l2pen:

  Small positive number for the amount of penalization for both the
  cells' and the genes' latent vectors as well as the coefficients.

- library_multipler:

  Vector of positive numerics of length \\n\\: the per-cell multiplier
  \\s_i\\ in the likelihood, entering by family (for example, the
  Poisson mean is \\s_i e^{\theta\_{ij}}\\, and the Gaussian has mean
  \\s_i\theta\_{ij}\\ and variance \\s_i\gamma_j^2\\); `"neg_binom2"`
  ignores it. This is used as an alternative interpretation of how
  library-size affects a cell's gene expression (instead of using the
  library size as a covariate to be regressed out).

- max_iter:

  Positive integer for number of iterations.

- nuisance_vec:

  Vector of positive numerics of length \\p\\, representing each gene's
  nuisance parameter when using an exponential-family distribution that
  requires one: the standard deviation for `"gaussian"`, the mean
  divided by the standard deviation (the inverse coefficient of
  variation) for `"curved_gaussian"`, and the size (number of failures)
  for `"neg_binom"` and `"neg_binom2"`. It is ignored by `"poisson"`,
  `"exponential"` and `"bernoulli"`. The default `NULL` uses `1` for
  every gene, which is a placeholder rather than an estimate; supply
  your own values for the families that use it.

- offset_variables:

  A vector of strings depicting which column names in
  `input_obj$covariates` be treated as an offset during the optimization
  (i.e., their coefficients will not change throughout the
  optimization).

- tol:

  Relative tolerance of the stopping rule: iteration stops once the
  objective changes by at most `tol` times
  `max(1, |previous objective|)`.

- verbose:

  Integer

- ...:

  Additional parameters

## Value

a `list` with elements `x_mat`, `y_mat`, `z_mat`, `covariates`,
`library_multipler`, `loss` (the objective after every iteration),
`nuisance_vec` and `param`. A warning is raised if any row or column
update ended in a failed line search.
