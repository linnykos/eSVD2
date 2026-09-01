# helper-fixtures.R
#
# Shared fixtures for the suite, per UNIT_TEST_PLAN.md sections 0 (C-03) and 1.
#
# The fixtures are BUILT here rather than stored as .rda. That is a deliberate
# departure from an earlier draft of the plan, and it buys three things: the
# tarball carries no fixture bytes at all (section 2.8 of CRAN_READINESS.md is
# about a 10.1 MB .RData), the provenance of every number is the code below
# rather than a binary blob, and a test that needs a *fitted* object exercises
# the pipeline instead of reading its own answer back. Building is cached in
# `.fixture_cache` so the cost is paid once per `devtools::test()` run.
#
# Section 1.1 of the plan asks that F-TINY carry two all-zero genes. It does not
# here: `gene_status` is not implemented yet, so all-zero genes in the shared
# fixture would make every unrelated test fail for the same reason and drown the
# signal. `.tiny_counts_with_zero_genes()` supplies them to section 2.16 only.

.fixture_cache <- new.env(parent = emptyenv())

# Builds the raw inputs for F-TINY: 120 cells, 20 genes, 6 individuals
# (3 case / 3 control), rank-2 latent structure plus a factor and a numeric
# covariate. Small enough that a full fit runs in well under a second and large
# enough that the k = 2 factorization is not degenerate (question Q-FIX-1).
.build_tiny_data <- function(num_cells_per_individual = 20,
                             num_genes = 20,
                             num_individuals = 6,
                             seed_number = 10){
  set.seed(seed_number)

  n <- num_cells_per_individual * num_individuals
  p <- num_genes
  k <- 2

  individual_vec <- factor(rep(paste0("indiv_", seq_len(num_individuals)),
                               each = num_cells_per_individual))
  # First half of the individuals are controls, second half cases. Kept
  # contiguous so a test that needs "an individual" can take a block of rows.
  cc_by_individual <- rep(c(0, 1), each = num_individuals / 2)
  cc_vec <- rep(cc_by_individual, each = num_cells_per_individual)

  sex_by_individual <- rep(c("F", "M"), times = num_individuals / 2)
  sex_vec <- factor(rep(sex_by_individual, each = num_cells_per_individual))
  age_by_individual <- stats::rnorm(num_individuals, mean = 40, sd = 8)
  age_vec <- rep(age_by_individual, each = num_cells_per_individual)

  covariate_df <- data.frame(CC = factor(cc_vec, levels = c(0, 1)),
                             Sex = sex_vec,
                             Age = age_vec)
  rownames(covariate_df) <- paste0("cell_", seq_len(n))

  # The latent structure. Kept small in magnitude so exp() of the natural
  # parameter stays in a range where Poisson counts are neither all zero nor
  # astronomically large.
  x_mat <- matrix(stats::rnorm(n * k, sd = 0.4), nrow = n, ncol = k)
  y_mat <- matrix(stats::rnorm(p * k, sd = 0.4), nrow = p, ncol = k)

  # A modest per-gene intercept plus a case-control effect on the first five
  # genes only, so the fixture has both signal and null genes.
  gene_intercept_vec <- stats::runif(p, min = 1, max = 2.5)
  cc_effect_vec <- c(rep(0.8, 5), rep(0, p - 5))

  nat_mat <- tcrossprod(x_mat, y_mat)
  nat_mat <- nat_mat + rep(gene_intercept_vec, each = n)
  nat_mat <- nat_mat + outer(cc_vec, cc_effect_vec)

  rownames(nat_mat) <- rownames(covariate_df)
  # No underscore in the gene names: `SeuratObject::CreateSeuratObject()`
  # rewrites "gene_1" to "gene-1", and every name-based assertion in the
  # gene-status tests compares Seurat-derived output against `colnames(dat)`.
  colnames(nat_mat) <- paste0("gene", seq_len(p))

  # `nuisance_true_vec` is the Gamma RATE beta_j = 1/gamma_j, matching what the
  # package's `nuisance_vec` holds -- NOT the paper's over-dispersion gamma_j.
  # Getting this backwards would make T-NUIS-05 assert the reciprocal of the
  # truth and still look plausible.
  nuisance_true_vec <- stats::runif(p, min = 2, max = 8)
  library_size_vec <- stats::runif(n, min = 0.8, max = 1.2)

  # Generated from the eSVD hierarchical model directly rather than through
  # `generate_data`. `generate_data(family = "neg_binom2")` draws with a
  # CONSTANT `size`, whereas eSVD assumes lambda ~ Gamma(mean = mu, var =
  # gamma*mu), whose shape mu*beta varies cell by cell. Only the latter gives
  # `estimate_nuisance` a truth it is actually estimating.
  set.seed(seed_number + 1)
  mean_mat <- exp(nat_mat)
  shape_mat <- .mult_mat_vec(mean_mat, nuisance_true_vec)
  lambda_mat <- matrix(stats::rgamma(n * p,
                                     shape = as.numeric(shape_mat),
                                     rate = rep(nuisance_true_vec, each = n)),
                       nrow = n, ncol = p)
  dat <- matrix(stats::rpois(n * p,
                             lambda = as.numeric(.mult_vec_mat(library_size_vec,
                                                               lambda_mat))),
                nrow = n, ncol = p)
  rownames(dat) <- rownames(nat_mat)
  colnames(dat) <- colnames(nat_mat)

  # A gene that came out all-zero by chance would silently change what several
  # tests are testing, so force a floor rather than leaving it to the seed.
  zero_gene_idx <- which(Matrix::colSums(dat) == 0)
  if(length(zero_gene_idx) > 0){
    dat[1, zero_gene_idx] <- 1
  }

  list(cc_vec = cc_vec,
       covariate_df = covariate_df,
       dat = dat,
       individual_vec = individual_vec,
       library_size_vec = library_size_vec,
       nat_mat = nat_mat,
       nuisance_true_vec = nuisance_true_vec)
}

# F-TINY. Raw inputs only; anything fitted is built by `.small_esvd_obj()`.
.tiny_data <- function(){
  if(is.null(.fixture_cache$tiny)){
    .fixture_cache$tiny <- .build_tiny_data()
  }
  .fixture_cache$tiny
}

# F-SMALL: 400 cells, 40 genes, 8 individuals (4/4). Used where a stable Welch
# degree of freedom or a defensible recovery tolerance is needed (T-NUIS-05,
# T-PROP-07), which 120 cells cannot support.
.small_data <- function(){
  if(is.null(.fixture_cache$small)){
    .fixture_cache$small <- .build_tiny_data(num_cells_per_individual = 50,
                                             num_genes = 40,
                                             num_individuals = 8,
                                             seed_number = 20)
  }
  .fixture_cache$small
}

# The count matrix of F-TINY, as a dense matrix.
.tiny_counts <- function(){
  .tiny_data()$dat
}

# The count matrix of F-TINY as a dgCMatrix, for the sparse code paths.
.tiny_counts_sparse <- function(){
  methods::as(methods::as(.tiny_counts(), "dMatrix"), "CsparseMatrix")
}

# F-TINY with two genes forced to all-zero, for section 2.16. The zeros are put
# at positions 3 and 12 rather than at an end, so a reinsertion that appends
# instead of restoring position is caught.
.tiny_counts_with_zero_genes <- function(zero_idx = c(3, 12)){
  dat <- .tiny_counts()
  dat[, zero_idx] <- 0
  attr(dat, "zero_idx") <- zero_idx
  dat
}

# The formatted covariate matrix for a given count matrix and F-TINY metadata.
.tiny_covariates <- function(dat = .tiny_counts(),
                             covariate_df = .tiny_data()$covariate_df){
  format_covariates(dat = dat,
                    covariate_df = covariate_df,
                    rescale_numeric_variables = "Age")
}

# A fitted eSVD object built from F-TINY, run through the whole pipeline up to
# and including `compute_pvalue`. This is what roughly forty of the tests below
# need, and building it once is what keeps the suite inside the 90 s budget of
# convention C-08.
#
# `k = 2` matches the rank used to generate the data. `max_iter` is deliberately
# small: these tests assert structure and invariants, not convergence quality.
.small_esvd_obj <- function(max_iter = 10){
  cache_key <- paste0("esvd_obj_", max_iter)
  if(!is.null(.fixture_cache[[cache_key]])) return(.fixture_cache[[cache_key]])

  dat_list <- .tiny_data()
  covariates <- .tiny_covariates()

  esvd_obj <- suppressWarnings(
    initialize_esvd(dat = dat_list$dat,
                    covariates = covariates,
                    metadata_individual = dat_list$individual_vec,
                    bool_intercept = TRUE,
                    case_control_variable = "CC_1",
                    k = 2,
                    lambda = 0.1,
                    metadata_case_control = covariates[, "CC_1"],
                    verbose = 0)
  )

  esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                                 fit_name = "fit_Init",
                                                 omitted_variables = "Log_UMI")

  esvd_obj <- suppressWarnings(
    opt_esvd(input_obj = esvd_obj,
             l2pen = 0.1,
             max_iter = max_iter,
             offset_variables = setdiff(colnames(esvd_obj$covariates), "CC_1"),
             tol = 1e-6,
             fit_name = "fit_First",
             fit_previous = "fit_Init",
             verbose = 0)
  )

  esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                                 fit_name = "fit_First",
                                                 omitted_variables = "Log_UMI")

  esvd_obj <- suppressWarnings(
    estimate_nuisance(input_obj = esvd_obj,
                      bool_covariates_as_library = TRUE,
                      verbose = 0)
  )

  esvd_obj <- compute_posterior(input_obj = esvd_obj,
                                alpha_max = 2 * max(dat_list$dat),
                                bool_covariates_as_library = TRUE,
                                library_min = 0.1)

  esvd_obj <- compute_test_statistic(input_obj = esvd_obj, verbose = 0)
  esvd_obj <- suppressWarnings(compute_pvalue(input_obj = esvd_obj))

  .fixture_cache[[cache_key]] <- esvd_obj
  esvd_obj
}

# A feasible parameter point for each of the seven families, for section 3.2.
# The families have different natural-parameter domains -- `curved_gaussian`
# needs theta > 0, `exponential` and `neg_binom` need theta < 0, the rest are
# unconstrained -- so a single random point does not work for all seven. This is
# F-DERIV of the plan, built rather than stored.
.feasible_point <- function(family,
                            num_cells = 5,
                            num_genes = 4,
                            num_factors = 2,
                            seed_number = 10){
  set.seed(seed_number)

  n <- num_cells
  p <- num_genes
  k <- num_factors

  x_mat <- matrix(stats::runif(n * k, min = 0.3, max = 0.9), nrow = n, ncol = k)
  y_mat <- matrix(stats::runif(p * k, min = 0.3, max = 0.9), nrow = p, ncol = k)

  # Flip the sign of the gene loadings for the families whose natural parameter
  # must be negative. `tcrossprod` of two positive matrices is positive.
  if(family %in% c("exponential", "neg_binom")){
    y_mat <- -y_mat
  }

  nat_mat <- tcrossprod(x_mat, y_mat)

  # `gamma` is consumed by four of the seven families. Supplying it for all
  # seven costs nothing and avoids the NA objective of T-OPT-05.
  gamma_vec <- stats::runif(p, min = 1, max = 3)
  s_vec <- rep(1, n)

  set.seed(seed_number + 1)
  dat <- switch(family,
                gaussian = matrix(stats::rnorm(n * p), nrow = n, ncol = p),
                curved_gaussian = matrix(stats::runif(n * p, 0.5, 2),
                                         nrow = n, ncol = p),
                exponential = matrix(stats::rexp(n * p), nrow = n, ncol = p),
                poisson = matrix(stats::rpois(n * p, lambda = 3),
                                 nrow = n, ncol = p),
                neg_binom = matrix(stats::rpois(n * p, lambda = 3),
                                   nrow = n, ncol = p),
                neg_binom2 = matrix(stats::rpois(n * p, lambda = 3),
                                    nrow = n, ncol = p),
                bernoulli = matrix(stats::rbinom(n * p, size = 1, prob = 0.5),
                                   nrow = n, ncol = p))
  dat <- dat * 1.0

  list(dat = dat,
       gamma_vec = gamma_vec,
       nat_mat = nat_mat,
       s_vec = s_vec,
       x_mat = x_mat,
       y_mat = y_mat)
}

# The seven families, in one place so a sweep cannot silently miss one.
.all_families <- function(){
  c("gaussian", "curved_gaussian", "exponential", "poisson", "neg_binom",
    "neg_binom2", "bernoulli")
}

# Builds a minimal SeuratObject around a count matrix and the F-TINY metadata,
# for the tests in sections 2.12 and 2.17 that go through `eSVD()` or
# `eSVD_helper()` -- both of which take `seurat_obj` rather than a matrix.
#
# `dat` here is cells-by-genes, matching what the rest of the suite uses;
# Seurat wants genes-by-cells, hence the transpose.
.tiny_seurat <- function(dat = .tiny_counts(),
                         covariate_df = .tiny_data()$covariate_df,
                         individual_vec = .tiny_data()$individual_vec){
  skip_if_not_installed("SeuratObject")

  count_mat <- Matrix::t(methods::as(methods::as(dat * 1.0, "dMatrix"),
                                     "CsparseMatrix"))
  colnames(count_mat) <- rownames(dat)
  rownames(count_mat) <- colnames(dat)

  meta_df <- covariate_df
  meta_df$Individual <- individual_vec
  rownames(meta_df) <- rownames(dat)

  seurat_obj <- SeuratObject::CreateSeuratObject(counts = count_mat,
                                                 meta.data = meta_df,
                                                 min.cells = 0,
                                                 min.features = 0)
  seurat_obj
}

# Thin wrapper so the section 2.16 and 2.17 tests read as one call. It routes
# through `eSVD_helper()`, which under question Q-COH-7 is where all the
# filtering, gene-status labelling and reinsertion live.
.helper_run <- function(dat, k = 2, ...){
  seurat_obj <- .tiny_seurat(dat = dat)

  eSVD_helper(batch_var_prefix = NULL,
              case_control_levels = c("0", "1"),
              case_control_var = "CC",
              categorical_vars = c("Sex"),
              id_var = "Individual",
              numerical_vars = "Age",
              seurat_obj = seurat_obj,
              k = k,
              ...)
}

# The same, through `eSVD()` directly, for the tests that must distinguish the
# wrapper's behaviour from the pipeline's.
.esvd_run <- function(dat, k = 2, ...){
  seurat_obj <- .tiny_seurat(dat = dat)

  eSVD(batch_var_prefix = NULL,
       case_control_levels = c("0", "1"),
       case_control_var = "CC",
       categorical_vars = c("Sex"),
       id_var = "Individual",
       numerical_vars = "Age",
       seurat_obj = seurat_obj,
       k = k,
       ...)
}
