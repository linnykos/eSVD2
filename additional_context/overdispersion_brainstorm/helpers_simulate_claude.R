# A generator whose overdispersion depends on the gene's expression
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Sourced by 00_simulate_extra_claude.R. `.simulate_cohort()` in
# version_comparison/helpers_claude.R draws each gene's rate independently of
# its expression, so a trend between the two has nothing to find there. This
# is the same generator with one change: the rate follows DESeq2's parametric
# trend, written for the eSVD model.
#
# In the generator's units the library size averages 1, the count of gene j in
# cell i has mean m_ji = exp(nat_ji) and variance m_ji * (1 + 1 / beta_j), and
# the negative binomial size in the usual Var = m + alpha * m^2 form is
# 1 / alpha_ji = m_ji * beta_j. DESeq2's trend is
#   alpha_j = asymptotic + extra_poisson / m_j,
# which at the gene's typical mean m_j = exp(intercept_j) gives
#   beta_j = 1 / (asymptotic * m_j + extra_poisson).
# Highly expressed genes are therefore MORE overdispersed relative to Poisson
# (smaller rate), and lowly expressed genes are nearly Poisson. Each gene then
# departs from the trend by a log-normal factor.
.simulate_cohort_trend <- function(asymptotic_dispersion = 0.1,
                                   bool_covariates = TRUE,
                                   cc_effect_vec = NULL,
                                   extra_poisson = 0.02,
                                   gene_intercept_range = c(-2, 2),
                                   k = 3,
                                   num_cells_per_individual = 30,
                                   num_de = 30,
                                   num_genes = 300,
                                   num_individuals = 20,
                                   trend_noise_sd = 0.3,
                                   seed_number = 10){
  set.seed(seed_number)

  n <- num_cells_per_individual * num_individuals
  p <- num_genes

  individual_vec <- factor(rep(paste0("indiv", seq_len(num_individuals)),
                               each = num_cells_per_individual))
  cc_by_individual <- rep(c(0, 1), each = num_individuals / 2)
  cc_vec <- rep(cc_by_individual, each = num_cells_per_individual)

  covariate_df <- data.frame(CC = factor(cc_vec, levels = c(0, 1)))
  if(bool_covariates){
    sex_by_individual <- rep(c("F", "M"), times = num_individuals / 2)
    age_by_individual <- stats::rnorm(num_individuals, mean = 40, sd = 8)
    covariate_df$Sex <- factor(rep(sex_by_individual,
                                   each = num_cells_per_individual))
    covariate_df$Age <- rep(age_by_individual,
                            each = num_cells_per_individual)
  }
  rownames(covariate_df) <- paste0("cell", seq_len(n))

  x_mat <- matrix(stats::rnorm(n * k, sd = 0.4), nrow = n, ncol = k)
  y_mat <- matrix(stats::rnorm(p * k, sd = 0.4), nrow = p, ncol = k)
  gene_intercept_vec <- stats::runif(p,
                                     min = gene_intercept_range[1],
                                     max = gene_intercept_range[2])
  if(is.null(cc_effect_vec)){
    cc_effect_vec <- c(rep(c(0.5, -0.5), length.out = num_de),
                       rep(0, p - num_de))
  }
  stopifnot(length(cc_effect_vec) == p)

  indiv_effect_mat <- matrix(stats::rnorm(num_individuals * p, sd = 0.1),
                             nrow = num_individuals, ncol = p)
  nat_mat <- tcrossprod(x_mat, y_mat)
  nat_mat <- nat_mat + rep(gene_intercept_vec, each = n)
  nat_mat <- nat_mat + outer(cc_vec, cc_effect_vec)
  nat_mat <- nat_mat + indiv_effect_mat[as.integer(individual_vec), ,
                                        drop = FALSE]
  if(bool_covariates){
    sex_effect_vec <- stats::rnorm(p, sd = 0.2)
    age_effect_vec <- stats::rnorm(p, sd = 0.2)
    age_scaled_vec <- as.numeric(scale(covariate_df$Age))
    nat_mat <- nat_mat + outer(as.numeric(covariate_df$Sex == "M"),
                               sex_effect_vec)
    nat_mat <- nat_mat + outer(age_scaled_vec, age_effect_vec)
  }
  gene_vec <- paste0("gene", seq_len(p))
  dimnames(nat_mat) <- list(rownames(covariate_df), gene_vec)

  # The one change from `.simulate_cohort()`.
  trend_vec <- 1 / (asymptotic_dispersion * exp(gene_intercept_vec) +
                      extra_poisson)
  nuisance_true_vec <- trend_vec * exp(stats::rnorm(p, sd = trend_noise_sd))
  library_size_vec <- stats::runif(n, min = 0.5, max = 1.5)

  set.seed(seed_number + 1)
  mean_mat <- exp(nat_mat)
  lambda_mat <- matrix(stats::rgamma(n * p,
                                     shape = as.numeric(mean_mat) *
                                       rep(nuisance_true_vec, each = n),
                                     rate = rep(nuisance_true_vec, each = n)),
                       nrow = n, ncol = p)
  dat <- matrix(stats::rpois(n * p, lambda = library_size_vec * lambda_mat),
                nrow = n, ncol = p)
  dimnames(dat) <- dimnames(nat_mat)

  zero_gene_idx <- which(colSums(dat) == 0)
  if(length(zero_gene_idx) > 0) dat[1, zero_gene_idx] <- 1

  truth_log2fc_vec <- log2(colMeans(mean_mat[cc_vec == 1, , drop = FALSE]) /
                             colMeans(mean_mat[cc_vec == 0, , drop = FALSE]))

  list(cc_vec = cc_vec,
       covariate_df = covariate_df,
       dat = dat,
       gene_intercept_vec = stats::setNames(gene_intercept_vec, gene_vec),
       gene_vec = gene_vec,
       individual_vec = individual_vec,
       is_de_vec = cc_effect_vec != 0,
       nuisance_true_vec = stats::setNames(nuisance_true_vec, gene_vec),
       nuisance_trend_vec = stats::setNames(trend_vec, gene_vec),
       truth_log2fc_vec = truth_log2fc_vec)
}
