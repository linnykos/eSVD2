# Shared helpers for the 2026-11-07 standard-error exploration
# Drafted by Claude for Kevin Z. Lin, 2026-11-07
#
# Sourced by 01_real-data-checks_claude.R, 02_replicate-cohorts_claude.R and
# 03_figure4-pairs_claude.R. See README.md in this folder for the design.
# Every logFC and SE below is on the log2 scale.

# Loading eSVD2 1.2.0 ----------------------------------------------------------

# The system library holds eSVD2 1.0.1.2; the version Tati ran (1.2.0) is the
# private library the version comparison built from the working tree.
.load_esvd2_devel <- function(exploration_dir){
  lib_dir <- file.path(exploration_dir, "..", "version_comparison", "lib",
                       "devel")
  stopifnot(dir.exists(lib_dir))
  .libPaths(c(normalizePath(lib_dir), .libPaths()))
  suppressPackageStartupMessages(library(eSVD2))
  version_val <- as.character(utils::packageVersion("eSVD2"))
  if(version_val != "1.2.0"){
    stop("expected eSVD2 1.2.0 from `", lib_dir, "`, found ", version_val)
  }
  invisible(version_val)
}

# The generator ----------------------------------------------------------------

# Gene-level truth shared by every cohort drawn from the same population: the
# gene intercepts, the gene loadings of the latent factors, the Gamma rates
# (`nuisance_range` is the range of the RATE beta_j, the quantity eSVD2 stores
# in `nuisance_vec`; a large rate is a nearly Poisson gene), and the
# case-control effect on the natural-log scale.
#
# `cc_effect_vec` may be supplied (as by `.simulate_effect_pair()`); otherwise
# `num_de` genes carry an effect of +/- `de_effect` with signs balanced, or all
# positive when `bool_asymmetric` (the compositional scenario).
.simulate_population <- function(bool_asymmetric = FALSE,
                                 cc_effect_vec = NULL,
                                 de_effect = 0.8,
                                 gene_intercept_range = c(0, 1.5),
                                 k_true = 3,
                                 num_de = 30,
                                 num_genes = 300,
                                 nuisance_range = c(2, 8),
                                 seed_number = 10){
  set.seed(seed_number)
  p <- num_genes
  gene_vec <- paste0("gene", seq_len(p))
  gene_intercept_vec <- stats::runif(p,
                                     min = gene_intercept_range[1],
                                     max = gene_intercept_range[2])
  y_mat <- matrix(stats::rnorm(p * k_true, sd = 0.4), nrow = p, ncol = k_true)
  nuisance_true_vec <- exp(stats::runif(p,
                                        min = log(nuisance_range[1]),
                                        max = log(nuisance_range[2])))

  if(is.null(cc_effect_vec)){
    cc_effect_vec <- rep(0, p)
    if(num_de > 0){
      sign_vec <- if(bool_asymmetric){
        rep(1, num_de)
      } else {
        rep(c(1, -1), length.out = num_de)
      }
      cc_effect_vec[seq_len(num_de)] <- sign_vec * de_effect
    }
  }
  stopifnot(length(cc_effect_vec) == p)
  names(cc_effect_vec) <- gene_vec
  names(nuisance_true_vec) <- gene_vec
  names(gene_intercept_vec) <- gene_vec

  list(cc_effect_vec = cc_effect_vec,
       gene_intercept_vec = gene_intercept_vec,
       gene_vec = gene_vec,
       k_true = k_true,
       nuisance_true_vec = nuisance_true_vec,
       truth_log2fc_vec = cc_effect_vec / log(2),
       y_mat = y_mat)
}

# Two populations A and B whose effects have a known ratio of sizes `k_true`
# (S_A / S_B) and a known correlation `r_true`, on the DE genes. Both share
# the intercepts, loadings and rates; only the case-control effects differ.
# The population-level S_A^2 is mean_g(beta_A^2) over ALL genes, which the
# caller reads off `truth_log2fc_vec`.
.simulate_effect_pair <- function(effect_sd_a = 0.6,
                                  k_true = 2,
                                  num_de = 45,
                                  num_genes = 300,
                                  r_true = 0.8,
                                  seed_number = 10,
                                  ...){
  set.seed(seed_number + 1000)
  shared_vec <- stats::rnorm(num_de)
  own_a_vec <- stats::rnorm(num_de)
  own_b_vec <- stats::rnorm(num_de)
  effect_a_vec <- rep(0, num_genes)
  effect_b_vec <- rep(0, num_genes)
  effect_a_vec[seq_len(num_de)] <- effect_sd_a *
    (sqrt(r_true) * shared_vec + sqrt(1 - r_true) * own_a_vec)
  effect_b_vec[seq_len(num_de)] <- (effect_sd_a / k_true) *
    (sqrt(r_true) * shared_vec + sqrt(1 - r_true) * own_b_vec)

  population_a <- .simulate_population(cc_effect_vec = effect_a_vec,
                                       num_genes = num_genes,
                                       seed_number = seed_number,
                                       ...)
  population_b <- population_a
  population_b$cc_effect_vec <- stats::setNames(effect_b_vec,
                                                population_a$gene_vec)
  population_b$truth_log2fc_vec <- population_b$cc_effect_vec / log(2)

  list(population_a = population_a,
       population_b = population_b)
}

# One cohort from a population: `num_individuals` donors (half case, half
# control, so the count must be even), each with `num_cells_per_individual`
# cells. Counts follow the eSVD hierarchical model directly:
# lambda ~ Gamma(mean = mu, rate = beta_j), A ~ Poisson(library * lambda),
# with donor random effects (sd `indiv_sd`) so that cells of one donor are
# not independent replicates. The truth for a ratio of arithmetic means is
# `cc_effect_vec / log(2)` exactly, because the latent factors, the donor
# effects and the library sizes are independent of the arm.
.simulate_cohort <- function(population,
                             indiv_sd = 0.1,
                             library_range = c(0.5, 1.5),
                             num_cells_per_individual = 50,
                             num_individuals = 20,
                             seed_number = 10){
  stopifnot(num_individuals %% 2 == 0)
  set.seed(seed_number)
  p <- length(population$gene_vec)
  n <- num_cells_per_individual * num_individuals
  individual_vec <- factor(rep(paste0("indiv", seq_len(num_individuals)),
                               each = num_cells_per_individual))
  cc_by_individual <- rep(c(0, 1), each = num_individuals / 2)
  cc_vec <- rep(cc_by_individual, each = num_cells_per_individual)
  covariate_df <- data.frame(CC = factor(cc_vec, levels = c(0, 1)))
  rownames(covariate_df) <- paste0("cell", seq_len(n))

  x_mat <- matrix(stats::rnorm(n * population$k_true, sd = 0.4),
                  nrow = n, ncol = population$k_true)
  indiv_effect_mat <- matrix(stats::rnorm(num_individuals * p, sd = indiv_sd),
                             nrow = num_individuals, ncol = p)
  nat_mat <- tcrossprod(x_mat, population$y_mat)
  nat_mat <- nat_mat + rep(population$gene_intercept_vec, each = n)
  nat_mat <- nat_mat + outer(cc_vec, population$cc_effect_vec)
  nat_mat <- nat_mat + indiv_effect_mat[as.integer(individual_vec), ,
                                        drop = FALSE]
  dimnames(nat_mat) <- list(rownames(covariate_df), population$gene_vec)
  library_size_vec <- stats::runif(n, min = library_range[1],
                                   max = library_range[2])

  set.seed(seed_number + 1)
  mean_mat <- exp(nat_mat)
  rate_vec <- rep(population$nuisance_true_vec, each = n)
  lambda_mat <- matrix(stats::rgamma(n * p,
                                     shape = as.numeric(mean_mat) * rate_vec,
                                     rate = rate_vec),
                       nrow = n, ncol = p)
  dat <- matrix(stats::rpois(n * p, lambda = library_size_vec * lambda_mat),
                nrow = n, ncol = p)
  dimnames(dat) <- dimnames(nat_mat)
  # An all-zero gene is refused by 1.2.0; keep it out of the comparison.
  zero_gene_idx <- which(colSums(dat) == 0)
  if(length(zero_gene_idx) > 0) dat[1, zero_gene_idx] <- 1

  list(cc_vec = cc_vec,
       covariate_df = covariate_df,
       dat = dat,
       gene_vec = population$gene_vec,
       individual_vec = individual_vec,
       is_de_vec = population$cc_effect_vec != 0,
       truth_log2fc_vec = population$truth_log2fc_vec)
}

# The eSVD pipeline ------------------------------------------------------------

# The stages in the order `eSVD()` runs them under 1.2.0 (the matrix path,
# `bool_diet = FALSE`, so that the posterior matrices are available), with
# every warning recorded and muffled. `bool_return_components = TRUE` keeps
# the posterior denominator `library + nuisance`, from which the per-cell
# library size and the shrinkage weight are recovered.
.run_esvd_pipeline <- function(dat,
                               covariate_df,
                               individual_vec,
                               cap_multiplier = 10,
                               k = 3,
                               max_iter = 100){
  warning_vec <- character(0)
  record_warning <- function(w){
    warning_vec <<- c(warning_vec, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  stage_val <- "format_covariates"
  esvd_obj <- tryCatch(withCallingHandlers({
    covariates <- eSVD2::format_covariates(dat = dat,
                                           covariate_df = covariate_df)
    cc_var <- "CC_1"
    stage_val <- "initialize_esvd"
    obj <- eSVD2::initialize_esvd(dat = dat,
                                  covariates = covariates,
                                  metadata_individual = individual_vec,
                                  bool_intercept = TRUE,
                                  case_control_variable = cc_var,
                                  k = k,
                                  lambda = 0.1,
                                  metadata_case_control = covariates[, cc_var],
                                  verbose = 0)
    obj <- eSVD2::reparameterization_esvd_covariates(
      input_obj = obj,
      fit_name = "fit_Init",
      omitted_variables = "Log_UMI"
    )
    stage_val <- "opt_esvd (first)"
    obj <- eSVD2::opt_esvd(input_obj = obj,
                           l2pen = 0.1,
                           max_iter = max_iter,
                           offset_variables = setdiff(colnames(obj$covariates),
                                                      cc_var),
                           tol = 1e-6,
                           fit_name = "fit_First",
                           fit_previous = "fit_Init",
                           verbose = 0)
    obj <- eSVD2::reparameterization_esvd_covariates(
      input_obj = obj,
      fit_name = "fit_First",
      omitted_variables = "Log_UMI"
    )
    stage_val <- "opt_esvd (second)"
    obj <- eSVD2::opt_esvd(input_obj = obj,
                           l2pen = 0.1,
                           max_iter = max_iter,
                           offset_variables = NULL,
                           tol = 1e-6,
                           fit_name = "fit_Second",
                           fit_previous = "fit_First",
                           verbose = 0)
    obj <- eSVD2::reparameterization_esvd_covariates(
      input_obj = obj,
      fit_name = "fit_Second",
      omitted_variables = NULL
    )
    stage_val <- "estimate_nuisance"
    obj <- eSVD2::estimate_nuisance(input_obj = obj,
                                    bool_covariates_as_library = TRUE,
                                    cap_multiplier = cap_multiplier,
                                    verbose = 0)
    stage_val <- "compute_posterior"
    obj <- eSVD2::compute_posterior(input_obj = obj,
                                    alpha_max = 2 * max(dat),
                                    bool_covariates_as_library = TRUE,
                                    bool_return_components = TRUE,
                                    library_min = 0.1)
    stage_val <- "compute_test_statistic"
    obj <- eSVD2::compute_test_statistic(input_obj = obj, verbose = 0)
    stage_val <- "compute_pvalue"
    eSVD2::compute_pvalue(input_obj = obj)
  }, warning = record_warning),
  error = function(e){e})

  if(inherits(esvd_obj, "error")){
    return(list(error_message = conditionMessage(esvd_obj),
                esvd_obj = NULL,
                stage = stage_val,
                warning_vec = warning_vec))
  }
  list(error_message = NULL,
       esvd_obj = esvd_obj,
       stage = "done",
       warning_vec = warning_vec)
}

# The nuisance rate as the posterior used it: floored at its 1% quantile and
# recentred when the geometric mean exceeds 1, the two steps of
# `compute_posterior.default()` under its defaults.
.posterior_nuisance <- function(esvd_obj){
  latest_fit <- esvd_obj[["latest_Fit"]]
  nuisance_vec <- esvd_obj[[latest_fit]]$nuisance_vec
  nuisance_vec <- pmax(nuisance_vec, stats::quantile(nuisance_vec, 0.01))
  if(mean(log10(nuisance_vec)) > 0){
    nuisance_vec <- 10^(log10(nuisance_vec) - mean(log10(nuisance_vec)))
  }
  nuisance_vec
}

# Per-donor summaries of one fitted object, the intermediate every candidate
# needs: for each donor and gene the mean posterior mean `m_ig`, the mean
# posterior variance, the mean of the unshrunk `A / library`, and the mean
# shrinkage weight `library / (library + nuisance)`, which is the share of
# the posterior mean that comes from the donor's own counts.
.donor_summaries <- function(esvd_obj){
  latest_fit <- esvd_obj[["latest_Fit"]]
  fit <- esvd_obj[[latest_fit]]
  stopifnot(!is.null(fit$denominator_mat), !is.null(esvd_obj$dat))
  dat <- as.matrix(esvd_obj$dat)
  nuisance_vec <- .posterior_nuisance(esvd_obj)
  library_mat <- sweep(fit$denominator_mat, MARGIN = 2, STATS = nuisance_vec,
                       FUN = "-")
  stopifnot(all(library_mat > 0))
  raw_mat <- dat / library_mat
  weight_mat <- library_mat / fit$denominator_mat

  case_individuals <- esvd_obj$param$test_case_individuals
  control_individuals <- esvd_obj$param$test_control_individuals
  idx_list <- eSVD2:::.determine_individual_indices(
    case_individuals = case_individuals,
    control_individuals = control_individuals,
    individual_vec = esvd_obj$individual
  )
  all_idx_list <- c(idx_list$case_indiv_idx, idx_list$control_indiv_idx)
  avg_mat <- eSVD2:::.construct_averaging_matrix(idx_list = all_idx_list,
                                                 n = nrow(dat))
  donor_vec <- c(as.character(case_individuals),
                 as.character(control_individuals))
  arm_vec <- c(rep(1, length(case_individuals)),
               rep(0, length(control_individuals)))
  to_matrix <- function(mat){
    out <- as.matrix(avg_mat %*% mat)
    dimnames(out) <- list(donor_vec, colnames(dat))
    out
  }

  list(arm_vec = arm_vec,
       donor_vec = donor_vec,
       mean_mat = to_matrix(fit$posterior_mean_mat),
       num_cells_vec = stats::setNames(sapply(all_idx_list, length), donor_vec),
       raw_mat = to_matrix(raw_mat),
       var_mat = to_matrix(fit$posterior_var_mat),
       weight_mat = to_matrix(weight_mat))
}

# The log2 fold change of a ratio of arm means of per-donor values, with the
# delta-method SE from a per-arm variance `var_fn(values of the arm)` that
# is already the variance OF THE DONOR MEAN's parent (the Welch convention:
# divided by the number of donors here).
.lfc_from_donor_means <- function(mean_mat, arm_vec, var_fn){
  case_mat <- mean_mat[arm_vec == 1, , drop = FALSE]
  control_mat <- mean_mat[arm_vec == 0, , drop = FALSE]
  n1 <- nrow(case_mat)
  n0 <- nrow(control_mat)
  case_mean <- colMeans(case_mat)
  control_mean <- colMeans(control_mat)
  case_var <- var_fn(case_mat)
  control_var <- var_fn(control_mat)
  lfc_vec <- log2(case_mean / control_mean)
  se_vec <- sqrt(case_var / (n1 * case_mean^2) +
                   control_var / (n0 * control_mean^2)) / log(2)
  # The Welch statistic and its degrees of freedom on the same two variances.
  teststat_vec <- (case_mean - control_mean) /
    sqrt(case_var / n1 + control_var / n0)
  df_vec <- (case_var / n1 + control_var / n0)^2 /
    ((case_var / n1)^2 / (n1 - 1) + (control_var / n0)^2 / (n0 - 1))
  list(df_vec = df_vec,
       lfc_vec = lfc_vec,
       se_vec = se_vec,
       teststat_vec = teststat_vec)
}

.sample_var_by_column <- function(mat){
  apply(mat, 2, stats::var)
}

# Leave-one-donor-out SD of the log2 fold change, the fit held fixed.
.jackknife_se <- function(mean_mat, arm_vec){
  n <- nrow(mean_mat)
  lfc_mat <- t(sapply(seq_len(n), function(i){
    keep_mat <- mean_mat[-i, , drop = FALSE]
    keep_arm <- arm_vec[-i]
    log2(colMeans(keep_mat[keep_arm == 1, , drop = FALSE]) /
           colMeans(keep_mat[keep_arm == 0, , drop = FALSE]))
  }))
  centre_vec <- colMeans(lfc_mat)
  sqrt((n - 1) / n * colSums(sweep(lfc_mat, 2, centre_vec)^2))
}

# Every candidate triple from one fitted (non-diet) object, one row per gene.
# Columns are `<candidate>_lfc`, `<candidate>_se`, `<candidate>_pvalue`.
# The `stored` p-value is eSVD2's own (empirical null); the other p-values
# are Wald, `2 * pnorm(-|lfc / se|)`, so that the triple is self-consistent.
.compute_candidates <- function(esvd_obj){
  summaries <- .donor_summaries(esvd_obj)
  arm_vec <- summaries$arm_vec
  mean_mat <- summaries$mean_mat
  gene_vec <- colnames(mean_mat)
  pvalue_list <- esvd_obj$pvalue_list
  latest_fit <- esvd_obj[["latest_Fit"]]
  wald <- function(lfc_vec, se_vec){2 * stats::pnorm(-abs(lfc_vec / se_vec))}

  stored_lfc <- as.numeric(esvd_obj$log2fc_vec[gene_vec])
  stored_se <- as.numeric(esvd_obj$log2fc_se_vec[gene_vec])
  null_mean <- pvalue_list$null_mean
  null_sd <- pvalue_list$null_sd

  between <- .lfc_from_donor_means(mean_mat, arm_vec, .sample_var_by_column)
  # The within term as the uncertainty of a donor mean would be: the mean
  # per-cell posterior variance of the donor divided by its cell count.
  within_donor_mat <- summaries$var_mat / summaries$num_cells_vec
  between_plus_var <- function(mat){
    rows <- match(rownames(mat), rownames(within_donor_mat))
    .sample_var_by_column(mat) +
      colMeans(within_donor_mat[rows, , drop = FALSE])
  }
  between_plus <- .lfc_from_donor_means(mean_mat, arm_vec, between_plus_var)
  raw <- .lfc_from_donor_means(summaries$raw_mat, arm_vec,
                               .sample_var_by_column)
  jackknife_se <- .jackknife_se(mean_mat, arm_vec)

  candidate_df <- data.frame(
    gene = gene_vec,
    nuisance = as.numeric(esvd_obj[[latest_fit]]$nuisance_vec[gene_vec]),
    nuisance_status = as.character(
      esvd_obj[[latest_fit]]$nuisance_status[gene_vec]),
    shrink_weight = colMeans(summaries$weight_mat),
    mean_count = colMeans(as.matrix(esvd_obj$dat))[gene_vec],
    welch_df = as.numeric(pvalue_list$df_vec[gene_vec]),
    gaussian_z = as.numeric(pvalue_list$gaussian_teststat[gene_vec]),
    null_mean = null_mean,
    null_sd = null_sd,
    stored_lfc = stored_lfc,
    stored_se = stored_se,
    stored_pvalue = 10^(-as.numeric(pvalue_list$log10pvalue[gene_vec])),
    null_rescaled_lfc = stored_lfc - null_mean * stored_se,
    null_rescaled_se = null_sd * stored_se,
    null_rescaled_pvalue = wald(stored_lfc - null_mean * stored_se,
                                null_sd * stored_se),
    between_lfc = between$lfc_vec,
    between_se = between$se_vec,
    between_pvalue = wald(between$lfc_vec, between$se_vec),
    between_welch_pvalue = 2 * stats::pt(-abs(between$teststat_vec),
                                         df = between$df_vec),
    between_plus_lfc = between_plus$lfc_vec,
    between_plus_se = between_plus$se_vec,
    between_plus_pvalue = wald(between_plus$lfc_vec, between_plus$se_vec),
    jackknife_lfc = stored_lfc,
    jackknife_se = jackknife_se,
    jackknife_pvalue = wald(stored_lfc, jackknife_se),
    raw_lfc = raw$lfc_vec,
    raw_se = raw$se_vec,
    raw_pvalue = wald(raw$lfc_vec, raw$se_vec),
    stringsAsFactors = FALSE
  )
  .add_derived_candidates(candidate_df)
}

# Candidates that are functions of the others, so that a cached table from
# an earlier run can be completed without a refit.
#   deshrunk: the between-donor SE divided by the gene's mean shrinkage
#   weight. The posterior pulls every donor mean toward the fit by the
#   factor `1 - w`, which shrinks the between-donor spread by about `w` but
#   leaves the arm contrast (and so the logFC's sampling variance) where the
#   unshrunk data put it; dividing by `w` undoes the first without touching
#   the second. Computable from a diet object (the weight needs only the
#   fitted library and the nuisance rate).
#   raw_null: option A applied on top of `raw`. The empirical null of
#   `multtest()` is fitted to the Wald statistics `raw_lfc / raw_se` of the
#   cohort (one fit per cohort, so `candidate_df` must hold one cohort), and
#   the triple is (lfc - null_mean * se, null_sd * se). If `raw` is
#   calibrated the null SD comes out near 1 and the null mean carries any
#   shared shift of the log fold changes, such as the depth adjustment's.
.add_derived_candidates <- function(candidate_df){
  if(is.null(candidate_df$deshrunk_se)){
    candidate_df$deshrunk_lfc <- candidate_df$stored_lfc
    candidate_df$deshrunk_se <- candidate_df$between_se /
      candidate_df$shrink_weight
    candidate_df$deshrunk_pvalue <- 2 * stats::pnorm(
      -abs(candidate_df$deshrunk_lfc / candidate_df$deshrunk_se))
  }
  if(is.null(candidate_df$raw_null_se)){
    z_vec <- candidate_df$raw_lfc / candidate_df$raw_se
    bool_ok <- is.finite(z_vec)
    # `multtest()` warns when locfdr falls back on its cruder estimators,
    # which it does on a few hundred genes; the estimate is still used.
    null_res <- suppressWarnings(eSVD2::multtest(z_vec[bool_ok]))
    candidate_df$raw_null_mean <- null_res$null_mean
    candidate_df$raw_null_sd <- null_res$null_sd
    candidate_df$raw_null_lfc <- candidate_df$raw_lfc -
      null_res$null_mean * candidate_df$raw_se
    candidate_df$raw_null_se <- null_res$null_sd * candidate_df$raw_se
    candidate_df$raw_null_pvalue <- 2 * stats::pnorm(
      -abs(candidate_df$raw_null_lfc / candidate_df$raw_null_se))
  }
  candidate_df
}

# The DESeq2 reference -----------------------------------------------------------

# DESeq2 on pseudobulk counts (one column per donor), design `~ CC`, the
# Wald triple Figure 4 was built for. Returns one row per gene with NA where
# DESeq2 gives none.
.run_deseq2_pseudobulk <- function(dat, cc_vec, individual_vec){
  individual_chr <- as.character(individual_vec)
  pseudobulk_mat <- t(rowsum(dat, group = individual_chr))
  individual_df <- unique(data.frame(Individual = individual_chr,
                                     CC = cc_vec,
                                     stringsAsFactors = FALSE))
  rownames(individual_df) <- individual_df$Individual
  individual_df <- individual_df[colnames(pseudobulk_mat), , drop = FALSE]
  individual_df$CC <- factor(individual_df$CC, levels = c(0, 1))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_mat,
                                        colData = individual_df,
                                        design = ~ CC)
  dds <- suppressMessages(suppressWarnings(DESeq2::DESeq(dds, quiet = TRUE)))
  res <- DESeq2::results(dds, name = "CC_1_vs_0")
  data.frame(gene = rownames(res),
             deseq2_lfc = res$log2FoldChange,
             deseq2_se = res$lfcSE,
             deseq2_pvalue = res$pvalue,
             stringsAsFactors = FALSE)
}

# The Figure 4 moments ---------------------------------------------------------

# Copied from `Writeup50_figure4-from-summary-statistics_claude.Rmd` so the
# numbers here are the numbers that Rmd would compute:
#   S_A^2 = mean(lfc_A^2) - mean(SE_A^2), rho_A = S_A^2 / mean(lfc_A^2),
#   k = S_A / S_B, r = mean(lfc_A lfc_B) / (S_A S_B),
# and slope = (1/k) rho_A r on the same genes.
.replication_slope <- function(lfc_a_vec, lfc_b_vec){
  sum(lfc_a_vec * lfc_b_vec) / sum(lfc_a_vec^2)
}

.decompose <- function(lfc_a_vec, lfc_b_vec, se_a_vec, se_b_vec){
  total_ms_a <- mean(lfc_a_vec^2)
  total_ms_b <- mean(lfc_b_vec^2)
  signal_ms_a <- total_ms_a - mean(se_a_vec^2)
  signal_ms_b <- total_ms_b - mean(se_b_vec^2)
  bool_valid <- signal_ms_a > 0 & signal_ms_b > 0
  signal_rms_a <- if(signal_ms_a > 0) sqrt(signal_ms_a) else NA_real_
  signal_rms_b <- if(signal_ms_b > 0) sqrt(signal_ms_b) else NA_real_
  cross_mean <- mean(lfc_a_vec * lfc_b_vec)

  data.frame(
    num_genes = length(lfc_a_vec),
    slope = .replication_slope(lfc_a_vec, lfc_b_vec),
    noise_rms_a = sqrt(mean(se_a_vec^2)),
    signal_rms_a = signal_rms_a,
    signal_rms_b = signal_rms_b,
    rho_a = signal_ms_a / total_ms_a,
    rho_b = signal_ms_b / total_ms_b,
    k = if(bool_valid) signal_rms_a / signal_rms_b else NA_real_,
    r = if(bool_valid) cross_mean / (signal_rms_a * signal_rms_b) else NA_real_,
    bool_valid = bool_valid
  )
}

# The Rmd's universe: finite in both cohorts, SE at or below the given
# quantile of each cohort's SEs.
.moments_universe <- function(lfc_a_vec, lfc_b_vec, se_a_vec, se_b_vec,
                              quantile_level = 0.95){
  bool_vec <- is.finite(lfc_a_vec) & is.finite(lfc_b_vec) &
    is.finite(se_a_vec) & is.finite(se_b_vec)
  if(quantile_level < 1){
    cutoff_a <- stats::quantile(se_a_vec[bool_vec], quantile_level)
    cutoff_b <- stats::quantile(se_b_vec[bool_vec], quantile_level)
    bool_vec <- bool_vec & se_a_vec <= cutoff_a & se_b_vec <= cutoff_b
  }
  bool_vec
}

# The truth the decomposition should recover, from the population effects
# (log2 scale) and the realized estimates of A, over the same genes.
.decompose_truth <- function(truth_a_vec, truth_b_vec, lfc_a_vec){
  signal_rms_a <- sqrt(mean(truth_a_vec^2))
  signal_rms_b <- sqrt(mean(truth_b_vec^2))
  data.frame(
    signal_rms_a = signal_rms_a,
    signal_rms_b = signal_rms_b,
    rho_a = mean(truth_a_vec^2) / mean(lfc_a_vec^2),
    k = signal_rms_a / signal_rms_b,
    r = mean(truth_a_vec * truth_b_vec) / (signal_rms_a * signal_rms_b),
    slope_truth = .replication_slope(truth_a_vec, truth_b_vec)
  )
}

# Bookkeeping ------------------------------------------------------------------

.candidate_vec <- function(){
  c("stored", "null_rescaled", "between", "between_plus", "jackknife",
    "deshrunk", "raw", "raw_null", "deseq2")
}

.append_progress <- function(progress_file, text){
  cat(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", text, "\n"),
      file = progress_file, append = TRUE)
}
