# Shared helpers for the master (3d5f7bf) vs devel (1.1.0) comparison
# Drafted by Claude for Kevin Z. Lin, 2026-09-29
#
# Sourced by 01_simulate_data_claude.R, 02_run_regimes_claude.R and
# 03_corner_cases_claude.R. The simulators need no eSVD2 except
# `.simulate_generate_null()`; the pipeline helpers call whichever eSVD2 the
# calling script attached.

# The generator ----------------------------------------------------------------

# The generator of `.build_tiny_data()` in tests/testthat/helper-fixtures.R,
# with the knobs this comparison sweeps exposed as arguments. Counts follow the
# eSVD hierarchical model directly: lambda ~ Gamma(mean = mu, rate = beta_j),
# then A ~ Poisson(library * lambda). `nuisance_range` is the range of that
# RATE beta_j, the quantity the package stores in `nuisance_vec`; a large rate
# is a nearly Poisson gene.
.simulate_cohort <- function(bool_covariates = TRUE,
                             cc_effect_vec = NULL,
                             gene_intercept_range = c(0, 1.5),
                             k = 3,
                             num_cells_per_individual = 40,
                             num_de = 15,
                             num_genes = 150,
                             num_individuals = 20,
                             nuisance_range = c(2, 8),
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
    cc_effect_vec <- c(rep(0.8, num_de), rep(0, p - num_de))
  }
  stopifnot(length(cc_effect_vec) == p)

  # Individual-level random effects, so that the cells of one donor are not
  # independent replicates and a test at the individual level is the right one.
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

  nuisance_true_vec <- exp(stats::runif(p,
                                        min = log(nuisance_range[1]),
                                        max = log(nuisance_range[2])))
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

  # An all-zero gene is refused by devel and is its own corner case (script
  # 03); keep it out of the head-to-head regimes.
  zero_gene_idx <- which(colSums(dat) == 0)
  if(length(zero_gene_idx) > 0) dat[1, zero_gene_idx] <- 1

  # The truth is a ratio of arithmetic means of exp(nat), which is what both
  # the logFC and the Welch statistic contrast (see CLAUDE_kevin.md).
  truth_log2fc_vec <- log2(colMeans(mean_mat[cc_vec == 1, , drop = FALSE]) /
                             colMeans(mean_mat[cc_vec == 0, , drop = FALSE]))

  list(cc_vec = cc_vec,
       covariate_df = covariate_df,
       dat = dat,
       gene_vec = gene_vec,
       individual_vec = individual_vec,
       is_de_vec = cc_effect_vec != 0,
       nuisance_true_vec = stats::setNames(nuisance_true_vec, gene_vec),
       truth_log2fc_vec = truth_log2fc_vec)
}

# The package's own generator, reshaped into the list `.simulate_cohort()`
# returns. Its covariates arrive already formatted (Log_UMI = log1p of the cell
# total) and its genes 1-10 carry the planted effect. It has no per-gene rate
# truth: its dispersion enters differently from the eSVD model.
.simulate_generate_null <- function(num_cells_per_individual,
                                    num_genes,
                                    num_individuals,
                                    seed_number){
  set.seed(seed_number)
  sim <- eSVD2::generate_null(cell_per_person = num_cells_per_individual,
                              num_genes = num_genes,
                              num_individuals = num_individuals)
  count_mat <- as.matrix(sim$obs_mat)
  gene_vec <- paste0("gene", seq_len(ncol(count_mat)))
  colnames(count_mat) <- gene_vec
  cc_vec <- sim$covariates[, "CC"]
  mean_mat <- exp(sim$nat_mat)

  list(cc_vec = cc_vec,
       covariates = sim$covariates,
       dat = count_mat,
       gene_vec = gene_vec,
       individual_vec = sim$metadata_individual,
       is_de_vec = seq_len(num_genes) <= 10,
       nuisance_true_vec = NULL,
       truth_log2fc_vec = stats::setNames(
         log2(colMeans(mean_mat[cc_vec == 1, , drop = FALSE]) /
                colMeans(mean_mat[cc_vec == 0, , drop = FALSE])),
         gene_vec))
}

# The pipeline -----------------------------------------------------------------

# master's compute_pvalue() computes the log p-values with Rmpfr::pnorm() inside
# sapply(), which can return a list of mpfr scalars rather than a numeric vector.
.as_numeric <- function(x){
  if(is.list(x)) return(vapply(x, function(v){as.numeric(v)}, numeric(1)))
  as.numeric(x)
}

# The pipeline in the order `eSVD()` runs it (both versions), with every warning
# recorded and muffled so one noisy stage does not hide the next. An error is
# returned, not thrown, together with the stage that raised it, because where a
# version fails is itself a result of this comparison.
.run_pipeline <- function(dat,
                          covariates,
                          cc_var,
                          individual_vec,
                          k = 3,
                          max_iter = 50,
                          nuisance_override_vec = NULL){
  stage_val <- "format"
  warning_vec <- character(0)
  record_warning <- function(w){
    warning_vec <<- c(warning_vec, conditionMessage(w))
    invokeRestart("muffleWarning")
  }

  esvd_obj <- tryCatch(withCallingHandlers({
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
    obj <- eSVD2::reparameterization_esvd_covariates(input_obj = obj,
                                                     fit_name = "fit_Init",
                                                     omitted_variables = "Log_UMI")
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
    obj <- eSVD2::reparameterization_esvd_covariates(input_obj = obj,
                                                     fit_name = "fit_First",
                                                     omitted_variables = "Log_UMI")
    stage_val <- "opt_esvd (second)"
    obj <- eSVD2::opt_esvd(input_obj = obj,
                           l2pen = 0.1,
                           max_iter = max_iter,
                           offset_variables = NULL,
                           tol = 1e-6,
                           fit_name = "fit_Second",
                           fit_previous = "fit_First",
                           verbose = 0)
    obj <- eSVD2::reparameterization_esvd_covariates(input_obj = obj,
                                                     fit_name = "fit_Second",
                                                     omitted_variables = NULL)
    stage_val <- "estimate_nuisance"
    obj <- eSVD2::estimate_nuisance(input_obj = obj,
                                    bool_covariates_as_library = TRUE,
                                    verbose = 0)
    if(!is.null(nuisance_override_vec)){
      stopifnot(length(nuisance_override_vec) == ncol(dat))
      obj[[obj[["latest_Fit"]]]]$nuisance_vec[] <- nuisance_override_vec
    }
    stage_val <- "compute_posterior"
    obj <- eSVD2::compute_posterior(input_obj = obj,
                                    alpha_max = 2 * max(dat),
                                    bool_covariates_as_library = TRUE,
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

# Pulls the per-gene results out of a fitted object. Only the fields both
# versions store are required; the logFC standard error exists only in devel.
.extract_gene_df <- function(esvd_obj,
                             gene_vec){
  latest_fit <- esvd_obj[["latest_Fit"]]
  pvalue_list <- esvd_obj$pvalue_list
  logfc_se_vec <- esvd_obj$log2fc_se_vec
  if(is.null(logfc_se_vec)) logfc_se_vec <- rep(NA, length(gene_vec))

  data.frame(gene = gene_vec,
             case_mean = as.numeric(esvd_obj$case_mean),
             control_mean = as.numeric(esvd_obj$control_mean),
             logFC = log2(as.numeric(esvd_obj$case_mean) /
                            as.numeric(esvd_obj$control_mean)),
             logFC_se = as.numeric(logfc_se_vec),
             teststat = as.numeric(esvd_obj$teststat_vec),
             df = as.numeric(pvalue_list$df_vec),
             gaussian_teststat = .as_numeric(pvalue_list$gaussian_teststat),
             log10p = .as_numeric(pvalue_list$log10pvalue),
             fdr = as.numeric(pvalue_list$fdr_vec),
             nuisance = as.numeric(esvd_obj[[latest_fit]]$nuisance_vec))
}

