# Shared machinery of the second lfcse loop: bank location, design,
# generator, the eSVD-DE fit, and the public view. Drafted by Claude for
# Kevin Z. Lin, 2026-11-11, from autoresearch_lfcse/harness/generator_claude.R.
#
# Sourced by prepare_claude.R (the bank builder) and evaluate_claude.R (the
# scorer). FIXED during a run: a change to anything in this file invalidates
# the bank (rebuild, re-baseline, rewrite notes/BANK.md) and needs a new
# `.bank_version()`. The design is notes/MEMO_claude.md, section 3.
# Every logFC and SE is on the log2 scale; every generator effect is on the
# natural-log scale.

# Location and version ---------------------------------------------------------

.bank_version <- function(){
  "lfcse-bank-2"
}

# The bank lives outside the repository (which sits in Dropbox): the root is
# `ESVD2_BANK_DIR`, defaulting to R's per-user cache directory for eSVD2.
.bank_dir <- function(){
  root_dir <- Sys.getenv("ESVD2_BANK_DIR", unset = "")
  if(root_dir == "") root_dir <- tools::R_user_dir("eSVD2", which = "cache")
  file.path(root_dir, .bank_version())
}

# The system library holds eSVD2 1.0.1.2; the loop uses 1.2.0, the private
# library the version comparison built from the working tree.
.load_esvd2 <- function(repo_dir){
  lib_dir <- file.path(repo_dir, "additional_context", "version_comparison",
                       "lib", "devel")
  if(!dir.exists(lib_dir)) stop("no eSVD2 library at `", lib_dir, "`")
  .libPaths(c(normalizePath(lib_dir), .libPaths()))
  suppressPackageStartupMessages(library(eSVD2))
  version_val <- as.character(utils::packageVersion("eSVD2"))
  if(version_val != "1.2.0"){
    stop("expected eSVD2 1.2.0 from `", lib_dir, "`, found ", version_val)
  }
  invisible(version_val)
}

# Design -----------------------------------------------------------------------

# The rank every cohort is fitted at. The generator has 3 cell-state
# factors and, at `depth = error`, a fourth donor-level direction; 5 is
# above both, as a user's rank is above the unknown truth.
.fit_rank <- function(){
  5
}

.factor_names <- function(){
  c("donors", "cells", "expression", "effects", "covariates", "donor_sd",
    "state", "arm_sd", "depth")
}

# The nine factors of the memo. The full factorial has 1728 cells; a split
# holds a 1/6 fraction of 288: a third of the 27 combinations of the three
# 3-level factors (the sum of their level indices modulo 3 is 0 on dev and
# 1 on test) crossed with half of the 64 combinations of the six 2-level
# factors (the sum of their level indices is even on dev and odd on test).
# Every level of every factor, and every pair of levels of two factors,
# appears equally often within a split.
.bank_design <- function(split){
  stopifnot(split %in% c("dev", "test"))
  level_list <- list(donors = c(10, 20, 40),
                     effects = c("sparse", "sparse_up", "dense"),
                     state = c("none", "donor", "donor_arm"),
                     cells = c("equal", "unequal"),
                     expression = c("high", "low"),
                     covariates = c("none", "sexage"),
                     donor_sd = c(0.2, 0.5),
                     arm_sd = c("equal", "unequal"),
                     depth = c("none", "error"))
  index_df <- expand.grid(lapply(level_list, seq_along),
                          KEEP.OUT.ATTRS = FALSE)
  index_df$cell_idx <- seq_len(nrow(index_df))
  three_sum <- rowSums(index_df[, c("donors", "effects", "state")] - 1)
  two_sum <- rowSums(index_df[, c("cells", "expression", "covariates",
                                  "donor_sd", "arm_sd", "depth")] - 1)
  remainder <- if(split == "dev") 0 else 1
  index_df <- index_df[three_sum %% 3 == remainder &
                         two_sum %% 2 == remainder, ]
  design_df <- data.frame(cell_idx = index_df$cell_idx,
                          stringsAsFactors = FALSE)
  for(name in names(level_list)){
    design_df[[name]] <- level_list[[name]][index_df[[name]]]
  }
  design_df <- design_df[, c(.factor_names(), "cell_idx")]
  rownames(design_df) <- NULL
  design_df
}

# One seed per example from its cell and split. Within an example the
# population uses `seed` and `seed + 1` and the cohort `seed + 10` and
# `seed + 11`; cells are 100 apart so no two examples share a stream.
.example_seed <- function(cell_idx, split){
  split_idx <- match(split, c("dev", "test"))
  stopifnot(!is.na(split_idx))
  1000000 * split_idx + 100 * cell_idx
}

.example_id <- function(cell_idx, split){
  paste0(split, "_", formatC(cell_idx, width = 4, flag = "0"))
}

# The generator ----------------------------------------------------------------

# Gene-level truth of an example.
#
# Intercepts: uniform on `c(0, 1.5)` at `expression = high` (the first
#   bank's range) and on `c(-3.5, 0)` at `low`, where most genes have
#   under one count per cell and the cap on the nuisance rate binds for most
#   genes. A block of 15% of the genes sits one unit higher; at
#   `depth = error` it carries a donor-level shift (see the cohort).
# Loadings: `k_true` columns of SD 0.4, centred over the genes, so that a
#   shift of the cell states moves no gene on average.
# Rates: log-uniform on 0.5 to 200 (the first bank's `wide`).
# Case-control effects (the marginal effect, natural log):
#   sparse:    10% of genes with SD 0.6, the others exactly 0.
#   sparse_up: the same sizes, every effect positive (absolute values).
#   dense:     every gene nonzero; 85% with SD 0.1 and 15% with SD 0.3.
# Arm shift of the cell states (`state = donor_arm`): a vector of SD
#   `state_shift_sd` per factor, added to every case cell's state, so the
#   arms differ in cell state and every gene has a nonzero effect.
# Donor SD per gene: `donor_sd` times a log-normal factor of log SD 0.3
#   (the Velmeshev estimate varies across genes by about that much within a
#   quartile). At `state = none` all of it is idiosyncratic to the gene. At
#   `state = donor` and `donor_arm`, `state_share` of the variance is carried
#   on average by donor-level cell states (low rank across genes) and the
#   rest stays idiosyncratic.
.simulate_population <- function(effects,
                                 expression,
                                 state,
                                 donor_sd,
                                 bool_covariates,
                                 block_fraction = 0.15,
                                 covariate_effect_sd = 0.2,
                                 k_true = 3,
                                 num_genes = 600,
                                 state_share = 0.75,
                                 state_shift_sd = 0.25,
                                 seed_number = 10){
  stopifnot(effects %in% c("sparse", "sparse_up", "dense"),
            expression %in% c("high", "low"),
            state %in% c("none", "donor", "donor_arm"))
  p <- num_genes
  gene_vec <- paste0("gene", seq_len(p))
  gene_intercept_range <- if(expression == "high") c(0, 1.5) else c(-3.5, 0)

  set.seed(seed_number)
  gene_intercept_vec <- stats::runif(p,
                                     min = gene_intercept_range[1],
                                     max = gene_intercept_range[2])
  bool_block_vec <- seq_len(p) %in% sample(p)[seq_len(round(block_fraction * p))]
  gene_intercept_vec[bool_block_vec] <- gene_intercept_vec[bool_block_vec] + 1
  y_mat <- matrix(stats::rnorm(p * k_true, sd = 0.4), nrow = p, ncol = k_true)
  y_mat <- scale(y_mat, center = TRUE, scale = FALSE)
  attr(y_mat, "scaled:center") <- NULL
  nuisance_true_vec <- exp(stats::runif(p, min = log(0.5), max = log(200)))
  # Drawn whether or not the example uses them, so that the stream below
  # does not depend on the levels of `covariates` and `state`.
  sex_effect_vec <- stats::rnorm(p, sd = covariate_effect_sd)
  age_effect_vec <- stats::rnorm(p, sd = covariate_effect_sd)
  if(!bool_covariates){
    sex_effect_vec <- rep(0, p)
    age_effect_vec <- rep(0, p)
  }
  state_shift_vec <- stats::rnorm(k_true, sd = state_shift_sd)
  if(state != "donor_arm") state_shift_vec <- rep(0, k_true)
  donor_sd_vec <- donor_sd * exp(stats::rnorm(p, sd = 0.3))
  if(state == "none"){
    state_sd <- 0
    idio_sd_vec <- donor_sd_vec
  } else {
    # The states' variance for gene g is `state_sd^2 * sum(y_g^2)`; its
    # mean over the genes is set to `state_share` of `donor_sd^2`.
    state_sd <- sqrt(state_share * donor_sd^2 / mean(rowSums(y_mat^2)))
    idio_sd_vec <- sqrt(1 - state_share) * donor_sd_vec
  }

  set.seed(seed_number + 1)
  base_vec <- stats::rnorm(p)
  # Positions permuted by position, not by value (`sample(p)`).
  order_vec <- sample(p)
  if(effects == "dense"){
    effect_sd_vec <- rep(0.1, p)
    effect_sd_vec[order_vec[seq_len(round(0.15 * p))]] <- 0.3
  } else {
    effect_sd_vec <- rep(0, p)
    effect_sd_vec[order_vec[seq_len(round(0.10 * p))]] <- 0.6
  }
  cc_effect_vec <- effect_sd_vec * base_vec
  if(effects == "sparse_up") cc_effect_vec <- abs(cc_effect_vec)

  name_it <- function(vec){stats::setNames(vec, gene_vec)}
  rownames(y_mat) <- gene_vec
  list(age_effect_vec = name_it(age_effect_vec),
       bool_block_vec = name_it(bool_block_vec),
       cc_effect_vec = name_it(cc_effect_vec),
       gene_intercept_vec = name_it(gene_intercept_vec),
       gene_vec = gene_vec,
       idio_sd_vec = name_it(idio_sd_vec),
       k_true = k_true,
       nuisance_true_vec = name_it(nuisance_true_vec),
       sex_effect_vec = name_it(sex_effect_vec),
       state_sd = state_sd,
       state_shift_vec = state_shift_vec,
       y_mat = y_mat)
}

# Cells per donor. `equal` is 50 for every donor. `unequal` is log-normal
# with median 150; a log SD of 0.7 puts 30 and 600 about two SDs out
# (log(600 / 150) = 1.39), and draws beyond them are clipped there.
.simulate_cells_per_donor <- function(cells, num_donors){
  if(cells == "equal") return(rep(50L, num_donors))
  num_vec <- round(exp(stats::rnorm(num_donors, mean = log(150), sd = 0.7)))
  as.integer(pmin(pmax(num_vec, 30), 600))
}

# One cohort: `num_donors` donors, the first half controls and the second
# half cases. Counts follow the eSVD model,
# lambda ~ Gamma(mean = mu, rate = beta_j), A ~ Poisson(library * lambda),
# with log mu = intercept + x_c' y_g + (case-control term) + (donor terms) +
# (covariate terms). The state of a cell is
#   x_c = (its donor's state) + (the arm shift, for a case) + N(0, 0.4^2),
# and a donor acts on gene g through its state (low rank across genes),
# through an idiosyncratic effect of SD `idio_sd_vec[g]`, and, at
# `depth = error`, through one shift of SD `block_sd` common to the block
# genes, which moves the donor's total count and so every other gene's
# count relative to the total.
#
# At `arm_sd = unequal` the donor states and the idiosyncratic effects of
# the cases have 1.5 times, and those of the controls 0.5 times, the
# variance they have at `equal`.
#
# Every donor term is centred so that its exponential has mean 1 in each
# arm (the gene's log-normal mean correction is subtracted). The truth
#   (cc_effect_vec + y_mat %*% state_shift_vec) / log(2)
# is therefore the log2 ratio of mean expression between the populations of
# case and control donors, at equal depth and equal covariates, and a
# contrast of mean log expression differs from it wherever the arms differ
# in donor variance.
.simulate_cohort <- function(population,
                             cells,
                             num_donors,
                             bool_covariates,
                             bool_arm_sd_unequal,
                             bool_depth_error,
                             block_sd = 0.5,
                             library_range = c(0.5, 1.5),
                             seed_number = 10){
  stopifnot(num_donors %% 2 == 0)
  p <- length(population$gene_vec)
  k_true <- population$k_true
  num_arm <- num_donors / 2

  set.seed(seed_number)
  num_cells_vec <- .simulate_cells_per_donor(cells = cells,
                                             num_donors = num_donors)
  n <- sum(num_cells_vec)
  donor_vec <- paste0("indiv", seq_len(num_donors))
  cc_donor_vec <- rep(c(0, 1), each = num_arm)
  # Sex: round(0.4 * num_arm) of the controls and round(0.6 * num_arm) of
  # the cases, the donors chosen at random within the arm.
  sex_donor_vec <- rep(0, num_donors)
  control_pos <- sample(num_arm)[seq_len(round(0.4 * num_arm))]
  case_pos <- num_arm + sample(num_arm)[seq_len(round(0.6 * num_arm))]
  sex_donor_vec[c(control_pos, case_pos)] <- 1
  age_donor_vec <- as.numeric(scale(stats::rnorm(num_donors)))

  donor_idx_vec <- rep(seq_len(num_donors), times = num_cells_vec)
  cc_vec <- cc_donor_vec[donor_idx_vec]
  individual_vec <- factor(donor_vec[donor_idx_vec], levels = donor_vec)
  cell_vec <- paste0("cell", seq_len(n))
  covariate_df <- data.frame(CC = factor(cc_vec, levels = c(0, 1)))
  if(bool_covariates){
    covariate_df$Sex <- factor(sex_donor_vec[donor_idx_vec], levels = c(0, 1))
    covariate_df$Age <- age_donor_vec[donor_idx_vec]
  }
  rownames(covariate_df) <- cell_vec

  # The SD multiplier of each donor's effects.
  arm_multiplier_vec <- rep(1, num_donors)
  if(bool_arm_sd_unequal){
    arm_multiplier_vec <- ifelse(cc_donor_vec == 1, sqrt(1.5), sqrt(0.5))
  }
  donor_state_mat <- matrix(stats::rnorm(num_donors * k_true), nrow = num_donors,
                            ncol = k_true) *
    (population$state_sd * arm_multiplier_vec)
  idio_mat <- matrix(stats::rnorm(num_donors * p), nrow = num_donors,
                     ncol = p)
  idio_sd_mat <- outer(arm_multiplier_vec, as.numeric(population$idio_sd_vec))
  block_shift_vec <- stats::rnorm(num_donors, sd = block_sd)
  if(!bool_depth_error) block_shift_vec <- rep(0, num_donors)

  # The donor's whole effect on each gene, each part less its log-normal
  # mean correction.
  y_square_vec <- rowSums(population$y_mat^2)
  donor_effect_mat <- tcrossprod(donor_state_mat, population$y_mat) -
    outer((population$state_sd * arm_multiplier_vec)^2, y_square_vec) / 2 +
    idio_mat * idio_sd_mat - idio_sd_mat^2 / 2
  if(bool_depth_error){
    donor_effect_mat <- donor_effect_mat +
      outer(block_shift_vec - block_sd^2 / 2,
            as.numeric(population$bool_block_vec))
  }
  dimnames(donor_effect_mat) <- list(donor_vec, population$gene_vec)

  cell_state_mat <- matrix(stats::rnorm(n * k_true, sd = 0.4), nrow = n,
                           ncol = k_true)
  marginal_effect_vec <- as.numeric(population$cc_effect_vec) +
    as.numeric(population$y_mat %*% population$state_shift_vec)
  nat_mat <- tcrossprod(cell_state_mat, population$y_mat)
  nat_mat <- nat_mat + rep(population$gene_intercept_vec, each = n)
  nat_mat <- nat_mat + outer(cc_vec, marginal_effect_vec)
  nat_mat <- nat_mat + donor_effect_mat[donor_idx_vec, , drop = FALSE]
  if(bool_covariates){
    nat_mat <- nat_mat + outer(sex_donor_vec[donor_idx_vec],
                               as.numeric(population$sex_effect_vec))
    nat_mat <- nat_mat + outer(age_donor_vec[donor_idx_vec],
                               as.numeric(population$age_effect_vec))
  }
  library_size_vec <- stats::runif(n, min = library_range[1],
                                   max = library_range[2])

  set.seed(seed_number + 1)
  rate_vec <- rep(population$nuisance_true_vec, each = n)
  lambda_mat <- matrix(stats::rgamma(n * p,
                                     shape = as.numeric(exp(nat_mat)) * rate_vec,
                                     rate = rate_vec),
                       nrow = n, ncol = p)
  dat <- matrix(stats::rpois(n * p, lambda = library_size_vec * lambda_mat),
                nrow = n, ncol = p)
  dimnames(dat) <- list(cell_vec, population$gene_vec)
  # An all-zero gene is refused by 1.2.0; one count keeps it in the cohort.
  zero_gene_idx <- which(colSums(dat) == 0)
  if(length(zero_gene_idx) > 0) dat[1, zero_gene_idx] <- 1L
  storage.mode(dat) <- "integer"

  donor_df <- data.frame(donor = donor_vec,
                         cc = cc_donor_vec,
                         num_cells = num_cells_vec,
                         sex = sex_donor_vec,
                         age = age_donor_vec,
                         block_shift = block_shift_vec,
                         stringsAsFactors = FALSE)
  list(cc_vec = cc_vec,
       covariate_df = covariate_df,
       dat = dat,
       donor_df = donor_df,
       donor_effect_mat = donor_effect_mat,
       donor_state_mat = donor_state_mat,
       individual_vec = individual_vec,
       num_zero_genes = length(zero_gene_idx),
       truth_log2fc_vec = stats::setNames(marginal_effect_vec / log(2),
                                          population$gene_vec))
}

# The eSVD pipeline ------------------------------------------------------------

# The stages in the order `eSVD()` runs them under 1.2.0 (matrix path,
# nothing dropped), with every warning recorded and muffled. Copied from
# the SE exploration's `helpers_claude.R` so that the bank does not depend
# on a file outside the loop folder.
.run_esvd_pipeline <- function(dat,
                               covariate_df,
                               individual_vec,
                               cap_multiplier = 10,
                               k = 5,
                               max_iter = 100){
  warning_vec <- character(0)
  record_warning <- function(w){
    warning_vec <<- c(warning_vec, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  storage.mode(dat) <- "double"
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

# The public view --------------------------------------------------------------

# What the bank stores of one fitted cohort: everything of the memo's view
# except the three cells x genes matrices of doubles (library, posterior
# mean, posterior variance), which `.build_view()` recomputes from the
# stored fit with the same call the fit made. Storing them would put the
# bank near 20 GB. Nothing here is derived from the truth.
.slim_cohort <- function(esvd_obj, cc_vec){
  latest_fit <- esvd_obj[["latest_Fit"]]
  fit <- esvd_obj[[latest_fit]]
  gene_vec <- colnames(esvd_obj$dat)
  pvalue_list <- esvd_obj$pvalue_list
  stored_df <- data.frame(
    gene = gene_vec,
    log2fc = as.numeric(esvd_obj$log2fc_vec[gene_vec]),
    log2fc_se = as.numeric(esvd_obj$log2fc_se_vec[gene_vec]),
    teststat = as.numeric(esvd_obj$teststat_vec[gene_vec]),
    df = as.numeric(pvalue_list$df_vec[gene_vec]),
    gaussian_z = as.numeric(pvalue_list$gaussian_teststat[gene_vec]),
    pvalue = 10^(-as.numeric(pvalue_list$log10pvalue[gene_vec])),
    stringsAsFactors = FALSE
  )
  attr(stored_df, "null_mean") <- pvalue_list$null_mean
  attr(stored_df, "null_sd") <- pvalue_list$null_sd

  dat <- as.matrix(esvd_obj$dat)
  alpha_max <- 2 * max(dat)
  storage.mode(dat) <- "integer"
  list(alpha_max = alpha_max,
       case_control_vec = as.numeric(cc_vec),
       covariates = esvd_obj$covariates,
       dat = dat,
       individual_vec = esvd_obj$individual,
       nuisance_fit_vec = fit$nuisance_vec,
       nuisance_status_vec = fit$nuisance_status,
       stored_df = stored_df,
       x_mat = fit$x_mat,
       y_mat = fit$y_mat,
       z_mat = fit$z_mat)
}

# The view handed to the trial function, rebuilt from a slim cohort. The
# posterior call repeats the fit's own (`compute_posterior.default()` with
# the library formed from every covariate but the case-control one, the
# intercept included, floored at 0.1). `nuisance_vec` is the rate as the
# posterior used it: floored at its 1% quantile and recentred when the
# geometric mean exceeds 1; the library is the posterior's denominator less
# that rate.
.build_view <- function(slim){
  dat <- slim$dat
  storage.mode(dat) <- "double"
  nuisance_vec <- slim$nuisance_fit_vec
  nuisance_vec <- pmax(nuisance_vec, stats::quantile(nuisance_vec, 0.01))
  if(mean(log10(nuisance_vec)) > 0){
    nuisance_vec <- 10^(log10(nuisance_vec) - mean(log10(nuisance_vec)))
  }
  res <- eSVD2::compute_posterior(
    input_obj = dat,
    case_control_variable = "CC_1",
    covariates = slim$covariates,
    esvd_res = list(x_mat = slim$x_mat, y_mat = slim$y_mat,
                    z_mat = slim$z_mat),
    library_size_variable = "Log_UMI",
    nuisance_vec = slim$nuisance_fit_vec,
    alpha_max = slim$alpha_max,
    bool_covariates_as_library = TRUE,
    bool_library_includes_interept = TRUE,
    bool_return_components = TRUE,
    library_min = 0.1
  )
  library_mat <- sweep(res$denominator_mat, MARGIN = 2, STATS = nuisance_vec,
                       FUN = "-")
  dimnames(library_mat) <- dimnames(dat)
  stopifnot(all(library_mat > 0))

  list(case_control_vec = slim$case_control_vec,
       covariates = slim$covariates,
       dat = dat,
       individual_vec = slim$individual_vec,
       library_mat = library_mat,
       nuisance_status_vec = slim$nuisance_status_vec,
       nuisance_vec = nuisance_vec,
       posterior_mean_mat = res$posterior_mean_mat,
       posterior_var_mat = res$posterior_var_mat,
       stored_df = slim$stored_df,
       x_mat = slim$x_mat,
       y_mat = slim$y_mat,
       z_mat = slim$z_mat)
}

# The DESeq2 reference ---------------------------------------------------------

# DESeq2 on pseudobulk counts (one column per donor), design `~ CC`, or
# `~ Sex + Age + CC` when the cohort has covariates. One row per gene, NA
# where DESeq2 gives none.
.run_deseq2_pseudobulk <- function(dat, donor_df, individual_vec,
                                   bool_covariates){
  pseudobulk_mat <- t(rowsum(dat, group = as.character(individual_vec)))
  col_df <- donor_df[match(colnames(pseudobulk_mat), donor_df$donor), ,
                     drop = FALSE]
  rownames(col_df) <- col_df$donor
  col_df$CC <- factor(col_df$cc, levels = c(0, 1))
  col_df$Sex <- factor(col_df$sex, levels = c(0, 1))
  col_df$Age <- col_df$age
  design_formula <- if(bool_covariates) ~ Sex + Age + CC else ~ CC
  res <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_mat,
                                          colData = col_df,
                                          design = design_formula)
    dds <- suppressMessages(suppressWarnings(DESeq2::DESeq(dds,
                                                           quiet = TRUE)))
    DESeq2::results(dds, name = "CC_1_vs_0")
  }, error = function(e){e})
  gene_vec <- colnames(dat)
  if(inherits(res, "error")){
    return(data.frame(gene = gene_vec, deseq2_lfc = NA_real_,
                      deseq2_se = NA_real_, deseq2_pvalue = NA_real_,
                      stringsAsFactors = FALSE))
  }
  idx_vec <- match(gene_vec, rownames(res))
  data.frame(gene = gene_vec,
             deseq2_lfc = res$log2FoldChange[idx_vec],
             deseq2_se = res$lfcSE[idx_vec],
             deseq2_pvalue = res$pvalue[idx_vec],
             stringsAsFactors = FALSE)
}

# The donor_wls reference ------------------------------------------------------

# The rule the first loop kept (`candidates/donor_wls_claude.R`, a copy of
# autoresearch_lfcse/notes/RULE_nov11_claude.R), which reads the counts,
# the design and the fitted library and nothing of the embedding. It is out
# of bounds for this loop and is stored with the truth as a reference line.
.run_donor_wls <- function(view, loop_dir){
  rule_env <- new.env()
  sys.source(file.path(loop_dir, "candidates", "donor_wls_claude.R"),
             envir = rule_env)
  res <- tryCatch(rule_env$report_lfc_se(view), error = function(e){e})
  gene_vec <- colnames(view$dat)
  if(inherits(res, "error")){
    return(data.frame(gene = gene_vec, wls_lfc = NA_real_, wls_se = NA_real_,
                      stringsAsFactors = FALSE))
  }
  data.frame(gene = gene_vec, wls_lfc = res$log2fc, wls_se = res$se,
             stringsAsFactors = FALSE)
}

# One example ------------------------------------------------------------------

# Gene class by the true log2 fold change: exactly null, weak (below 0.5 in
# absolute value, the genes the validity gate regresses on), strong.
.gene_class <- function(truth_log2fc_vec){
  class_vec <- rep("strong", length(truth_log2fc_vec))
  class_vec[abs(truth_log2fc_vec) < 0.5] <- "weak"
  class_vec[truth_log2fc_vec == 0] <- "null"
  stats::setNames(class_vec, names(truth_log2fc_vec))
}

# Generates, fits and summarizes one example (one cohort). Returns `public`
# (the slim cohort and nothing else), `truth`, and a one-row `info_df` of
# realized values for the manifest. A failed fit gives `error_message` and
# no `public`.
.make_example <- function(design_row, split, loop_dir){
  seed_val <- .example_seed(cell_idx = design_row$cell_idx, split = split)
  bool_covariates <- design_row$covariates == "sexage"
  population <- .simulate_population(effects = design_row$effects,
                                     expression = design_row$expression,
                                     state = design_row$state,
                                     donor_sd = design_row$donor_sd,
                                     bool_covariates = bool_covariates,
                                     seed_number = seed_val)
  sim <- .simulate_cohort(population = population,
                          cells = design_row$cells,
                          num_donors = design_row$donors,
                          bool_covariates = bool_covariates,
                          bool_arm_sd_unequal = design_row$arm_sd == "unequal",
                          bool_depth_error = design_row$depth == "error",
                          seed_number = seed_val + 10)
  time_start <- proc.time()[["elapsed"]]
  res <- .run_esvd_pipeline(dat = sim$dat,
                            covariate_df = sim$covariate_df,
                            individual_vec = sim$individual_vec,
                            k = .fit_rank())
  fit_seconds <- proc.time()[["elapsed"]] - time_start
  if(!is.null(res$error_message)){
    return(list(error_message = paste0(res$stage, ": ", res$error_message)))
  }
  slim <- .slim_cohort(esvd_obj = res$esvd_obj, cc_vec = sim$cc_vec)

  # The rebuilt view must be the fit's own posterior.
  fit <- res$esvd_obj[[res$esvd_obj[["latest_Fit"]]]]
  view <- .build_view(slim)
  view_diff <- max(abs(view$posterior_mean_mat - fit$posterior_mean_mat),
                   abs(view$posterior_var_mat - fit$posterior_var_mat))

  deseq2_df <- .run_deseq2_pseudobulk(dat = sim$dat,
                                      donor_df = sim$donor_df,
                                      individual_vec = sim$individual_vec,
                                      bool_covariates = bool_covariates)
  wls_df <- .run_donor_wls(view = view, loop_dir = loop_dir)
  status_vec <- as.character(slim$nuisance_status_vec)
  info_df <- data.frame(fit_seconds = fit_seconds,
                        frac_boundary = mean(status_vec == "boundary"),
                        frac_capped = mean(status_vec == "capped"),
                        max_cells = max(sim$donor_df$num_cells),
                        mean_count = mean(sim$dat),
                        min_cells = min(sim$donor_df$num_cells),
                        num_cells = nrow(sim$dat),
                        num_exact_null = sum(sim$truth_log2fc_vec == 0),
                        num_warnings = length(res$warning_vec),
                        num_zero_genes = sim$num_zero_genes,
                        view_diff = view_diff)

  list(error_message = NULL,
       info_df = info_df,
       public = slim,
       truth = list(deseq2_df = deseq2_df,
                    donor_df = sim$donor_df,
                    donor_effect_mat = sim$donor_effect_mat,
                    donor_state_mat = sim$donor_state_mat,
                    gene_class_vec = .gene_class(sim$truth_log2fc_vec),
                    population = population,
                    truth_log2fc_vec = sim$truth_log2fc_vec,
                    warning_vec = res$warning_vec,
                    wls_df = wls_df))
}
