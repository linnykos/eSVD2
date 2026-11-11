# Shared machinery of the lfcse loop: bank location, design, generator, the
# eSVD-DE fit, and the public view. Drafted by Claude for Kevin Z. Lin,
# 2026-11-11.
#
# Sourced by prepare_claude.R (the bank builder) and evaluate_claude.R (the
# scorer). FIXED during a run: a change to anything in this file invalidates
# the bank (rebuild, re-baseline, rewrite notes/BANK.md) and needs a new
# `.bank_version()`. The design is notes/MEMO_claude.md, section 3.
# Every logFC and SE is on the log2 scale; every generator effect is on the
# natural-log scale.

# Location and version ---------------------------------------------------------

.bank_version <- function(){
  "lfcse-bank-1"
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

# The two levels of the donor SD (log scale): the lower and upper quartiles
# across genes of the SD estimated on the fitted Velmeshev L2/3 cohort
# (0.209 and 0.518; diagnostics/donor_sd_asd_claude.R), rounded.
.donor_sd_levels <- function(){
  c(0.2, 0.5)
}

# The six factors of the memo.
.bank_design <- function(donor_sd_vec = .donor_sd_levels()){
  design_df <- expand.grid(donors = c(10, 20, 40),
                           cells = c("equal", "unequal"),
                           rates = c("narrow", "wide"),
                           effects = c("sparse", "sparse_up", "dense"),
                           covariates = c("none", "sexage"),
                           donor_sd = donor_sd_vec,
                           KEEP.OUT.ATTRS = FALSE,
                           stringsAsFactors = FALSE)
  design_df$cell_idx <- seq_len(nrow(design_df))
  design_df
}

# One seed per example from its cell and split. Within an example the
# population uses `seed`, the effects `seed + 1`, cohort A `seed + 10` and
# `seed + 11`, cohort B `seed + 20` and `seed + 21`; cells are 100 apart so
# no two examples share a stream.
.example_seed <- function(cell_idx, split){
  split_idx <- match(split, c("dev", "test"))
  stopifnot(!is.na(split_idx))
  1000000 * split_idx + 100 * cell_idx
}

.example_id <- function(cell_idx, split){
  paste0(split, "_", formatC(cell_idx, width = 3, flag = "0"))
}

# The generator ----------------------------------------------------------------

# Gene-level truth shared by the two cohorts of an example: intercepts,
# loadings, Gamma rates (a large rate is a nearly Poisson gene), the
# covariate effects, and the two correlated case-control effect vectors.
#
# The effects of A and B are `sd_g * (sqrt(r) * shared + sqrt(1 - r) * own)`
# with B's SD divided by `k_ratio`, so on the genes that carry an effect the
# population correlation is `r_true` and the size ratio S_A / S_B is
# `k_ratio`. (Under `sparse_up` the absolute values raise the uncentred
# correlation to about 0.85; the scorer uses each example's realized k and
# r, stored with its truth.)
#   sparse:    10% of genes with SD 0.6, the others exactly 0.
#   sparse_up: the same sizes, every effect positive (absolute values).
#   dense:     every gene nonzero; 85% with SD 0.1 and 15% with SD 0.3.
.simulate_population <- function(effects,
                                 rates,
                                 bool_covariates,
                                 covariate_effect_sd = 0.2,
                                 gene_intercept_range = c(0, 1.5),
                                 k_ratio = 2,
                                 k_true = 3,
                                 num_genes = 600,
                                 r_true = 0.8,
                                 seed_number = 10){
  stopifnot(effects %in% c("sparse", "sparse_up", "dense"),
            rates %in% c("narrow", "wide"))
  p <- num_genes
  gene_vec <- paste0("gene", seq_len(p))
  nuisance_range <- if(rates == "narrow") c(2, 8) else c(0.5, 200)

  set.seed(seed_number)
  gene_intercept_vec <- stats::runif(p,
                                     min = gene_intercept_range[1],
                                     max = gene_intercept_range[2])
  y_mat <- matrix(stats::rnorm(p * k_true, sd = 0.4), nrow = p, ncol = k_true)
  nuisance_true_vec <- exp(stats::runif(p,
                                        min = log(nuisance_range[1]),
                                        max = log(nuisance_range[2])))
  # Drawn whether or not the example has covariates, so that the stream
  # below does not depend on the level of `covariates`.
  sex_effect_vec <- stats::rnorm(p, sd = covariate_effect_sd)
  age_effect_vec <- stats::rnorm(p, sd = covariate_effect_sd)
  if(!bool_covariates){
    sex_effect_vec <- rep(0, p)
    age_effect_vec <- rep(0, p)
  }

  set.seed(seed_number + 1)
  shared_vec <- stats::rnorm(p)
  own_a_vec <- stats::rnorm(p)
  own_b_vec <- stats::rnorm(p)
  base_a_vec <- sqrt(r_true) * shared_vec + sqrt(1 - r_true) * own_a_vec
  base_b_vec <- sqrt(r_true) * shared_vec + sqrt(1 - r_true) * own_b_vec
  # Positions permuted by position, not by value (`sample(p)`).
  order_vec <- sample(p)
  if(effects == "dense"){
    effect_sd_vec <- rep(0.1, p)
    effect_sd_vec[order_vec[seq_len(round(0.15 * p))]] <- 0.3
  } else {
    effect_sd_vec <- rep(0, p)
    effect_sd_vec[order_vec[seq_len(round(0.10 * p))]] <- 0.6
  }
  effect_a_vec <- effect_sd_vec * base_a_vec
  effect_b_vec <- (effect_sd_vec / k_ratio) * base_b_vec
  if(effects == "sparse_up"){
    effect_a_vec <- abs(effect_a_vec)
    effect_b_vec <- abs(effect_b_vec)
  }

  name_it <- function(vec){stats::setNames(vec, gene_vec)}
  rownames(y_mat) <- gene_vec
  list(age_effect_vec = name_it(age_effect_vec),
       cc_effect_a_vec = name_it(effect_a_vec),
       cc_effect_b_vec = name_it(effect_b_vec),
       effect_sd_vec = name_it(effect_sd_vec),
       gene_intercept_vec = name_it(gene_intercept_vec),
       gene_vec = gene_vec,
       k_true = k_true,
       nuisance_true_vec = name_it(nuisance_true_vec),
       sex_effect_vec = name_it(sex_effect_vec),
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
# with a donor random effect of SD `donor_sd` on the log scale. With
# covariates, each donor has a sex (60% of cases and 40% of controls are 1)
# and a standardized age, each acting on every gene through the
# population's effect vectors.
#
# The truth `cc_effect_vec / log(2)` is the log2 ratio of mean expression
# at equal depth and equal covariates: the latent factors, the donor
# effects and the library sizes are independent of the arm.
.simulate_cohort <- function(population,
                             cc_effect_vec,
                             cells,
                             donor_sd,
                             num_donors,
                             bool_covariates,
                             library_range = c(0.5, 1.5),
                             seed_number = 10){
  stopifnot(num_donors %% 2 == 0, length(cc_effect_vec) ==
              length(population$gene_vec))
  p <- length(population$gene_vec)
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

  x_mat <- matrix(stats::rnorm(n * population$k_true, sd = 0.4),
                  nrow = n, ncol = population$k_true)
  donor_effect_mat <- matrix(stats::rnorm(num_donors * p, sd = donor_sd),
                             nrow = num_donors, ncol = p)
  nat_mat <- tcrossprod(x_mat, population$y_mat)
  nat_mat <- nat_mat + rep(population$gene_intercept_vec, each = n)
  nat_mat <- nat_mat + outer(cc_vec, as.numeric(cc_effect_vec))
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
                         stringsAsFactors = FALSE)
  list(cc_vec = cc_vec,
       covariate_df = covariate_df,
       dat = dat,
       donor_df = donor_df,
       individual_vec = individual_vec,
       num_zero_genes = length(zero_gene_idx),
       truth_log2fc_vec = cc_effect_vec / log(2))
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
                               k = 3,
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

# One example ------------------------------------------------------------------

# Gene class by the true log2 fold change: exactly null, weak (below 0.5 in
# absolute value, the genes the validity gate regresses on), strong.
.gene_class <- function(truth_log2fc_vec){
  class_vec <- rep("strong", length(truth_log2fc_vec))
  class_vec[abs(truth_log2fc_vec) < 0.5] <- "weak"
  class_vec[truth_log2fc_vec == 0] <- "null"
  stats::setNames(class_vec, names(truth_log2fc_vec))
}

# Generates, fits and summarizes one example (a pair of cohorts A and B).
# Returns `public` (the two slim cohorts and nothing else), `truth`, and a
# one-row `info_df` of realized values for the manifest. A failed fit gives
# `error_message` and no `public`.
.make_example <- function(design_row, split){
  seed_val <- .example_seed(cell_idx = design_row$cell_idx, split = split)
  bool_covariates <- design_row$covariates == "sexage"
  population <- .simulate_population(effects = design_row$effects,
                                     rates = design_row$rates,
                                     bool_covariates = bool_covariates,
                                     seed_number = seed_val)

  cohort_list <- list()
  for(cohort in c("a", "b")){
    cc_effect_vec <- population[[paste0("cc_effect_", cohort, "_vec")]]
    sim <- .simulate_cohort(population = population,
                            cc_effect_vec = cc_effect_vec,
                            cells = design_row$cells,
                            donor_sd = design_row$donor_sd,
                            num_donors = design_row$donors,
                            bool_covariates = bool_covariates,
                            seed_number = seed_val +
                              if(cohort == "a") 10 else 20)
    time_start <- proc.time()[["elapsed"]]
    res <- .run_esvd_pipeline(dat = sim$dat,
                              covariate_df = sim$covariate_df,
                              individual_vec = sim$individual_vec,
                              k = population$k_true)
    fit_seconds <- proc.time()[["elapsed"]] - time_start
    if(!is.null(res$error_message)){
      return(list(error_message = paste0("cohort ", cohort, ", ", res$stage,
                                         ": ", res$error_message)))
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
    status_vec <- as.character(slim$nuisance_status_vec)
    cohort_list[[cohort]] <- list(
      slim = slim,
      truth = list(deseq2_df = deseq2_df,
                   donor_df = sim$donor_df,
                   gene_class_vec = .gene_class(sim$truth_log2fc_vec),
                   truth_log2fc_vec = sim$truth_log2fc_vec,
                   warning_vec = res$warning_vec),
      info = list(fit_seconds = fit_seconds,
                  frac_boundary = mean(status_vec == "boundary"),
                  frac_capped = mean(status_vec == "capped"),
                  max_cells = max(sim$donor_df$num_cells),
                  min_cells = min(sim$donor_df$num_cells),
                  num_cells = nrow(sim$dat),
                  num_warnings = length(res$warning_vec),
                  num_zero_genes = sim$num_zero_genes,
                  sex_case = mean(sim$donor_df$sex[sim$donor_df$cc == 1]),
                  sex_control = mean(sim$donor_df$sex[sim$donor_df$cc == 0]),
                  view_diff = view_diff)
    )
  }

  truth_a_vec <- cohort_list$a$truth$truth_log2fc_vec
  truth_b_vec <- cohort_list$b$truth$truth_log2fc_vec
  signal_rms_a <- sqrt(mean(truth_a_vec^2))
  signal_rms_b <- sqrt(mean(truth_b_vec^2))
  pair <- list(k = signal_rms_a / signal_rms_b,
               r = mean(truth_a_vec * truth_b_vec) /
                 (signal_rms_a * signal_rms_b))

  info_a <- cohort_list$a$info
  info_b <- cohort_list$b$info
  names(info_a) <- paste0(names(info_a), "_a")
  names(info_b) <- paste0(names(info_b), "_b")
  info_df <- as.data.frame(c(info_a, info_b, list(truth_k = pair$k,
                                                  truth_r = pair$r)))

  list(error_message = NULL,
       info_df = info_df,
       public = list(a = cohort_list$a$slim, b = cohort_list$b$slim),
       truth = list(a = cohort_list$a$truth,
                    b = cohort_list$b$truth,
                    pair = pair,
                    population = population))
}
