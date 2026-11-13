# Shared machinery of the beyond-score study: the bank `lfcse-bank-3`
# (location, design, generator) and what both later steps need.
# Drafted by Claude for Kevin Z. Lin, 2026-11-12, from
# autoresearch_lfcse2/harness/generator_claude.R.
#
# Sourced by 01_build-bank_claude.R and 02_apply-rules_claude.R after the
# second loop's generator, whose `.load_esvd2()`, `.run_esvd_pipeline()`,
# `.slim_cohort()`, `.build_view()` and `.gene_class()` are used unchanged.
# The functions below replace that file's bank version, design, population
# and cohort. Every logFC and SE is on the log2 scale; every generator
# effect is on the natural-log scale.
#
# What the third bank adds to the second (Kevin's note of 2026-11-12):
#   confound   covariates that move expression more than the arm does, at
#              the donor level (sex, age) or the cell level (Pct, Phase),
#              and that are unbalanced between the arms
#   sparsity   every cohort holds genes from about 10% to 99% zeros
#   deg        few (5%) or many (30%, mostly up) true DEGs among exact
#              nulls, or every gene weakly nonzero
#   replicates every example is two cohorts drawn from one population

# Location and version ---------------------------------------------------------

.bank_version <- function(){
  "lfcse-bank-3"
}

.bank_dir <- function(){
  root_dir <- Sys.getenv("ESVD2_BANK_DIR", unset = "")
  if(root_dir == "") root_dir <- tools::R_user_dir("eSVD2", which = "cache")
  file.path(root_dir, .bank_version())
}

.replicate_names <- function(){
  c("a", "b")
}

# Design -----------------------------------------------------------------------

.factor_names <- function(){
  c("donors", "confound", "deg", "state", "cells", "donor_sd", "arm_sd",
    "depth")
}

# Eight factors; the full factorial has 3^4 * 2^4 = 1296 cells. The bank
# holds a 1/6 fraction of 216: a third of the 81 combinations of the four
# 3-level factors (the sum of their level indices is 0 modulo 3) crossed
# with half of the 16 combinations of the four 2-level factors (the sum of
# their level indices is even). Every level of every factor, and every pair
# of levels of two factors, appears equally often. There is one split: no
# rule is tuned on this bank.
.bank_design <- function(){
  level_list <- list(donors = c(10, 20, 40),
                     confound = c("mild", "donor", "cell"),
                     deg = c("few", "many", "dense"),
                     state = c("none", "donor", "donor_arm"),
                     cells = c("equal", "unequal"),
                     donor_sd = c(0.2, 0.5),
                     arm_sd = c("equal", "unequal"),
                     depth = c("none", "error"))
  index_df <- expand.grid(lapply(level_list, seq_along),
                          KEEP.OUT.ATTRS = FALSE)
  index_df$cell_idx <- seq_len(nrow(index_df))
  three_sum <- rowSums(index_df[, c("donors", "confound", "deg", "state")] - 1)
  two_sum <- rowSums(index_df[, c("cells", "donor_sd", "arm_sd", "depth")] - 1)
  index_df <- index_df[three_sum %% 3 == 0 & two_sum %% 2 == 0, ]
  design_df <- data.frame(cell_idx = index_df$cell_idx,
                          stringsAsFactors = FALSE)
  for(name in names(level_list)){
    design_df[[name]] <- level_list[[name]][index_df[[name]]]
  }
  design_df <- design_df[, c(.factor_names(), "cell_idx")]
  rownames(design_df) <- NULL
  design_df
}

# One seed per example. The population uses `seed` and `seed + 1`,
# replicate a `seed + 10` and `seed + 11`, replicate b `seed + 20` and
# `seed + 21`; examples are 100 apart so no two share a stream.
.example_seed <- function(cell_idx){
  3000000 + 100 * cell_idx
}

.example_id <- function(cell_idx){
  paste0("ex_", formatC(cell_idx, width = 4, flag = "0"))
}

# The generator ----------------------------------------------------------------

# Gene-level truth of an example. As the second bank's population, except:
#
# Intercepts: uniform on c(-4, 1.5) in every cohort, so that the mean count
#   per cell runs from about 0.02 to 4.5 (and e times that in the block)
#   and a cohort holds genes that are nearly all zeros beside genes that
#   are rarely zero.
# Case-control effects (natural log):
#   few:   5% of genes with SD 0.6, both signs, the others exactly 0.
#   many:  30% of genes with SD 0.6, three in four positive, the others
#          exactly 0. The positive excess moves the total count of a case
#          cell, which every library estimate built on totals inherits.
#   dense: every gene nonzero; 85% with SD 0.1 and 15% with SD 0.3.
# Covariate effects, per gene:
#   sex, age (donor level): SD 0.2 at `confound = mild` and `cell`, 0.5 at
#          `donor`, against the 0.6 of a DEG that 5% or 30% of genes have.
#   Pct (cell level, continuous, SD about 1 within a donor): SD 0.3.
#   Phase (cell level, binary): SD 0.5.
#   The two cell-level effects are zero unless `confound = cell`.
.simulate_population <- function(deg,
                                 state,
                                 confound,
                                 donor_sd,
                                 block_fraction = 0.15,
                                 gene_intercept_range = c(-4, 1.5),
                                 k_true = 3,
                                 num_genes = 600,
                                 state_share = 0.75,
                                 state_shift_sd = 0.25,
                                 seed_number = 10){
  stopifnot(deg %in% c("few", "many", "dense"),
            state %in% c("none", "donor", "donor_arm"),
            confound %in% c("mild", "donor", "cell"))
  p <- num_genes
  gene_vec <- paste0("gene", seq_len(p))

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
  # Standard normal draws, scaled below, so that the stream does not depend
  # on the level of `confound`.
  sex_base_vec <- stats::rnorm(p)
  age_base_vec <- stats::rnorm(p)
  pct_base_vec <- stats::rnorm(p)
  phase_base_vec <- stats::rnorm(p)
  donor_covariate_sd <- if(confound == "donor") 0.5 else 0.2
  bool_cell <- confound == "cell"
  state_shift_vec <- stats::rnorm(k_true, sd = state_shift_sd)
  if(state != "donor_arm") state_shift_vec <- rep(0, k_true)
  donor_sd_vec <- donor_sd * exp(stats::rnorm(p, sd = 0.3))
  if(state == "none"){
    state_sd <- 0
    idio_sd_vec <- donor_sd_vec
  } else {
    state_sd <- sqrt(state_share * donor_sd^2 / mean(rowSums(y_mat^2)))
    idio_sd_vec <- sqrt(1 - state_share) * donor_sd_vec
  }

  set.seed(seed_number + 1)
  base_vec <- stats::rnorm(p)
  # Positions permuted by position, not by value (`sample(p)`).
  order_vec <- sample(p)
  sign_vec <- ifelse(stats::runif(p) < 0.75, 1, -1)
  if(deg == "dense"){
    effect_sd_vec <- rep(0.1, p)
    effect_sd_vec[order_vec[seq_len(round(0.15 * p))]] <- 0.3
    cc_effect_vec <- effect_sd_vec * base_vec
  } else {
    deg_fraction <- if(deg == "few") 0.05 else 0.30
    effect_sd_vec <- rep(0, p)
    effect_sd_vec[order_vec[seq_len(round(deg_fraction * p))]] <- 0.6
    cc_effect_vec <- effect_sd_vec * base_vec
    if(deg == "many") cc_effect_vec <- sign_vec * abs(cc_effect_vec)
  }

  name_it <- function(vec){stats::setNames(vec, gene_vec)}
  rownames(y_mat) <- gene_vec
  list(age_effect_vec = name_it(donor_covariate_sd * age_base_vec),
       bool_block_vec = name_it(bool_block_vec),
       cc_effect_vec = name_it(cc_effect_vec),
       gene_intercept_vec = name_it(gene_intercept_vec),
       gene_vec = gene_vec,
       idio_sd_vec = name_it(idio_sd_vec),
       k_true = k_true,
       nuisance_true_vec = name_it(nuisance_true_vec),
       pct_effect_vec = name_it(bool_cell * 0.3 * pct_base_vec),
       phase_effect_vec = name_it(bool_cell * 0.5 * phase_base_vec),
       sex_effect_vec = name_it(donor_covariate_sd * sex_base_vec),
       state_sd = state_sd,
       state_shift_vec = state_shift_vec,
       y_mat = y_mat)
}

.simulate_cells_per_donor <- function(cells, num_donors){
  if(cells == "equal") return(rep(50L, num_donors))
  num_vec <- round(exp(stats::rnorm(num_donors, mean = log(150), sd = 0.7)))
  as.integer(pmin(pmax(num_vec, 30), 600))
}

# One cohort. As the second bank's cohort (donor states, idiosyncratic
# donor effects, unequal arm variance, the depth-error block), except for
# the covariates:
#
# Sex and age are in every cohort. At `confound = mild` and `cell` they are
#   nearly balanced (sex 40% of controls and 60% of cases; age independent
#   of the arm). At `confound = donor` sex is 20% of controls and 80% of
#   cases and the cases are older by one SD.
# At `confound = cell` each cell also has
#   Pct   = (donor shift, SD 0.5) + 0.8 for a case + N(0, 1), and
#   Phase = 1 with probability plogis(qlogis(0.2 or 0.5 by arm) + donor
#           shift of SD 0.5),
#   so the arms differ in both and a donor-level summary cannot remove
#   their effect without the per-gene coefficients.
#
# The truth is the case-control effect at equal covariates and equal depth:
# the marginal one, (cc_effect_vec + y_mat %*% state_shift_vec) / log(2),
# includes the arms' difference in cell state; the direct one is
# cc_effect_vec / log(2). They differ only at `state = donor_arm`.
.simulate_cohort <- function(population,
                             cells,
                             num_donors,
                             confound,
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
  sex_fraction_vec <- if(confound == "donor") c(0.2, 0.8) else c(0.4, 0.6)
  sex_donor_vec <- rep(0, num_donors)
  control_pos <- sample(num_arm)[seq_len(round(sex_fraction_vec[1] * num_arm))]
  case_pos <- num_arm +
    sample(num_arm)[seq_len(round(sex_fraction_vec[2] * num_arm))]
  sex_donor_vec[c(control_pos, case_pos)] <- 1
  age_shift <- if(confound == "donor") 1 else 0
  age_donor_vec <- as.numeric(scale(stats::rnorm(num_donors) +
                                      age_shift * cc_donor_vec))

  donor_idx_vec <- rep(seq_len(num_donors), times = num_cells_vec)
  cc_vec <- cc_donor_vec[donor_idx_vec]
  individual_vec <- factor(donor_vec[donor_idx_vec], levels = donor_vec)
  cell_vec <- paste0("cell", seq_len(n))

  # Drawn at every level of `confound` so that the stream below does not
  # depend on it.
  pct_donor_vec <- stats::rnorm(num_donors, sd = 0.5)
  phase_donor_vec <- stats::rnorm(num_donors, sd = 0.5)
  pct_vec <- pct_donor_vec[donor_idx_vec] + 0.8 * cc_vec + stats::rnorm(n)
  phase_prob_vec <- stats::plogis(stats::qlogis(ifelse(cc_vec == 1, 0.5, 0.2)) +
                                    phase_donor_vec[donor_idx_vec])
  phase_vec <- as.numeric(stats::runif(n) < phase_prob_vec)

  covariate_df <- data.frame(CC = factor(cc_vec, levels = c(0, 1)),
                             Sex = factor(sex_donor_vec[donor_idx_vec],
                                          levels = c(0, 1)),
                             Age = age_donor_vec[donor_idx_vec])
  if(confound == "cell"){
    covariate_df$Pct <- pct_vec
    covariate_df$Phase <- factor(phase_vec, levels = c(0, 1))
  }
  rownames(covariate_df) <- cell_vec

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
  # mean correction, so that its exponential has mean 1 in each arm.
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
  # The cell-level covariates' part of the log mean; zero unless
  # `confound = cell`.
  cell_covariate_mat <- outer(pct_vec, as.numeric(population$pct_effect_vec)) +
    outer(phase_vec, as.numeric(population$phase_effect_vec))
  nat_mat <- tcrossprod(cell_state_mat, population$y_mat)
  nat_mat <- nat_mat + rep(population$gene_intercept_vec, each = n)
  nat_mat <- nat_mat + outer(cc_vec, marginal_effect_vec)
  nat_mat <- nat_mat + donor_effect_mat[donor_idx_vec, , drop = FALSE]
  nat_mat <- nat_mat + outer(sex_donor_vec[donor_idx_vec],
                             as.numeric(population$sex_effect_vec))
  nat_mat <- nat_mat + outer(age_donor_vec[donor_idx_vec],
                             as.numeric(population$age_effect_vec))
  nat_mat <- nat_mat + cell_covariate_mat
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
       library_true_vec = stats::setNames(library_size_vec, cell_vec),
       num_zero_genes = length(zero_gene_idx),
       phase_vec = phase_vec,
       pct_vec = pct_vec,
       truth_log2fc_vec = stats::setNames(marginal_effect_vec / log(2),
                                          population$gene_vec))
}

# One example ------------------------------------------------------------------

# Generates one population and its two replicate cohorts, and fits each
# cohort end to end. Returns `public_list` (one slim cohort per replicate,
# nothing of the truth), `truth`, and a one-row `info_df` for the manifest.
.make_example <- function(design_row){
  seed_val <- .example_seed(cell_idx = design_row$cell_idx)
  population <- .simulate_population(deg = design_row$deg,
                                     state = design_row$state,
                                     confound = design_row$confound,
                                     donor_sd = design_row$donor_sd,
                                     seed_number = seed_val)
  replicate_vec <- .replicate_names()
  public_list <- list()
  cohort_list <- list()
  info_list <- list()
  for(replicate_idx in seq_along(replicate_vec)){
    replicate <- replicate_vec[replicate_idx]
    sim <- .simulate_cohort(population = population,
                            cells = design_row$cells,
                            num_donors = design_row$donors,
                            confound = design_row$confound,
                            bool_arm_sd_unequal = design_row$arm_sd == "unequal",
                            bool_depth_error = design_row$depth == "error",
                            seed_number = seed_val + 10 * replicate_idx)
    time_start <- proc.time()[["elapsed"]]
    res <- .run_esvd_pipeline(dat = sim$dat,
                              covariate_df = sim$covariate_df,
                              individual_vec = sim$individual_vec,
                              k = .fit_rank())
    fit_seconds <- proc.time()[["elapsed"]] - time_start
    if(!is.null(res$error_message)){
      return(list(error_message = paste0("replicate ", replicate, ", ",
                                         res$stage, ": ", res$error_message)))
    }
    slim <- .slim_cohort(esvd_obj = res$esvd_obj, cc_vec = sim$cc_vec)
    # The package's own multiple-testing call, kept for the comparison of
    # significance with the reported interval.
    slim$stored_df$fdr <- as.numeric(
      res$esvd_obj$pvalue_list$fdr_vec[slim$stored_df$gene]
    )
    fit <- res$esvd_obj[[res$esvd_obj[["latest_Fit"]]]]
    view <- .build_view(slim)
    view_diff <- max(abs(view$posterior_mean_mat - fit$posterior_mean_mat),
                     abs(view$posterior_var_mat - fit$posterior_var_mat))
    rm(view, fit)
    status_vec <- as.character(slim$nuisance_status_vec)
    public_list[[replicate]] <- slim
    cohort_list[[replicate]] <- list(donor_df = sim$donor_df,
                                     donor_effect_mat = sim$donor_effect_mat,
                                     donor_state_mat = sim$donor_state_mat,
                                     library_true_vec = sim$library_true_vec,
                                     phase_vec = sim$phase_vec,
                                     pct_vec = sim$pct_vec,
                                     warning_vec = res$warning_vec)
    info_list[[replicate]] <- c(fit_seconds = fit_seconds,
                                frac_boundary = mean(status_vec == "boundary"),
                                frac_capped = mean(status_vec == "capped"),
                                frac_zero = mean(sim$dat == 0),
                                mean_count = mean(sim$dat),
                                num_cells = nrow(sim$dat),
                                num_warnings = length(res$warning_vec),
                                num_zero_genes = sim$num_zero_genes,
                                view_diff = view_diff)
  }
  info_mat <- do.call(rbind, info_list)
  info_df <- as.data.frame(as.list(c(
    stats::setNames(info_mat["a", ], paste0(colnames(info_mat), "_a")),
    stats::setNames(info_mat["b", ], paste0(colnames(info_mat), "_b"))
  )))
  marginal_vec <- (as.numeric(population$cc_effect_vec) +
                     as.numeric(population$y_mat %*%
                                  population$state_shift_vec)) / log(2)
  info_df$num_exact_null <- sum(population$cc_effect_vec == 0)

  list(error_message = NULL,
       info_df = info_df,
       public_list = public_list,
       truth = list(cohort_list = cohort_list,
                    population = population,
                    truth_direct_vec = population$cc_effect_vec / log(2),
                    truth_marginal_vec = stats::setNames(marginal_vec,
                                                         population$gene_vec)))
}
