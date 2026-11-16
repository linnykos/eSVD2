# Shared machinery of the second beyond-score study: the bank
# `lfcse-bank-4` (location, design, generator).
# Drafted by Claude for Kevin Z. Lin, 2026-11-16, from
# ../2026-11-12_beyond-score/generator_claude.R.
#
# Sourced by 01_build-bank_claude.R and 02_apply-rules_claude.R after the
# second loop's generator, whose `.load_esvd2()`, `.run_esvd_pipeline()`,
# `.slim_cohort()` and `.build_view()` are used unchanged. Every logFC and
# SE is on the log2 scale; every generator effect is on the natural-log
# scale.
#
# What the fourth bank adds to the third: the four things the Gabitto
# split-half study listed under "What the autoresearch bank should add"
# (../2026-11-15_gabitto/2026-11-15_gabitto-split-half_claude.Rmd), on
# cohorts shaped like the Gabitto and Green halves.
#   rank       every cohort is fitted twice, at rank 5 and at rank 30 (the
#              default of `eSVD_helper()`, near or above 40 donors)
#   balance    half of the donors are controls, or one third of them
#   deg        weak effects on every gene, or on half of the genes with the
#              other half exactly null, so that each rule's calls can be
#              counted, with a capped share near the real 17%
#   shape      40 donors of many cells (a Gabitto half), 80 donors (the
#              Gabitto cohort), or 200 donors of few cells (a Green half),
#              at one total number of cells
# Each example is still two replicate cohorts of one population, which is
# what the gate on retained signal needs.
#
# What the Green split-half study added (../2026-11-16_green/), in a second
# block of examples:
#   tail       a few donors lie far below the others on the log scale for
#              a gene, which is what broke `donor_wls` on the Green halves
#   capped     a share of capped genes near the Green halves' 43%

# Location and version ---------------------------------------------------------

.bank_version <- function(){
  "lfcse-bank-4"
}

.bank_dir <- function(){
  root_dir <- Sys.getenv("ESVD2_BANK_DIR", unset = "")
  if(root_dir == "") root_dir <- tools::R_user_dir("eSVD2", which = "cache")
  file.path(root_dir, .bank_version())
}

.replicate_names <- function(){
  c("a", "b")
}

# A trial run (`ESVD2_BANK4_TRIAL` set, with `ESVD2_BANK_DIR` pointing at a
# scratch folder) builds two small examples to exercise the three steps.
.bool_trial <- function(){
  Sys.getenv("ESVD2_BANK4_TRIAL", unset = "") != ""
}

# Cells per cohort, in expectation.
.cells_total <- function(){
  if(.bool_trial()) 1500 else 10000
}

# The ranks every cohort is fitted at.
.fit_ranks <- function(){
  c(5, 30)
}

# Design -----------------------------------------------------------------------

.factor_names <- function(){
  c("donors", "balance", "deg", "state", "donor_sd", "draw", "tail", "capped")
}

# The first block: five factors crossed in full (3 * 2^4 = 48
# combinations), and two populations (`draw`) per combination at 40 donors,
# the shape the recommendation is weighted towards: 64 examples, with
# Gaussian donor effects (`tail = gaussian`) and the lower capped share
# (`capped = low`).
# The second block, from the Green study, at 40 and at 200 donors with a
# third of them controls: 16 examples with outlying donors (`tail =
# outlier`; both levels of `deg`, `state` and `donor_sd`), and 8 with the
# higher capped share (`capped = high`; `state = donor`, both levels of
# `deg` and `donor_sd`).
# 88 examples, 176 cohorts, 352 fits. There is one split: no rule is tuned
# on this bank. The first block's rows keep their `cell_idx`, and so their
# seeds, whatever is appended.
.bank_design <- function(){
  level_list <- list(donors = c(40, 80, 200),
                     balance = c("half", "third"),
                     deg = c("dense", "half_null"),
                     state = c("none", "donor"),
                     donor_sd = c(0.2, 0.4),
                     draw = c(1, 2))
  design_df <- expand.grid(level_list, KEEP.OUT.ATTRS = FALSE,
                           stringsAsFactors = FALSE)
  design_df <- design_df[design_df$draw == 1 | design_df$donors == 40, ]
  design_df$tail <- "gaussian"
  design_df$capped <- "low"

  green_list <- list(donors = c(40, 200),
                     balance = "third",
                     deg = c("dense", "half_null"),
                     state = c("none", "donor"),
                     donor_sd = c(0.2, 0.4),
                     draw = 1)
  outlier_df <- expand.grid(green_list, KEEP.OUT.ATTRS = FALSE,
                            stringsAsFactors = FALSE)
  outlier_df$tail <- "outlier"
  outlier_df$capped <- "low"
  capped_df <- outlier_df[outlier_df$state == "donor", ]
  capped_df$tail <- "gaussian"
  capped_df$capped <- "high"

  design_df$cell_idx <- seq_len(nrow(design_df))
  # The outlier examples were first built with a milder mechanism under
  # the indices 65 to 80 (2026-11-16, superseded and no longer in the
  # design); the recalibrated ones take 101 to 116, so that their ids and
  # seeds are new, and the `capped = high` examples keep 81 to 88.
  outlier_df$cell_idx <- 100 + seq_len(nrow(outlier_df))
  capped_df$cell_idx <- 80 + seq_len(nrow(capped_df))
  design_df <- rbind(design_df, outlier_df, capped_df)
  rownames(design_df) <- NULL
  # A trial run of the three steps on a scratch bank keeps three examples.
  if(.bool_trial()) design_df <- design_df[c(4, 70, 86), ]
  design_df
}

# One seed per example. The population uses `seed` and `seed + 1`,
# replicate a `seed + 10` and `seed + 11`, replicate b `seed + 20` and
# `seed + 21`; examples are 100 apart so no two share a stream.
.example_seed <- function(cell_idx){
  4000000 + 100 * cell_idx
}

.example_id <- function(cell_idx){
  paste0("ex_", formatC(cell_idx, width = 4, flag = "0"))
}

# The generator ----------------------------------------------------------------

# Gene-level truth of an example. As the third bank's population, except:
#
# Intercepts: uniform on c(-3.2, 1), the 5% to 95% range of the mean count
#   per cell of a Gabitto half (0.055 to 2.5).
# Case-control effects (natural log), sized to the Gabitto halves, whose
#   `donor_wls` logFC covary between halves as true effects of SD 0.11
#   would:
#   dense:     every gene nonzero; 90% with SD 0.08 and 10% with SD 0.25.
#   half_null: half of the genes exactly 0; 40% with SD 0.1 and 10% with
#              SD 0.3.
# Over-dispersion: the Gamma rate is log-uniform on `nuisance_range`,
#   chosen by a pilot so that the capped and boundary genes of a 40-donor
#   cohort fitted at rank 30 are 19% (the Gabitto halves have 17%; the
#   third bank's range of 0.5 to 200 would give over 40%). At rank 5 the
#   same cohort has fewer.
#   At `capped = high` the range is 0.5 to 50, for a capped share near the
#   Green halves' 43%.
# Covariates: four at the donor level, as in Tati's fits (sex, APOE4,
#   PMI, age), each with a per-gene effect of SD 0.2. None varies within a
#   donor, and there is no depth-error block.
.simulate_population <- function(deg,
                                 state,
                                 donor_sd,
                                 capped = "low",
                                 covariate_sd = 0.2,
                                 gene_intercept_range = c(-3.2, 1),
                                 k_true = 3,
                                 nuisance_range = c(0.15, 15),
                                 num_genes = 600,
                                 state_share = 0.75,
                                 seed_number = 10){
  stopifnot(deg %in% c("dense", "half_null"),
            state %in% c("none", "donor"),
            capped %in% c("low", "high"))
  if(capped == "high") nuisance_range <- c(0.5, 50)
  p <- num_genes
  gene_vec <- paste0("gene", seq_len(p))

  set.seed(seed_number)
  gene_intercept_vec <- stats::runif(p,
                                     min = gene_intercept_range[1],
                                     max = gene_intercept_range[2])
  y_mat <- matrix(stats::rnorm(p * k_true, sd = 0.4), nrow = p, ncol = k_true)
  y_mat <- scale(y_mat, center = TRUE, scale = FALSE)
  attr(y_mat, "scaled:center") <- NULL
  nuisance_true_vec <- exp(stats::runif(p, min = log(nuisance_range[1]),
                                        max = log(nuisance_range[2])))
  covariate_effect_mat <- matrix(stats::rnorm(p * 4, sd = covariate_sd),
                                 nrow = p, ncol = 4)
  colnames(covariate_effect_mat) <- c("sex", "apoe", "pmi", "age")
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
  num_strong <- round(0.10 * p)
  if(deg == "dense"){
    effect_sd_vec <- rep(0.08, p)
    effect_sd_vec[order_vec[seq_len(num_strong)]] <- 0.25
  } else {
    effect_sd_vec <- rep(0, p)
    effect_sd_vec[order_vec[seq_len(round(0.5 * p))]] <- 0.1
    effect_sd_vec[order_vec[seq_len(num_strong)]] <- 0.3
  }
  cc_effect_vec <- effect_sd_vec * base_vec

  name_it <- function(vec){stats::setNames(vec, gene_vec)}
  rownames(y_mat) <- gene_vec
  rownames(covariate_effect_mat) <- gene_vec
  list(cc_effect_vec = name_it(cc_effect_vec),
       covariate_effect_mat = covariate_effect_mat,
       gene_intercept_vec = name_it(gene_intercept_vec),
       gene_vec = gene_vec,
       idio_sd_vec = name_it(idio_sd_vec),
       k_true = k_true,
       nuisance_true_vec = name_it(nuisance_true_vec),
       state_sd = state_sd,
       y_mat = y_mat)
}

# Cells per donor: log-normal around `num_cells_total / num_donors` with
# the log-scale SD of the Gabitto donors (0.46), at least 20.
.simulate_cells_per_donor <- function(num_donors, num_cells_total,
                                      log_sd = 0.46){
  mean_val <- num_cells_total / num_donors
  num_vec <- round(mean_val * exp(stats::rnorm(num_donors, sd = log_sd) -
                                    log_sd^2 / 2))
  as.integer(pmax(num_vec, 20))
}

# One cohort. Donor states and idiosyncratic donor effects are the third
# bank's (equal donor variance in the arms, no depth error).
#
# Arms: half of the donors are controls (`balance = half`) or a third of
#   them (`third`; 13 of 40, as in a Gabitto half).
# Covariates, all at the donor level: sex (half of each arm), APOE4 (12% of
#   controls and 37% of cases, the Gabitto shares), PMI and age
#   (standardized, independent of the arm).
#
# Outlying donors (`tail = outlier`), calibrated to the Green cohort. Of
#   the log donor means of its genes with under 80% zeros, 2.3% lie more
#   than 3 MADs below the gene's median (0.8% in Gabitto; 0.1% under a
#   Gaussian), by 2.1 on the natural-log scale on average (median 1.7),
#   and they sit in a tenth of the donors, each with 12 to 25% of its
#   genes affected. Here a tenth of the donors are outlying, 15% of an
#   outlying donor's genes (0.5% of any other donor's) are affected, and
#   an affected donor-gene pair has its log expression lowered by 1 plus
#   an exponential of mean 1 (mean 2, median 1.7). Winsorizing each gene's
#   log donor means restored the agreement of `donor_wls` between the
#   Green halves (../2026-11-16_green/05_diagnose-donor-wls_claude.R): an
#   outlier on the low side moves a mean of logarithms far more than an
#   arithmetic mean. The outliers are drawn alike in both arms, so no
#   truth changes.
#
# The truth is the case-control effect at equal covariates,
# cc_effect_vec / log(2). No arm differs in cell state, so the marginal
# and the direct truths coincide.
.simulate_cohort <- function(population,
                             num_donors,
                             balance,
                             tail = "gaussian",
                             library_range = c(0.5, 1.5),
                             num_cells_total = .cells_total(),
                             seed_number = 10){
  stopifnot(balance %in% c("half", "third"),
            tail %in% c("gaussian", "outlier"))
  p <- length(population$gene_vec)
  k_true <- population$k_true
  num_control <- if(balance == "half") num_donors / 2 else round(num_donors / 3)
  num_case <- num_donors - num_control

  set.seed(seed_number)
  num_cells_vec <- .simulate_cells_per_donor(num_donors = num_donors,
                                             num_cells_total = num_cells_total)
  n <- sum(num_cells_vec)
  donor_vec <- paste0("indiv", seq_len(num_donors))
  cc_donor_vec <- rep(c(0, 1), times = c(num_control, num_case))
  .draw_binary <- function(fraction_vec){
    value_vec <- rep(0, num_donors)
    control_pos <- sample(num_control)[seq_len(round(fraction_vec[1] *
                                                       num_control))]
    case_pos <- num_control +
      sample(num_case)[seq_len(round(fraction_vec[2] * num_case))]
    value_vec[c(control_pos, case_pos)] <- 1
    value_vec
  }
  sex_donor_vec <- .draw_binary(c(0.5, 0.5))
  apoe_donor_vec <- .draw_binary(c(0.12, 0.37))
  pmi_donor_vec <- as.numeric(scale(stats::rnorm(num_donors)))
  age_donor_vec <- as.numeric(scale(stats::rnorm(num_donors)))
  donor_covariate_mat <- cbind(sex = sex_donor_vec, apoe = apoe_donor_vec,
                               pmi = pmi_donor_vec, age = age_donor_vec)

  donor_idx_vec <- rep(seq_len(num_donors), times = num_cells_vec)
  cc_vec <- cc_donor_vec[donor_idx_vec]
  individual_vec <- factor(donor_vec[donor_idx_vec], levels = donor_vec)
  cell_vec <- paste0("cell", seq_len(n))
  covariate_df <- data.frame(CC = factor(cc_vec, levels = c(0, 1)),
                             Sex = factor(sex_donor_vec[donor_idx_vec],
                                          levels = c(0, 1)),
                             Apoe = factor(apoe_donor_vec[donor_idx_vec],
                                           levels = c(0, 1)),
                             Pmi = pmi_donor_vec[donor_idx_vec],
                             Age = age_donor_vec[donor_idx_vec])
  rownames(covariate_df) <- cell_vec

  donor_state_mat <- matrix(stats::rnorm(num_donors * k_true),
                            nrow = num_donors, ncol = k_true) *
    population$state_sd
  idio_mat <- matrix(stats::rnorm(num_donors * p), nrow = num_donors,
                     ncol = p)
  idio_sd_mat <- matrix(as.numeric(population$idio_sd_vec), nrow = num_donors,
                        ncol = p, byrow = TRUE)
  # The donor's whole effect on each gene, each part less its log-normal
  # mean correction, so that its exponential has mean 1 in each arm.
  y_square_vec <- rowSums(population$y_mat^2)
  donor_effect_mat <- tcrossprod(donor_state_mat, population$y_mat) -
    matrix(population$state_sd^2 * y_square_vec / 2, nrow = num_donors,
           ncol = p, byrow = TRUE) +
    idio_mat * idio_sd_mat - idio_sd_mat^2 / 2
  donor_effect_mat <- donor_effect_mat +
    tcrossprod(donor_covariate_mat, population$covariate_effect_mat)
  donor_effect_mat <- donor_effect_mat +
    outer(cc_donor_vec, as.numeric(population$cc_effect_vec))
  dimnames(donor_effect_mat) <- list(donor_vec, population$gene_vec)

  cell_state_mat <- matrix(stats::rnorm(n * k_true, sd = 0.4), nrow = n,
                           ncol = k_true)
  if(tail == "outlier"){
    # Drawn after everything the Gaussian cohorts draw at this seed. The
    # probability is recycled down the rows, one value per donor.
    outlier_prob_vec <- ifelse(stats::runif(num_donors) < 0.10, 0.15, 0.005)
    bool_outlier_mat <- matrix(stats::runif(num_donors * p) < outlier_prob_vec,
                               nrow = num_donors, ncol = p)
    num_outliers <- sum(bool_outlier_mat)
    donor_effect_mat[bool_outlier_mat] <- donor_effect_mat[bool_outlier_mat] -
      (1 + stats::rexp(num_outliers, rate = 1))
  }
  nat_mat <- tcrossprod(cell_state_mat, population$y_mat)
  nat_mat <- nat_mat + rep(population$gene_intercept_vec, each = n)
  nat_mat <- nat_mat + donor_effect_mat[donor_idx_vec, , drop = FALSE]
  library_size_vec <- stats::runif(n, min = library_range[1],
                                   max = library_range[2])

  set.seed(seed_number + 1)
  rate_vec <- rep(population$nuisance_true_vec, each = n)
  lambda_mat <- matrix(stats::rgamma(n * p,
                                     shape = as.numeric(exp(nat_mat)) * rate_vec,
                                     rate = rate_vec),
                       nrow = n, ncol = p)
  rm(nat_mat, rate_vec)
  dat <- matrix(stats::rpois(n * p, lambda = library_size_vec * lambda_mat),
                nrow = n, ncol = p)
  rm(lambda_mat)
  dimnames(dat) <- list(cell_vec, population$gene_vec)
  # An all-zero gene is refused by 1.2.0; one count keeps it in the cohort.
  zero_gene_idx <- which(colSums(dat) == 0)
  if(length(zero_gene_idx) > 0) dat[1, zero_gene_idx] <- 1L
  storage.mode(dat) <- "integer"

  donor_df <- data.frame(donor = donor_vec,
                         cc = cc_donor_vec,
                         num_cells = num_cells_vec,
                         donor_covariate_mat,
                         stringsAsFactors = FALSE)
  list(cc_vec = cc_vec,
       covariate_df = covariate_df,
       dat = dat,
       donor_df = donor_df,
       individual_vec = individual_vec,
       num_zero_genes = length(zero_gene_idx))
}

# One example ------------------------------------------------------------------

# Generates one population and its two replicate cohorts, and fits each
# cohort end to end at every rank of `.fit_ranks()`. Returns `public_list`
# (one slim cohort per replicate and rank, named `<replicate>_k<rank>`,
# nothing of the truth), `truth`, and a one-row `info_df` for the manifest.
.make_example <- function(design_row){
  seed_val <- .example_seed(cell_idx = design_row$cell_idx)
  population <- .simulate_population(deg = design_row$deg,
                                     state = design_row$state,
                                     donor_sd = design_row$donor_sd,
                                     capped = design_row$capped,
                                     seed_number = seed_val)
  replicate_vec <- .replicate_names()
  public_list <- list()
  cohort_list <- list()
  info_list <- list()
  for(replicate_idx in seq_along(replicate_vec)){
    replicate <- replicate_vec[replicate_idx]
    sim <- .simulate_cohort(population = population,
                            num_donors = design_row$donors,
                            balance = design_row$balance,
                            tail = design_row$tail,
                            seed_number = seed_val + 10 * replicate_idx)
    cohort_list[[replicate]] <- list(donor_df = sim$donor_df)
    for(rank in .fit_ranks()){
      label <- paste0(replicate, "_k", rank)
      time_start <- proc.time()[["elapsed"]]
      res <- .run_esvd_pipeline(dat = sim$dat,
                                covariate_df = sim$covariate_df,
                                individual_vec = sim$individual_vec,
                                k = rank)
      fit_seconds <- proc.time()[["elapsed"]] - time_start
      if(!is.null(res$error_message)){
        return(list(error_message = paste0(label, ", ", res$stage, ": ",
                                           res$error_message)))
      }
      slim <- .slim_cohort(esvd_obj = res$esvd_obj, cc_vec = sim$cc_vec)
      # The package's own multiple-testing call, kept for the comparison of
      # its calls with each rule's.
      slim$stored_df$fdr <- as.numeric(
        res$esvd_obj$pvalue_list$fdr_vec[slim$stored_df$gene]
      )
      status_vec <- as.character(slim$nuisance_status_vec)
      public_list[[label]] <- slim
      info_list[[label]] <- c(fit_seconds = fit_seconds,
                              frac_boundary = mean(status_vec == "boundary"),
                              frac_capped = mean(status_vec == "capped"),
                              frac_zero = mean(sim$dat == 0),
                              num_cells = nrow(sim$dat),
                              num_warnings = length(res$warning_vec),
                              num_zero_genes = sim$num_zero_genes)
      rm(res, slim)
    }
  }
  info_df <- as.data.frame(as.list(unlist(lapply(names(info_list),
                                                 function(label){
    stats::setNames(info_list[[label]],
                    paste0(names(info_list[[label]]), "_", label))
  }))))
  info_df$num_exact_null <- sum(population$cc_effect_vec == 0)

  truth_vec <- population$cc_effect_vec / log(2)
  list(error_message = NULL,
       info_df = info_df,
       public_list = public_list,
       truth = list(cohort_list = cohort_list,
                    population = population,
                    truth_direct_vec = truth_vec,
                    truth_marginal_vec = truth_vec))
}
