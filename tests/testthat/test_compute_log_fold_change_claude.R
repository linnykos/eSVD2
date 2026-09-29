# UNIT_TEST_PLAN.md section 2.18 -- T-LFC-01 .. T-LFC-19.
#
# The log2 fold change and its standard error:
#
#   log2fc    = log2(case_mean / control_mean)
#   log2fc_se = (1 / ln 2) * sqrt(case_var / (n1 * case_mean^2) +
#                                 control_var / (n0 * control_mean^2))
#
# where `case_mean` / `control_mean` are the averages over INDIVIDUALS of their
# mean posterior expression, `case_var` / `control_var` are the mixture
# variances the Welch statistic already uses (population variance, no Bessel
# correction; Kevin, 2026-08-29), and n1 / n0 count individuals, never cells.
#
# Three sections:
#   A. the formula, on hand-built posterior matrices (exact oracles);
#   B. the plumbing, on fitted toy data (both pipelines, the wrappers);
#   C. whether the number behaves as a standard error (bootstrap, truth,
#      repeated sampling).
#
# Every fixture is built here in code from a fixed seed. The fitted ones are
# cached in `.fixture_cache` (helper-fixtures.R) so each is fitted once per run.

# ---- fixtures ---------------------------------------------------------------

# Hand-built posterior matrices. Deliberately NOT aligned or balanced: 3 case
# and 5 control individuals, a different number of cells per individual, the
# cells shuffled so no individual is a contiguous block of rows, and the
# individuals passed to the function in an order that is neither their
# appearance order nor their factor-level order. Positional indexing that is
# right by accident on a tidy fixture is wrong here.
.lfc_matrix_fixture <- function(seed_number = 10){
  set.seed(seed_number)

  case_individuals <- c("indiv_7", "indiv_2", "indiv_5")
  control_individuals <- c("indiv_4", "indiv_8", "indiv_1", "indiv_6",
                           "indiv_3")
  num_cells_vec <- c(indiv_1 = 6, indiv_2 = 3, indiv_3 = 4, indiv_4 = 7,
                     indiv_5 = 5, indiv_6 = 3, indiv_7 = 4, indiv_8 = 3)
  individual_vec <- rep(names(num_cells_vec), times = num_cells_vec)
  individual_vec <- individual_vec[sample(length(individual_vec))]

  n <- length(individual_vec)
  p <- 6
  # Gene-specific scales spanning two orders of magnitude, so a formula that
  # forgets to divide by the mean is off by a different factor in every gene.
  gene_scale_vec <- c(0.2, 1, 5, 20, 0.5, 2)
  posterior_mean_mat <- matrix(stats::rgamma(n * p, shape = 4, rate = 1),
                               nrow = n, ncol = p)
  posterior_mean_mat <- posterior_mean_mat * rep(gene_scale_vec, each = n)
  posterior_var_mat <- matrix(stats::runif(n * p, min = 0.05, max = 0.5),
                              nrow = n, ncol = p)
  posterior_var_mat <- posterior_var_mat * rep(gene_scale_vec^2, each = n)
  dimnames(posterior_mean_mat) <- list(paste0("cell_", seq_len(n)),
                                       paste0("gene", seq_len(p)))
  dimnames(posterior_var_mat) <- dimnames(posterior_mean_mat)

  list(case_individuals = case_individuals,
       control_individuals = control_individuals,
       individual_vec = factor(individual_vec),
       posterior_mean_mat = posterior_mean_mat,
       posterior_var_mat = posterior_var_mat)
}

.lfc_new_fields <- function(){
  c("case_var", "control_var", "log2fc_se_vec", "log2fc_vec")
}

# Every section-A test goes through here, and the length check is what keeps
# them from passing vacuously: `expect_equal(NULL, NULL)` succeeds, so before
# the fields existed a comparison of two absent vectors was green.
.lfc_run_default <- function(fixture){
  res <- compute_test_statistic(
    input_obj = fixture$posterior_mean_mat,
    posterior_var_mat = fixture$posterior_var_mat,
    case_individuals = fixture$case_individuals,
    control_individuals = fixture$control_individuals,
    individual_vec = fixture$individual_vec,
    verbose = 0
  )

  for(element_name in .lfc_new_fields()){
    expect_length(res[[element_name]], ncol(fixture$posterior_mean_mat))
  }
  res
}

# The oracle for section A: the same quantities by explicit loops over genes,
# arms and individuals. It shares no code with the package (no averaging
# matrix, no `colMeans`), and it writes the between-individual variance in its
# centred form, mean((m - mean(m))^2), where the package uses
# mean(m^2) - mean(m)^2.
.lfc_oracle <- function(posterior_mean_mat,
                        posterior_var_mat,
                        case_individuals,
                        control_individuals,
                        individual_vec){
  individual_vec <- as.character(individual_vec)
  arm_list <- list(case = as.character(case_individuals),
                   control = as.character(control_individuals))
  p <- ncol(posterior_mean_mat)

  res <- list(case_mean = numeric(p), case_var = numeric(p),
              control_mean = numeric(p), control_var = numeric(p),
              log2fc_se_vec = numeric(p), log2fc_vec = numeric(p),
              within_se_vec = numeric(p))
  for(j in seq_len(p)){
    arm_mean <- c(case = NA, control = NA)
    arm_within <- c(case = NA, control = NA)
    arm_between <- c(case = NA, control = NA)
    for(arm in names(arm_list)){
      indiv_mean_vec <- numeric(0)
      indiv_var_vec <- numeric(0)
      for(indiv in arm_list[[arm]]){
        cell_idx <- which(individual_vec == indiv)
        indiv_mean_vec <- c(indiv_mean_vec,
                            mean(posterior_mean_mat[cell_idx, j]))
        indiv_var_vec <- c(indiv_var_vec,
                           mean(posterior_var_mat[cell_idx, j]))
      }
      arm_mean[arm] <- mean(indiv_mean_vec)
      arm_within[arm] <- mean(indiv_var_vec)
      arm_between[arm] <- mean((indiv_mean_vec - mean(indiv_mean_vec))^2)
    }
    num_case <- length(arm_list$case)
    num_control <- length(arm_list$control)
    arm_var <- arm_within + arm_between

    res$case_mean[j] <- arm_mean["case"]
    res$control_mean[j] <- arm_mean["control"]
    res$case_var[j] <- arm_var["case"]
    res$control_var[j] <- arm_var["control"]
    res$log2fc_vec[j] <- log(arm_mean["case"] / arm_mean["control"]) / log(2)
    res$log2fc_se_vec[j] <- sqrt(
      arm_var["case"] / (num_case * arm_mean["case"]^2) +
        arm_var["control"] / (num_control * arm_mean["control"]^2)
    ) / log(2)
    # The part of the SE that comes from the per-cell posterior variance
    # alone; section C subtracts it to isolate the between-individual part.
    res$within_se_vec[j] <- sqrt(
      arm_within["case"] / (num_case * arm_mean["case"]^2) +
        arm_within["control"] / (num_control * arm_mean["control"]^2)
    ) / log(2)
  }

  res
}

# Identical cells within each individual and zero posterior variance, so the
# within-individual average is exact and the mixture variance reduces to the
# population variance of the individual means. That reduction is what lets
# `stats::var` and `stats::t.test` act as external oracles. Three cells per
# individual because `compute_test_statistic` refuses fewer.
.lfc_donor_level_fixture <- function(num_case, num_control,
                                     cells_per_donor = 3, num_genes = 5,
                                     seed_number = 10){
  set.seed(seed_number)

  num_individuals <- num_case + num_control
  donor_mean_mat <- matrix(stats::runif(num_individuals * num_genes,
                                        min = 2, max = 9),
                           nrow = num_individuals, ncol = num_genes)
  individual_names <- paste0("indiv_", seq_len(num_individuals))
  rownames(donor_mean_mat) <- individual_names
  colnames(donor_mean_mat) <- paste0("gene", seq_len(num_genes))

  row_idx <- rep(seq_len(num_individuals), each = cells_per_donor)
  posterior_mean_mat <- donor_mean_mat[row_idx, , drop = FALSE]
  rownames(posterior_mean_mat) <- paste0("cell_", seq_along(row_idx))
  posterior_var_mat <- posterior_mean_mat * 0

  list(case_individuals = individual_names[seq_len(num_case)],
       control_individuals = individual_names[num_case + seq_len(num_control)],
       donor_mean_mat = donor_mean_mat,
       individual_vec = factor(individual_names[row_idx],
                               levels = individual_names),
       posterior_mean_mat = posterior_mean_mat,
       posterior_var_mat = posterior_var_mat)
}

# The whole pipeline in the order `eSVD()` runs it, on the matrix path (the
# posterior matrices are kept, which section C needs): two fits, the
# case-control coefficient held fixed and then freed. `.small_esvd_obj()` in
# helper-fixtures.R stops after the first fit; a fold change should be judged
# on the fit a user actually gets.
.lfc_fit_pipeline <- function(dat,
                              covariates,
                              individual_vec,
                              case_control_variable,
                              max_iter = 10){
  esvd_obj <- suppressWarnings(
    initialize_esvd(dat = dat,
                    covariates = covariates,
                    metadata_individual = individual_vec,
                    bool_intercept = TRUE,
                    case_control_variable = case_control_variable,
                    k = 2,
                    lambda = 0.1,
                    metadata_case_control = covariates[, case_control_variable],
                    verbose = 0)
  )
  esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                                 fit_name = "fit_Init",
                                                 omitted_variables = "Log_UMI")

  esvd_obj <- suppressWarnings(
    opt_esvd(input_obj = esvd_obj,
             l2pen = 0.1,
             max_iter = max_iter,
             offset_variables = setdiff(colnames(esvd_obj$covariates),
                                        case_control_variable),
             tol = 1e-6,
             fit_name = "fit_First",
             fit_previous = "fit_Init",
             verbose = 0)
  )
  esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                                 fit_name = "fit_First",
                                                 omitted_variables = "Log_UMI")

  esvd_obj <- suppressWarnings(
    opt_esvd(input_obj = esvd_obj,
             l2pen = 0.1,
             max_iter = max_iter,
             offset_variables = NULL,
             tol = 1e-6,
             fit_name = "fit_Second",
             fit_previous = "fit_First",
             verbose = 0)
  )
  esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
                                                 fit_name = "fit_Second",
                                                 omitted_variables = NULL)

  esvd_obj <- suppressWarnings(
    estimate_nuisance(input_obj = esvd_obj,
                      bool_covariates_as_library = TRUE,
                      verbose = 0)
  )
  esvd_obj <- compute_posterior(input_obj = esvd_obj,
                                alpha_max = 2 * max(dat),
                                bool_covariates_as_library = TRUE,
                                library_min = 0.1)
  compute_test_statistic(input_obj = esvd_obj, verbose = 0)
}

# F-SMALL (400 cells, 40 genes, 4 case / 4 control individuals), fitted. Genes
# 1-5 carry a true case-control effect of +0.8 on the natural-log scale. The
# raw inputs are `.small_data()`.
.lfc_fitted_small <- function(){
  if(!is.null(.fixture_cache$lfc_fitted_small)){
    return(.fixture_cache$lfc_fitted_small)
  }

  dat_list <- .small_data()
  covariates <- .tiny_covariates(dat = dat_list$dat,
                                 covariate_df = dat_list$covariate_df)
  esvd_obj <- .lfc_fit_pipeline(dat = dat_list$dat,
                                covariates = covariates,
                                individual_vec = dat_list$individual_vec,
                                case_control_variable = "CC_1")

  .fixture_cache$lfc_fitted_small <- esvd_obj
  esvd_obj
}

# One cohort from `generate_null()`, fitted. Unlike F-SMALL this generator
# draws a random effect per individual, so the individuals genuinely differ
# and the between-individual variance is not just cell-sampling noise.
.lfc_fitted_null <- function(seed_number,
                             cell_per_person = 20,
                             num_genes = 40,
                             num_individuals = 8){
  cache_key <- paste0("lfc_fitted_null_", seed_number, "_", cell_per_person,
                      "_", num_genes, "_", num_individuals)
  if(!is.null(.fixture_cache[[cache_key]])) return(.fixture_cache[[cache_key]])

  set.seed(seed_number)
  null_list <- generate_null(cell_per_person = cell_per_person,
                             num_genes = num_genes,
                             num_individuals = num_individuals)
  esvd_obj <- .lfc_fit_pipeline(
    dat = as.matrix(null_list$obs_mat),
    covariates = null_list$covariates,
    individual_vec = factor(null_list$metadata_individual),
    case_control_variable = "CC"
  )

  res <- list(esvd_obj = esvd_obj, nat_mat = null_list$nat_mat)
  .fixture_cache[[cache_key]] <- res
  res
}

# The log2 fold change the generator planted, for every gene: the log2 ratio
# of the arms' mean expression, exp(nat_mat) being the mean of the Gamma the
# counts are drawn through. It is read off the generator's `nat_mat` and never
# touches a fitted object.
.lfc_planted_log2fc <- function(nat_mat, cc_vec){
  case_idx <- which(cc_vec == 1)
  control_idx <- which(cc_vec == 0)
  mean_mat <- exp(nat_mat)

  log(apply(mean_mat[case_idx, , drop = FALSE], 2, mean) /
        apply(mean_mat[control_idx, , drop = FALSE], 2, mean)) / log(2)
}

# Which genes have a nuisance estimate that has diverged. `nuisance_vec` is
# the Gamma RATE, and both generators draw it from a known range: 2 to 8 in
# F-SMALL, 0.1 to 10 in `generate_null()`. An estimate above 1e4 is a
# thousand times the largest true value. In practice the estimates split
# cleanly, the bulk below 1e3 and the diverged ones near 1e7, so the exact
# threshold does not matter.
.lfc_nuisance_diverged <- function(esvd_obj, threshold = 1e4){
  latest_fit <- esvd_obj[["latest_Fit"]]
  esvd_obj[[latest_fit]]$nuisance_vec > threshold
}

# The larger of the two arms' between-individual coefficients of variation of
# the individual means, per gene: the quantity the accuracy of the delta
# method depends on.
.lfc_between_cv <- function(avg_list){
  cv_fun <- function(mat){
    apply(mat, 2, function(x){sqrt(mean((x - mean(x))^2)) / mean(x)})
  }
  pmax(cv_fun(avg_list$donor_mean_mat[avg_list$case_individuals, ,
                                      drop = FALSE]),
       cv_fun(avg_list$donor_mean_mat[avg_list$control_individuals, ,
                                      drop = FALSE]))
}

# The per-individual averages of the two posterior matrices, by `rowsum()`
# rather than the package's averaging matrix.
.lfc_donor_averages <- function(esvd_obj){
  latest_fit <- esvd_obj[["latest_Fit"]]
  individual_vec <- as.character(esvd_obj[["individual"]])
  num_cells_vec <- table(individual_vec)

  mean_sum_mat <- rowsum(esvd_obj[[latest_fit]]$posterior_mean_mat,
                         group = individual_vec)
  var_sum_mat <- rowsum(esvd_obj[[latest_fit]]$posterior_var_mat,
                        group = individual_vec)
  num_cells_vec <- as.numeric(num_cells_vec[rownames(mean_sum_mat)])

  cc_by_individual <- tapply(esvd_obj[["case_control"]], individual_vec,
                             function(x){unique(x)})
  stopifnot(all(sapply(cc_by_individual, length) == 1))

  list(case_individuals = names(cc_by_individual)[cc_by_individual == 1],
       control_individuals = names(cc_by_individual)[cc_by_individual == 0],
       donor_mean_mat = mean_sum_mat / num_cells_vec,
       donor_var_mat = var_sum_mat / num_cells_vec)
}

# The same decomposition as `.lfc_oracle`, from per-individual averages.
.lfc_from_donor_averages <- function(avg_list){
  population_var <- function(mat){
    apply(mat, 2, function(x){mean((x - mean(x))^2)})
  }
  case_mean_mat <- avg_list$donor_mean_mat[avg_list$case_individuals, ,
                                           drop = FALSE]
  control_mean_mat <- avg_list$donor_mean_mat[avg_list$control_individuals, ,
                                              drop = FALSE]
  case_var_mat <- avg_list$donor_var_mat[avg_list$case_individuals, ,
                                         drop = FALSE]
  control_var_mat <- avg_list$donor_var_mat[avg_list$control_individuals, ,
                                            drop = FALSE]
  num_case <- nrow(case_mean_mat)
  num_control <- nrow(control_mean_mat)

  case_mean <- apply(case_mean_mat, 2, mean)
  control_mean <- apply(control_mean_mat, 2, mean)
  case_within <- apply(case_var_mat, 2, mean)
  control_within <- apply(control_var_mat, 2, mean)
  case_between <- population_var(case_mean_mat)
  control_between <- population_var(control_mean_mat)

  list(between_se_vec = sqrt(case_between / (num_case * case_mean^2) +
                               control_between /
                               (num_control * control_mean^2)) / log(2),
       case_mean = case_mean,
       case_var = case_within + case_between,
       control_mean = control_mean,
       control_var = control_within + control_between,
       log2fc_se_vec = sqrt((case_within + case_between) /
                              (num_case * case_mean^2) +
                              (control_within + control_between) /
                              (num_control * control_mean^2)) / log(2),
       log2fc_vec = log(case_mean / control_mean) / log(2),
       within_se_vec = sqrt(case_within / (num_case * case_mean^2) +
                              control_within /
                              (num_control * control_mean^2)) / log(2))
}

# Resamples INDIVIDUALS with replacement within each arm, holding the fit and
# every cell's posterior fixed, and returns the bootstrap SD of the log2 fold
# change for every gene.
.lfc_donor_bootstrap <- function(avg_list, num_bootstrap = 2000,
                                 seed_number = 10){
  case_mat <- avg_list$donor_mean_mat[avg_list$case_individuals, ,
                                      drop = FALSE]
  control_mat <- avg_list$donor_mean_mat[avg_list$control_individuals, ,
                                         drop = FALSE]
  num_case <- nrow(case_mat)
  num_control <- nrow(control_mat)

  set.seed(seed_number)
  boot_mat <- matrix(NA_real_, nrow = num_bootstrap, ncol = ncol(case_mat))
  for(b in seq_len(num_bootstrap)){
    # `sample.int()` rather than `sample()`: the latter on a length-1 vector
    # permutes 1:x.
    case_idx <- sample.int(num_case, size = num_case, replace = TRUE)
    control_idx <- sample.int(num_control, size = num_control, replace = TRUE)
    boot_mat[b, ] <- log(
      apply(case_mat[case_idx, , drop = FALSE], 2, mean) /
        apply(control_mat[control_idx, , drop = FALSE], 2, mean)
    ) / log(2)
  }

  apply(boot_mat, 2, stats::sd)
}

# ---- A. the formula ---------------------------------------------------------

test_that("T-LFC-01: the four new vectors equal an explicit per-individual loop", {
  fixture <- .lfc_matrix_fixture()
  res <- .lfc_run_default(fixture)
  oracle <- .lfc_oracle(
    posterior_mean_mat = fixture$posterior_mean_mat,
    posterior_var_mat = fixture$posterior_var_mat,
    case_individuals = fixture$case_individuals,
    control_individuals = fixture$control_individuals,
    individual_vec = fixture$individual_vec
  )

  expect_true(all(.lfc_new_fields() %in% names(res)))
  for(element_name in c("case_mean", "control_mean", .lfc_new_fields())){
    expect_equal(unname(res[[element_name]]), oracle[[element_name]],
                 tolerance = 1e-10, info = element_name)
  }

  # The natural-log SE, which is the scale NEBULA reports on, is the log2 SE
  # times ln 2. Written out so the scale is pinned by a test and not only by
  # the documentation.
  natural_log_se_vec <- sqrt(
    oracle$case_var / (3 * oracle$case_mean^2) +
      oracle$control_var / (5 * oracle$control_mean^2)
  )
  expect_equal(unname(res$log2fc_se_vec) * log(2), natural_log_se_vec,
               tolerance = 1e-10)
})

## The delta method, checked without using its closed form. The SE of g(a, b)
## for independent a and b is sqrt(grad' Sigma grad); here the gradient of
## g = log2(a / b) comes from numerical differentiation, so an algebra slip in
## the package's closed form (a missing 1 / ln 2, a mean that is not squared,
## the two arms' means swapped in the denominators) cannot also be in the
## oracle.
test_that("T-LFC-02: the SE is the delta method with a numerical gradient", {
  skip_if_not_installed("numDeriv")

  fixture <- .lfc_matrix_fixture()
  res <- .lfc_run_default(fixture)
  num_case <- length(fixture$case_individuals)
  num_control <- length(fixture$control_individuals)

  for(gene_idx in seq_len(ncol(fixture$posterior_mean_mat))){
    label <- paste0("gene ", gene_idx)
    gradient_vec <- numDeriv::grad(
      func = function(x){log2(x[1] / x[2])},
      x = c(res$case_mean[gene_idx], res$control_mean[gene_idx])
    )
    sigma_mat <- diag(c(res$case_var[gene_idx] / num_case,
                        res$control_var[gene_idx] / num_control))
    delta_se <- sqrt(as.numeric(
      crossprod(gradient_vec, sigma_mat %*% gradient_vec)
    ))

    expect_equal(as.numeric(res$log2fc_se_vec[gene_idx]), delta_se,
                 tolerance = 1e-7, info = label)
  }
})

## The degenerate case that makes an external oracle. With identical cells and
## zero posterior variance the mixture variance is the POPULATION variance of
## the individual means, so
##
##   case_var = (n1 - 1) / n1 * stats::var(case individual means)
##
## and the linear-scale SE is Welch's standard error with each arm's sample
## variance shrunk by (n - 1) / n. The factor is deliberate (no Bessel
## correction; Kevin, 2026-08-29) and is the same one T-TSTAT-01a pins for the
## test statistic.
test_that("T-LFC-03: on donor-level data the SE reduces to the textbook form", {
  grid <- expand.grid(num_case = c(3, 5, 12), num_control = c(3, 7))

  for(i in seq_len(nrow(grid))){
    num_case <- grid$num_case[i]
    num_control <- grid$num_control[i]
    label <- paste0("num_case = ", num_case, ", num_control = ", num_control)

    fixture <- .lfc_donor_level_fixture(num_case = num_case,
                                        num_control = num_control,
                                        seed_number = 10 + i)
    res <- .lfc_run_default(fixture)

    case_mat <- fixture$donor_mean_mat[fixture$case_individuals, ,
                                       drop = FALSE]
    control_mat <- fixture$donor_mean_mat[fixture$control_individuals, ,
                                          drop = FALSE]
    case_sample_var <- apply(case_mat, 2, stats::var)
    control_sample_var <- apply(control_mat, 2, stats::var)
    case_factor <- (num_case - 1) / num_case
    control_factor <- (num_control - 1) / num_control

    expect_equal(unname(res$case_var), unname(case_factor * case_sample_var),
                 tolerance = 1e-10, info = label)
    expect_equal(unname(res$control_var),
                 unname(control_factor * control_sample_var),
                 tolerance = 1e-10, info = label)

    textbook_se_vec <- sqrt(
      case_factor * case_sample_var / (num_case * apply(case_mat, 2, mean)^2) +
        control_factor * control_sample_var /
        (num_control * apply(control_mat, 2, mean)^2)
    ) / log(2)
    expect_equal(unname(res$log2fc_se_vec), unname(textbook_se_vec),
                 tolerance = 1e-10, info = label)

    # With equal arms the factor is common to both, and `stats::t.test`
    # becomes the oracle for the linear-scale SE directly.
    if(num_case == num_control){
      linear_se_vec <- sqrt(res$case_var / num_case +
                              res$control_var / num_control)
      for(gene_idx in seq_len(ncol(case_mat))){
        t_res <- stats::t.test(case_mat[, gene_idx], control_mat[, gene_idx],
                               var.equal = FALSE)
        expect_equal(as.numeric(linear_se_vec[gene_idx]),
                     as.numeric(t_res$stderr) * sqrt(case_factor),
                     tolerance = 1e-10,
                     info = paste0(label, ", gene ", gene_idx))
      }
    }
  }
})

## Ties the two stored variances to the statistic that is already tested
## (T-TSTAT-01a and onward): they must be the variances the Welch statistic
## divided by, not a second computation that could drift from it.
test_that("T-LFC-04: the stored variances reproduce the test statistic", {
  fixture <- .lfc_matrix_fixture()
  res <- .lfc_run_default(fixture)
  num_case <- length(fixture$case_individuals)
  num_control <- length(fixture$control_individuals)

  linear_se_vec <- sqrt(res$case_var / num_case +
                          res$control_var / num_control)
  expect_equal(unname((res$case_mean - res$control_mean) / linear_se_vec),
               unname(res$teststat_vec), tolerance = 1e-10)
})

test_that("T-LFC-05: swapping the arms negates the fold change and keeps the SE", {
  fixture <- .lfc_matrix_fixture()
  res <- .lfc_run_default(fixture)

  swapped <- fixture
  swapped$case_individuals <- fixture$control_individuals
  swapped$control_individuals <- fixture$case_individuals
  res_swapped <- .lfc_run_default(swapped)

  expect_equal(unname(res_swapped$log2fc_vec), -unname(res$log2fc_vec),
               tolerance = 1e-10)
  expect_equal(unname(res_swapped$log2fc_se_vec), unname(res$log2fc_se_vec),
               tolerance = 1e-10)
  expect_equal(unname(res_swapped$case_var), unname(res$control_var),
               tolerance = 1e-10)
  expect_equal(unname(res_swapped$control_var), unname(res$case_var),
               tolerance = 1e-10)
})

## A log ratio has no units. Rescaling expression by c multiplies every
## posterior mean by c and every posterior variance by c^2, and must change
## neither the fold change nor its SE. (This is a statement about the two
## matrices handed to the function, not about rescaling the COUNTS and
## refitting, which the model is not invariant to.)
test_that("T-LFC-06: rescaling expression leaves the fold change and SE alone", {
  fixture <- .lfc_matrix_fixture()
  res <- .lfc_run_default(fixture)

  for(scale_val in c(0.01, 0.1, 10, 100)){
    label <- paste0("scale = ", scale_val)
    rescaled <- fixture
    rescaled$posterior_mean_mat <- fixture$posterior_mean_mat * scale_val
    rescaled$posterior_var_mat <- fixture$posterior_var_mat * scale_val^2
    res_rescaled <- .lfc_run_default(rescaled)

    expect_equal(res_rescaled$log2fc_vec, res$log2fc_vec,
                 tolerance = 1e-10, info = label)
    expect_equal(res_rescaled$log2fc_se_vec, res$log2fc_se_vec,
                 tolerance = 1e-10, info = label)
    # Not vacuous: the variances themselves did move.
    expect_equal(res_rescaled$case_var, res$case_var * scale_val^2,
                 tolerance = 1e-10, info = label)
  }
})

## The property that makes this SE comparable to a pseudobulk one: the unit of
## replication is the individual. Sequencing k times as many cells from the
## same individuals must not shrink the SE at all, and recruiting k times as
## many individuals must shrink it by exactly sqrt(k). A cell-level SE (the
## pseudoreplication the method exists to avoid) fails the first half.
test_that("T-LFC-07: the SE scales with individuals, not with cells", {
  fixture <- .lfc_matrix_fixture()
  res <- .lfc_run_default(fixture)

  for(num_copies in c(2, 5, 20)){
    label <- paste0("num_copies = ", num_copies)
    n <- nrow(fixture$posterior_mean_mat)
    row_idx <- rep(seq_len(n), times = num_copies)

    # k times the cells, same individuals.
    more_cells <- fixture
    more_cells$posterior_mean_mat <- fixture$posterior_mean_mat[row_idx, ,
                                                                drop = FALSE]
    more_cells$posterior_var_mat <- fixture$posterior_var_mat[row_idx, ,
                                                              drop = FALSE]
    rownames(more_cells$posterior_mean_mat) <- paste0("cell_",
                                                      seq_along(row_idx))
    rownames(more_cells$posterior_var_mat) <- paste0("cell_",
                                                     seq_along(row_idx))
    more_cells$individual_vec <- factor(
      as.character(fixture$individual_vec)[row_idx]
    )
    res_cells <- .lfc_run_default(more_cells)

    expect_equal(res_cells$log2fc_vec, res$log2fc_vec,
                 tolerance = 1e-10, info = label)
    expect_equal(res_cells$log2fc_se_vec, res$log2fc_se_vec,
                 tolerance = 1e-10, info = label)

    # k times the individuals, each a copy of an original one.
    copy_vec <- rep(seq_len(num_copies), each = n)
    more_donors <- more_cells
    more_donors$individual_vec <- factor(
      paste0(as.character(fixture$individual_vec)[row_idx], "_copy", copy_vec)
    )
    more_donors$case_individuals <- as.vector(
      outer(fixture$case_individuals, seq_len(num_copies),
            FUN = function(x, y){paste0(x, "_copy", y)})
    )
    more_donors$control_individuals <- as.vector(
      outer(fixture$control_individuals, seq_len(num_copies),
            FUN = function(x, y){paste0(x, "_copy", y)})
    )
    res_donors <- .lfc_run_default(more_donors)

    expect_equal(res_donors$log2fc_vec, res$log2fc_vec,
                 tolerance = 1e-10, info = label)
    expect_equal(res_donors$log2fc_se_vec,
                 res$log2fc_se_vec / sqrt(num_copies),
                 tolerance = 1e-10, info = label)
  }
})

## Decision D4. A posterior mean from `compute_posterior` is always positive,
## so this is only reachable through the matrix method with a user's own
## matrices. `log2()` of a non-positive ratio is NaN or -Inf without a word;
## the contract is NA for both numbers, one warning that says how many genes,
## and every other gene untouched.
test_that("T-LFC-08: a non-positive arm mean gives NA with a warning", {
  fixture <- .lfc_matrix_fixture()
  res <- .lfc_run_default(fixture)

  bad_gene_idx <- c(2, 5)
  control_cell_idx <- which(as.character(fixture$individual_vec) %in%
                              fixture$control_individuals)
  broken <- fixture
  # Gene 2: the control mean is negative. Gene 5: it is exactly zero.
  broken$posterior_mean_mat[control_cell_idx, 2] <- -1
  broken$posterior_mean_mat[control_cell_idx, 5] <- 0

  expect_warning(res_broken <- .lfc_run_default(broken),
                 regexp = "2 gene\\(s\\).*not positive")

  expect_true(all(is.na(res_broken$log2fc_vec[bad_gene_idx])))
  expect_true(all(is.na(res_broken$log2fc_se_vec[bad_gene_idx])))
  expect_equal(res_broken$log2fc_vec[-bad_gene_idx],
               res$log2fc_vec[-bad_gene_idx], tolerance = 1e-10)
  expect_equal(res_broken$log2fc_se_vec[-bad_gene_idx],
               res$log2fc_se_vec[-bad_gene_idx], tolerance = 1e-10)
  expect_true(all(is.finite(res_broken$log2fc_se_vec[-bad_gene_idx])))

  # And the ordinary fixture is silent.
  expect_silent(.lfc_run_default(fixture))
})

## The same contract at the helper, for the inputs the matrix method cannot
## easily be made to produce. Found by code review of the first draft, which
## (a) returned an SE of NaN, not NA, for a gene with a NaN variance, because
## it blanked the inputs instead of the outputs and NaN / NA is NaN; and
## (b) treated a gene with ONE NA input as padded, so a negative mean beside
## an NA variance skipped the guard and went through `log2()`.
test_that("T-LFC-08b: every invalid input gives NA in both outputs, never NaN or Inf", {
  gene_vec <- paste0("gene", 1:8)
  case_mean <- stats::setNames(c(2, 2, -1, 2, NA, 2, 0, Inf), gene_vec)
  control_mean <- stats::setNames(c(1, 1, 1, 1, NA, 1, 1, 1), gene_vec)
  case_var <- stats::setNames(c(0.5, NaN, NA, -0.1, NA, NA, 0.5, 0.5),
                              gene_vec)
  control_var <- stats::setNames(c(0.5, 0.5, 0.5, 0.5, NA, 0.5, 0.5, 0.5),
                                 gene_vec)
  # gene1 valid; gene2 NaN variance; gene3 negative mean beside an NA
  # variance; gene4 negative variance; gene5 padded (all four NA); gene6 one
  # NA among valid inputs; gene7 zero mean; gene8 infinite mean.

  expect_warning(
    res <- .compute_log2_fold_change(case_mean = case_mean,
                                     control_mean = control_mean,
                                     case_var = case_var,
                                     control_var = control_var,
                                     num_case = 4,
                                     num_control = 4),
    regexp = "^6 gene\\(s\\)"
  )

  for(element_name in c("log2fc_vec", "log2fc_se_vec")){
    vec <- res[[element_name]]
    expect_identical(names(vec), gene_vec, info = element_name)
    expect_true(is.finite(vec[1]), info = element_name)
    expect_true(all(is.na(vec[-1])), info = element_name)
    # NA_real_ and not NaN: `is.na()` is TRUE for both.
    expect_false(any(is.nan(vec)), info = element_name)
  }
  expect_equal(unname(res$log2fc_vec[1]), 1)
  expect_equal(unname(res$log2fc_se_vec[1]),
               sqrt(0.5 / (4 * 2^2) + 0.5 / (4 * 1^2)) / log(2))

  # The padded gene alone raises nothing.
  expect_silent(
    .compute_log2_fold_change(case_mean = case_mean[c(1, 5)],
                              control_mean = control_mean[c(1, 5)],
                              case_var = case_var[c(1, 5)],
                              control_var = control_var[c(1, 5)],
                              num_case = 4,
                              num_control = 4)
  )
})

test_that("T-LFC-09: every new vector carries the gene names", {
  fixture <- .lfc_matrix_fixture()
  res <- .lfc_run_default(fixture)
  gene_vec <- colnames(fixture$posterior_mean_mat)

  for(element_name in .lfc_new_fields()){
    expect_identical(names(res[[element_name]]), gene_vec,
                     info = element_name)
    expect_true(is.numeric(res[[element_name]]), info = element_name)
    expect_null(dim(res[[element_name]]), info = element_name)
  }

  esvd_obj <- .small_esvd_obj()
  for(element_name in .lfc_new_fields()){
    expect_identical(names(esvd_obj[[element_name]]), colnames(esvd_obj$dat),
                     info = element_name)
  }
})

# ---- B. the plumbing --------------------------------------------------------

test_that("T-LFC-10: the per-gene pipeline agrees with the matrix pipeline", {
  esvd_obj <- .small_esvd_obj()

  # Strip the matrix pipeline's results so the per-gene path starts clean,
  # including the recorded individuals, which the per-gene path must record
  # for itself.
  per_gene_input <- esvd_obj
  latest_fit <- per_gene_input[["latest_Fit"]]
  per_gene_input[[latest_fit]]$posterior_mean_mat <- NULL
  per_gene_input[[latest_fit]]$posterior_var_mat <- NULL
  for(element_name in c("teststat_vec", "case_mean", "control_mean",
                        "pvalue_list", .lfc_new_fields())){
    per_gene_input[[element_name]] <- NULL
  }
  per_gene_input$param$test_case_individuals <- NULL
  per_gene_input$param$test_control_individuals <- NULL

  res_per_gene <- .muffle_locfdr_fallback(
    compute_test_per_gene(input_obj = per_gene_input,
                          alpha_max = 2 * max(esvd_obj$dat),
                          bool_covariates_as_library = TRUE,
                          library_min = 0.1,
                          verbose = 0)
  )

  for(element_name in .lfc_new_fields()){
    expect_false(is.null(res_per_gene[[element_name]]), info = element_name)
    expect_equal(res_per_gene[[element_name]], esvd_obj[[element_name]],
                 tolerance = 1e-8, info = element_name)
  }
  expect_setequal(as.character(res_per_gene$param$test_case_individuals),
                  as.character(esvd_obj$param$test_case_individuals))
  expect_setequal(as.character(res_per_gene$param$test_control_individuals),
                  as.character(esvd_obj$param$test_control_individuals))
})

test_that("T-LFC-11: compute_log_fold_change matches the posterior matrices", {
  esvd_obj <- .lfc_fitted_small()
  oracle <- .lfc_from_donor_averages(.lfc_donor_averages(esvd_obj))

  stripped <- esvd_obj
  stripped$log2fc_vec <- NULL
  stripped$log2fc_se_vec <- NULL
  res <- compute_log_fold_change(input_obj = stripped)

  expect_true(inherits(res, "eSVD"))
  for(element_name in .lfc_new_fields()){
    expect_equal(res[[element_name]], oracle[[element_name]],
                 tolerance = 1e-10, info = element_name)
  }

  # What the test-statistic step stored is what the standalone call returns,
  # and a second call changes nothing.
  expect_equal(res$log2fc_vec, esvd_obj$log2fc_vec, tolerance = 1e-12)
  expect_equal(res$log2fc_se_vec, esvd_obj$log2fc_se_vec, tolerance = 1e-12)
  expect_identical(compute_log_fold_change(input_obj = res), res)
})

test_that("T-LFC-12: report_results carries the fold change and its SE", {
  esvd_obj <- .small_esvd_obj()
  res <- report_results(esvd_obj)

  expect_true(all(c("logFC", "logFC_se") %in% colnames(res)))
  expect_identical(which(colnames(res) == "logFC_se"),
                   which(colnames(res) == "logFC") + 1L)
  expect_equal(res$logFC_se, unname(esvd_obj$log2fc_se_vec),
               tolerance = 1e-12)
  expect_equal(res$logFC, unname(esvd_obj$log2fc_vec), tolerance = 1e-12)
  expect_identical(res$genes, names(esvd_obj$log2fc_se_vec))

  expect_true(all(is.finite(res$logFC_se)))
  expect_true(all(res$logFC_se > 0))

  # Both columns are read from the stored vectors, so a gene that is NA there
  # is NA in both. Recomputing `logFC` from the arm means would put a finite
  # number, or Inf, beside an NA standard error.
  blanked_obj <- esvd_obj
  blanked_obj$log2fc_vec[4] <- NA_real_
  blanked_obj$log2fc_se_vec[4] <- NA_real_
  res_blanked <- report_results(blanked_obj)
  expect_true(is.na(res_blanked$logFC[4]))
  expect_true(is.na(res_blanked$logFC_se[4]))
  expect_equal(res_blanked$logFC[-4], res$logFC[-4], tolerance = 1e-12)
})

test_that("T-LFC-13: both eSVD() paths and eSVD_helper() carry the new vectors", {
  skip_if_not_installed("SeuratObject")

  res_diet <- .esvd_run(.tiny_counts(), bool_diet = TRUE, max_iter = 10)
  res_full <- .esvd_run(.tiny_counts(), bool_diet = FALSE, max_iter = 10)
  for(element_name in .lfc_new_fields()){
    expect_false(is.null(res_diet[[element_name]]), info = element_name)
    expect_equal(res_diet[[element_name]], res_full[[element_name]],
                 tolerance = 1e-8, info = element_name)
  }

  # The diet object has dropped `dat` and the posterior matrices, and the
  # standalone function still works on it: that is why the variances are
  # stored rather than recomputed.
  expect_null(res_diet[["dat"]])
  res_again <- compute_log_fold_change(input_obj = res_diet)
  expect_equal(res_again$log2fc_se_vec, res_diet$log2fc_se_vec,
               tolerance = 1e-12)

  # All-zero genes sit mid-vector (positions 3 and 12), so a reinsertion that
  # appends instead of restoring position is caught.
  dat <- .tiny_counts_with_zero_genes()
  res_helper <- .helper_run(dat, max_iter = 10)
  zero_idx <- which(res_helper[["gene_status"]] == "all_zero")
  expect_identical(unname(zero_idx), c(3L, 12L))

  for(element_name in .lfc_new_fields()){
    vec <- res_helper[[element_name]]
    expect_identical(names(vec), colnames(dat), info = element_name)
    expect_true(all(is.na(vec[zero_idx])), info = element_name)
    expect_true(all(is.finite(vec[-zero_idx])), info = element_name)
  }

  report_df <- report_results(res_helper)
  expect_true(all(is.na(report_df$logFC_se[zero_idx])))
  expect_true(all(is.finite(report_df$logFC_se[-zero_idx])))

  # Recomputing on the padded object is silent: NA in, NA out, which is not
  # the non-positive-mean case of T-LFC-08.
  expect_silent(res_padded <- compute_log_fold_change(input_obj = res_helper))
  expect_equal(res_padded$log2fc_se_vec, res_helper$log2fc_se_vec,
               tolerance = 1e-12)
})

## Decision D5. An object saved by eSVD2 < 1.1.0 has the means but not the
## variances, and for a `bool_diet = TRUE` object nothing is left to recompute
## them from. `compute_log_fold_change` must say so; `report_results` must
## keep returning the data frame it always returned, with the SE column
## missing loudly rather than absent.
test_that("T-LFC-14: an object without the variances is refused by name", {
  esvd_obj <- .small_esvd_obj()
  old_obj <- esvd_obj
  for(element_name in .lfc_new_fields()){
    old_obj[[element_name]] <- NULL
  }

  expect_error(compute_log_fold_change(input_obj = old_obj),
               regexp = "case_var")

  expect_warning(res <- report_results(old_obj), regexp = "logFC_se")
  expect_true(is.data.frame(res))
  expect_true(all(is.na(res$logFC_se)))
  expect_equal(res$logFC,
               unname(log2(esvd_obj$case_mean / esvd_obj$control_mean)),
               tolerance = 1e-12)
  expect_equal(res$log10pvalue, unname(esvd_obj$pvalue_list$log10pvalue),
               tolerance = 1e-12)
})

## Found by code review. `compute_log_fold_change` reads the number of
## individuals per arm from `param`, and `.combine_two_named_lists` never
## overwrites an entry that is already there (T-UTIL-03 pins that). So after a
## rerun on a changed cohort `param` would still list the OLD individuals, and
## a later `compute_log_fold_change` would divide by the wrong number. Both
## test functions now overwrite the two entries.
##
## The stale state is planted directly: two extra names are appended to the
## recorded case arm, as if an earlier run had had five case individuals.
test_that("T-LFC-19: rerunning the test step refreshes the recorded individuals", {
  esvd_obj <- .small_esvd_obj()
  true_case_vec <- as.character(esvd_obj$param$test_case_individuals)
  true_control_vec <- as.character(esvd_obj$param$test_control_individuals)

  stale_obj <- esvd_obj
  stale_obj$param$test_case_individuals <- c(true_case_vec, "indiv_old_1",
                                             "indiv_old_2")
  # The stale count does change the answer, so the test below has teeth.
  res_stale <- compute_log_fold_change(input_obj = stale_obj)
  expect_true(all(res_stale$log2fc_se_vec < esvd_obj$log2fc_se_vec))

  res_matrix <- compute_test_statistic(input_obj = stale_obj, verbose = 0)
  res_per_gene <- .muffle_locfdr_fallback(
    compute_test_per_gene(input_obj = stale_obj,
                          alpha_max = 2 * max(esvd_obj$dat),
                          bool_covariates_as_library = TRUE,
                          library_min = 0.1,
                          verbose = 0)
  )

  res_list <- list("compute_test_statistic" = res_matrix,
                   "compute_test_per_gene" = res_per_gene)
  for(path_name in names(res_list)){
    res <- res_list[[path_name]]
    expect_setequal(as.character(res$param$test_case_individuals),
                    true_case_vec)
    expect_setequal(as.character(res$param$test_control_individuals),
                    true_control_vec)

    res_again <- compute_log_fold_change(input_obj = res)
    expect_equal(res_again$log2fc_se_vec, esvd_obj$log2fc_se_vec,
                 tolerance = 1e-8, info = path_name)
  }
})

# ---- C. does it behave as a standard error ----------------------------------

## The donor bootstrap: hold the fit and every cell's posterior fixed,
## resample INDIVIDUALS with replacement within each arm, recompute the log2
## fold change. Its SD is a direct, formula-free measurement of how much the
## fold change moves when the individuals are redrawn and nothing else is.
##
## What it can and cannot check. The reported SE has two parts,
##
##   log2fc_se^2 = between^2 + within^2,
##
## `between` from the spread of the individual means and `within` from the
## average per-cell posterior variance. The bootstrap only sees individuals
## move, so it measures `between`, and the assertion is that
## sqrt(log2fc_se^2 - within^2) matches the bootstrap SD. That checks the
## delta-method linearization and the 1/n scaling on real fitted values, at 4
## individuals per arm, where a linearization is least safe.
##
## It says nothing about `within`, which is most of the SE: a median of 98% of
## log2fc_se^2 on F-SMALL and 95% on `generate_null()`. The reported SE is
## therefore several times the bootstrap SD (median 7.8 and 4.3 times). What
## that larger number is worth is T-LFC-17's question, not this test's.
##
## The regime. A delta method is a first-order expansion, and its error grows
## with the coefficient of variation of what is being averaged. A plain-R
## probe (gamma individual means, 3 to 15 per arm, no package code) put the
## bootstrap / delta ratio inside [0.97, 1.07] when the between-individual CV
## was at most 0.2, and as high as 1.26 at a CV of 0.5. So the assertion is
## made on the genes with CV at most 0.2.
##
## Tolerance. The bootstrap SD of B = 2000 draws has a relative standard error
## of 1 / sqrt(2B) = 1.6%, and the largest of 40 of those is about 3 of them,
## 5%. Together with the probe's 7% that is the 10% used here.
##
## Correction to the first draft, recorded because it changes what the test
## claims. The first draft applied the 10% to every gene and failed on
## `generate_null()`, whose CV reaches 1.0; the ratio there reached 1.20. It
## also asserted that the reported SE is never below the bootstrap SD. That
## is not a property of the statistic: where `within` is near zero and the CV
## is high, the reported SE is the understated `between` and falls below the
## bootstrap. On these two fits the smallest reported SE / bootstrap SD is
## 0.96, so the assertion passed only through its 5% allowance, and it was
## removed. Measured outside the regime, for the record and not asserted:
## ratios of 0.97 to 1.20 over the 18 genes of `generate_null()` with CV
## above 0.2. Inside it: 0.98 to 1.03 on both fits.
test_that("T-LFC-15: the between part of the SE agrees with a bootstrap over individuals", {
  # Thresholds on the output of an iterative fit. The pipeline amplifies
  # last-digit differences (CLAUDE_kevin.md), so a different BLAS could move
  # one gene across a threshold; run these with NOT_CRAN=true.
  skip_on_cran()

  fitted_list <- list(
    "F-SMALL" = .lfc_fitted_small(),
    "generate_null" = .lfc_fitted_null(seed_number = 10)$esvd_obj
  )

  for(fixture_name in names(fitted_list)){
    esvd_obj <- fitted_list[[fixture_name]]
    avg_list <- .lfc_donor_averages(esvd_obj)
    oracle <- .lfc_from_donor_averages(avg_list)
    boot_sd_vec <- .lfc_donor_bootstrap(avg_list, num_bootstrap = 2000)

    regime_idx <- which(.lfc_between_cv(avg_list) <= 0.2)
    # Not vacuous: the regime holds a substantial share of the genes.
    expect_true(length(regime_idx) >= 15,
                info = paste0(fixture_name, ": ", length(regime_idx),
                              " genes have CV <= 0.2"))

    between_se_vec <- sqrt(esvd_obj$log2fc_se_vec^2 - oracle$within_se_vec^2)
    ratio_vec <- (boot_sd_vec / between_se_vec)[regime_idx]
    label <- paste0(fixture_name, ": bootstrap / between-SE ranges over [",
                    signif(min(ratio_vec), 3), ", ",
                    signif(max(ratio_vec), 3), "] on ", length(regime_idx),
                    " genes")

    expect_true(all(is.finite(ratio_vec)), info = label)
    expect_true(all(abs(ratio_vec - 1) <= 0.10), info = label)
    expect_true(abs(stats::median(ratio_vec) - 1) <= 0.03, info = label)
  }
})

## Recovery of a known effect. Both generators plant one:
##
##   F-SMALL        genes 1-5, +0.8 on the natural-log scale, a log2 fold
##                  change of 0.8 / ln 2 = 1.154;
##   generate_null  genes 1-10, control = case -/+ 1 on the natural-log scale
##                  with a random sign per gene, plus the "null_large_var"
##                  genes (odd positions from 11), which are not null for
##                  this estimand (see T-LFC-17).
##
## The truth for every gene is `.lfc_planted_log2fc()` of the generator's own
## `nat_mat`.
##
## Assertions: on the planted genes the direction is right; and the truth is
## inside estimate +/- 3 SE. The 3 is the conventional bound, not a number
## tuned to these fixtures.
##
## Two things a reader should know before trusting the second assertion.
##
## 1. It is made on the genes whose nuisance estimate has not diverged, and
##    that restriction is a finding, not a convenience. Where the estimated
##    rate is near 1e7 the posterior collapses onto the fitted prediction,
##    `within` vanishes, and the SE is a few hundredths while the error in the
##    fold change is a few tenths. On F-SMALL that is 7 of 40 genes, 5 of
##    which miss their truth by more than 3 SE, among them the planted gene 2
##    (estimate 0.90, SE 0.018, truth 1.13: 12.5 SE out). On
##    `generate_null()` it is 4 of 40, all 4 more than 3 SE out and the worst
##    62. On the genes that are kept, the largest miss is 1.4 SE on F-SMALL
##    and 1.6 SE on `generate_null()`. T-LFC-18 pins the mechanism. Whether
##    the divergence itself should be repaired is question Q10 of
##    CRAN_READINESS.md.
##
## 2. The SE is wide enough to hide a bias. On F-SMALL the 35 null genes have
##    a median estimate of -0.21 against a median truth of -0.01, and the
##    planted genes come out near 0.9 against 1.15. The cause is the depth
##    adjustment: eSVD adjusts each cell for its observed total count, the
##    five planted genes lift a case cell's total by 2^0.23, and every gene's
##    fold change is lowered by about that much. It is the composition effect
##    of any total-count normalization. With standard errors of 0.3 to 0.6 it
##    sits well inside 3 SE, so this test passes with it present and does not
##    detect it.
test_that("T-LFC-16: a planted effect is recovered within 3 SE", {
  # See T-LFC-15 for why this does not run on CRAN.
  skip_on_cran()

  small_data <- .small_data()
  null_list <- .lfc_fitted_null(seed_number = 10)
  fixture_list <- list(
    "F-SMALL" = list(de_idx = 1:5,
                     esvd_obj = .lfc_fitted_small(),
                     truth_vec = .lfc_planted_log2fc(
                       nat_mat = small_data$nat_mat,
                       cc_vec = small_data$cc_vec
                     )),
    "generate_null" = list(de_idx = 1:10,
                           esvd_obj = null_list$esvd_obj,
                           truth_vec = .lfc_planted_log2fc(
                             nat_mat = null_list$nat_mat,
                             cc_vec = null_list$esvd_obj[["case_control"]]
                           ))
  )

  for(fixture_name in names(fixture_list)){
    fixture <- fixture_list[[fixture_name]]
    esvd_obj <- fixture$esvd_obj
    truth_vec <- fixture$truth_vec
    de_idx <- fixture$de_idx
    z_vec <- (esvd_obj$log2fc_vec - truth_vec) / esvd_obj$log2fc_se_vec

    # The planted effects are large: at least 1 in absolute log2 terms.
    expect_true(all(abs(truth_vec[de_idx]) > 1), info = fixture_name)
    expect_identical(unname(sign(esvd_obj$log2fc_vec[de_idx])),
                     unname(sign(truth_vec[de_idx])), info = fixture_name)

    keep_idx <- which(!.lfc_nuisance_diverged(esvd_obj))
    expect_true(length(keep_idx) >= 30,
                info = paste0(fixture_name, ": ", length(keep_idx),
                              " genes have an ordinary nuisance estimate"))
    label <- paste0(fixture_name, ": largest |estimate - truth| / SE is ",
                    signif(max(abs(z_vec[keep_idx])), 3), " at ",
                    names(z_vec)[keep_idx][which.max(abs(z_vec[keep_idx]))])
    expect_true(all(abs(z_vec[keep_idx]) <= 3), info = label)
  }
})

## Repeated sampling: the one check that does not hold the fit fixed. Thirty
## cohorts are drawn from `generate_null()` and each is refitted from scratch,
## so the spread of the fold change across cohorts includes everything the
## plug-in SE ignores (the factorization, the nuisance estimate). This is the
## test of whether the number is a standard error in the sense DESeq2's,
## dreamlet's and NEBULA's are: the size of the error one makes across
## cohorts.
##
## Only the "null_interleaved" genes (even positions from 12) are used. They
## are the genes whose two arms are drawn from the same distribution, so the
## true log2 fold change is exactly 0 whatever the depth adjustment does. The
## "null_large_var" genes (odd positions from 11) are NOT null for this
## estimand: their arms share a mean on the log scale but differ in
## within-individual SD (0.1 against 0.75), and E[exp(N(m, s^2))] =
## exp(m + s^2 / 2), so the ratio of arithmetic means is exp(0.276), a log2
## fold change of 0.40. They are null for a difference of mean logs and not
## for a log of a ratio of means.
##
## Two quantities:
##
##   calibration ratio = sqrt(mean(log2fc^2) / mean(log2fc_se^2)),
##
## the realized root mean squared error over the root mean claimed variance
## (1 is exact, below 1 is conservative); and the share of genes whose
## estimate +/- 2 SE covers 0.
##
## Bounds, fixed before the first run. Genes within a cohort share a fit, so
## the 30 cohorts are the independent units. The mean of 30 chi-square-like
## terms has relative SE sqrt(2 / 30) = 0.26, which is 0.13 on the
## square-root scale; one and a half of those gives a ratio of at most 1.2.
## Coverage is held to 95% less 3 binomial standard errors at n = 30, which is
## 0.831.
##
## Measured on the first run (seeds 101 to 130):
##
##                                  n    ratio   coverage
##   all genes                     450   0.805    0.884
##   nuisance not diverged         394   0.774    0.959
##   nuisance diverged              56   1.538    0.357
##
## So on the 88% of genes with an ordinary nuisance estimate the SE is
## slightly conservative and its coverage is nominal, even though it is
## several times the fit-fixed bootstrap SD of T-LFC-15: refitting moves the
## fold change far more than resampling individuals does, and the `within`
## term is what accounts for it. On the 12% with a diverged estimate the SE
## is anti-conservative. The assertions are made on all genes, as planned,
## and on the not-diverged genes; the diverged genes are reported above and
## not asserted, since asserting their failure would pin a defect.
##
## There is deliberately no lower bound on the ratio. How conservative the SE
## is depends on the data, so a lower bound would be a snapshot of this
## generator and not a property.
test_that("T-LFC-17: across refitted cohorts the SE is not anti-conservative", {
  skip_on_cran()

  num_cohorts <- 30
  null_idx <- seq(12, 40, by = 2)
  log2fc_mat <- matrix(NA_real_, nrow = num_cohorts, ncol = length(null_idx))
  log2fc_se_mat <- log2fc_mat
  diverged_mat <- matrix(NA, nrow = num_cohorts, ncol = length(null_idx))

  for(cohort_idx in seq_len(num_cohorts)){
    esvd_obj <- .lfc_fitted_null(seed_number = 100 + cohort_idx)$esvd_obj
    log2fc_mat[cohort_idx, ] <- esvd_obj$log2fc_vec[null_idx]
    log2fc_se_mat[cohort_idx, ] <- esvd_obj$log2fc_se_vec[null_idx]
    diverged_mat[cohort_idx, ] <- .lfc_nuisance_diverged(esvd_obj)[null_idx]
  }

  expect_true(all(is.finite(log2fc_mat)))
  expect_true(all(is.finite(log2fc_se_mat)))
  expect_true(all(log2fc_se_mat > 0))

  coverage_bound <- 0.95 - 3 * sqrt(0.95 * 0.05 / num_cohorts)
  subset_list <- list("all genes" = rep(TRUE, length(log2fc_mat)),
                      "nuisance not diverged" = as.vector(!diverged_mat))
  for(subset_name in names(subset_list)){
    keep_vec <- subset_list[[subset_name]]
    log2fc_vec <- as.vector(log2fc_mat)[keep_vec]
    log2fc_se_vec <- as.vector(log2fc_se_mat)[keep_vec]

    calibration_ratio <- sqrt(mean(log2fc_vec^2) / mean(log2fc_se_vec^2))
    coverage_val <- mean(abs(log2fc_vec) <= 2 * log2fc_se_vec)
    label <- paste0(subset_name, " (n = ", sum(keep_vec),
                    "): calibration ratio = ", signif(calibration_ratio, 3),
                    ", coverage of 0 by +/- 2 SE = ", signif(coverage_val, 3))

    # The subset is most of the genes, so this is not a test of a remnant.
    expect_true(mean(keep_vec) >= 0.75, info = label)
    expect_true(calibration_ratio <= 1.2, info = label)
    expect_true(coverage_val >= coverage_bound, info = label)
  }
})

## The mechanism behind the restriction in T-LFC-16 and T-LFC-17, as a
## property of the model and not as a pinned defect.
##
## A cell's posterior is Gamma(y + r * mu, s + r), with r the nuisance rate,
## mu the fitted prediction and s the adjusted depth. Its mean is
## (y + r * mu) / (s + r) and its variance (y + r * mu) / (s + r)^2. As r
## grows the mean tends to mu and the variance to mu / r, that is, to 0: the
## data stop mattering, every cell is assigned its fitted value, and the
## `within` part of the SE disappears. What is left is the spread of the
## fitted predictions across individuals, which carries no information about
## how wrong the fit is.
##
## The test keeps the F-SMALL fit and replaces the rates: every gene at 1,
## then every gene at 1e6. Setting them by hand, and turning off both the
## stabilization and the lower-quantile floor, means the rates handed in are
## the rates used and the test does not depend on what `estimate_nuisance`
## returned. It checks both halves: `within` falls by about 1e3 (the square
## root of 1e6), and the SE falls onto `between`.
test_that("T-LFC-18: a very large nuisance rate removes the within part of the SE", {
  esvd_obj <- .lfc_fitted_small()
  latest_fit <- esvd_obj[["latest_Fit"]]
  small_data <- .small_data()

  recompute <- function(rate_val){
    obj <- esvd_obj
    obj[[latest_fit]]$nuisance_vec[] <- rate_val
    obj <- compute_posterior(input_obj = obj,
                             alpha_max = 2 * max(small_data$dat),
                             bool_covariates_as_library = TRUE,
                             bool_stabilize_underdispersion = FALSE,
                             library_min = 0.1,
                             nuisance_lower_quantile = 0)
    obj <- compute_test_statistic(input_obj = obj, verbose = 0)
    oracle <- .lfc_from_donor_averages(.lfc_donor_averages(obj))

    list(between_se_vec = oracle$between_se_vec,
         log2fc_se_vec = obj$log2fc_se_vec,
         within_se_vec = oracle$within_se_vec)
  }
  res_base <- recompute(rate_val = 1)
  res_large <- recompute(rate_val = 1e6)

  # At a rate of 1, `within` is most of the SE ...
  within_share_vec <- res_base$within_se_vec^2 / res_base$log2fc_se_vec^2
  expect_true(all(within_share_vec > 0.5),
              info = paste0("smallest within share is ",
                            signif(min(within_share_vec), 3)))
  # ... and at 1e6 it has fallen by three orders of magnitude.
  shrink_vec <- res_large$within_se_vec / res_base$within_se_vec
  expect_true(all(shrink_vec < 1e-2),
              info = paste0("largest within ratio is ",
                            signif(max(shrink_vec), 3)))

  # The SE that is left is the between part alone.
  expect_equal(unname(res_large$log2fc_se_vec),
               unname(res_large$between_se_vec), tolerance = 1e-2)
  expect_true(all(res_large$log2fc_se_vec < res_base$log2fc_se_vec))
})
