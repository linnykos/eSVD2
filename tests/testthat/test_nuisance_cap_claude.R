# The cap on the nuisance rate -- T-CAP-01 .. T-CAP-09.
#
# A boundary gene has no finite MLE, so it is set to the cap itself.
#
# `estimate_nuisance()` returns, for gene j,
#
#   nuisance_vec[j] = max( min( MLE_j, cap_multiplier * m_j ), min_val * m_j ),
#   m_j = median_i s_ji,
#
# where MLE_j is the Gamma RATE that `gamma_rate()` returns (large = little
# over-dispersion) and s_ji is `library_mat`. Beside it, each gene gets a
# status:
#
#   failed     both `gamma_rate` and `log_gamma_rate` failed
#   boundary   D_j <= 0, with D_j = sum_i [ (A_ji - m_ji)^2 - A_ji ] / mu_ji and
#              m = mu * s. Around the Poisson limit the log-likelihood is
#              l_Poisson + D_j / (2 * beta), so D_j <= 0 means it rises all the
#              way to beta = Inf and no finite MLE exists
#   capped     D_j > 0 and MLE_j above the cap
#   estimated  the gene keeps its own MLE
#
# Oracles are tagged per test: [oracle] is an independent computation,
# [invariant] an identity the rule implies.
#
# Every fixture is built here from a fixed seed.

# ---- fixtures ---------------------------------------------------------------

# Twelve genes in three blocks of four, 300 cells. The blocks are interleaved
# (gene 1 over-dispersed, gene 2 Poisson, gene 3 under-dispersed, ...) so an
# assignment by position that is right on a sorted fixture is wrong here.
#
# The library size differs by GENE as well as by cell, as it does in a fit
# (the library includes the gene's intercept), so the cap is a different
# number for every gene. A cap that used one median for all genes would pass
# on a fixture whose library is the same in every column.
.build_cap_data <- function(num_cells = 300,
                            seed_number = 10){
  n <- num_cells
  p <- 12
  block_vec <- rep(c("overdispersed", "poisson", "underdispersed"), times = 4)

  set.seed(seed_number)
  cell_size_vec <- stats::runif(n, min = 0.5, max = 2)
  gene_size_vec <- stats::runif(p, min = 0.3, max = 4)
  library_mat <- outer(cell_size_vec, gene_size_vec)
  mean_mat <- matrix(stats::runif(n * p, min = 2, max = 6), nrow = n, ncol = p)
  fitted_mat <- mean_mat * library_mat

  set.seed(seed_number + 1)
  dat <- matrix(0, nrow = n, ncol = p)
  for(j in seq_len(p)){
    if(block_vec[j] == "overdispersed"){
      # The model itself, with a rate of 2: A ~ NB(size = mu * beta, mean = mu * s)
      dat[, j] <- stats::rnbinom(n, size = mean_mat[, j] * 2,
                                 mu = fitted_mat[, j])
    } else if(block_vec[j] == "poisson"){
      dat[, j] <- stats::rpois(n, lambda = fitted_mat[, j])
    } else {
      # Binomial with the model's mean and half the Poisson variance.
      dat[, j] <- stats::rbinom(n, size = ceiling(2 * fitted_mat[, j]),
                                prob = fitted_mat[, j] / ceiling(2 * fitted_mat[, j]))
    }
  }

  rownames(dat) <- paste0("cell", seq_len(n))
  colnames(dat) <- paste0("gene", seq_len(p))
  dimnames(mean_mat) <- dimnames(dat)
  dimnames(library_mat) <- dimnames(dat)

  list(block_vec = block_vec,
       dat = dat,
       library_mat = library_mat,
       mean_mat = mean_mat)
}

.cap_data <- function(){
  if(is.null(.fixture_cache$cap_data)){
    .fixture_cache$cap_data <- .build_cap_data()
  }
  .fixture_cache$cap_data
}

# The median of each gene's library size, computed without `matrixStats`.
.median_library <- function(library_mat){
  apply(library_mat, 2, stats::median)
}

# The marginal log-likelihood of one gene at rate `beta`, from R's own
# negative binomial rather than the formula in `src/gamma_rate.cpp`.
.loglik_at_rate <- function(x_vec, mu_vec, s_vec, beta){
  sum(stats::dnbinom(x_vec, size = mu_vec * beta, mu = mu_vec * s_vec,
                     log = TRUE))
}

# ---- the rule ---------------------------------------------------------------

## [oracle] the median is recomputed with `stats::median`, and the uncapped
## value comes from the same function at `cap_multiplier = Inf`.
##
## The fixture must hold genes on BOTH sides of the cap, or "capped genes equal
## the cap" and "the others keep their estimate" are each vacuous. That is
## asserted first.
test_that("T-CAP-01: a rate above the cap is set to cap_multiplier times the median library size, and no other rate moves", {
  dat_list <- .cap_data()
  median_vec <- .median_library(dat_list$library_mat)

  uncapped_vec <- estimate_nuisance(input_obj = dat_list$dat,
                                    mean_mat = dat_list$mean_mat,
                                    library_mat = dat_list$library_mat,
                                    cap_multiplier = Inf)

  for(cap_multiplier in c(0.5, 1, 10)){
    label <- paste0("cap_multiplier = ", cap_multiplier)
    capped_vec <- estimate_nuisance(input_obj = dat_list$dat,
                                    mean_mat = dat_list$mean_mat,
                                    library_mat = dat_list$library_mat,
                                    cap_multiplier = cap_multiplier)

    above_idx <- which(uncapped_vec > cap_multiplier * median_vec)
    below_idx <- setdiff(seq_along(uncapped_vec), above_idx)
    expect_true(length(above_idx) > 0, info = label)
    expect_true(length(below_idx) > 0, info = label)

    expect_equal(as.numeric(capped_vec[above_idx]),
                 as.numeric(cap_multiplier * median_vec[above_idx]),
                 tolerance = 1e-12, info = label)
    expect_identical(capped_vec[below_idx], uncapped_vec[below_idx],
                     info = label)
    expect_identical(names(capped_vec), colnames(dat_list$dat), info = label)
  }
})

## [oracle] the two C++ routes called directly, gene by gene, in the order the
## package tries them. `Inf` is the documented way to get the behaviour of
## version 1.1.0.
##
## `gamma_rate()` does not always return for a gene with no finite MLE: on
## this fixture its Newton iteration raises "no root to be found" (near 2.4e7)
## for two of the Poisson genes, 8 and 11, and the package then takes
## `exp(log_gamma_rate())`, which stops at its upper bracket
## exp(10) = 22026.5. So without a cap the rate of such a gene is either
## about 1e7 or exactly 22026.5, depending on which route answered.
test_that("T-CAP-02: cap_multiplier = Inf returns the maximum-likelihood rate", {
  dat_list <- .cap_data()

  res <- estimate_nuisance(input_obj = dat_list$dat,
                           mean_mat = dat_list$mean_mat,
                           library_mat = dat_list$library_mat,
                           cap_multiplier = Inf)

  num_direct <- 0
  for(j in seq_len(ncol(dat_list$dat))){
    reference <- tryCatch(
      gamma_rate(x = dat_list$dat[, j],
                 mu = dat_list$mean_mat[, j],
                 s = dat_list$library_mat[, j]),
      error = function(e){NULL}
    )
    if(is.null(reference)){
      reference <- exp(log_gamma_rate(x = dat_list$dat[, j],
                                      mu = dat_list$mean_mat[, j],
                                      s = dat_list$library_mat[, j]))
    } else {
      num_direct <- num_direct + 1
    }
    expect_equal(as.numeric(res[j]), reference, tolerance = 1e-12,
                 info = paste0("gene ", j))
  }
  # Most genes are answered by the first route, or the loop above compares
  # the fallback with itself.
  expect_true(num_direct >= 10)

  # The under-dispersed genes are the ones the cap exists for.
  expect_true(all(res[dat_list$block_vec == "underdispersed"] > 1e4))
})

test_that("T-CAP-02b: the default cap_multiplier is 10", {
  dat_list <- .cap_data()

  res_default <- estimate_nuisance(input_obj = dat_list$dat,
                                   mean_mat = dat_list$mean_mat,
                                   library_mat = dat_list$library_mat)
  res_ten <- estimate_nuisance(input_obj = dat_list$dat,
                               mean_mat = dat_list$mean_mat,
                               library_mat = dat_list$library_mat,
                               cap_multiplier = 10)
  res_inf <- estimate_nuisance(input_obj = dat_list$dat,
                               mean_mat = dat_list$mean_mat,
                               library_mat = dat_list$library_mat,
                               cap_multiplier = Inf)

  expect_identical(res_default, res_ten)
  expect_false(identical(res_default, res_inf))
})

## [invariant] a looser cap can only raise a rate and can only release genes.
test_that("T-CAP-03: rates rise and the number of capped genes falls as cap_multiplier grows", {
  dat_list <- .cap_data()
  cap_vec <- c(0.5, 1, 2, 10, 100, Inf)

  res_list <- lapply(cap_vec, function(cap_multiplier){
    .estimate_nuisance_matrix(dat = dat_list$dat,
                              mean_mat = dat_list$mean_mat,
                              library_mat = dat_list$library_mat,
                              bool_use_log = FALSE,
                              cap_multiplier = cap_multiplier,
                              min_val = 1e-4,
                              verbose = 0)
  })
  num_capped_vec <- sapply(res_list, function(res){res$num_capped})

  for(i in seq_len(length(cap_vec) - 1)){
    label <- paste0("cap_multiplier ", cap_vec[i], " against ", cap_vec[i + 1])
    expect_true(all(res_list[[i]]$nuisance_vec <= res_list[[i + 1]]$nuisance_vec),
                info = label)
    expect_true(num_capped_vec[i] >= num_capped_vec[i + 1], info = label)
  }

  # The sweep moves: without this, every comparison above holds with equality.
  expect_true(num_capped_vec[1] > num_capped_vec[4])
  expect_true(num_capped_vec[4] > 0)
  expect_equal(num_capped_vec[length(cap_vec)], 0)
})

## [invariant] The floor shares the units of the cap (Kevin, 2026-09-29), so
## on the unit-free scale every rate lies in [min_val, cap_multiplier] and no
## gene is above its cap. Under the earlier absolute floor, a gene whose
## median library size was below `min_val / cap_multiplier` sat above its cap.
test_that("T-CAP-03b: the floor is min_val times the median library size, and no rate is above its cap", {
  dat_list <- .cap_data()
  median_vec <- .median_library(dat_list$library_mat)

  grid <- expand.grid(cap_multiplier = c(1, 10),
                      min_val = c(1e-4, 0.5))
  for(i in seq_len(nrow(grid))){
    cap_multiplier <- grid$cap_multiplier[i]
    min_val <- grid$min_val[i]
    label <- paste0("cap_multiplier = ", cap_multiplier,
                    ", min_val = ", min_val)
    res <- estimate_nuisance(input_obj = dat_list$dat,
                             mean_mat = dat_list$mean_mat,
                             library_mat = dat_list$library_mat,
                             cap_multiplier = cap_multiplier,
                             min_val = min_val)
    unit_free_vec <- as.numeric(res) / median_vec

    expect_true(all(unit_free_vec >= min_val * (1 - 1e-8)), info = label)
    expect_true(all(unit_free_vec <= cap_multiplier * (1 + 1e-8)),
                info = label)
  }
})

## [oracle] the rule applied by hand to three genes whose median library
## sizes differ by a factor of 16, so that a floor or a cap that used one
## number for all genes gives the wrong answer on two of them.
test_that("T-CAP-03c: a rate below the floor is lifted to min_val times the gene's own median library size", {
  median_vec <- c(0.5, 2, 8)
  mle_vec <- c(1e-9, 1, 1e9)

  res <- .apply_nuisance_cap(nuisance_mle_vec = mle_vec,
                             library_median_vec = median_vec,
                             bool_boundary_vec = rep(FALSE, 3),
                             bool_failed_vec = rep(FALSE, 3),
                             cap_multiplier = 10,
                             min_val = 1e-2)

  expect_equal(as.numeric(res$nuisance_vec),
               c(1e-2 * 0.5, 1, 10 * 8),
               tolerance = 1e-12)
  expect_equal(as.character(res$nuisance_status),
               c("estimated", "estimated", "capped"))
  expect_equal(res$num_capped, 1)
})

## The floor and the cap share units, so a floor at or above the cap would
## put every gene above its cap; it is refused, on both methods.
test_that("T-CAP-03d: min_val must be one positive finite number below cap_multiplier", {
  dat_list <- .cap_data()
  bad_list <- list(zero = 0,
                   negative = -1,
                   at_cap = 10,
                   above_cap = 20,
                   infinite = Inf,
                   na = NA_real_,
                   two = c(1e-4, 1e-3),
                   character = "a",
                   null = NULL)

  for(name in names(bad_list)){
    expect_error(estimate_nuisance(input_obj = dat_list$dat,
                                   mean_mat = dat_list$mean_mat,
                                   library_mat = dat_list$library_mat,
                                   cap_multiplier = 10,
                                   min_val = bad_list[[name]]),
                 regexp = "`min_val`", info = name)
  }

  esvd_obj <- .small_esvd_obj()
  expect_error(estimate_nuisance(input_obj = esvd_obj,
                                 cap_multiplier = 1,
                                 min_val = 1),
               regexp = "`min_val`")
  # The default floor is below any ordinary cap.
  expect_error(estimate_nuisance(input_obj = esvd_obj,
                                 cap_multiplier = 1e-3),
               regexp = NA)
})

# ---- the status -------------------------------------------------------------

## [oracle] R's `dnbinom` and `dpois`. For a boundary gene the likelihood at a
## very large rate is still below the Poisson likelihood and still rising; for
## a gene with D > 0 the likelihood at a large rate is above the Poisson one,
## so the maximum is interior.
##
## The rates compared are 1e3 and 1e4, where D / (2 * beta) is about 1e-2 and
## 1e-3: large enough to be far outside the rounding error of a sum of 300
## log-densities, and large enough that the expansion holds.
test_that("T-CAP-04: a gene is at the boundary exactly when the Poisson likelihood is approached from below", {
  dat_list <- .cap_data()

  res <- .estimate_nuisance_matrix(dat = dat_list$dat,
                                   mean_mat = dat_list$mean_mat,
                                   library_mat = dat_list$library_mat,
                                   bool_use_log = FALSE,
                                   cap_multiplier = 10,
                                   min_val = 1e-4,
                                   verbose = 0)
  bool_boundary_vec <- res$nuisance_status == "boundary"

  # Both kinds are present.
  expect_true(sum(bool_boundary_vec) >= 4)
  expect_true(sum(!bool_boundary_vec) >= 4)
  expect_true(all(bool_boundary_vec[dat_list$block_vec == "underdispersed"]))
  expect_true(all(!bool_boundary_vec[dat_list$block_vec == "overdispersed"]))

  for(j in seq_len(ncol(dat_list$dat))){
    label <- paste0("gene ", j, " (", dat_list$block_vec[j], ")")
    x_vec <- dat_list$dat[, j]
    mu_vec <- dat_list$mean_mat[, j]
    s_vec <- dat_list$library_mat[, j]

    loglik_poisson <- sum(stats::dpois(x_vec, lambda = mu_vec * s_vec,
                                       log = TRUE))
    loglik_large <- .loglik_at_rate(x_vec, mu_vec, s_vec, beta = 1e3)
    loglik_larger <- .loglik_at_rate(x_vec, mu_vec, s_vec, beta = 1e4)

    if(bool_boundary_vec[j]){
      expect_true(loglik_large < loglik_larger, info = label)
      expect_true(loglik_larger < loglik_poisson, info = label)
    } else {
      expect_true(loglik_larger > loglik_poisson, info = label)
    }
  }
})

## Finding of the code review (2026-09-29). A boundary gene has no finite
## maximum-likelihood rate, so `min(MLE, cap)` is the cap. The rule was
## written as `pmin(value the optimizer stopped at, cap)`, and that value is
## arbitrary: about 1e7 on the first route and exactly exp(10) = 22026 on the
## log route. Where the library size is in the thousands (the library without
## the intercept is on the scale of the counts per cell) the cap is above
## 22026, and the boundary genes kept 22026 and were not counted as capped.
##
## [oracle] constructed input, the expected values written out by hand.
test_that("T-CAP-04c: a boundary gene is at the cap, whatever the optimizer returned", {
  nuisance_mle_vec <- c(gene1 = 5, gene2 = 22026, gene3 = 22026, gene4 = 50)
  library_median_vec <- c(gene1 = 1, gene2 = 3000, gene3 = 1, gene4 = 1)
  bool_boundary_vec <- c(FALSE, TRUE, TRUE, FALSE)

  res <- .apply_nuisance_cap(nuisance_mle_vec = nuisance_mle_vec,
                             library_median_vec = library_median_vec,
                             bool_boundary_vec = bool_boundary_vec,
                             bool_failed_vec = rep(FALSE, 4),
                             cap_multiplier = 10,
                             min_val = 1e-4)

  expect_equal(res$nuisance_vec,
               c(gene1 = 5, gene2 = 30000, gene3 = 10, gene4 = 10))
  expect_equal(as.character(res$nuisance_status),
               c("estimated", "boundary", "boundary", "capped"))
  expect_equal(res$num_boundary, 2)
  expect_equal(res$num_capped, 3)

  # Without a cap there is nothing to set the gene to, and it keeps what the
  # optimizer returned, as in version 1.1.0.
  res_inf <- .apply_nuisance_cap(nuisance_mle_vec = nuisance_mle_vec,
                                 library_median_vec = library_median_vec,
                                 bool_boundary_vec = bool_boundary_vec,
                                 bool_failed_vec = rep(FALSE, 4),
                                 cap_multiplier = Inf,
                                 min_val = 1e-4)
  expect_equal(res_inf$nuisance_vec, nuisance_mle_vec)
  expect_equal(res_inf$num_capped, 0)
  expect_equal(res_inf$num_boundary, 2)
})

## The same, end to end, in the regime where it matters: the fixture with its
## library multiplied by 5000 and its mean divided by 5000, which leaves every
## fitted count, and so every gene's status, as it was.
test_that("T-CAP-04d: on a library in the thousands every boundary gene is at the cap", {
  dat_list <- .cap_data()
  scale_val <- 5000
  library_mat <- dat_list$library_mat * scale_val
  mean_mat <- dat_list$mean_mat / scale_val
  cap_vec <- 10 * .median_library(library_mat)

  for(bool_use_log in c(FALSE, TRUE)){
    label <- paste0("bool_use_log = ", bool_use_log)
    res <- suppressWarnings(
      .estimate_nuisance_matrix(dat = dat_list$dat,
                                mean_mat = mean_mat,
                                library_mat = library_mat,
                                bool_use_log = bool_use_log,
                                cap_multiplier = 10,
                                min_val = 1e-4,
                                verbose = 0)
    )
    bool_boundary_vec <- res$nuisance_status == "boundary"

    expect_true(sum(bool_boundary_vec) >= 4, info = label)
    expect_equal(as.numeric(res$nuisance_vec[bool_boundary_vec]),
                 as.numeric(cap_vec[bool_boundary_vec]),
                 tolerance = 1e-12, info = label)
    expect_true(res$num_capped >= sum(bool_boundary_vec), info = label)
  }
})

## [oracle] the statistic, written out from its definition on the whole
## matrix, against the per-gene helper called column by column (the helper
## went per-gene when it moved into the loop of `.estimate_nuisance_matrix`;
## T-CAP-11 covers dense against sparse counts through that loop).
test_that("T-CAP-04b: the boundary statistic is sum of ((A - m)^2 - A) / mu", {
  dat_list <- .cap_data()
  p <- ncol(dat_list$dat)

  res <- sapply(seq_len(p), function(j){
    .compute_boundary_statistic(x_vec = dat_list$dat[, j],
                                mu_vec = dat_list$mean_mat[, j],
                                s_vec = dat_list$library_mat[, j])
  })

  fitted_mat <- dat_list$mean_mat * dat_list$library_mat
  reference <- colSums(((dat_list$dat - fitted_mat)^2 - dat_list$dat) /
                         dat_list$mean_mat)

  expect_equal(as.numeric(res), as.numeric(reference), tolerance = 1e-10)
  # Both signs are present, or the boundary tests are vacuous.
  expect_true(any(res > 0) && any(res <= 0))
})

## [invariant] the three records of what the cap did agree with each other:
## the per-gene status, the counts, and the comparison of the two rate vectors.
test_that("T-CAP-05: the counts in param agree with the status and with the rates", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  for(cap_multiplier in c(1, 10, Inf)){
    label <- paste0("cap_multiplier = ", cap_multiplier)
    res <- estimate_nuisance(input_obj = esvd_obj,
                             bool_covariates_as_library = TRUE,
                             cap_multiplier = cap_multiplier)
    fit <- res[[latest_fit]]
    gene_vec <- colnames(esvd_obj$dat)

    for(element_name in c("nuisance_vec", "nuisance_mle_vec",
                          "nuisance_library_median_vec", "nuisance_status",
                          "gene_mean_count_vec", "gene_sparsity_vec")){
      expect_identical(names(fit[[element_name]]), gene_vec,
                       info = paste0(label, ", ", element_name))
    }
    expect_true(is.factor(fit$nuisance_status), info = label)
    expect_identical(levels(fit$nuisance_status),
                     c("estimated", "capped", "boundary", "failed"),
                     info = label)

    cap_vec <- cap_multiplier * fit$nuisance_library_median_vec
    # A boundary gene is set to the cap whatever value it holds (T-CAP-04c).
    bool_at_cap_vec <- fit$nuisance_status == "boundary" &
      is.finite(cap_multiplier)
    bool_replaced_vec <- fit$nuisance_mle_vec > cap_vec | bool_at_cap_vec
    expect_equal(res$param$nuisance_num_capped, sum(bool_replaced_vec),
                 info = label)
    expect_equal(res$param$nuisance_num_boundary,
                 sum(fit$nuisance_status == "boundary"), info = label)
    expect_equal(res$param$nuisance_cap_multiplier, cap_multiplier,
                 info = label)

    # A gene is `capped` when the cap replaced its rate and it is not at the
    # boundary; `estimated` genes kept their rate.
    expect_true(all(bool_replaced_vec[fit$nuisance_status == "capped"]),
                info = label)
    expect_true(all(!bool_replaced_vec[fit$nuisance_status == "estimated"]),
                info = label)
    expected_vec <- pmin(fit$nuisance_mle_vec, cap_vec)
    expected_vec[bool_at_cap_vec] <- cap_vec[bool_at_cap_vec]
    expect_equal(as.numeric(fit$nuisance_vec), as.numeric(expected_vec),
                 tolerance = 1e-12, info = label)
  }
})

## [oracle] the median library size and the two gene summaries, recomputed
## from the matrices the function is documented to use.
test_that("T-CAP-05b: the stored median library size and gene summaries are what they say", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  res <- estimate_nuisance(input_obj = esvd_obj,
                           bool_covariates_as_library = TRUE)
  fit <- res[[latest_fit]]

  # `bool_covariates_as_library = TRUE` with the intercept: every covariate
  # but the case-control one is part of the library.
  covariates <- esvd_obj$covariates
  library_idx <- which(colnames(covariates) != "CC_1")
  library_mat <- exp(tcrossprod(covariates[, library_idx],
                                fit$z_mat[, library_idx]))

  expect_equal(as.numeric(fit$nuisance_library_median_vec),
               as.numeric(.median_library(library_mat)), tolerance = 1e-10)
  expect_equal(as.numeric(fit$gene_mean_count_vec),
               as.numeric(colMeans(as.matrix(esvd_obj$dat))),
               tolerance = 1e-12)
  expect_equal(as.numeric(fit$gene_sparsity_vec),
               as.numeric(colMeans(as.matrix(esvd_obj$dat) == 0)),
               tolerance = 1e-12)
})

## A gene that fails both routes used to be indistinguishable from a gene
## estimated at `min_val`. It has its own status, and it is not counted as
## capped: the cap did not act on it.
test_that("T-CAP-06: a gene whose estimation fails has status failed", {
  dat_list <- .cap_data()
  mean_mat <- dat_list$mean_mat
  mean_mat[, 5] <- NA_real_

  expect_warning(
    res <- .estimate_nuisance_matrix(dat = dat_list$dat,
                                     mean_mat = mean_mat,
                                     library_mat = dat_list$library_mat,
                                     bool_use_log = FALSE,
                                     cap_multiplier = 10,
                                     min_val = 1e-4,
                                     verbose = 0),
    regexp = "nuisance estimation failed for 1 of 12"
  )

  expect_equal(as.character(res$nuisance_status[5]), "failed")
  expect_equal(as.numeric(res$nuisance_vec[5]),
               1e-4 * as.numeric(.median_library(dat_list$library_mat)[5]),
               tolerance = 1e-12)
  expect_equal(res$num_failed, 1)
  expect_true(all(res$nuisance_status[-5] != "failed"))

  reference <- suppressWarnings(
    .estimate_nuisance_matrix(dat = dat_list$dat,
                              mean_mat = dat_list$mean_mat,
                              library_mat = dat_list$library_mat,
                              bool_use_log = FALSE,
                              cap_multiplier = 10,
                              min_val = 1e-4,
                              verbose = 0)
  )
  expect_identical(res$nuisance_status[-5], reference$nuisance_status[-5])
})

# ---- the argument -----------------------------------------------------------

test_that("T-CAP-07: cap_multiplier must be one positive number", {
  dat_list <- .cap_data()
  esvd_obj <- .small_esvd_obj()

  bad_list <- list(zero = 0,
                   negative = -1,
                   missing_value = NA_real_,
                   two_values = c(1, 10),
                   a_string = "10",
                   nothing = NULL)

  for(i in seq_along(bad_list)){
    label <- names(bad_list)[i]
    expect_error(estimate_nuisance(input_obj = dat_list$dat,
                                   mean_mat = dat_list$mean_mat,
                                   library_mat = dat_list$library_mat,
                                   cap_multiplier = bad_list[[i]]),
                 regexp = "cap_multiplier", info = label)
    expect_error(estimate_nuisance(input_obj = esvd_obj,
                                   cap_multiplier = bad_list[[i]]),
                 regexp = "cap_multiplier", info = label)
  }
})

## `.combine_two_named_lists()` keeps an entry that is already there
## (T-UTIL-03), so a second call at another cap would leave the first call's
## multiplier and counts in `param`, beside rates they do not describe.
test_that("T-CAP-08: estimating again at another cap overwrites what param recorded", {
  esvd_obj <- .small_esvd_obj()

  res_one <- estimate_nuisance(input_obj = esvd_obj,
                               bool_covariates_as_library = TRUE,
                               cap_multiplier = 1)
  res_two <- estimate_nuisance(input_obj = res_one,
                               bool_covariates_as_library = FALSE,
                               cap_multiplier = Inf)

  expect_equal(res_one$param$nuisance_cap_multiplier, 1)
  expect_true(res_one$param$nuisance_num_capped > 0)

  expect_equal(res_two$param$nuisance_cap_multiplier, Inf)
  expect_equal(res_two$param$nuisance_num_capped, 0)
  expect_false(res_two$param$nuisance_bool_covariates_as_library)
})

# ---- the wrappers -----------------------------------------------------------

## The cap acts after the fit, so two runs that differ only in the cap share
## their factorization and their maximum-likelihood rates exactly, and differ
## in the capped rates.
test_that("T-CAP-09: eSVD() passes cap_multiplier to the nuisance estimate, under either bool_diet", {
  skip_on_cran()
  dat <- .tiny_counts()

  for(bool_diet in c(TRUE, FALSE)){
    label <- paste0("bool_diet = ", bool_diet)
    res_default <- .esvd_run(dat = dat, bool_diet = bool_diet)
    res_ten <- .esvd_run(dat = dat, bool_diet = bool_diet, cap_multiplier = 10)
    res_one <- .esvd_run(dat = dat, bool_diet = bool_diet, cap_multiplier = 1)

    fit_default <- res_default[[res_default$latest_Fit]]
    fit_ten <- res_ten[[res_ten$latest_Fit]]
    fit_one <- res_one[[res_one$latest_Fit]]

    expect_identical(fit_default$nuisance_vec, fit_ten$nuisance_vec,
                     info = label)
    expect_equal(res_one$param$nuisance_cap_multiplier, 1, info = label)

    expect_identical(fit_one$y_mat, fit_ten$y_mat, info = label)
    expect_identical(fit_one$nuisance_mle_vec, fit_ten$nuisance_mle_vec,
                     info = label)
    expect_true(all(fit_one$nuisance_vec <= fit_ten$nuisance_vec),
                info = label)
    expect_true(res_one$param$nuisance_num_capped >
                  res_ten$param$nuisance_num_capped, info = label)
    expect_equal(as.numeric(fit_one$nuisance_vec),
                 as.numeric(pmin(fit_one$nuisance_mle_vec,
                                 fit_one$nuisance_library_median_vec)),
                 tolerance = 1e-12, info = label)

    expect_false(isTRUE(all.equal(res_one$teststat_vec, res_ten$teststat_vec)),
                 info = label)
  }
})

test_that("T-CAP-09b: eSVD_helper() passes cap_multiplier and pads the new per-gene elements", {
  skip_on_cran()
  dat <- .tiny_counts_with_zero_genes()
  zero_idx <- attr(dat, "zero_idx")
  attr(dat, "zero_idx") <- NULL

  res <- .helper_run(dat = dat, cap_multiplier = 1)
  fit <- res[[res$latest_Fit]]

  expect_equal(res$param$nuisance_cap_multiplier, 1)
  for(element_name in c("nuisance_vec", "nuisance_mle_vec",
                        "nuisance_library_median_vec", "nuisance_status",
                        "gene_mean_count_vec", "gene_sparsity_vec")){
    expect_identical(names(fit[[element_name]]), colnames(dat),
                     info = element_name)
    expect_true(all(is.na(fit[[element_name]][zero_idx])),
                info = element_name)
    expect_true(all(!is.na(fit[[element_name]][-zero_idx])),
                info = element_name)
  }
  expect_true(is.factor(fit$nuisance_status))
  expect_identical(levels(fit$nuisance_status),
                   c("estimated", "capped", "boundary", "failed"))
})

test_that("T-CAP-09c: report_results() reports each gene's nuisance status", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  result_df <- report_results(esvd_obj)

  expect_identical(colnames(result_df),
                   c("genes", "logFC", "logFC_se", "log10pvalue", "pvalue",
                     "pvalue_adj", "nuisance_status"))
  expect_identical(as.character(result_df$nuisance_status),
                   as.character(esvd_obj[[latest_fit]]$nuisance_status))
  expect_identical(result_df$genes,
                   names(esvd_obj[[latest_fit]]$nuisance_status))
})

## An object saved before the cap existed has no status. It still gets its
## data frame, with the column present and empty.
test_that("T-CAP-09d: report_results() on an object without a status gives NA", {
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]
  esvd_obj[[latest_fit]]$nuisance_status <- NULL

  result_df <- report_results(esvd_obj)

  expect_equal(ncol(result_df), 7)
  expect_true(all(is.na(result_df$nuisance_status)))
})

# ---- the library columns ----------------------------------------------------

## [oracle] the column sets of `.library_column_oracle()`, written by hand.
## Pinned before the rule's two inline copies (`compute_posterior.default`,
## `compute_test_per_gene`) were routed through this helper, so that the
## refactor is checked against the sets and not against itself. The columns
## of F-TINY are Intercept, Log_UMI, Age, CC_1, Sex_M, in that order, so the
## case-control column sits between two library columns and an index that is
## off by one lands on a real column.
test_that("T-CAP-10: .nuisance_library_idx() picks the columns written out by hand, on every setting", {
  covariates <- .tiny_covariates()
  expect_equal(colnames(covariates),
               c("Intercept", "Log_UMI", "Age", "CC_1", "Sex_M"))

  for(case in .library_column_oracle()){
    label <- paste0("cov_lib = ", case$cov_lib, ", incl_int = ",
                    case$incl_int, ", cc = ",
                    if(is.null(case$cc)) "NULL" else case$cc)
    idx_vec <- .nuisance_library_idx(
      covariates = covariates,
      case_control_variable = case$cc,
      library_size_variable = "Log_UMI",
      bool_covariates_as_library = case$cov_lib,
      bool_library_includes_interept = case$incl_int
    )

    expect_true(is.integer(idx_vec), info = label)
    expect_setequal(colnames(covariates)[idx_vec], case$columns)
    expect_equal(colnames(covariates)[idx_vec],
                 intersect(colnames(covariates), case$columns),
                 info = label)
    # Increasing, so that `covariates[, idx]` and `z_mat[, idx]` line up
    # column for column.
    expect_true(all(diff(idx_vec) > 0), info = label)
  }

  # An empty case-control vector means the same as NULL.
  expect_equal(
    .nuisance_library_idx(covariates = covariates,
                          case_control_variable = character(0),
                          library_size_variable = "Log_UMI",
                          bool_covariates_as_library = TRUE,
                          bool_library_includes_interept = TRUE),
    seq_len(5)
  )
})

## [invariant] dense and sparse counts, and the two estimation routes, go
## through one code path; the whole output list must agree, not only the
## rates. Pinned before the boundary statistic moved into the per-gene loop.
test_that("T-CAP-11: .estimate_nuisance_matrix() returns the same list on dense and sparse counts", {
  dat_list <- .cap_data()
  dat_sparse <- methods::as(methods::as(dat_list$dat * 1.0, "dMatrix"),
                            "CsparseMatrix")

  for(bool_use_log in c(FALSE, TRUE)){
    label <- paste0("bool_use_log = ", bool_use_log)
    res_dense <- .estimate_nuisance_matrix(dat = dat_list$dat,
                                           mean_mat = dat_list$mean_mat,
                                           library_mat = dat_list$library_mat,
                                           bool_use_log = bool_use_log,
                                           cap_multiplier = 10,
                                           min_val = 1e-4,
                                           verbose = 0)
    res_sparse <- .estimate_nuisance_matrix(dat = dat_sparse,
                                            mean_mat = dat_list$mean_mat,
                                            library_mat = dat_list$library_mat,
                                            bool_use_log = bool_use_log,
                                            cap_multiplier = 10,
                                            min_val = 1e-4,
                                            verbose = 0)

    expect_equal(names(res_dense), names(res_sparse), info = label)
    expect_equal(res_dense, res_sparse, tolerance = 1e-12, info = label)
    expect_true(sum(res_dense$nuisance_status == "boundary") >= 4,
                info = label)
    expect_equal(names(res_dense$nuisance_vec), colnames(dat_list$dat),
                 info = label)
  }
})
