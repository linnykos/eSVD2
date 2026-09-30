# The diagnostic plots -- T-DIAG-01 .. T-DIAG-13.
#
# `plot_nuisance()` draws each gene's unit-free rate,
#
#   nuisance_vec[j] / median_i s_ji,
#
# against the gene's mean count, its -log10 p-value or its sparsity. The cap is
# the horizontal line at `cap_multiplier`.
#
# `plot_fitted_vs_observed()` draws, for a cell and a gene, the fitted count
# m = mu * s against the observed count A, both as log(1 + .), and marks the
# pairs with |A - m| > num_sd * SD, where under the model
#
#   A ~ NB(size = mu * beta, mean = m),   Var(A) = m + m^2 / (mu * beta)
#                                                = m * (1 + s / beta).
#
# The tests assert on the data and the layers of the returned plot, and that
# the plot builds. They do not compare images.
#
# Every fixture is built here from a fixed seed.

# ---- fixtures ---------------------------------------------------------------

# An eSVD object whose "fit" is the truth the counts were drawn from, built by
# hand so that the fitted count and its SD are known without going through any
# function of the package. 400 cells, 30 genes, 8 individuals.
#
# The library is the intercept plus `Log_UMI`, so it differs by gene; the
# rates are set in units of each gene's median library size, between 0.2 and
# 1 times it. That is a variance of 2 to 6 times the Poisson one, enough
# over-dispersion that getting the rate wrong shows.
.build_truth_obj <- function(num_cells_per_individual = 50,
                             num_genes = 30,
                             num_individuals = 8,
                             seed_number = 10){
  n <- num_cells_per_individual * num_individuals
  p <- num_genes
  k <- 2

  set.seed(seed_number)
  individual_vec <- factor(rep(paste0("indiv", seq_len(num_individuals)),
                               each = num_cells_per_individual))
  cc_vec <- rep(rep(c(0, 1), each = num_individuals / 2),
                each = num_cells_per_individual)
  covariates <- cbind(Intercept = 1,
                      CC_1 = cc_vec,
                      Log_UMI = stats::rnorm(n, mean = 0, sd = 0.3))
  rownames(covariates) <- paste0("cell", seq_len(n))

  x_mat <- matrix(stats::rnorm(n * k, sd = 0.3), nrow = n, ncol = k)
  y_mat <- matrix(stats::rnorm(p * k, sd = 0.3), nrow = p, ncol = k)
  z_mat <- cbind(Intercept = stats::runif(p, min = 0, max = 2.5),
                 CC_1 = stats::rnorm(p, sd = 0.2),
                 Log_UMI = 1)
  rownames(x_mat) <- rownames(covariates)
  rownames(y_mat) <- paste0("gene", seq_len(p))
  rownames(z_mat) <- rownames(y_mat)

  library_mat <- exp(tcrossprod(covariates[, c("Intercept", "Log_UMI")],
                                z_mat[, c("Intercept", "Log_UMI")]))
  mean_mat <- exp(tcrossprod(x_mat, y_mat) +
                    tcrossprod(covariates[, "CC_1", drop = FALSE],
                               z_mat[, "CC_1", drop = FALSE]))
  library_median_vec <- apply(library_mat, 2, stats::median)
  nuisance_vec <- stats::runif(p, min = 0.2, max = 1) * library_median_vec
  names(nuisance_vec) <- rownames(y_mat)

  # The hierarchical model itself: lambda ~ Gamma(mu * beta, beta), then
  # A ~ Poisson(s * lambda).
  set.seed(seed_number + 1)
  lambda_mat <- matrix(stats::rgamma(n * p,
                                     shape = as.numeric(.mult_mat_vec(mean_mat, nuisance_vec)),
                                     rate = rep(nuisance_vec, each = n)),
                       nrow = n, ncol = p)
  dat <- matrix(stats::rpois(n * p, lambda = as.numeric(library_mat * lambda_mat)),
                nrow = n, ncol = p)
  dimnames(dat) <- list(rownames(covariates), rownames(y_mat))
  dimnames(mean_mat) <- dimnames(dat)
  dimnames(library_mat) <- dimnames(dat)

  fit <- .form_esvd_fit(x_mat = x_mat, y_mat = y_mat, z_mat = z_mat)
  fit$nuisance_vec <- nuisance_vec
  fit$nuisance_mle_vec <- nuisance_vec
  fit$nuisance_library_median_vec <- stats::setNames(library_median_vec,
                                                     rownames(y_mat))
  fit$nuisance_status <- stats::setNames(
    factor(rep("estimated", p),
           levels = c("estimated", "capped", "boundary", "failed")),
    rownames(y_mat)
  )
  fit$gene_mean_count_vec <- colMeans(dat)
  fit$gene_sparsity_vec <- colMeans(dat == 0)

  esvd_obj <- structure(
    list(dat = dat,
         covariates = covariates,
         latest_Fit = "fit_Second",
         fit_Second = fit,
         param = list(init_case_control_variable = "CC_1",
                      init_library_size_variable = "Log_UMI",
                      nuisance_bool_covariates_as_library = FALSE,
                      nuisance_bool_library_includes_interept = TRUE,
                      nuisance_cap_multiplier = 10,
                      nuisance_min_val = 1e-4),
         case_control = cc_vec,
         individual = individual_vec),
    class = "eSVD"
  )

  list(esvd_obj = esvd_obj,
       library_mat = library_mat,
       mean_mat = mean_mat)
}

.truth_obj <- function(){
  if(is.null(.fixture_cache$truth_obj)){
    .fixture_cache$truth_obj <- .build_truth_obj()
  }
  .fixture_cache$truth_obj
}

.geom_names <- function(plot_obj){
  sapply(plot_obj$layers, function(layer){class(layer$geom)[1]})
}

# ---- plot_nuisance ----------------------------------------------------------

test_that("T-DIAG-01: plot_nuisance returns a plot with one row per analyzed gene, for each x_axis", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  esvd_obj <- .small_esvd_obj()
  gene_vec <- colnames(esvd_obj$dat)

  for(x_axis in c("mean_expression", "log10pvalue", "sparsity")){
    plot1 <- plot_nuisance(input_obj = esvd_obj, x_axis = x_axis)

    expect_true(inherits(plot1, "ggplot"), info = x_axis)
    expect_identical(as.character(plot1$data$gene), gene_vec, info = x_axis)
    expect_true(all(is.finite(plot1$data$x_value)), info = x_axis)
    expect_true(all(is.finite(plot1$data$unit_free_rate)), info = x_axis)
    expect_error(ggplot2::ggplot_build(plot1), regexp = NA, info = x_axis)
  }
})

## [oracle] each axis, from the quantity it is documented to show.
test_that("T-DIAG-02: plot_nuisance draws the rate over the median library size against the stated quantity", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  esvd_obj <- .small_esvd_obj()
  fit <- esvd_obj[[esvd_obj$latest_Fit]]
  dat <- as.matrix(esvd_obj$dat)

  # The fixture has genes on both sides of the cap of 10. (At a cap of 1
  # every gene of this fixture is above it.)
  expect_equal(esvd_obj$param$nuisance_cap_multiplier, 10)
  expect_true(sum(fit$nuisance_mle_vec > 10 * fit$nuisance_library_median_vec) >= 3)
  expect_true(sum(fit$nuisance_mle_vec < 10 * fit$nuisance_library_median_vec) >= 3)

  plot1 <- plot_nuisance(input_obj = esvd_obj, x_axis = "mean_expression")
  expect_equal(plot1$data$unit_free_rate,
               as.numeric(fit$nuisance_vec / fit$nuisance_library_median_vec),
               tolerance = 1e-12)
  expect_true(all(plot1$data$unit_free_rate <= 10 + 1e-10))
  expect_true(any(plot1$data$unit_free_rate < 9))
  expect_equal(plot1$data$x_value, as.numeric(colMeans(dat)),
               tolerance = 1e-12)
  expect_identical(as.character(plot1$data$nuisance_status),
                   as.character(fit$nuisance_status))

  plot2 <- plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity")
  expect_equal(plot2$data$x_value, as.numeric(colMeans(dat == 0)),
               tolerance = 1e-12)

  plot3 <- plot_nuisance(input_obj = esvd_obj, x_axis = "log10pvalue")
  expect_equal(plot3$data$x_value,
               as.numeric(esvd_obj$pvalue_list$log10pvalue),
               tolerance = 1e-12)

  # Before the cap: the same genes, their maximum-likelihood rates.
  plot4 <- plot_nuisance(input_obj = esvd_obj, x_axis = "mean_expression",
                         bool_uncapped = TRUE)
  expect_equal(plot4$data$unit_free_rate,
               as.numeric(fit$nuisance_mle_vec / fit$nuisance_library_median_vec),
               tolerance = 1e-12)
  expect_true(any(plot4$data$unit_free_rate > 100))
})

test_that("T-DIAG-03: plot_nuisance draws the cap as a line, and no line when there is no cap", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  esvd_obj <- .small_esvd_obj()

  plot1 <- plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity")
  expect_true("GeomHline" %in% .geom_names(plot1))
  hline_idx <- which(.geom_names(plot1) == "GeomHline")
  expect_equal(plot1$layers[[hline_idx]]$data$yintercept, 10)

  uncapped_obj <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = esvd_obj, cap_multiplier = Inf)
  )
  plot2 <- plot_nuisance(input_obj = uncapped_obj, x_axis = "sparsity")
  expect_false("GeomHline" %in% .geom_names(plot2))
  expect_error(ggplot2::ggplot_build(plot2), regexp = NA)
})

## Finding of the by-hand run (2026-09-29): under the absolute floor of the
## first draft, a gene whose fit is degenerate (median library size near
## 1e-8, cap 1e-7) was held at `min_val` (1e-4), thousands of times above its
## cap, and the plot drew it as "capped" far above the line. The floor is now
## a multiple of the library size (Kevin, 2026-09-29), so the plot must show
## no gene above the line, whatever the library size.
##
## The gene is built through `.apply_nuisance_cap()` and not by hand, so that
## the test follows the rule if the rule changes.
test_that("T-DIAG-03b: plot_nuisance draws no gene above the cap, whatever its library size", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  esvd_obj <- .truth_obj()$esvd_obj

  fit <- esvd_obj$fit_Second
  p <- length(fit$nuisance_vec)
  fit$nuisance_library_median_vec[c(3, 7)] <- 1e-8
  cap_res <- .apply_nuisance_cap(
    nuisance_mle_vec = fit$nuisance_mle_vec,
    library_median_vec = fit$nuisance_library_median_vec,
    bool_boundary_vec = rep(FALSE, p),
    bool_failed_vec = rep(FALSE, p),
    cap_multiplier = 10,
    min_val = 1e-4
  )
  fit$nuisance_vec <- cap_res$nuisance_vec
  fit$nuisance_status <- cap_res$nuisance_status
  esvd_obj$fit_Second <- fit

  plot2 <- plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity")
  expect_true(all(plot2$data$unit_free_rate <= 10 * (1 + 1e-8)))
  expect_equal(as.numeric(plot2$data$unit_free_rate[c(3, 7)]), c(10, 10),
               tolerance = 1e-8)
  expect_equal(as.character(plot2$data$nuisance_status[c(3, 7)]),
               c("capped", "capped"))
  expect_false(grepl("above the cap", plot2$labels$subtitle, fixed = TRUE))
  expect_error(ggplot2::ggplot_build(plot2), regexp = NA)

  # Drawn before the cap, a rate above the line is what the plot is for.
  plot3 <- plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity",
                         bool_uncapped = TRUE)
  expect_true(any(plot3$data$unit_free_rate[c(3, 7)] > 10))
})

test_that("T-DIAG-04: plot_nuisance labels the genes asked for, and by default the capped genes with the smallest p-values", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  esvd_obj <- .muffle_locfdr_fallback(
    recompute_pvalue(input_obj = .small_esvd_obj(), cap_multiplier = 1)
  )
  fit <- esvd_obj[[esvd_obj$latest_Fit]]
  gene_vec <- colnames(esvd_obj$dat)

  # Reversed, and not the first genes: the labels follow the names.
  asked_vec <- gene_vec[c(17, 9, 4)]
  plot1 <- plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity",
                         genes = asked_vec)
  expect_setequal(as.character(plot1$data$gene[plot1$data$bool_label]),
                  asked_vec)
  expect_true("GeomTextRepel" %in% .geom_names(plot1))

  # The default: among the genes that are not `estimated`, the `num_label`
  # with the largest -log10 p-value.
  flagged_vec <- gene_vec[as.character(fit$nuisance_status) != "estimated"]
  expect_true(length(flagged_vec) > 3)
  log10pvalue_vec <- esvd_obj$pvalue_list$log10pvalue[flagged_vec]
  expected_vec <- names(sort(log10pvalue_vec, decreasing = TRUE))[1:3]

  plot2 <- plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity",
                         num_label = 3)
  expect_setequal(as.character(plot2$data$gene[plot2$data$bool_label]),
                  expected_vec)

  plot3 <- plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity",
                         num_label = 0)
  expect_false(any(plot3$data$bool_label))
  expect_error(ggplot2::ggplot_build(plot3), regexp = NA)
})

test_that("T-DIAG-05: plot_nuisance refuses what it cannot draw", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  esvd_obj <- .small_esvd_obj()
  latest_fit <- esvd_obj[["latest_Fit"]]

  expect_error(plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity",
                             genes = c("gene1", "not_a_gene")),
               regexp = "not_a_gene")
  expect_error(plot_nuisance(input_obj = esvd_obj, x_axis = "expression"),
               regexp = "x_axis")

  untested_obj <- esvd_obj
  untested_obj$pvalue_list <- NULL
  expect_error(plot_nuisance(input_obj = untested_obj,
                             x_axis = "log10pvalue"),
               regexp = "pvalue_list")
  # The other two axes do not need the test.
  expect_error(plot_nuisance(input_obj = untested_obj, x_axis = "sparsity"),
               regexp = NA)

  old_obj <- esvd_obj
  old_obj[[latest_fit]]$nuisance_library_median_vec <- NULL
  expect_error(plot_nuisance(input_obj = old_obj, x_axis = "sparsity"),
               regexp = "nuisance_library_median_vec")

  unnamed_obj <- esvd_obj
  rownames(unnamed_obj[[latest_fit]]$y_mat) <- NULL
  expect_error(plot_nuisance(input_obj = unnamed_obj, x_axis = "sparsity"),
               regexp = "gene names")
})

## The three scatter plots need no counts, so they work on the object `eSVD()`
## returns by default. The reinserted all-zero genes are not drawn.
test_that("T-DIAG-06: plot_nuisance works on a diet object and leaves out the genes that were not analyzed", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  dat <- .tiny_counts_with_zero_genes()
  zero_idx <- attr(dat, "zero_idx")
  attr(dat, "zero_idx") <- NULL
  esvd_obj <- .helper_run(dat = dat, bool_diet = TRUE)
  expect_null(esvd_obj$dat)

  for(x_axis in c("mean_expression", "log10pvalue", "sparsity")){
    plot1 <- plot_nuisance(input_obj = esvd_obj, x_axis = x_axis)

    expect_identical(as.character(plot1$data$gene),
                     colnames(dat)[-zero_idx], info = x_axis)
    expect_true(all(is.finite(plot1$data$unit_free_rate)), info = x_axis)
    expect_error(ggplot2::ggplot_build(plot1), regexp = NA, info = x_axis)
  }
})

# ---- plot_fitted_vs_observed ------------------------------------------------

## [oracle] the fitted count and its variance from the generator's own
## matrices, the variance in the negative-binomial form m + m^2 / size.
test_that("T-DIAG-07: plot_fitted_vs_observed marks the pairs beyond num_sd standard deviations", {
  skip_if_not_installed("ggplot2")
  truth_list <- .truth_obj()
  esvd_obj <- truth_list$esvd_obj
  nuisance_vec <- esvd_obj$fit_Second$nuisance_vec
  # Not the first genes, and not in the order of the fit.
  gene_vec <- colnames(esvd_obj$dat)[c(22, 3, 11)]

  for(num_sd in c(1, 3)){
    label <- paste0("num_sd = ", num_sd)
    plot1 <- plot_fitted_vs_observed(input_obj = esvd_obj,
                                     genes = gene_vec,
                                     num_sd = num_sd)
    plot_df <- plot1$data

    expect_true(inherits(plot1, "ggplot"), info = label)
    expect_equal(nrow(plot_df), length(gene_vec) * nrow(esvd_obj$dat),
                 info = label)
    expect_setequal(as.character(plot_df$gene), gene_vec)

    cell_vec <- as.character(plot_df$cell)
    pair_mat <- cbind(cell_vec, as.character(plot_df$gene))
    observed_vec <- esvd_obj$dat[pair_mat]
    mean_vec <- truth_list$mean_mat[pair_mat]
    fitted_vec <- mean_vec * truth_list$library_mat[pair_mat]
    size_vec <- mean_vec * nuisance_vec[as.character(plot_df$gene)]
    sd_vec <- sqrt(fitted_vec + fitted_vec^2 / size_vec)

    expect_equal(plot_df$observed, as.numeric(observed_vec), info = label)
    expect_equal(plot_df$fitted, as.numeric(fitted_vec), tolerance = 1e-10,
                 info = label)
    expect_equal(plot_df$sd, as.numeric(sd_vec), tolerance = 1e-10,
                 info = label)
    expect_identical(plot_df$bool_outside,
                     as.logical(abs(observed_vec - fitted_vec) > num_sd * sd_vec),
                     info = label)
    expect_equal(plot_df$log_observed, log1p(as.numeric(observed_vec)),
                 info = label)
    expect_equal(plot_df$log_fitted, log1p(as.numeric(fitted_vec)),
                 tolerance = 1e-10, info = label)

    # Both kinds of point are present, or the comparison of flags is empty.
    expect_true(any(plot_df$bool_outside), info = label)
    expect_true(any(!plot_df$bool_outside), info = label)
    expect_error(ggplot2::ggplot_build(plot1), regexp = NA, info = label)
  }
})

## [invariant] the bar is the interval of counts within num_sd standard
## deviations of the fit, drawn on the y-axis at the observed count. So a pair
## is marked exactly when its bar does not reach the line y = x.
test_that("T-DIAG-08: a pair is marked exactly when its bar misses the diagonal", {
  skip_if_not_installed("ggplot2")
  esvd_obj <- .truth_obj()$esvd_obj

  plot1 <- plot_fitted_vs_observed(input_obj = esvd_obj,
                                   bool_draw_bars = TRUE,
                                   max_points = 5000)
  plot_df <- plot1$data

  expect_true(all(plot_df$log_lower <= plot_df$log_fitted))
  expect_true(all(plot_df$log_upper >= plot_df$log_fitted))
  expect_true(all(plot_df$log_lower >= 0))

  bool_crosses_vec <- plot_df$log_lower <= plot_df$log_observed &
    plot_df$log_observed <= plot_df$log_upper
  expect_identical(plot_df$bool_outside, !bool_crosses_vec)
  expect_true(any(plot_df$bool_outside))
})

## Found by breaking the code on purpose (2026-09-29): with the lower end of
## the bar halved, every test passed. On counts drawn from the model the fit
## minus 3 SD is almost always negative, so the lower end is 0 and no count
## can lie below it: T-DIAG-08 never saw a pair marked from BELOW.
##
## Here one gene is given a rate a thousand times its library size, so that
## its SD is the Poisson one, sqrt(m), and the first ten cells a count of 0.
## For a fitted count above 9 the lower end, m - 3 * sqrt(m), is positive.
##
## [oracle] the ends of the bar written out from their definition.
test_that("T-DIAG-08b: a count far below the fit is marked, and the bar stops at fit minus num_sd SD", {
  skip_if_not_installed("ggplot2")
  truth_list <- .truth_obj()
  esvd_obj <- truth_list$esvd_obj

  # The gene with the largest fitted counts, and its ten cells with the
  # largest of them.
  fitted_mat <- truth_list$mean_mat * truth_list$library_mat
  gene_idx <- which.max(apply(fitted_mat, 2, stats::median))
  gene_name <- colnames(fitted_mat)[gene_idx]
  cell_idx_vec <- order(fitted_mat[, gene_idx], decreasing = TRUE)[seq_len(10)]
  expect_true(all(fitted_mat[cell_idx_vec, gene_idx] > 9))

  fit <- esvd_obj$fit_Second
  fit$nuisance_vec[gene_name] <- 1000 * fit$nuisance_library_median_vec[gene_name]
  esvd_obj$fit_Second <- fit
  esvd_obj$dat[cell_idx_vec, gene_name] <- 0

  plot1 <- plot_fitted_vs_observed(input_obj = esvd_obj,
                                   bool_draw_bars = TRUE,
                                   genes = gene_name)
  plot_df <- plot1$data
  plot_df <- plot_df[match(rownames(esvd_obj$dat), plot_df$cell), ]

  m_vec <- fitted_mat[, gene_idx]
  sd_vec <- sqrt(m_vec * (1 + truth_list$library_mat[, gene_idx] /
                            fit$nuisance_vec[gene_name]))
  expect_equal(plot_df$log_lower, as.numeric(log1p(pmax(m_vec - 3 * sd_vec, 0))),
               tolerance = 1e-10)
  expect_equal(plot_df$log_upper, as.numeric(log1p(m_vec + 3 * sd_vec)),
               tolerance = 1e-10)

  expect_true(all(plot_df$log_lower[cell_idx_vec] > 0))
  expect_true(all(plot_df$bool_outside[cell_idx_vec]))
  expect_true(all(plot_df$log_observed[cell_idx_vec] <
                    plot_df$log_lower[cell_idx_vec]))

  bool_crosses_vec <- plot_df$log_lower <= plot_df$log_observed &
    plot_df$log_observed <= plot_df$log_upper
  expect_identical(plot_df$bool_outside, !bool_crosses_vec)
})

test_that("T-DIAG-09: the bars and the diagonal are drawn when they should be", {
  skip_if_not_installed("ggplot2")
  esvd_obj <- .truth_obj()$esvd_obj

  plot1 <- plot_fitted_vs_observed(input_obj = esvd_obj, max_points = 2000)
  expect_true("GeomPoint" %in% .geom_names(plot1))
  expect_true("GeomAbline" %in% .geom_names(plot1))
  expect_false("GeomLinerange" %in% .geom_names(plot1))

  plot2 <- plot_fitted_vs_observed(input_obj = esvd_obj, max_points = 2000,
                                   bool_draw_bars = TRUE)
  expect_true("GeomLinerange" %in% .geom_names(plot2))
  expect_error(ggplot2::ggplot_build(plot2), regexp = NA)
})

## [oracle] simulation, the counts drawn from the model and the fit set to the
## truth. Measured share of pairs beyond 3 SD, three seeds, 12000 pairs each:
##
##     rate used / true rate     0.1       1        10
##     seed 10                  0.0002   0.0141   0.0467
##     seed 20                  0.0001   0.0163   0.0534
##     seed 30                  0.0003   0.0129   0.0491
##
## At the truth the share is above the Gaussian 0.3% because a count is skewed
## to the right. A rate that is too large is an over-dispersion that is too
## small: the bars are too short and the share more than triples. A rate that
## is too small makes the bars too long and nearly nothing is marked.
##
## This is the test that the plot shows what it is for. The thresholds are
## set at about half the measured contrast.
test_that("T-DIAG-10: the marked share is small at the true rates, grows when the over-dispersion is understated and vanishes when it is overstated", {
  skip_if_not_installed("ggplot2")

  for(seed_number in c(10, 20, 30)){
    label <- paste0("seed ", seed_number)
    esvd_obj <- .build_truth_obj(seed_number = seed_number)$esvd_obj
    num_pairs <- nrow(esvd_obj$dat) * ncol(esvd_obj$dat)

    fraction_vec <- sapply(c(0.1, 1, 10), function(multiplier){
      obj <- esvd_obj
      obj$fit_Second$nuisance_vec <- multiplier * esvd_obj$fit_Second$nuisance_vec
      plot1 <- plot_fitted_vs_observed(input_obj = obj,
                                       max_points = num_pairs)
      stopifnot(nrow(plot1$data) == num_pairs)
      mean(plot1$data$bool_outside)
    })
    label <- paste0(label, ": shares ",
                    paste0(signif(fraction_vec, 3), collapse = ", "))

    expect_true(fraction_vec[2] > 0.003 && fraction_vec[2] < 0.025,
                info = label)
    expect_true(fraction_vec[3] > 2 * fraction_vec[2], info = label)
    expect_true(fraction_vec[1] < fraction_vec[2] / 5, info = label)
  }
})

test_that("T-DIAG-11: max_points bounds the pairs drawn, and seed_number decides which", {
  skip_if_not_installed("ggplot2")
  esvd_obj <- .truth_obj()$esvd_obj
  num_pairs <- nrow(esvd_obj$dat) * ncol(esvd_obj$dat)

  plot1 <- plot_fitted_vs_observed(input_obj = esvd_obj, max_points = 500)
  plot2 <- plot_fitted_vs_observed(input_obj = esvd_obj, max_points = 500)
  plot3 <- plot_fitted_vs_observed(input_obj = esvd_obj, max_points = 500,
                                   seed_number = 11)
  pair_vec1 <- paste0(plot1$data$cell, "-", plot1$data$gene)
  pair_vec3 <- paste0(plot3$data$cell, "-", plot3$data$gene)

  expect_equal(nrow(plot1$data), 500)
  expect_equal(anyDuplicated(pair_vec1), 0)
  expect_identical(plot1$data, plot2$data)
  expect_false(setequal(pair_vec1, pair_vec3))
  # The sample spans the genes and is not the first 500 pairs.
  expect_true(length(unique(plot1$data$gene)) > 20)

  # More than there are: every pair, once.
  plot4 <- plot_fitted_vs_observed(input_obj = esvd_obj,
                                   max_points = 10 * num_pairs)
  expect_equal(nrow(plot4$data), num_pairs)
  expect_equal(anyDuplicated(paste0(plot4$data$cell, "-", plot4$data$gene)), 0)
})

## On a fitted object the fit and the library are the ones `estimate_nuisance`
## used, whatever the posterior did with the rates afterwards. [oracle] the
## matrices rebuilt as `estimate_nuisance.eSVD` documents them.
test_that("T-DIAG-12: on a fitted object the fitted count is exp of the natural parameter, from counts or from the Seurat object", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  seurat_obj <- .tiny_seurat()
  full_obj <- .esvd_run(dat = .tiny_counts(), bool_diet = FALSE)
  diet_obj <- .esvd_run(dat = .tiny_counts(), bool_diet = TRUE)
  fit <- full_obj[[full_obj$latest_Fit]]
  num_pairs <- nrow(full_obj$dat) * ncol(full_obj$dat)

  plot_full <- plot_fitted_vs_observed(input_obj = full_obj,
                                       max_points = num_pairs)
  plot_diet <- plot_fitted_vs_observed(input_obj = diet_obj,
                                       max_points = num_pairs,
                                       seurat_obj = seurat_obj)
  expect_equal(plot_diet$data, plot_full$data, tolerance = 1e-10)
  expect_error(plot_fitted_vs_observed(input_obj = diet_obj),
               regexp = "seurat_obj")

  plot_df <- plot_full$data
  pair_mat <- cbind(as.character(plot_df$cell), as.character(plot_df$gene))
  fitted_mat <- exp(tcrossprod(fit$x_mat, fit$y_mat) +
                      tcrossprod(full_obj$covariates, fit$z_mat))
  # `eSVD()` estimates the rates with every covariate but the case-control
  # one in the library.
  library_idx <- which(colnames(full_obj$covariates) != "CC_1")
  library_mat <- exp(tcrossprod(full_obj$covariates[, library_idx],
                                fit$z_mat[, library_idx]))
  dimnames(fitted_mat) <- dimnames(full_obj$dat)
  dimnames(library_mat) <- dimnames(full_obj$dat)
  fitted_vec <- fitted_mat[pair_mat]
  sd_vec <- sqrt(fitted_vec * (1 + library_mat[pair_mat] /
                                 fit$nuisance_vec[pair_mat[, 2]]))

  expect_equal(plot_df$observed, as.numeric(full_obj$dat[pair_mat]))
  expect_equal(plot_df$fitted, as.numeric(fitted_vec), tolerance = 1e-10)
  expect_equal(plot_df$sd, as.numeric(sd_vec), tolerance = 1e-10)
})

test_that("T-DIAG-13: plot_fitted_vs_observed refuses what it cannot draw", {
  skip_if_not_installed("ggplot2")
  esvd_obj <- .truth_obj()$esvd_obj

  expect_error(plot_fitted_vs_observed(input_obj = esvd_obj,
                                       genes = c("gene1", "not_a_gene")),
               regexp = "not_a_gene")
  expect_error(plot_fitted_vs_observed(input_obj = esvd_obj, num_sd = 0),
               regexp = "num_sd")
  expect_error(plot_fitted_vs_observed(input_obj = esvd_obj, max_points = 0),
               regexp = "max_points")

  no_rate_obj <- esvd_obj
  no_rate_obj$fit_Second$nuisance_vec <- NULL
  expect_error(plot_fitted_vs_observed(input_obj = no_rate_obj),
               regexp = "nuisance_vec")

  unnamed_obj <- esvd_obj
  rownames(unnamed_obj$fit_Second$y_mat) <- NULL
  expect_error(plot_fitted_vs_observed(input_obj = unnamed_obj),
               regexp = "gene names")
})
