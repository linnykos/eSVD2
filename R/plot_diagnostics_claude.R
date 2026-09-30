# `.data` is the pronoun of tidy evaluation, used inside `ggplot2::aes()`.
# 'ggplot2' is in Suggests, so it cannot be imported, and `ggplot2::.data$x`
# does not evaluate; declaring the name is what keeps `R CMD check` from
# reporting it as an undefined global.
utils::globalVariables(".data")

#' Plot each gene's nuisance rate against a summary of the gene
#'
#' A diagnostic of the over-dispersion. Each point is a gene. The y-axis is the
#' gene's \strong{unit-free rate}, its Gamma rate divided by the median over
#' cells of its library size, \eqn{\hat\beta_j / \mathrm{median}_i\,
#' \ell_{ji}}, on a \eqn{\log_{10}} scale. A \strong{larger} value means
#' \strong{less} over-dispersion. The cap of \code{estimate_nuisance} is the
#' same number for every gene on this scale, and is drawn as a dashed line at
#' \code{cap_multiplier}; the genes it acted on sit on the line.
#'
#' What to look for. Most genes should be \code{estimated}, below the line.
#' A large share at the line means the cap, and not the data, is setting the
#' over-dispersion. Against \code{"mean_expression"} and \code{"sparsity"}, the
#' genes at the line are expected among the lowly expressed and the sparse,
#' whose counts carry little information about the rate. Against
#' \code{"log10pvalue"}, the genes with the smallest p-values should not be
#' predominantly the ones at the line: a gene whose rate is very large has a
#' posterior that follows the fit, and an inflated statistic. Drawing the
#' plot with \code{bool_uncapped = TRUE} shows how far above the line the
#' rates were.
#'
#' The function needs neither the counts nor the covariates, so it works on
#' an object built with \code{bool_diet = TRUE}. Genes that were not analyzed
#' (the all-zero genes \code{eSVD_helper} reinserts) are not drawn.
#'
#' @param input_obj      \code{eSVD} object after \code{estimate_nuisance}
#'                       (version 1.2.0 or later), and after the test when
#'                       \code{x_axis = "log10pvalue"}.
#' @param x_axis         One of \code{"mean_expression"} (the gene's mean count
#'                       over cells, \eqn{\log_{10}} scale),
#'                       \code{"log10pvalue"} (\eqn{-\log_{10}} of its
#'                       p-value) or \code{"sparsity"} (the fraction of cells
#'                       in which its count is zero).
#' @param bool_uncapped  If \code{TRUE}, draw the rate before the cap
#'                       (\code{nuisance_mle_vec}). Default \code{FALSE}, the
#'                       rate the analysis used.
#' @param genes          \code{NULL}, or a character vector of genes to label.
#' @param num_label      When \code{genes} is \code{NULL}, how many genes to
#'                       label: the genes whose status is not
#'                       \code{estimated}, in decreasing order of
#'                       \eqn{-\log_{10}} p-value (or of the rate before the
#'                       cap, if the object has not been tested). \code{0}
#'                       labels none.
#'
#' @returns A \code{ggplot} object. Its \code{data} has one row per analyzed
#' gene, in the order of the fit, and columns \code{gene}, \code{x_value},
#' \code{unit_free_rate}, \code{nuisance_status} and \code{bool_label}.
#' @examples
#' set.seed(10)
#' sim <- generate_null(cell_per_person = 15, num_genes = 40,
#'                      num_individuals = 8)
#' esvd_obj <- initialize_esvd(dat = sim$obs_mat,
#'                             covariates = sim$covariates,
#'                             metadata_individual = sim$metadata_individual,
#'                             case_control_variable = "CC",
#'                             bool_intercept = TRUE,
#'                             k = 2,
#'                             lambda = 0.1)
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_Init",
#'                                                omitted_variables = "Log_UMI")
#' esvd_obj <- opt_esvd(input_obj = esvd_obj,
#'                      max_iter = 5,
#'                      offset_variables = setdiff(colnames(esvd_obj$covariates), "CC"),
#'                      fit_name = "fit_First",
#'                      fit_previous = "fit_Init")
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_First",
#'                                                omitted_variables = "Log_UMI")
#' esvd_obj <- estimate_nuisance(input_obj = esvd_obj)
#' if(requireNamespace("ggplot2", quietly = TRUE) &&
#'    requireNamespace("ggrepel", quietly = TRUE)){
#'   plot1 <- plot_nuisance(input_obj = esvd_obj, x_axis = "sparsity")
#'   plot(plot1)
#' }
#' @export
plot_nuisance <- function(input_obj,
                          x_axis,
                          bool_uncapped = FALSE,
                          genes = NULL,
                          num_label = 10){
  .check_plot_packages(c("ggplot2", "ggrepel"), "plot_nuisance")
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            is.logical(bool_uncapped), length(bool_uncapped) == 1,
            is.null(genes) || is.character(genes),
            is.numeric(num_label), length(num_label) == 1, num_label >= 0)
  .check_fit_is_named(input_obj, bool_cells = FALSE)
  x_axis_vec <- c("mean_expression", "log10pvalue", "sparsity")
  if(!is.character(x_axis) || length(x_axis) != 1 || !x_axis %in% x_axis_vec){
    stop("`x_axis` must be one of \"", paste0(x_axis_vec, collapse = "\", \""),
         "\"")
  }

  plot_df <- .form_nuisance_df(input_obj = input_obj,
                               x_axis = x_axis,
                               bool_uncapped = bool_uncapped)
  plot_df$bool_label <- .choose_labeled_genes(input_obj = input_obj,
                                              plot_df = plot_df,
                                              genes = genes,
                                              num_label = num_label)
  cap_multiplier <- .get_param_or_default(input_obj,
                                          "nuisance_cap_multiplier", Inf)

  xlab_val <- c(mean_expression = "Mean count of the gene",
                log10pvalue = "-log10 p-value",
                sparsity = "Fraction of cells with a zero count")[[x_axis]]
  status_table <- table(plot_df$nuisance_status)
  subtitle_val <- paste0(
    nrow(plot_df), " genes: ",
    paste0(as.integer(status_table), " ", names(status_table), collapse = ", "),
    if(is.finite(cap_multiplier)) paste0(". Dashed line: the cap, ",
                                        cap_multiplier)
  )

  plot1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x_value,
                                                 y = .data$unit_free_rate))
  if(is.finite(cap_multiplier)){
    plot1 <- plot1 + ggplot2::geom_hline(
      data = data.frame(yintercept = cap_multiplier),
      mapping = ggplot2::aes(yintercept = .data$yintercept),
      linetype = "dashed", color = "#52514e"
    )
  }
  # `show.legend = TRUE`: a status no gene has is still drawn in the legend.
  plot1 <- plot1 + ggplot2::geom_point(
    ggplot2::aes(color = .data$nuisance_status,
                 shape = .data$nuisance_status),
    size = 2, show.legend = TRUE
  )
  # Every level is kept in the legend, so a color means the same status in
  # every plot, whichever statuses the data hold.
  plot1 <- plot1 + ggplot2::scale_color_manual(
    values = .nuisance_status_colors(), drop = FALSE, name = "Status"
  )
  plot1 <- plot1 + ggplot2::scale_shape_manual(
    values = .nuisance_status_shapes(), drop = FALSE, name = "Status"
  )
  if(any(plot_df$bool_label)){
    plot1 <- plot1 + ggrepel::geom_text_repel(
      data = plot_df[plot_df$bool_label, , drop = FALSE],
      mapping = ggplot2::aes(label = .data$gene),
      color = "#0b0b0b", size = 3, min.segment.length = 0,
      max.overlaps = Inf
    )
  }
  plot1 <- plot1 + ggplot2::scale_y_log10()
  if(x_axis == "mean_expression") plot1 <- plot1 + ggplot2::scale_x_log10()
  plot1 <- plot1 + ggplot2::labs(
    x = xlab_val,
    y = paste0(if(bool_uncapped) "Rate before the cap" else "Rate",
               " / median library size\n(larger = less over-dispersion)"),
    title = paste0("Nuisance rate against ",
                   c(mean_expression = "mean expression",
                     log10pvalue = "significance",
                     sparsity = "sparsity")[[x_axis]]),
    subtitle = subtitle_val
  )
  plot1 <- plot1 + ggplot2::theme_bw()

  plot1
}

#' Plot the fitted count against the observed count
#'
#' A diagnostic of whether the over-dispersion accounts for how far the counts
#' lie from the fit. Each point is one cell and one gene. The x-axis is the
#' observed count \eqn{A_{ji}} and the y-axis the fitted count
#' \eqn{m_{ji} = \ell_{ji}\mu_{ji}}, both drawn as \eqn{\log(1 + \cdot)}, so a
#' count that equals its fit is on the line \eqn{y = x}. Under the model
#' \eqn{A_{ji}} is negative binomial with mean \eqn{m_{ji}} and variance
#' \eqn{m_{ji}(1 + \ell_{ji}/\beta_j)}, and a point is drawn in red when its
#' count is more than \code{num_sd} standard deviations from its fit.
#'
#' With \code{bool_draw_bars = TRUE} each point carries a vertical bar, the
#' interval \eqn{m_{ji} \pm} \code{num_sd} standard deviations (truncated at
#' zero), drawn along the y-axis at the observed count. A point is red exactly
#' when its bar does not reach the line \eqn{y = x}.
#'
#' What to look for. The points should follow the line, and only a small
#' share should be red. With \code{num_sd = 3}, counts simulated from the
#' model and drawn with their true rates gave a share of 1 to 2 percent, more
#' than the 0.3 percent of a Gaussian because a count is skewed to the right.
#' A rate ten times too large gave about 5 percent, and a rate ten times too
#' small under 0.1 percent. So a larger share, or bars that are
#' systematically too short to reach the line, means the over-dispersion is
#' \emph{under}estimated (the rate is too large), and bars that reach far past
#' the line with almost no red point mean it is \emph{over}estimated. Points
#' that leave the line together, on one side, are a problem of the fit and
#' not of the over-dispersion.
#'
#' The rate is \code{nuisance_vec}, the rate after the cap and before the
#' adjustments \code{compute_posterior} makes to it, and the library size is
#' the one \code{estimate_nuisance} used.
#'
#' @param input_obj       \code{eSVD} object after \code{estimate_nuisance}.
#' @param bool_draw_bars  If \code{TRUE}, draw each point's bar. Default
#'                        \code{FALSE}.
#' @param genes           \code{NULL}, for all the genes, or a character
#'                        vector of genes, each drawn in its own panel with
#'                        all of its cells.
#' @param max_points      When \code{genes} is \code{NULL}, the largest number
#'                        of points to draw. If there are more pairs of a cell
#'                        and a gene than this, that many are sampled; the
#'                        fitted count is computed for the sampled pairs only,
#'                        so no cell-by-gene matrix is allocated.
#' @param num_sd          One positive number, the number of standard
#'                        deviations beyond which a point is red.
#' @param seurat_obj      \code{NULL}, or the \code{Seurat} object the analysis
#'                        was run on. Required when \code{input_obj} has no
#'                        counts (it was built with \code{bool_diet = TRUE});
#'                        see \code{recompute_pvalue}.
#' @param seed_number     Seed for the sample of pairs; \code{NULL} leaves the
#'                        random number stream alone.
#'
#' @returns A \code{ggplot} object. Its \code{data} has one row per point and
#' columns \code{cell}, \code{gene}, \code{observed}, \code{fitted},
#' \code{sd}, \code{bool_outside} (is the point red?), \code{status} (the
#' same, as the factor the legend is drawn from), and the coordinates
#' drawn: \code{log_observed}, \code{log_fitted}, \code{log_lower} and
#' \code{log_upper} (the ends of the bar). The share of red points is
#' \code{mean(plot$data$bool_outside)}.
#' @examples
#' set.seed(10)
#' sim <- generate_null(cell_per_person = 15, num_genes = 40,
#'                      num_individuals = 8)
#' esvd_obj <- initialize_esvd(dat = sim$obs_mat,
#'                             covariates = sim$covariates,
#'                             metadata_individual = sim$metadata_individual,
#'                             case_control_variable = "CC",
#'                             bool_intercept = TRUE,
#'                             k = 2,
#'                             lambda = 0.1)
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_Init",
#'                                                omitted_variables = "Log_UMI")
#' esvd_obj <- opt_esvd(input_obj = esvd_obj,
#'                      max_iter = 5,
#'                      offset_variables = setdiff(colnames(esvd_obj$covariates), "CC"),
#'                      fit_name = "fit_First",
#'                      fit_previous = "fit_Init")
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_First",
#'                                                omitted_variables = "Log_UMI")
#' esvd_obj <- estimate_nuisance(input_obj = esvd_obj)
#' if(requireNamespace("ggplot2", quietly = TRUE)){
#'   plot1 <- plot_fitted_vs_observed(input_obj = esvd_obj,
#'                                    genes = colnames(sim$obs_mat)[1:2],
#'                                    bool_draw_bars = TRUE)
#'   plot(plot1)
#'   mean(plot1$data$bool_outside)
#' }
#' @export
plot_fitted_vs_observed <- function(input_obj,
                                    bool_draw_bars = FALSE,
                                    genes = NULL,
                                    max_points = 50000,
                                    num_sd = 3,
                                    seurat_obj = NULL,
                                    seed_number = 10){
  .check_plot_packages("ggplot2", "plot_fitted_vs_observed")
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            is.logical(bool_draw_bars), length(bool_draw_bars) == 1,
            is.null(genes) || is.character(genes))
  if(!is.numeric(num_sd) || length(num_sd) != 1 || is.na(num_sd) || num_sd <= 0){
    stop("`num_sd` must be one positive number")
  }
  if(!is.numeric(max_points) || length(max_points) != 1 || is.na(max_points) ||
     max_points < 1){
    stop("`max_points` must be one number that is at least 1")
  }
  .check_fit_is_named(input_obj)
  if(is.null(input_obj[[input_obj[["latest_Fit"]]]]$nuisance_vec)){
    stop("`input_obj$", input_obj[["latest_Fit"]], "` has no `nuisance_vec`; ",
         "run `estimate_nuisance` first")
  }

  input_list <- .get_counts_and_covariates(input_obj = input_obj,
                                           seurat_obj = seurat_obj)
  dat <- input_list$dat
  if(!is.null(genes)){
    missing_genes <- setdiff(genes, colnames(dat))
    if(length(missing_genes) > 0){
      stop("gene(s) `", paste0(missing_genes, collapse = "`, `"),
           "` in `genes` are not among the genes that were analyzed")
    }
  }

  pair_list <- .choose_pairs(n = nrow(dat),
                             gene_idx_vec = if(is.null(genes)) seq_len(ncol(dat)) else match(unique(genes), colnames(dat)),
                             bool_sample = is.null(genes),
                             max_points = max_points,
                             seed_number = seed_number)
  plot_df <- .form_fitted_df(input_obj = input_obj,
                             covariates = input_list$covariates,
                             dat = dat,
                             cell_idx_vec = pair_list$cell_idx_vec,
                             gene_idx_vec = pair_list$gene_idx_vec,
                             num_sd = num_sd)
  # The red points are drawn last, on top.
  plot_df <- plot_df[order(plot_df$bool_outside), , drop = FALSE]
  rownames(plot_df) <- NULL

  num_outside <- sum(plot_df$bool_outside)
  subtitle_val <- paste0(
    num_outside, " of ", nrow(plot_df), " points (",
    signif(100 * num_outside / nrow(plot_df), 2), "%) are more than ",
    num_sd, " SD from the fit",
    if(pair_list$bool_sampled) paste0("; a sample of the ",
                                     pair_list$num_total, " pairs")
  )
  color_vec <- c("FALSE" = "#52514e", "TRUE" = "#e34948")
  label_vec <- c("FALSE" = paste0("within ", num_sd, " SD"),
                 "TRUE" = paste0("beyond ", num_sd, " SD"))
  plot_df$status <- factor(as.character(plot_df$bool_outside),
                           levels = c("FALSE", "TRUE"))

  plot1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$log_observed,
                                                 y = .data$log_fitted))
  if(bool_draw_bars){
    plot1 <- plot1 + ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$log_lower, ymax = .data$log_upper,
                   color = .data$status),
      alpha = 0.25, linewidth = 0.3, show.legend = FALSE
    )
  }
  plot1 <- plot1 + ggplot2::geom_abline(slope = 1, intercept = 0,
                                        linetype = "dashed",
                                        color = "#0b0b0b")
  plot1 <- plot1 + ggplot2::geom_point(
    ggplot2::aes(color = .data$status, shape = .data$status,
                 alpha = .data$status),
    size = 1.2, show.legend = TRUE
  )
  plot1 <- plot1 + ggplot2::scale_color_manual(values = color_vec,
                                               labels = label_vec,
                                               drop = FALSE, name = NULL)
  # The many points within the interval are drawn faintly, the few beyond
  # it solid.
  plot1 <- plot1 + ggplot2::scale_alpha_manual(values = c("FALSE" = 0.4,
                                                          "TRUE" = 1),
                                               labels = label_vec,
                                               drop = FALSE, name = NULL)
  plot1 <- plot1 + ggplot2::scale_shape_manual(values = c("FALSE" = 16,
                                                          "TRUE" = 17),
                                               labels = label_vec,
                                               drop = FALSE, name = NULL)
  if(!is.null(genes)){
    plot1 <- plot1 + ggplot2::facet_wrap(ggplot2::vars(.data$gene))
  }
  # The same range on both axes, so that the line y = x is the diagonal of
  # the panel and a departure from it reads the same in either direction.
  limit_vec <- c(0, max(plot_df$log_observed, plot_df$log_fitted,
                        if(bool_draw_bars) plot_df$log_upper))
  plot1 <- plot1 + ggplot2::coord_cartesian(xlim = limit_vec,
                                            ylim = limit_vec)
  plot1 <- plot1 + ggplot2::labs(
    x = "Observed count, log(1 + count)",
    y = "Fitted count, log(1 + fit)",
    title = "Fitted against observed counts",
    subtitle = subtitle_val
  )
  plot1 <- plot1 + ggplot2::theme_bw()

  plot1
}

# ---- helpers ----------------------------------------------------------------

.check_plot_packages <- function(package_vec, function_name){
  for(package_name in package_vec){
    if(!requireNamespace(package_name, quietly = TRUE)){
      stop("package `", package_name, "` is required by `", function_name,
           "()`; install it with install.packages(\"", package_name, "\")")
    }
  }

  invisible()
}

# Neutral gray for the ordinary genes and three hues for the rest. The hues
# are separated for color-vision deficiency (checked for all pairs); the
# shape repeats the status, so it is never carried by color alone.
.nuisance_status_colors <- function(){
  c(estimated = "#8a8985", capped = "#2a78d6", boundary = "#eb6834",
    failed = "#1baf7a")
}

.nuisance_status_shapes <- function(){
  c(estimated = 16, capped = 17, boundary = 15, failed = 18)
}

#' The data frame that \code{plot_nuisance} draws
#'
#' @inheritParams plot_nuisance
#'
#' @returns Data frame with one row per analyzed gene, in the order of the
#' fit, and columns \code{gene}, \code{x_value}, \code{unit_free_rate} and
#' \code{nuisance_status}.
#' @noRd
.form_nuisance_df <- function(input_obj, x_axis, bool_uncapped){
  latest_Fit <- input_obj[["latest_Fit"]]
  fit <- input_obj[[latest_Fit]]

  rate_name <- if(bool_uncapped) "nuisance_mle_vec" else "nuisance_vec"
  x_name <- c(mean_expression = "gene_mean_count_vec",
              log10pvalue = NA,
              sparsity = "gene_sparsity_vec")[[x_axis]]
  needed_vec <- stats::na.omit(c(rate_name, "nuisance_library_median_vec",
                                 "nuisance_status", x_name))
  missing_elements <- setdiff(needed_vec, names(fit))
  if(length(missing_elements) > 0){
    stop("`input_obj$", latest_Fit, "` has no `",
         paste0(missing_elements, collapse = "`, `"),
         "`; run `estimate_nuisance` of eSVD2 1.2.0 or later on the object")
  }

  gene_vec <- .get_analyzed_genes(input_obj)
  if(x_axis == "log10pvalue"){
    if(is.null(input_obj[["pvalue_list"]])){
      stop("`input_obj` has no `pvalue_list`, which `x_axis = ",
           "\"log10pvalue\"` draws; run the test first, or choose another ",
           "`x_axis`")
    }
    x_vec <- input_obj[["pvalue_list"]]$log10pvalue
  } else {
    x_vec <- fit[[x_name]]
  }
  stopifnot(all(gene_vec %in% names(x_vec)),
            all(gene_vec %in% names(fit[[rate_name]])))

  data.frame(
    gene = gene_vec,
    x_value = as.numeric(x_vec[gene_vec]),
    unit_free_rate = as.numeric(fit[[rate_name]][gene_vec] /
                                  fit$nuisance_library_median_vec[gene_vec]),
    nuisance_status = factor(as.character(fit$nuisance_status[gene_vec]),
                             levels = .nuisance_status_levels()),
    stringsAsFactors = FALSE
  )
}

#' Which genes \code{plot_nuisance} labels
#'
#' @inheritParams plot_nuisance
#' @param plot_df  Output of \code{.form_nuisance_df}.
#'
#' @returns Logical vector, one entry per row of \code{plot_df}.
#' @noRd
.choose_labeled_genes <- function(input_obj, plot_df, genes, num_label){
  if(!is.null(genes)){
    missing_genes <- setdiff(genes, plot_df$gene)
    if(length(missing_genes) > 0){
      stop("gene(s) `", paste0(missing_genes, collapse = "`, `"),
           "` in `genes` are not among the genes that were analyzed")
    }
    return(plot_df$gene %in% genes)
  }

  flagged_idx <- which(as.character(plot_df$nuisance_status) != "estimated")
  num_label <- min(floor(num_label), length(flagged_idx))
  if(num_label == 0) return(rep(FALSE, nrow(plot_df)))

  # The flagged genes most likely to be reported: the most significant. An
  # object that has not been tested is ordered by how far the rate ran.
  if(!is.null(input_obj[["pvalue_list"]])){
    order_vec <- input_obj[["pvalue_list"]]$log10pvalue[plot_df$gene]
  } else {
    order_vec <- input_obj[[input_obj[["latest_Fit"]]]]$nuisance_mle_vec[plot_df$gene]
  }
  order_vec <- as.numeric(order_vec)[flagged_idx]
  label_idx <- flagged_idx[order(order_vec, decreasing = TRUE)[seq_len(num_label)]]

  seq_len(nrow(plot_df)) %in% label_idx
}

#' Choose the pairs of a cell and a gene to draw
#'
#' @param n             Number of cells.
#' @param gene_idx_vec  Indices of the genes to draw from.
#' @param bool_sample   Whether to sample when there are more than
#'                      \code{max_points} pairs.
#' @param max_points    Largest number of pairs when sampling.
#' @param seed_number   Seed, or \code{NULL}.
#'
#' @returns List with \code{bool_sampled}, \code{cell_idx_vec},
#' \code{gene_idx_vec} (the two of the same length, one entry per pair) and
#' \code{num_total}.
#' @noRd
.choose_pairs <- function(n, gene_idx_vec, bool_sample, max_points,
                          seed_number){
  # As doubles: the number of pairs of a large data set exceeds the largest
  # integer.
  num_genes <- length(gene_idx_vec)
  num_total <- as.numeric(n) * as.numeric(num_genes)
  bool_sampled <- bool_sample && num_total > max_points

  if(bool_sampled){
    if(!is.null(seed_number)) set.seed(seed_number)
    # `sample.int()` and not `sample()`, which on one number x samples 1:x.
    pair_idx_vec <- sort(sample.int(num_total, size = floor(max_points)))
  } else {
    pair_idx_vec <- seq_len(num_total)
  }

  # Pairs are numbered down the cells of the first gene, then the second.
  list(bool_sampled = bool_sampled,
       cell_idx_vec = ((pair_idx_vec - 1) %% n) + 1,
       gene_idx_vec = gene_idx_vec[((pair_idx_vec - 1) %/% n) + 1],
       num_total = num_total)
}

#' The data frame that \code{plot_fitted_vs_observed} draws
#'
#' The fitted count and the library size are computed for the given pairs
#' only. The library is the one \code{estimate_nuisance.eSVD} used, from the
#' settings it recorded in \code{param}.
#'
#' @param input_obj     \code{eSVD} object.
#' @param covariates    Covariate matrix of the fit, rows are cells.
#' @param dat           Count matrix, rows are cells and columns are the
#'                      analyzed genes.
#' @param cell_idx_vec  Row of \code{dat} of each pair.
#' @param gene_idx_vec  Column of \code{dat} of each pair.
#' @param num_sd        Number of standard deviations.
#'
#' @returns Data frame with one row per pair; see
#' \code{plot_fitted_vs_observed}.
#' @noRd
.form_fitted_df <- function(input_obj, covariates, dat, cell_idx_vec,
                            gene_idx_vec, num_sd){
  fit <- input_obj[[input_obj[["latest_Fit"]]]]
  gene_vec <- colnames(dat)
  stopifnot(identical(rownames(fit$x_mat), rownames(dat)),
            all(gene_vec %in% rownames(fit$y_mat)),
            identical(colnames(covariates), colnames(fit$z_mat)))
  y_mat <- fit$y_mat[gene_vec, , drop = FALSE]
  z_mat <- fit$z_mat[gene_vec, , drop = FALSE]
  nuisance_vec <- fit$nuisance_vec[gene_vec]

  library_idx <- .nuisance_library_idx(
    covariates = covariates,
    case_control_variable = .get_param_or_default(input_obj, "init_case_control_variable", NULL),
    library_size_variable = input_obj$param[["init_library_size_variable"]],
    bool_covariates_as_library = input_obj$param[["nuisance_bool_covariates_as_library"]],
    bool_library_includes_interept = input_obj$param[["nuisance_bool_library_includes_interept"]]
  )

  nat_vec <- rowSums(fit$x_mat[cell_idx_vec, , drop = FALSE] *
                       y_mat[gene_idx_vec, , drop = FALSE]) +
    rowSums(covariates[cell_idx_vec, , drop = FALSE] *
              z_mat[gene_idx_vec, , drop = FALSE])
  library_vec <- exp(rowSums(covariates[cell_idx_vec, library_idx, drop = FALSE] *
                               z_mat[gene_idx_vec, library_idx, drop = FALSE]))
  fitted_vec <- exp(nat_vec)
  sd_vec <- sqrt(fitted_vec * (1 + library_vec / nuisance_vec[gene_idx_vec]))
  observed_vec <- as.numeric(dat[cbind(cell_idx_vec, gene_idx_vec)])

  data.frame(
    cell = rownames(dat)[cell_idx_vec],
    gene = factor(gene_vec[gene_idx_vec],
                  levels = gene_vec[unique(gene_idx_vec)]),
    observed = observed_vec,
    fitted = as.numeric(fitted_vec),
    sd = as.numeric(sd_vec),
    bool_outside = as.logical(abs(observed_vec - fitted_vec) > num_sd * sd_vec),
    log_observed = log1p(observed_vec),
    log_fitted = log1p(as.numeric(fitted_vec)),
    log_lower = log1p(pmax(as.numeric(fitted_vec - num_sd * sd_vec), 0)),
    log_upper = log1p(as.numeric(fitted_vec + num_sd * sd_vec)),
    stringsAsFactors = FALSE
  )
}
