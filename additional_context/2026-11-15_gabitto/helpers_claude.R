# Gabitto split-half study: the replication statistics and the two figures
# Drafted by Claude for Kevin Z. Lin, 2026-11-15

# Sourced by step 4 and by the writeup, so that a number in a table and the
# same number in a figure title come from one function. Defines functions
# only.

# Tests ------------------------------------------------------------------------

# Two-sided Wald p-value of a logFC and its SE against the standard normal.
# Every rule is tested the same way so that the rules' own nominations
# differ only through the pair they report.
.wald_pvalue <- function(lfc_vec, se_vec){
  2 * stats::pnorm(-abs(lfc_vec / se_vec))
}

# Statistics of one gene set ---------------------------------------------------

# The Gaussian negative log-likelihood of an observed difference (nats).
.nll_vec <- function(diff_vec, var_vec){
  0.5 * log(2 * pi * var_vec) + diff_vec^2 / (2 * var_vec)
}

# Agreement between the two halves over the genes `idx`. The replicate
# score is the autoresearch score with the unknown truth replaced by the
# other half: if both halves estimate the same effect and both SEs are
# honest, lfc1 - lfc2 is N(0, se1^2 + se2^2).
.agreement_row <- function(pair_df, idx, min_genes = 5){
  sub_df <- pair_df[idx, , drop = FALSE]
  sub_df <- sub_df[is.finite(sub_df$lfc1) & is.finite(sub_df$lfc2) &
                     is.finite(sub_df$se1) & is.finite(sub_df$se2) &
                     sub_df$se1 > 0 & sub_df$se2 > 0, , drop = FALSE]
  num_genes <- nrow(sub_df)
  if(num_genes < min_genes){
    return(data.frame(num_genes = num_genes, score = NA_real_,
                      z_rms = NA_real_, coverage = NA_real_,
                      pearson = NA_real_, spearman = NA_real_,
                      sign_agree = NA_real_, rmsd = NA_real_,
                      sd_lfc = NA_real_))
  }
  diff_vec <- sub_df$lfc1 - sub_df$lfc2
  var_vec <- sub_df$se1^2 + sub_df$se2^2
  data.frame(num_genes = num_genes,
             score = mean(.nll_vec(diff_vec, var_vec)),
             # 1 when the SEs explain the disagreement between the halves;
             # above 1 when they are too small.
             z_rms = sqrt(mean(diff_vec^2 / var_vec)),
             coverage = mean(abs(diff_vec) < stats::qnorm(0.975) *
                               sqrt(var_vec)),
             pearson = stats::cor(sub_df$lfc1, sub_df$lfc2),
             spearman = stats::cor(sub_df$lfc1, sub_df$lfc2,
                                   method = "spearman"),
             sign_agree = mean(sign(sub_df$lfc1) == sign(sub_df$lfc2)),
             rmsd = sqrt(mean(diff_vec^2)),
             sd_lfc = stats::sd(c(sub_df$lfc1, sub_df$lfc2)))
}

# Discovery in one half, replication in the other, pooled over the two
# directions. `nominated1` and `nominated2` flag the genes each half
# nominates. The slope is the least-squares line through the origin of the
# replication logFC on the discovery logFC (Writeup49's convention): 1
# when the discovered effects reappear at full size, below 1 under the
# winner's curse.
.discovery_row <- function(pair_df, nominated1, nominated2, min_genes = 5){
  bool_ok <- is.finite(pair_df$lfc1) & is.finite(pair_df$lfc2) &
    is.finite(pair_df$se1) & is.finite(pair_df$se2)
  nominated1 <- nominated1 & bool_ok
  nominated2 <- nominated2 & bool_ok
  disc_vec <- c(pair_df$lfc1[nominated1], pair_df$lfc2[nominated2])
  rep_vec <- c(pair_df$lfc2[nominated1], pair_df$lfc1[nominated2])
  rep_p_vec <- c(pair_df$wald_p2[nominated1], pair_df$wald_p1[nominated2])
  num_disc <- length(disc_vec)
  if(num_disc < min_genes){
    return(data.frame(num_half1 = sum(nominated1), num_half2 = sum(nominated2),
                      num_both = sum(nominated1 & nominated2),
                      sign_rep = NA_real_, slope_rep = NA_real_,
                      rep_p05 = NA_real_))
  }
  bool_sign <- sign(disc_vec) == sign(rep_vec)
  data.frame(num_half1 = sum(nominated1), num_half2 = sum(nominated2),
             num_both = sum(nominated1 & nominated2),
             sign_rep = mean(bool_sign),
             slope_rep = sum(disc_vec * rep_vec) / sum(disc_vec^2),
             # Same sign and nominal Wald p < 0.05 in the replication half.
             rep_p05 = mean(bool_sign & rep_p_vec < 0.05))
}

# Figures ----------------------------------------------------------------------

.half_palette <- function(){
  c("neither" = "grey80", "half 1 only" = "#0072B2", "half 2 only" = "#E69F00",
    "both halves" = "#D55E00")
}

# Which half nominates each gene, as an ordered factor.
.nomination_factor <- function(nominated1, nominated2){
  label_vec <- rep("neither", length(nominated1))
  label_vec[nominated1 & !nominated2] <- "half 1 only"
  label_vec[!nominated1 & nominated2] <- "half 2 only"
  label_vec[nominated1 & nominated2] <- "both halves"
  factor(label_vec, levels = names(.half_palette()))
}

# Half 2's logFC against half 1's for one rule, one panel per split, the
# nominated genes coloured by the half that nominates them. `plot_df` has
# the columns split, lfc1, lfc2 and nomination (see `.nomination_factor`);
# `label_df` has one row per split with the text printed in the panel.
.plot_replication <- function(plot_df, label_df, title, bool_only_nominated,
                              quantile_limit = 0.999){
  plot_df <- plot_df[is.finite(plot_df$lfc1) & is.finite(plot_df$lfc2), ]
  if(bool_only_nominated){
    plot_df <- plot_df[plot_df$nomination != "neither", ]
  }
  # A few genes with a logFC in the tens would flatten every other gene;
  # the axes stop at a high quantile and those genes are drawn on the edge.
  limit <- stats::quantile(abs(c(plot_df$lfc1, plot_df$lfc2)), quantile_limit,
                           na.rm = TRUE)
  limit <- max(limit, 0.1)
  plot_df$lfc1 <- pmax(pmin(plot_df$lfc1, limit), -limit)
  plot_df$lfc2 <- pmax(pmin(plot_df$lfc2, limit), -limit)
  plot_df <- plot_df[order(plot_df$nomination), ]
  label_df$x <- -limit
  label_df$y <- limit

  plot1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = lfc1, y = lfc2,
                                                 color = nomination))
  plot1 <- plot1 + ggplot2::geom_hline(yintercept = 0, color = "grey60",
                                       linewidth = 0.3)
  plot1 <- plot1 + ggplot2::geom_vline(xintercept = 0, color = "grey60",
                                       linewidth = 0.3)
  plot1 <- plot1 + ggplot2::geom_abline(slope = 1, intercept = 0,
                                        color = "grey40", linetype = "dashed",
                                        linewidth = 0.3)
  plot1 <- plot1 + ggplot2::geom_point(size = 0.7, alpha = 0.7)
  plot1 <- plot1 + ggplot2::geom_text(data = label_df,
                                      ggplot2::aes(x = x, y = y, label = label),
                                      inherit.aes = FALSE, hjust = 0, vjust = 1,
                                      size = 2.9)
  plot1 <- plot1 + ggplot2::scale_color_manual(values = .half_palette(),
                                               drop = FALSE,
                                               name = "Nominated by")
  plot1 <- plot1 + ggplot2::facet_wrap(~ split)
  plot1 <- plot1 + ggplot2::coord_equal(xlim = c(-limit, limit),
                                        ylim = c(-limit, limit))
  plot1 <- plot1 + ggplot2::labs(x = "log2 fold change, half 1",
                                 y = "log2 fold change, half 2",
                                 title = title)
  plot1 <- plot1 + ggplot2::theme_bw(base_size = 10)
  plot1 <- plot1 + ggplot2::theme(legend.position = "bottom",
                                  plot.title = ggplot2::element_text(size = 10))
  plot1
}

# One statistic for every method, a point per split.
.plot_summary <- function(summary_df, column, xlab, method_vec,
                          reference = NULL){
  plot_df <- summary_df[, c("method", "split", column)]
  colnames(plot_df)[3] <- "value"
  plot_df$method <- factor(plot_df$method, levels = rev(method_vec))
  plot1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = value, y = method,
                                                 color = split, shape = split))
  if(!is.null(reference)){
    plot1 <- plot1 + ggplot2::geom_vline(xintercept = reference,
                                         color = "grey50", linetype = "dashed")
  }
  plot1 <- plot1 + ggplot2::geom_point(size = 2.2)
  plot1 <- plot1 + ggplot2::scale_color_manual(values = c(split1 = "#0072B2",
                                                          split2 = "#D55E00"))
  plot1 <- plot1 + ggplot2::labs(x = xlab, y = NULL)
  plot1 <- plot1 + ggplot2::theme_bw(base_size = 10)
  plot1 <- plot1 + ggplot2::theme(legend.position = "bottom")
  plot1
}
