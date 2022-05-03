# this file is for more experimental plotting functions


plot_volcano <- function(mat, pval_vec, de_idx, bool_iqr = T,
                         xlab = ifelse(bool_iqr, "IQR (signed by selection)",
                                       "Std (signed by selection)"),
                         ylab = "-log10(pval)",
                         xlim = NA,
                         xgrid = 5, ygrid = 5, cex = 1,
                         color_palette = rev(colorspace::sequential_hcl(10, palette = "Terrain 2")),
                         col_range = c(0,1),
                         ...){
  log_vec <- -log10(pval_vec)
  tmp <- max(log_vec[!is.infinite(log_vec)])
  log_vec[is.infinite(log_vec)] <- tmp

  if(bool_iqr){
    sd_vec <- matrixStats::rowDiffs(matrixStats::colQuantiles(mat, probs = c(0.25, 0.75)))
    mean_vec <- matrixStats::rowMedians(mat)
  } else {
    sd_vec <- matrixStats::colSds(mat)
    mean_vec <- matrixStats::rowMeans2(mat)
  }

  if(all(is.na(xlim))) xlim <- c(-1,1)*max(sd_vec)
  sd_vec[de_idx] <- -sd_vec[de_idx]
  ylim <- range(c(0, log_vec))

  # assign colors
  tmp <- seq(col_range[1], col_range[2], length.out = length(color_palette))
  col_vec <- sapply(mean_vec, function(x){
    color_palette[which.min(abs(x-tmp))]
  })

  graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)

  # draw grid
  for(i in seq(xlim[1], xlim[2], length.out = xgrid)){
    graphics::lines(rep(i, 2), c(-10,10)*max(abs(ylim)), lty = 3,
                    lwd = 0.5, col = "gray")
  }
  for(i in seq(ylim[1], ylim[2], length.out = ygrid)){
    graphics::lines(c(-10,10)*max(abs(xlim)), rep(i, 2), lty = 3,
                    lwd = 0.5, col = "gray")
  }

  graphics::points(x = sd_vec, y = log_vec, pch = 21, cex = cex,
                   bg = col_vec, col = "black")

  invisible()
}

plot_sd_scatter <- function(mat, membership_vec, de_idx,
                            bool_iqr = T,
                            xlab = ifelse(bool_iqr, "Between-celltype IQR", "Between-celltype Median"),
                            ylab = ifelse(bool_iqr, "Within-celltype IQR", "Within-celltype Median"),
                            xlim = NA,
                            gridsize = 5,
                            ...){
  uniq_val <- sort(unique(membership_vec))

  mean_mat <- sapply(uniq_val, function(celltype){
    idx <- which(membership_vec == celltype)
    if(bool_iqr) {
      matrixStats::colMedians(mat[idx,,drop = F])
    } else {
      matrixStats::colMeans2(mat[idx,,drop = F])
    }
  })

  if(bool_iqr){
    between_sd_vec <- matrixStats::rowDiffs(matrixStats::rowQuantiles(mean_mat, probs = c(0.25, 0.75)))
  } else {
    between_sd_vec <- matrixStats::rowSds(mean_mat)
  }

  sd_mat <- sapply(uniq_val, function(celltype){
    idx <- which(membership_vec == celltype)
    if(bool_iqr){
      matrixStats::rowDiffs(matrixStats::colQuantiles(mat[idx,,drop = F], probs = c(0.25, 0.75)))
    } else {
      matrixStats::colSds(mat[idx,,drop = F])
    }
  })
  if(bool_iqr){
    within_sd_vec <- matrixStats::rowMedians(sd_mat)
  } else {
    within_sd_vec <- matrixStats::rowSds(sd_mat)
  }

  if(all(is.na(xlim))) xlim <- range(c(between_sd_vec, within_sd_vec))
  graphics::plot(NA, xlim = xlim, ylim = xlim, xlab = xlab, ylab = ylab, ...)

  # draw grid
  for(i in seq(xlim[1], xlim[2], length.out = gridsize)){
    graphics::lines(rep(i, 2), c(-10,10)*max(abs(xlim)), lty = 3,
                    lwd = 0.5, col = "gray")
  }
  for(i in seq(xlim[1], xlim[2], length.out = gridsize)){
    graphics::lines(c(-10,10)*max(abs(xlim)), rep(i, 2), lty = 3,
                    lwd = 0.5, col = "gray")
  }

  col_vec <- rep("black", length(between_sd_vec))
  col_vec[de_idx] <- "red"
  graphics::points(x = between_sd_vec, y = within_sd_vec,
                   pch = 21, bg = col_vec, col = "black")

  invisible()
}
