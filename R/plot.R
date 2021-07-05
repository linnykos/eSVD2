plot_scores_heatmap <- function(score_mat, membership_vec = NA, num_col = 10,
                                bool_center = T, bool_scale = T,
                                bool_log = F, scaling_power = 1, luminosity = F,
                                ...){

  n <- nrow(score_mat)
  if(bool_center | bool_scale){
    score_mat <- scale(score_mat, center = bool_center, scale = bool_scale)
  }
  if(bool_log){
    score_mat <- log(abs(score_mat)+1)*sign(score_mat)
  }
  zlim <- range(score_mat)

  # construct colors. green is negative
  max_val <- max(abs(zlim)); min_val <- max(min(abs(zlim)), 1e-3)
  col_vec_neg <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(1,1,1), num_col,
                                   luminosity = luminosity)
  break_vec_neg <- -max_val*seq(1, 0, length.out = num_col+2)^scaling_power
  break_vec_neg <- break_vec_neg[-length(break_vec_neg)]

  # red is positive
  col_vec_pos <- .colorRamp_custom(c(1,1,1), c(0.803, 0.156, 0.211), num_col,
                                   luminosity = luminosity)
  break_vec_pos <- max_val*seq(0, 1, length.out = num_col+2)^scaling_power
  break_vec_pos <- break_vec_pos[-1]

  # combine the two
  break_vec <- c(break_vec_neg, break_vec_pos)
  col_vec <- c(col_vec_neg, "white", col_vec_pos)

  if(!all(is.na(membership_vec))){
    stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(score_mat))

    membership_vec <- as.numeric(membership_vec) ## convert to integers
    idx <- order(membership_vec, decreasing = F)
    breakpoints <- 1-which(abs(diff(sort(membership_vec, decreasing = F))) >= 1e-6)/n
  } else {
    idx <- 1:n
  }

  line_func <- function(){
    if(!all(is.na(membership_vec))){
      for(i in 1:length(breakpoints)){
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2.1, col = "white")
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2, lty = 2)
      }
    }
  }

  graphics::image(.rotate(score_mat[idx,,drop = F]), col = col_vec,
                  breaks = break_vec, ...)
  line_func()

  invisible()
}

######################

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
