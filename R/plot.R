plot_scores_heatmap <- function(score_mat, membership_vec = NA,
                                major_breakpoint = NA,
                                num_col = 10,
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
        if(!all(is.na(major_breakpoint)) && i %in% major_breakpoint){
          graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2.2, lty = 1)
        } else {
          graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2, lty = 2)
        }
      }
    }
  }

  graphics::image(.rotate(score_mat[idx,,drop = F]), col = col_vec,
                  breaks = break_vec, ...)
  line_func()

  invisible()
}

#######################

plot_scatterplot_poisson <- function(mat,
                                     mean_mat,
                                     xlim = NA,
                                     quantile_vec = c(0.1, 0.9),
                                     only_nonzero = T,
                                     max_num = 1e5,
                                     point_cex = 1,
                                     point_col = grDevices::rgb(0,0,0,0.1),
                                     mean_lwd = 2,
                                     polygon_density = 30,
                                     asp = T, ...){
  if(only_nonzero) {
    idx <- which(mat != 0)
  } else {
    idx <- 1:prod(dim(mat))
  }
  if(length(idx) > max_num){
    idx <- sample(idx, size = max_num)
  }

  if(all(is.na(xlim))) xlim <- range(c(0, mat[idx], mean_mat[idx]))
  mean_vec <- seq(0, xlim[2], length.out = 100)
  lower_vec <- stats::qpois(quantile_vec[1], lambda = mean_vec)
  upper_vec <- stats::qpois(quantile_vec[2], lambda = mean_vec)

  graphics::plot(NA,
                 xlim = xlim, ylim = xlim,
                 asp = asp,
                 ...)
  graphics::polygon(x = c(mean_vec, rev(mean_vec)),
                    y = c(lower_vec, rev(upper_vec)),
                    col = grDevices::rgb(1,0,0,0.2),
                    border = NA,
                    density = polygon_density,
                    angle = -45)
  graphics::lines(x = mean_vec, y = mean_vec,
                  col = "red",
                  lwd = mean_lwd)
  graphics::points(x = mean_mat[idx],
                   y = mat[idx],
                   pch = 16,
                   col = point_col,
                   cex = point_cex)

  invisible()
}

plot_scatterplot_nb <- function(mat,
                                mean_mat,
                                size_vec,
                                main = "",
                                include_percentage_in_main = T,
                                xlim = NA,
                                quantile_shoulder = 0.9,
                                only_nonzero = T,
                                max_num = 1e5,
                                included_col = "black",
                                excluded_col = "red",
                                cex = 1,
                                asp = T,
                                verbose = T,
                                ...){
  stopifnot(length(size_vec) == ncol(mat),
            all(dim(mat) == dim(mean_mat)))
  if(only_nonzero) {
    idx <- which(mat != 0)
  } else {
    idx <- 1:prod(dim(mat))
  }
  if(length(idx) > max_num){
    idx <- sample(idx, size = max_num)
  }

  size_mat <- matrix(size_vec,
                     nrow = nrow(mat),
                     ncol = ncol(mat),
                     byrow = T)

  x_vec <- mean_mat[idx]
  y_vec <- mat[idx]
  if(all(is.na(xlim))) xlim <- range(c(0, x_vec, y_vec))

  tmp_mat <- cbind(mat[idx], mean_mat[idx])
  bool_vec <- apply(tmp_mat, 1, function(x){all(x >= xlim[1]) & all(x <= xlim[2])})
  tmp_mat <- tmp_mat[which(bool_vec),]
  angle_val <- .compute_principal_angle(tmp_mat)
  print(angle_val)

  if(verbose) print("Starting to compute quantiles")
  tabulated_mat <- sapply(1:length(idx), function(counter){
    if(verbose & length(idx) > 10 & counter %% floor(length(idx)/10) == 0) cat('*')
    i <- idx[counter]

    quant <- stats::pnbinom(mean_mat[i],
                            size = size_mat[i],
                            mu = mean_mat[i])
    lower_quant <- (1-quantile_shoulder)*quant
    upper_quant <- quant + (1-quant)*quantile_shoulder

    lower_val <- stats::qnbinom(lower_quant,
                                size = size_mat[i],
                                mu = mean_mat[i])
    upper_val <- stats::qnbinom(upper_quant,
                                size = size_mat[i],
                                mu = mean_mat[i])

    bool <- mat[i] >= lower_val & mat[i] <= upper_val
    width <- upper_quant-lower_quant

    c(bool = as.numeric(bool), width = width)
  })

  col_vec <- c(excluded_col, included_col)[tabulated_mat["bool",]+1]
  observed_percentage <- round(length(which(col_vec == included_col))/length(col_vec), 2)
  expected_percentage <- round(mean(tabulated_mat["width",]), 2)

  if(include_percentage_in_main){
    main_modified = paste0(main, " (", 100*observed_percentage,  "% of ", 100*expected_percentage, "%)")
  } else {
    main_modified = main
  }
  graphics::plot(NA,
                 xlim = xlim,
                 ylim = xlim,
                 asp = asp,
                 main = main_modified,
                 ...)

  # plot diagonal
  seq_vec <- c(0, 2*max(abs(xlim)))
  graphics::lines(x = seq_vec,
                  y = seq_vec,
                  col = "black",
                  lwd = 2,
                  lty = 2)
  graphics::lines(x = seq_vec,
                  y = seq_vec * tan(angle_val*pi/180),
                  col = "green",
                  lwd = 2,
                  lty = 2)

  graphics::points(x = x_vec,
                   y = y_vec,
                   pch = 16,
                   cex = cex,
                   col = col_vec)

}

#############
#' Compute principal angle
#'
#' @param tmp_mat a matrix with \code{n} rows (for \code{n} samples) and \code{2} columns,
#' where the first column represents the observed data and the second column represents its
#' corresponding predicted values
#'
#' @return numeric
.compute_principal_angle <- function(tmp_mat){
  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  vec <- pca_res$rotation[,1]; vec <- vec/.l2norm(vec)
  if(sign(vec[1]) < 0)  vec <- -1*vec
  angle_val <- as.numeric(acos(as.numeric(c(0,1) %*% vec)))
  angle_val * 180/pi
}

