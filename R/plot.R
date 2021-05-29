plot_scores_heatmap <- function(score_mat, membership_vec = NA, num_col = 10,
                                bool_center = T, bool_scale = T,
                                bool_log = F, scaling_power = 1, luminosity = F,
                                ...){

  n <- nrow(obj$common_score)
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

plot_maplot <- function(){

}
