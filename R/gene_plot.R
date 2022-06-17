gene_plot <- function(input_obj,
                      what_1,
                      what_2 = NULL,
                      gene_list = NULL,
                      color_palette = NULL,
                      which_fit = "latest_Fit",
                      xlab = NULL,
                      ylab = NULL,
                      ...){
  what_vec <- c("nuisance", "teststat", colnames(input_obj$covariates),
                "sparsity")
  stopifnot(what_1 %in% what_vec,
            is.null(what_2) || what_2 %in% what_vec)
  if(!all(is.null(gene_list)) | !all(is.null(color_palette))){
    stopifnot(is.list(gene_list), length(gene_list) == length(color_palette),
              all(unlist(gene_list) %in% colnames(input_obj$dat)),
              !any(duplicated(unlist(gene_list))))
  }

  # first compute the relevant vectors
  vec_1 <- .compute_what_gene_vector(input_obj = input_obj,
                                     what = what_1,
                                     which_fit = which_fit)
  vec_2 <- .compute_what_gene_vector(input_obj = input_obj,
                                     what = what_2,
                                     which_fit = which_fit)
  stopifnot(length(names(vec_1)) > 0)
  if(!all(is.null(vec_2))){
    stopifnot(length(vec_1) == length(vec_2), all(names(vec_1) == names(vec_2)))
  }

  if(all(!is.null(vec_2))) {
    vec_all <- cbind(vec_1, vec_2)
  } else {
    vec_all <- matrix(vec_1, nrow = length(vec_1), ncol = 1)
  }

  # now plot
  if(!all(is.null(gene_list))){
    color_vec <- rep(NA, length(vec_1))
    for(i in 1:length(gene_list)){
      color_vec[which(names(vec_1) %in% gene_list[[i]])] <- color_palette[i]
    }
    shuff_idx <- sample(which(!is.na(color_vec)))
  }

  if(all(!is.null(vec_2))){
    if(is.null(xlab)) xlab <- what_1
    if(is.null(ylab)) ylab <- what_2

    graphics::plot(x = vec_1, y = vec_2,
                   cex = 1,
                   pch = 16,
                   col = rgb(0.5, 0.5, 0.5, 0.2),
                   xlab = xlab, ylab = ylab,
                   ...)

    if(!all(is.null(gene_list))){
      graphics::points(x = vec_1[shuff_idx],
                       y = vec_2[shuff_idx],
                       col = "white", pch = 16, cex = 3)
      graphics::points(x = vec_1[shuff_idx],
                       y = vec_2[shuff_idx],
                       col = color_vec[shuff_idx],
                       pch = 16, cex = 2.5)
    }

  } else {
    if(is.null(xlab)) xlab <- what_1
    graphics::hist(vec_1,
                   xlab = xlab,
                   ...)

    if(!all(is.null(gene_list))){
      for(idx in shuff_idx){
        graphics::rug(vec_1[idx], col = color_vec[idx],
                      side = 1, lwd = 2)
      }
    }
  }

  invisible()

}

.compute_what_gene_vector <- function(input_obj,
                                      what,
                                      which_fit){
  if(is.null(what)) return(NULL)
  if(which_fit == "latest_Fit"){
    latest_Fit <- .get_object(eSVD_obj = input_obj,
                              what_obj = "latest_Fit",
                              which_fit = NULL)
  } else {
    stopifnot(which_fit %in% names(input_obj))
    latest_Fit <- which_fit
  }

  if(what == "nuisance") {
    return(input_obj[[latest_Fit]]$nuisance_vec)

  } else if (what == "teststat"){
    return(input_obj$teststat_vec)

  } else if (what == "sparsity"){
    dat <- input_obj$dat
    n <- nrow(dat); p <- ncol(dat)
    vec <- sapply(1:p, function(j){
      (p - length(.nonzero_col(dat, col_idx = j, bool_value = F)))/p
    })
    return(vec)

  } else if (what %in% colnames(input_obj$covariates)){
    return(input_obj[[latest_Fit]]$z_mat[,what])

  } else {
    stop("what not found")
  }

}
