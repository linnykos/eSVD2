plot_scatterplot_mean <- function(mat,
                                  esvd_res,
                                  nuisance_vec,
                                  case_control_variable,
                                  alpha_max = 50,
                                  asp = T,
                                  cex = 1,
                                  max_num = 1e5,
                                  mean_type = "predicted",
                                  nuisance_lower_quantile = 0.01,
                                  only_nonzero = T,
                                  xlim = NA,
                                  verbose = T,
                                  ...){
  stopifnot(mean_type %in% c("predicted", "posterior"))

  if(mean_type == "predicted"){
    nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
    nat_mat2 <- tcrossprod(esvd_res$covariates,
                           esvd_res$b_mat)
    mean_mat <- exp(nat_mat1 + nat_mat2)
  } else {
    res <- compute_posterior(mat = mat,
                             esvd_res = esvd_res,
                             nuisance_vec = nuisance_vec,
                             case_control_variable = case_control_variable,
                             alpha_max = alpha_max,
                             nuisance_lower_quantile = nuisance_lower_quantile)
    mean_mat <- res$posterior_mean_mat
  }

  if(only_nonzero) {
    idx <- which(mat != 0)
  } else {
    idx <- 1:prod(dim(mat))
  }
  if(length(idx) > max_num){
    idx <- sample(idx, size = max_num)
  }

  x_vec <- mean_mat[idx]
  y_vec <- mat[idx]
  if(all(is.na(xlim))) xlim <- range(c(0, x_vec, y_vec))

  tmp_mat <- cbind(mat[idx], mean_mat[idx])
  bool_vec <- apply(tmp_mat, 1, function(x){all(x >= xlim[1]) & all(x <= xlim[2])})
  tmp_mat <- tmp_mat[which(bool_vec),]
  angle_val <- .compute_principal_angle(tmp_mat)
  print(angle_val)

  graphics::plot(NA,
                 xlim = xlim,
                 ylim = xlim,
                 asp = asp,
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
  invisible()
}

plot_gene_librarysize <- function(mat,
                                  col_individual = 1,
                                  col_mean1 = 2,
                                  col_mean2 = "white",
                                  library_size_vec = NULL,
                                  lwd_individual = 1,
                                  lwd_mean1 = 3,
                                  lwd_mean2 = 5,
                                  ngroups = 6,
                                  return_obj = F,
                                  verbose = T,
                                  ...){
  if(all(is.null(library_size_vec))) library_size_vec <- colSums(mat)
  avg_expression <- Matrix::colMeans(mat)
  res <- stats::kmeans(avg_expression, centers = ngroups)
  kmean_clust <- res$cluster
  clust_reorder <- sort(res$centers, decreasing = F)
  tmp <- kmean_clust
  for(i in 1:ngroups){
    tmp[kmean_clust == clust_reorder[i]] <- i
  }
  kmean_clust <- tmp

  p <- ncol(mat)
  pred_mat <- sapply(1:p, function(j){
    if(verbose && p > 10 && j %% floor(p/10) == 0) cat('*')
    bw_res <- np::npregbw(xdat = library_size_vec,
                          ydat = mat[,j])
    stats::predict(np::npreg(bw_res))
  })

  lis_obj <- lapply(1:ngroups, function(k){
    idx <- which(kmean_clust == k)
    pred_mat[,idx,drop = F]
  })
  names(lis_obj) <- paste0("cluster", 1:ngroups)
  lis_obj$library_size <- library_size_vec
  lis_obj$range <- range(pred_mat)

  ord_vec <- order(lis_obj$library_size)
  for(k in 1:ngroups){
    plot(NA, xlim = range(lis_obj$library_size), ylim = lis_obj$range,
         ...)
    for(j in ncol(lis_obj[[k]])){
      lines(x = lis_obj$library_size[ord_vec],
            y = lis_obj[[k]][ord_vec,j],
            col = col_individual,
            lwd = lwd_individual)
    }
    lines(x = lis_obj$library_size[ord_vec],
          y = Matrix::colMeans(lis_obj[[k]])[ord_vec],
          col = col_mean2,
          lwd = lwd_mean2)
    lines(x = lis_obj$library_size[ord_vec],
          y = Matrix::colMeans(lis_obj[[k]])[ord_vec],
          col = col_mean1,
          lwd = lwd_mean1)
  }

  if(return_obj) return(lis_obj) else invisible()
}

# https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html
plot_gene_variability <- function(mat,
                                  library_size_vec = NULL,
                                  ngroups_cells = 5,
                                  ngroups_genes = 6,
                                  verbose = T){
  if(all(is.null(library_size_vec))) library_size_vec <- colSums(mat)
  avg_expression <- Matrix::colMeans(mat)
  res <- stats::kmeans(avg_expression, centers = ngroups_genes)$cluster
  gene_clust <- res$cluster
  clust_reorder <- sort(res$centers, decreasing = F)
  tmp <- gene_clust
  for(i in 1:ngroups_genes){
    tmp[gene_clust == clust_reorder[i]] <- i
  }
  gene_clust <- tmp

  res <- stats::kmeans(library_size_vec, centers = ngroups_cells)$cluster
  cell_clust <- res$cluster
  clust_reorder <- sort(res$centers, decreasing = F)
  tmp <- cell_clust
  for(i in 1:ngroups_cells){
    tmp[cell_clust == clust_reorder[i]] <- i
  }
  cell_clust <- tmp

  df <- data.frame(gene_cluster = rep(1:ngroups_genes, times = ngroups_cells),
                   cell_cluster = rep(1:ngroups_cells, each = ngroups_genes))
  vec <- sapply(1:nrow(df), function(i){
    gene_idx <- which(gene_clust == df[i,"gene_cluster"])
    cell_idx <- which(cell_clust ==  df[i,"cell_cluster"])
    var_val <- mean(sapply(gene_idx, function(j){
      stats::var(mat[cell_idx,j])
    }))
  })
  df$variance <- vec

  df
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

