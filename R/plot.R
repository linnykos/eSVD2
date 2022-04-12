plot_scatterplot_mean <- function(mat,
                                  esvd_res,
                                  nuisance_vec,
                                  case_control_variable,
                                  alpha_max = 50,
                                  asp = T,
                                  bool_logscale = F,
                                  cex = 1,
                                  col_points = rgb(0.5, 0.5, 0.5, 0.1),
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

    offset_var <- setdiff(colnames(esvd_res_full$covariates), case_control_variable)
    library_mat <- exp(tcrossprod(
      esvd_res$covariates[,offset_var],
      esvd_res$b_mat[,offset_var]
    ))
    mat <- mat/library_mat
  }

  if(only_nonzero) {
    idx <- which(mat != 0)
  } else {
    idx <- 1:prod(dim(mat))
  }
  if(length(idx) > max_num){
    idx <- sample(idx, size = max_num)
  }

  tmp_mat <- cbind(mat[idx], mean_mat[idx])
  if(bool_logscale){
    tmp_mat <- log1p(tmp_mat)
  }
  x_vec <- tmp_mat[,2]
  y_vec <- tmp_mat[,1]
  if(all(is.na(xlim))) xlim <- range(c(0, x_vec, y_vec))

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
  graphics::points(x = x_vec,
                   y = y_vec,
                   pch = 16,
                   cex = cex,
                   col = col_points)
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
  invisible()
}

compute_gene_librarysize <- function(mat,
                                     avg_expression = NULL,
                                     library_size_vec = NULL,
                                     deg = 2,
                                     ngroups = 6,
                                     num_boot = 100,
                                     verbose = 1){
  if(all(is.null(library_size_vec))) library_size_vec <- colSums(mat)
  if(all(is.null(avg_expression))) avg_expression <- Matrix::colMeans(mat)
  res <- stats::kmeans(avg_expression, centers = ngroups)
  if(verbose > 0) print(cbind(res$size, res$centers))
  kmean_clust <- res$cluster
  clust_reorder <- order(res$centers, decreasing = F)
  tmp <- kmean_clust
  for(i in 1:ngroups){
    tmp[kmean_clust == clust_reorder[i]] <- i
  }
  kmean_clust <- tmp
  if(verbose > 0) print(table(kmean_clust))

  p <- ncol(mat)
  pred_mat <- sapply(1:p, function(j){
    if(verbose == 1 && p > 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose == 2) print(j)
    dat <- data.frame(y = mat[,j], x = library_size_vec)
    np_fit <- npregfast::frfast(y ~ x,
                                data = dat,
                                p = deg,
                                nboot = num_boot)

    np_fit$p[,1,1]
  })

  dat <- data.frame(y = mat[,1], x = library_size_vec)
  np_fit <- npregfast::frfast(y ~ x,
                              data = dat,
                              p = deg,
                              nboot = num_boot)
  fitted_library_vec <- np_fit$x

  lis_obj <- lapply(1:ngroups, function(k){
    idx <- which(kmean_clust == k)
    pred_mat[,idx,drop = F]
  })
  names(lis_obj) <- paste0("cluster", 1:ngroups)
  lis_obj$library_size <- fitted_library_vec
  lis_obj$range <- range(pred_mat)
  lis_obj
}

plot_gene_librarysize <- function(lis_obj,
                                  bool_fix_ylim = F,
                                  col_rgb_r = 0.5,
                                  col_rgb_g = 0.5,
                                  col_rgb_b = 0.5,
                                  col_mean1 = 2,
                                  col_mean2 = "white",
                                  lwd_individual = 1,
                                  lwd_mean1 = 3,
                                  lwd_mean2 = 5,
                                  ...){
  stopifnot(is.list(lis_obj), c("library_size", "range") %in% names(lis_obj))

  for(k in 1:length(grep("cluster", names(lis_obj)))){
    col_var <- rgb(col_rgb_r, col_rgb_g, col_rgb_b, pmax(pmin(1/sqrt(ncol(lis_obj[[k]])), 0.5), 0.1))

    if(bool_fix_ylim) ylim_tmp <- lis_obj$range else ylim_tmp <- pmax(quantile(lis_obj[[k]], probs = c(0, 0.95)), 0)
    plot(NA, xlim = range(lis_obj$library_size), ylim = ylim_tmp,
         main = paste0("Cluster ", k, "\n(", ncol(lis_obj[[k]]), " genes)"),
         ...)
    for(j in 1:ncol(lis_obj[[k]])){
      lines(x = lis_obj$library_size,
            y = lis_obj[[k]][,j],
            col = col_var,
            lwd = lwd_individual)
    }
    lines(x = lis_obj$library_size,
          y = Matrix::rowMeans(lis_obj[[k]]),
          col = col_mean2,
          lwd = lwd_mean2)
    lines(x = lis_obj$library_size,
          y = Matrix::rowMeans(lis_obj[[k]]),
          col = col_mean1,
          lwd = lwd_mean1)
  }

  invisible()
}

# https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html
plot_gene_variability <- function(mat,
                                  library_size_vec = NULL,
                                  ngroups_cells = 5,
                                  ngroups_genes = 6){
  if(all(is.null(library_size_vec))) library_size_vec <- colSums(mat)
  avg_expression <- Matrix::colMeans(mat)
  res <- stats::kmeans(avg_expression, centers = ngroups_genes)
  gene_clust <- res$cluster
  clust_reorder <- order(res$centers, decreasing = F)
  tmp <- gene_clust
  for(i in 1:ngroups_genes){
    tmp[gene_clust == clust_reorder[i]] <- i
  }
  gene_clust <- tmp

  res <- stats::kmeans(library_size_vec, centers = ngroups_cells)
  cell_clust <- res$cluster
  clust_reorder <- order(res$centers, decreasing = F)
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
  df$gene_cluster <- as.factor(df$gene_cluster)
  df$cell_cluster <- as.factor(df$cell_cluster)

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

