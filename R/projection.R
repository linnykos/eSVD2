.projection_kmeans <- function(mat, k, domain = NA, row = T){
 stopifnot(k <= min(dim(mat)), is.matrix(mat))

 if(!all(is.na(domain))){
  stopifnot(length(domain) == 2, domain[1] < domain[2])
  mat <- pmin(mat, domain[2]); mat <- pmax(mat, domain[1])
 }

 if(!row) mat <- t(mat)

 svd_res <- .svd_truncated(mat, K = k,
                           symmetric = F,
                           rescale = F,
                           mean_vec = NULL,
                           sd_vec = NULL,
                           K_full_rank = F)
 dr_mat <- .mult_mat_vec(svd_res$u, svd_res$d)

 kmean_res <- stats::kmeans(dr_mat, centers = k)
 row_mat <- t(sapply(1:k, function(i){
   idx <- which(kmean_res$cluster == i)
   matrixStats::colMeans2(mat[idx,,drop = F])
 }))
 new_mat <- row_mat[kmean_res$cluster,,drop = F]

 if(row) new_mat <- t(new_mat)
 new_mat
}
