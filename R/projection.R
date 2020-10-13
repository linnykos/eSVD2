.projection_kmeans <- function(mat, k, domain = NA, row = T){
 stopifnot(k <= min(dim(mat)))
 
 if(!all(is.na(domain))){
  stopifnot(length(domain) == 2, domain[1] < domain[2])
  mat <- pmin(mat, domain[2]); mat <- pmax(mat, domain[1])
 }

 if(!row) mat <- t(mat)
 
 kmean_res <- stats::kmeans(mat, centers = k)
 new_mat <- sapply(kmean_res$cluster, function(x){
  kmean_res$centers[x,]
 })
 
 if(row) new_mat <- t(new_mat)
 new_mat
}