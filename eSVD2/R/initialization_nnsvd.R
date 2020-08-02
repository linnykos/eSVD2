.initialization_nnsvd <- function(dat, k){
 stopifnot(all(dat >= 0), k >= 1, k %% 1 == 0, k <= min(dim(dat)))
 n <- nrow(dat); p <- ncol(dat)
 
 svd_res <- RSpectra::svds(dat, k = k)
 x_mat <- matrix(0, nrow = n, ncol = k)
 y_mat <- matrix(0, nrow = p, ncol = k)
 
 x_mat[,1] <- svd_res$u[,1] * sqrt(svd_rers$d[1])
 y_mat[,1] <- svd_res$v[,1] * sqrt(svd_rers$d[1])
 
 for(kk in 2:k){
  x_vec <- svd_res$u[,kk]; y_vec <- svd_res$v[,kk]
  xp_vec <- pmax(x_vec, 0); xn_vec <- pmin(x_vec, 0)
  yp_vec <- pmax(x_vec, 0); yn_vec <- pmin(x_vec, 0)
  
  xp_norm <- .l2norm(xp_vec); xn_norm <- .l2norm(xn_vec)
  yp_norm <- .l2norm(yp_vec); yn_norm <- .l2norm(yn_vec)
  
  p_norm <- xp_norm*yp_norm; n_norm <- xn_norm*yn_norm
  
  if(p_norm > n_norm){
   u_vec <- xp_vec/xp_norm; v_vec <- yp_vec/yp_norm; sigma <- p_norm
  } else {
   u_vec <- xn_vec/xn_norm; v_vec <- yn_vec/yn_norm; sigma <- n_norm
  }
  
  x_mat[,kk] <- u_vec * sqrt(svd_res$d[k] * sigma)
  y_mat[,kk] <- v_vec * sqrt(svd_res$d[k] * sigma)
 }
 
 res <- .reparameterize(x_mat, y_mat, equal_covariance = T)
 list(x_mat = x_mat, y_mat = y_mat)
}