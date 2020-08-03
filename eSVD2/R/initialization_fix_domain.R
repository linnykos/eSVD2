# assumes that domain (one of the values) is close to 0 and the other domain is a large magnitude
.fix_domain <- function(nat_mat, dat, family, domain, ...){
 if(abs(domain[2]) > abs(domain[1])) nat_mat <- -nat_mat
 max_val <- abs(max(domain))
 
 while(TRUE){
  bool_mat <- (nat_mat >= max_val)
  if(sum(bool_mat) == 0) break()
  
  idx <- which(nat_mat == which.max(nat_mat), arr.ind = T)
  idx_i <- idx[1]; idx_j <- idx[2]
  
  target_row <- dat[idx_i,]; target_column <- dat[,idx_j]
  
  old_row <- compute_mean(nat_mat[idx_i,], family = family, ...)
  old_column <- compute_mean(nat_mat[,idx_j], family = family, ...)
  
  rescale_factor <- max_val/nat_mat[idx_i, idx_j]
  new_row <- compute_mean(nat_mat[idx_i,]*rescale_factor, family = family, ...)
  new_column <- compute_mean(nat_mat[,idx_j]*rescale_factor, family = family, ...)
  
  row_error <- .l2norm(c(target_row, target_column) - c(new_row, old_column))
  column_error <- .l2norm(c(target_row, target_column) - c(old_row, new_column))
  
  if(row_error < column_error){
   nat_mat[idx_i,] <- nat_mat[idx_i,]*rescale_factor
  } else {
   nat_mat[,idx_j] <- nat_mat[,idx_j]*rescale_factor
  }
 }
 
 if(abs(domain[2]) > abs(domain[1])) nat_mat <- -nat_mat
}