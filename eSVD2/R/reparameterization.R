.identification <- function(cov_x, cov_y, check = F, tol = 1e-6){
 stopifnot(all(dim(cov_x) == dim(cov_y)), nrow(cov_x) == ncol(cov_x))
 if(nrow(cov_x) == 1){
  return(matrix((as.numeric(cov_y)/as.numeric(cov_x))^(1/4), 1, 1))
 }
 
 eigen_x <- eigen(cov_x)
 eigen_y <- eigen(cov_y)
 
 Vx <- eigen_x$vectors
 Vy <- eigen_y$vectors
 
 if(any(eigen_x$values <= tol) | any(eigen_y$values <= tol)) warning("Detecting rank defficiency in reparameterization step")
 
 Dx <- diag(eigen_x$values)
 Dy <- diag(eigen_y$values)
 
 tmp <- sqrt(Dy) %*% t(Vy) %*% Vx %*% sqrt(Dx)
 svd_tmp <- svd(tmp)
 R <- svd_tmp$u %*% t(svd_tmp$v)
 
 # run a check
 if(check){
  Q <- t(R) %*% tmp
  stopifnot(sum(abs(Q - t(Q))) <= 1e-6)
 }
 
 Dx_inv <- Dx; diag(Dx_inv) <- 1/diag(Dx)
 sym_prod <- Vx %*% sqrt(Dx_inv) %*% t(R) %*% sqrt(Dy) %*% t(Vy)
 sym_prod[which(abs(sym_prod) <= tol)] <- 0
 
 if(check){
  stopifnot(sum(abs(sym_prod - t(sym_prod))) <= 1e-6)
 }
 
 eigen_sym <- eigen(sym_prod)
 T_mat <- diag(sqrt(eigen_sym$values)) %*% t(eigen_sym$vectors)
 
 if(check){
  mat1 <- T_mat %*% cov_x %*% t(T_mat)
  
  T_mat_inv <- solve(T_mat)
  mat2 <- t(T_mat_inv) %*% cov_y %*% T_mat_inv
  stopifnot(sum(abs(mat1 - mat2)) <= 1e-6)
 }
 
 # adjust the transformation so it yields a diagonal matrix
 eig_res <- eigen(T_mat %*% cov_x %*% t(T_mat))
 
 t(eig_res$vectors) %*% T_mat
}


#' Function to reparameterize two matrices
#'
#' test
#'
#' Designed to output matrices of the same dimension as \code{x_mat}
#' and \code{y_mat}, but linearly transformed so \code{x_mat \%*\% t(y_mat)}
#' is preserved but either \code{x_mat \%*\% t(x_mat)} is diagonal and equal to
#' \code{y_mat \%*\% t(y_mat)} (if \code{equal_covariance} is \code{FALSE})
#' or \code{x_mat \%*\% t(x_mat)/nrow(x_mat)} is diagonal and equal to
#' \code{y_mat \%*\% t(y_mat)/nrow(y_mat)} (if \code{equal_covariance} is \code{TRUE})
#'
#' @param x_mat matrix of dimension \code{n} by \code{k}
#' @param y_mat matrix of dimension \code{p} by \code{k}
#' @param equal_covariance boolean
#'
#' @return list of two matrices
.reparameterize <- function(x_mat, y_mat, equal_covariance = F){
 stopifnot(ncol(x_mat) == ncol(y_mat))
 n <- nrow(x_mat); p <- nrow(y_mat)
 
 res <- .identification(t(x_mat) %*% x_mat, t(y_mat) %*% y_mat)
 
 if(equal_covariance){
  list(x_mat = (n/p)^(1/4)*x_mat %*% t(res), y_mat = (p/n)^(1/4)*y_mat %*% solve(res))
 } else {
  list(x_mat = x_mat %*% t(res), y_mat = y_mat %*% solve(res))
 }
}

.factorize_matrix <- function(mat, k, equal_covariance = T){
  stopifnot(k <= min(dim(mat)))
  
  svd_res <- RSpectra::svds(mat, k = k)
  x_mat <- svd_res$u[,1:k, drop = F] %*% .diag_matrix(sqrt(svd_res$d[1:k]))
  y_mat <- svd_res$v[,1:k, drop = F] %*% .diag_matrix(sqrt(svd_res$d[1:k]))
  
  .reparameterize(x_mat, y_mat, equal_covariance = equal_covariance)
}



