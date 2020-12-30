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

 Dx <- eigen_x$values
 Dy <- eigen_y$values

 # form R
 tmp <- crossprod(.mult_mat_vec(Vy, sqrt(Dy)), .mult_mat_vec(Vx, sqrt(Dx)))
 svd_tmp <- svd(tmp)
 Q <- tcrossprod(svd_tmp$u, svd_tmp$v)

 # run a check
 if(check){
  sym_mat <- crossprod(Q, tmp)
  stopifnot(sum(abs(sym_mat - t(sym_mat))) <= 1e-6)
 }

 # now form the symmetric matrix to later factorize
 sym_prod <- tcrossprod(tcrossprod(.mult_mat_vec(Vx, Dx^(-1/2)), Q), .mult_mat_vec(Vy, sqrt(Dy)))
 sym_prod[which(abs(sym_prod) <= tol)] <- 0

 if(check){
  stopifnot(sum(abs(sym_prod - t(sym_prod))) <= 1e-6)
 }

 eigen_sym <- eigen(sym_prod)
 W_mat <- .mult_vec_mat(sqrt(eigen_sym$values), t(eigen_sym$vectors))

 if(check){
  mat1 <- tcrossprod(W_mat %*% cov_x, W_mat)

  W_mat_inv <- solve(W_mat)
  mat2 <- crossprod(W_mat_inv, cov_y) %*% W_mat_inv
  stopifnot(sum(abs(mat1 - mat2)) <= 1e-6)
 }

 # adjust the transformation so it yields a diagonal matrix
 eig_res <- eigen(tcrossprod(W_mat %*% cov_x, W_mat))

 crossprod(eig_res$vectors, W_mat)
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

 res <- .identification(crossprod(x_mat), crossprod(y_mat))

 if(equal_covariance){
  list(x_mat = (n/p)^(1/4)*tcrossprod(x_mat, res), y_mat = (p/n)^(1/4)*y_mat %*% solve(res))
 } else {
  list(x_mat = tcrossprod(x_mat, res), y_mat = y_mat %*% solve(res))
 }
}

.factorize_matrix <- function(mat, k, equal_covariance = T){
  stopifnot(k <= min(dim(mat)))

  svd_res <- RSpectra::svds(mat, k = k)
  x_mat <- .mult_mat_vec(svd_res$u[,1:k, drop = F], sqrt(svd_res$d[1:k]))
  y_mat <- .mult_mat_vec(svd_res$v[,1:k, drop = F], sqrt(svd_res$d[1:k]))

  .reparameterize(x_mat, y_mat, equal_covariance = equal_covariance)
}



