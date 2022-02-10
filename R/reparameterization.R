.identification <- function(cov_x, cov_y, check = F, tol = 1e-6){
 stopifnot(all(dim(cov_x) == dim(cov_y)), nrow(cov_x) == ncol(cov_x))
 if(nrow(cov_x) == 1){
  return(matrix((as.numeric(cov_y)/as.numeric(cov_x))^(1/4), 1, 1))
 }

 eigen_x <- eigen(cov_x)
 eigen_y <- eigen(cov_y)

 Vx <- eigen_x$vectors
 Vy <- eigen_y$vectors

 if(any(eigen_x$values <= tol) | any(eigen_y$values <= tol))
   warning("Detecting rank defficiency in reparameterization step")

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
.reparameterize <- function(x_mat, y_mat, equal_covariance){
 stopifnot(ncol(x_mat) == ncol(y_mat))
 n <- nrow(x_mat); p <- nrow(y_mat)

 res <- .identification(crossprod(x_mat), crossprod(y_mat))

 if(equal_covariance){
  list(x_mat = (n/p)^(1/4)*tcrossprod(x_mat, res), y_mat = (p/n)^(1/4)*y_mat %*% solve(res))
 } else {
  list(x_mat = tcrossprod(x_mat, res), y_mat = y_mat %*% solve(res))
 }
}

.factorize_matrix <- function(mat, k, equal_covariance){
  stopifnot(k <= min(dim(mat)))

  svd_res <- .svd_truncated(mat, K = k, symmetric = F, rescale  = F,
                            mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
  x_mat <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_mat <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  .reparameterize(x_mat, y_mat, equal_covariance = equal_covariance)
}

#########

.svd_truncated <- function(mat, K, symmetric, rescale,
                           mean_vec, sd_vec,
                           K_full_rank){
  if(is.na(K)) K <- min(dim(mat))
  stopifnot(min(dim(mat)) >= K)
  if(K == min(dim(mat))) K_full_rank <- T

  ## [[note to self: probably factor this out and add sd_vec]]
  if(length(mean_vec) == 1 && !is.null(mean_vec)){
    if(mean_vec){
      if(inherits(x = mat, what = c('dgCMatrix', 'dgTMatrix'))){
        mean_vec <- Matrix::colMeans(mat)
      } else {
        mean_vec <- matrixStats::colMeans2(mat)
      }
    } else{
      mean_vec <- NULL
    }
  }

  if(min(dim(mat)) > 2*(K+2)){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      if(symmetric){
        tmp <- irlba::partial_eigen(mat, n = ifelse(K_full_rank, K, K+2),
                                    center = mean_vec, scale = sd_vec)
        list(u = tmp$vectors, d = tmp$values, v = tmp$vectors)
      } else {
        irlba::irlba(mat, nv = ifelse(K_full_rank, K, K+2),
                     center = mean_vec, scale = sd_vec)
      }
    }, warning = function(e){
      if(!all(is.null(mean_vec)) | !all(is.null(sd_vec))) print("mean_vec or sd_vec not used")
      RSpectra::svds(mat, k = ifelse(K_full_rank, K, K+2))
    }, error = function(e){
      if(!all(is.null(mean_vec)) | !all(is.null(sd_vec))) print("mean_vec or sd_vec not used")
      RSpectra::svds(mat, k = ifelse(K_full_rank, K, K+2))
    })
  } else {
    res <- svd(mat)
  }

  res$u <- res$u[,1:K, drop = F]; res$v <- res$v[,1:K, drop = F]; res$d <- res$d[1:K]

  # pass row-names and column-names
  if(length(rownames(mat)) != 0) rownames(res$u) <- rownames(mat)
  if(length(colnames(mat)) != 0) rownames(res$v) <- colnames(mat)

  # useful only if your application requires only the singular vectors
  # if the number of rows or columns is too large, the singular vectors themselves
  # are often a bit too small numerically
  if(rescale){
    n <- nrow(mat); p <- ncol(mat)
    res$u <- res$u * sqrt(n)
    res$v <- res$v * sqrt(p)
    res$d <- res$d / (sqrt(n)*sqrt(p))
  }

  res
}



