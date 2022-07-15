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

  svd_res <- .svd_safe(mat = mat,
                       check_stability = T,
                       K = k,
                       mean_vec = NULL,
                       rescale = F,
                       scale_max = NULL,
                       sd_vec = NULL)
  x_mat <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_mat <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  .reparameterize(x_mat, y_mat, equal_covariance = equal_covariance)
}

.reparameterization_esvd_covariates <- function(eSVD_obj,
                                                fit_name,
                                                omitted_variables,
                                                verbose = 0){

  x_mat <- eSVD_obj[[fit_name]]$x_mat
  y_mat <- eSVD_obj[[fit_name]]$y_mat
  z_mat <- eSVD_obj[[fit_name]]$z_mat
  k <- ncol(x_mat)
  p <- nrow(y_mat)

  covariate_mat <- eSVD_obj$covariates
  stopifnot("Intercept" %in% colnames(covariate_mat))
  covariate_mat <- covariate_mat[,which(!colnames(covariate_mat) %in% omitted_variables),drop=F]
  covariate_mat2 <- covariate_mat[,which(colnames(covariate_mat) != "Intercept")]

  for(ell in 1:k){
    tmp_df <- cbind(x_mat[,ell], covariate_mat2)
    colnames(tmp_df)[1] <- "x"
    tmp_df <- as.data.frame(tmp_df)

    lm_res <- stats::lm(x ~ . , data = tmp_df)
    coef_vec <- stats::coef(lm_res)
    names(coef_vec)[1] <- "Intercept"
    if(verbose > 0) print(paste0(ell, ": R2 of ", round(summary(lm_res)$r.squared, 2)))

    for(j in 1:p){
      z_mat[j,names(coef_vec)] <- z_mat[j,names(coef_vec)] + coef_vec*y_mat[j,ell]
    }

    x_mat[,ell] <- stats::residuals(lm_res)
  }

  res <- eSVD2:::.reparameterize(x_mat, y_mat, equal_covariance = T)
  x_mat <- res$x_mat; y_mat <- res$y_mat

  eSVD_obj[[fit_name]]$x_mat <- x_mat
  eSVD_obj[[fit_name]]$y_mat <- y_mat
  eSVD_obj[[fit_name]]$z_mat <- z_mat

  eSVD_obj
}


#########

.svd_safe <- function(mat,
                      check_stability, # boolean
                      K, # positive integer
                      mean_vec, # boolean, NULL or vector
                      rescale, # boolean
                      scale_max, # NULL or positive integer
                      sd_vec){ # boolean, NULL or vector
  if(is.na(K)) K <- min(dim(mat))
  stopifnot(min(dim(mat)) >= K)

  mean_vec <- .compute_matrix_mean(mat, mean_vec)
  sd_vec <- .compute_matrix_sd(mat, sd_vec)

  res <- .svd_in_sequence(check_stability = check_stability,
                          K = K,
                          mat = mat,
                          mean_vec = mean_vec,
                          scale_max = scale_max,
                          sd_vec = sd_vec)
  res <- list(d = res$d, u = res$u, v = res$v, method = res$method)
  class(res) <- "svd"

  # pass row-names and column-names
  if(length(rownames(mat)) > 0) rownames(res$u) <- rownames(mat)
  if(length(colnames(mat)) > 0) rownames(res$v) <- colnames(mat)

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

.svd_in_sequence <- function(check_stability,
                             K,
                             mat,
                             mean_vec,
                             scale_max,
                             sd_vec){
  res <- tryCatch({
    .irlba_custom(check_stability = check_stability,
                  K = K,
                  mat = mat,
                  mean_vec = mean_vec,
                  scale_max = scale_max,
                  sd_vec = sd_vec)
  },
  warning = function(e){NULL},
  error = function(e){NULL})
  if(!all(is.null(res))) {res$method <- "irlba"; return(res)}

  ##

  res <- tryCatch({
    .rpsectra_custom(check_stability = check_stability,
                     K = K,
                     mat = mat,
                     mean_vec = mean_vec,
                     scale_max = scale_max,
                     sd_vec = sd_vec)
  },
  warning = function(e){NULL},
  error = function(e){NULL})
  if(!all(is.null(res))) {res$method <- "RSpectra"; return(res)}

  ##

  if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
  if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
  if(!is.null(scale_max)){
    mat[mat > abs(scale_max)] <- abs(scale_max)
    mat[mat < -abs(scale_max)] <- -abs(scale_max)
  }
  res <- svd(mat)
  res$method <- "base"
  res
}

.irlba_custom <- function(check_stability,
                          K,
                          mat,
                          mean_vec,
                          scale_max,
                          sd_vec){
  if(inherits(mat, "dgCMatrix")){
    if(!all(is.null(scale_max))) warning("scale_max does not work with sparse matrices when using irlba")
    tmp <- irlba::irlba(A = mat,
                        nv = K,
                        work = min(c(K + 10, dim(mat))),
                        scale = sd_vec,
                        center = mean_vec)

    if(check_stability & K > 5) {
      tmp2 <- irlba::irlba(A = mat,
                           nv = 5,
                           scale = sd_vec,
                           center = mean_vec)
      ratio_vec <- tmp2$d/tmp$d[1:5]
      if(any(ratio_vec > 2) | any(ratio_vec < 1/2)) warning("irlba is potentially unstable")
    }

    return(tmp)

  } else {
    if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
    if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
    if(!is.null(scale_max)){
      mat[mat > abs(scale_max)] <- abs(scale_max)
      mat[mat < -abs(scale_max)] <- -abs(scale_max)
    }

    tmp <- irlba::irlba(A = mat, nv = K)

    if(check_stability & K > 5) {
      tmp2 <- irlba::irlba(A = mat, nv = 5)
      ratio_vec <- tmp2$d/tmp$d[1:5]
      if(any(ratio_vec > 2) | any(ratio_vec < 1/2)) warning("irlba is potentially unstable")
    }

    return(tmp)
  }
}

.rpsectra_custom <- function(check_stability,
                             K,
                             mat,
                             mean_vec,
                             scale_max,
                             sd_vec){

  if(inherits(mat, "dgCMatrix")){
    if(!all(is.null(mean_vec))) warning("mean_vec does not work with sparse matrices when using RSpectra")
    if(!all(is.null(sd_vec))) warning("sd_vec does not work with sparse matrices when using RSpectra")
    if(!all(is.null(scale_max))) warning("scale_max does not work with sparse matrices when using RSpectra")
  } else {
    if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
    if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
    if(!is.null(scale_max)){
      mat[mat > abs(scale_max)] <- abs(scale_max)
      mat[mat < -abs(scale_max)] <- -abs(scale_max)
    }
  }

  tmp <- RSpectra::svds(A = mat, k = K)

  if(check_stability & K > 5) {
    tmp2 <- RSpectra::svds(A = mat, k = 5)
    ratio_vec <- tmp2$d/tmp$d[1:5]
    if(any(ratio_vec > 2) | any(ratio_vec < 1/2)) warning("RSpectra is potentially unstable")
  }

  tmp
}

.compute_matrix_mean <- function(mat, mean_vec){
  if(length(mean_vec) == 1 && !is.null(mean_vec)){
    if(mean_vec){
      mean_vec <- Matrix::colMeans(mat)
    } else{
      mean_vec <- NULL
    }
  }

  mean_vec
}

.compute_matrix_sd <- function(mat, sd_vec){
  if(length(sd_vec) == 1 && !is.null(sd_vec)){
    if(sd_vec){
      if(inherits(x = mat, what = c('dgCMatrix', 'dgTMatrix'))){
        sd_vec <- sparseMatrixStats::colSds(mat)
      } else {
        sd_vec <- matrixStats::colSds(mat)
      }
    } else{
      sd_vec <- NULL
    }
  }

  sd_vec
}


