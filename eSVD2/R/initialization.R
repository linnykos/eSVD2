initialize_esvd <- function(dat, family, k, nuisance_param_vec = NA, library_size_vec = NA,
                            config = initalization_default()){
 stopifnot(all(dat[!is.na(dat)] >= 0))
 
 dat <- .matrix_completion(dat, k = k)
 init_res <- .determine_initial_matrix(dat, family = family, k = k, max_val = config$max_val,
                                       tol = config$tol)
 
 if(config$method == "kmean_row"){
  nat_mat <- .initialization_kmean(init_res$nat_mat, k = k, domain = init_res$domain, row = T)
 } else {
   stop("config method not found")
 }

 # reparameterize
 res <- .factorize_matrix(nat_mat, k = k, equal_covariance = T)
 .fix_rank_defficiency_initialization(res$x_mat, res$y_mat, domain = domain)
}

initalization_default <- function(method = "kmean_row", max_val = NA, tol = 1e-3){
 stopifnot(method %in% c("kmean_row"), tol > 0, tol <= 1)
 
 structure(list(method = method, max_val = max_val, tol = tol), class = "initialization_param")
}


################

.initialization_kmean <- function(mat, k, domain, row = T){
 .projection_kmeans(mat, k = k, domain = domain, row = row)
}

###############

#' Fill in missing values
#'
#' Uses \code{softImpute::softImpute} to fill in all the possible missing values.
#' This function enforces all the resulting entries to be non-negative.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param k positive integer less than \code{min(c(nrow(dat), ncol(dat)))}
#'
#' @return a \code{n} by \code{p} matrix
.matrix_completion <- function(dat, k){
 if(any(is.na(dat))){
  lambda0_val <- softImpute::lambda0(dat)
  res <- softImpute::softImpute(dat, rank.max = k, lambda = min(30, lambda0_val/100))
  diag_mat <- .diag_matrix(res$d[1:k])
  pred_naive <- res$u %*% diag_mat %*% t(res$v)
  dat[which(is.na(dat))] <- pred_naive[which(is.na(dat))]
 }
 
 pmax(dat, 0)
}

#' Initialize the matrix of natural parameters
#'
#' This function first transforms each entry in \code{dat} according to the inverse function that maps
#' natural parameters to their expectation (according to \code{eSVD:::.mean_transformation}) and then
#' uses \code{eSVD:::.project_rank_feasibility} to get a rank-\code{k} approximation of this matrix
#' that lies within the domain of \code{family}
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes.
#' @param k  positive integer less than \code{min(c(nrow(dat), ncol(dat)))}
#' @param family character (\code{"gaussian"}, \code{"exponential"}, \code{"poisson"}, \code{"neg_binom"},
#' or \code{"curved gaussian"})
#' @param max_val maximum magnitude of the inner product
#' @param tol numeric
#' @param ... extra arguments, such as nuisance parameters for \code{"neg_binom"}
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return \code{n} by \code{p} matrix
.determine_initial_matrix <- function(dat, family, k, max_val = NA, tol = 1e-3, ...){
 stopifnot((is.na(max_val) || max_val >= 0), all(dat >= 0))
 
 domain <- .determine_domain(family, tol)
 if(!is.na(max_val)) domain <- .intersect_intervals(domain, c(-max_val, max_val))
 
 dat[which(dat <= tol)] <- tol/2
 nat_mat <- .mean_transformation(dat, family, ...)
 nat_mat <- pmax(nat_mat, domain[1])
 nat_mat <- pmin(nat_mat, domain[2])
 
 list(nat_mat = nat_mat, domain = domain)
}


#' Fix rank defficiency among two matrices
#'
#' Given two matrices, \code{x_mat} and \code{y_mat} (both with the same number of columns),
#' adjust these two matrices so the \code{x_mat \%*\% t(y_mat)} is actually of the desired
#' rank. This function is needed since sometimes upstream, the matrices \code{x_mat} and \code{y_mat}
#' do not actually have the rank equal to the number of columns (i.e., empirically we have
#' observed that \code{x_mat} might have a column that is all constant).
#'
#' @param x_mat a numeric matrix
#' @param y_mat a numeric matrix with the same number of columns as \code{x_mat}
#' @param direction character either \code{"<="} or \code{">="} or \code{NA}
#'
#' @return a list of \code{x_mat} and \code{y_mat}
.fix_rank_defficiency_initialization <- function(x_mat, y_mat, domain){
 k <- ncol(x_mat)
 nat_mat <- x_mat %*% t(y_mat)
 k2 <- as.numeric(Matrix::rankMatrix(nat_mat))
 
 if(k != k2){
  stopifnot(k2 < k)
  sign_val <- ifelse(abs(domain[1]) < abs(domain[2]), 1, -1)
  
  sd_val <- mean(c(apply(x_mat[,1:k2, drop = F], 2, stats::sd),apply(y_mat[,1:k2, drop = F], 2, stats::sd)))
  for(i in (k2+1):k){
   x_mat[,i] <- abs(stats::rnorm(nrow(x_mat), sd = sd_val/10))
   y_mat[,i] <- sign_val*abs(stats::rnorm(nrow(y_mat), sd = sd_val/10))
  }
 }
 
 list(x_mat = x_mat, y_mat = y_mat)
}
