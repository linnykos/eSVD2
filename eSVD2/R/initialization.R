initialization <- function(dat, k, nuisance_param_vec = NA, library_size_vec = NA){
 
}

################

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
 dat[which(dat <= tol)] <- tol/2
 nat_mat <- .mean_transformation(dat, family, ...)
 direction <- .dictate_direction(family)
 
 ## INSERT HERE: many options to construct initial matrix
}