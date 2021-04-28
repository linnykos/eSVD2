#' Estimating the nuisance parameters given an initial estimate of natural matrix
#'
#' @param dat dataset where the \eqn{n} rows represent cells and \eqn{p} columns represent genes
#' @param init_nat_mat initial estimate of the natural parameter matrix, of also size \eqn{n \times p}
#' @param family A character string, one of \code{"gaussian"}, \code{"exponential"},
#'               \code{"poisson"}, \code{"neg_binom"}, \code{"curved_gaussian"},
#'               and \code{"bernoulli"}.
#' @param library_size_vec either \code{NA} or a single numeric (default is \code{1}) or
#' a length-\eqn{n} vector of numerics.
#' If \code{NA}, the library size will be estimated.
#'
#' @return vector of length \eqn{p}, representing the nuisance parameters, one for each
#' of the \eqn{p} genes
#' @export
initialize_nuisance_param <- function(dat, init_nat_mat, family,
                                      library_size_vec){
  stopifnot(all(dim(dat) == dim(init_nat_mat)), all(dat >= 0))
  stopifnot((family == "curved_gaussian" && all(init_nat_mat < 0)) ||
              (family == "neg_binom"))

  previous_family <- ifelse(family == "curved_gaussian", "exponential", "poisson")

  mean_mat <- compute_mean(init_nat_mat, family = previous_family,
                           library_size_vec = library_size_vec)

  param_vec <- sapply(1:ncol(dat), function(j){
    empirical_var_vec <- (dat[,j] - mean_mat[,j])^2

    if(family == "curved_gaussian"){
      tmp <- sqrt(mean_mat[,j]^2/(library_size_vec * empirical_var_vec))
      tmp <- 1/tmp^2
      min_val <- max(min(tmp), 0); max_val <- max(max(tmp), 10)

      val <- tryCatch(stats::uniroot(.root_curved_gaussian_closure(empirical_var_vec, mean_mat[,j], library_size_vec),
                              interval = c(min_val, max_val))$root,
               error = function(e){
                 1
               })
    } else {
      tmp <- mean_mat[,j]^2/(library_size_vec * (empirical_var_vec - mean_mat[,j]) -  mean_mat[,j])
      tmp <- 1/tmp
      min_val <- max(min(tmp), 0); max_val <- max(max(tmp), 10)

      tryCatch(stats::uniroot(.root_neg_binom_closure(empirical_var_vec, mean_mat[,j], library_size_vec),
                              interval = c(min_val, max_val))$root,
               error = function(e){
                 1/(10*max(dat[,j]))
               })
    }
  })

  if(family == "neg_binom") param_vec <- 1/param_vec else param_vec <- 1/sqrt(param_vec)

  param_vec
}

#########################

.root_curved_gaussian_closure <- function(empirical_var_vec, theoretical_mean_vec, library_size_vec){
  function(x){
    sum((empirical_var_vec - x*theoretical_mean_vec^2/library_size_vec) / (x^2*theoretical_mean_vec^2/library_size_vec))
  }
}

.root_neg_binom_closure <- function(empirical_var_vec, theoretical_mean_vec, library_size_vec){
  function(x){
    sum((empirical_var_vec - theoretical_mean_vec * (1 + x*theoretical_mean_vec/library_size_vec)) /
          2*library_size_vec*(1 + x*theoretical_mean_vec/library_size_vec)^2)
  }
}
