# natural parameter -> mean

#' Compute mean from natural parameter matrix
#'
#' @param nat_mat matrix
#' @param family character
#' @param nuisance_param_vec either \code{NA} or a single numeric or a length-\eqn{p}
#' vector of numerics representing nuisance parameters (for \code{family = "neg_binom"} or
#' \code{family = "curved_gausian"}).
#' @param library_size_vec either \code{NA} or a single numeric or
#' a length-\eqn{n} vector of numerics or \code{NULL}.
#'
#' @return matrix
#' @export
compute_mean <- function(nat_mat, family, nuisance_param_vec = NA,
                         library_size_vec = rep(1, nrow(nat_mat))){
  if(all(!is.na(nuisance_param_vec)) & length(nuisance_param_vec) == 1) {
    nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(nat_mat))
  }

  if(family == "gaussian") {
    nat_mat

  } else if(family == "poisson"){
    library_size_vec * exp(nat_mat)

  } else if(family == "neg_binom"){
    stopifnot(all(nat_mat < 0))
    .compute_mean_neg_binom(nat_mat, nuisance_param_vec = nuisance_param_vec,
                            library_size_vec = library_size_vec)

  } else if(family == "neg_binom2"){
    library_size_vec * exp(nat_mat)

  } else if(family == "exponential"){
    stopifnot(all(nat_mat < 0))
    library_size_vec * (-1/nat_mat)

  } else if(family == "curved_gaussian"){
    stopifnot(all(nat_mat > 0))
    library_size_vec * (1/nat_mat)

  } else if(family == "bernoulli"){
    stats::plogis(nat_mat)

  } else {
    stop("family not found")
  }
}

# Compute mean for negative binomial
.compute_mean_neg_binom <- function(nat_mat, nuisance_param_vec = NA,
                                    library_size_vec = rep(1, nrow(nat_mat))){
  stopifnot(length(nuisance_param_vec) == ncol(nat_mat))

  res <- exp(nat_mat)/(1-exp(nat_mat))
  res <- library_size_vec * (res * rep(nuisance_param_vec, each = nrow(nat_mat)))

  colnames(res) <- colnames(nat_mat)
  rownames(res) <- rownames(nat_mat)

  res
}

# Properly handle the library size parameter
.parse_library_size <- function(dat, library_size_vec) {
  n <- nrow(dat)
  stopifnot(length(library_size_vec) %in% c(1, n))

  if(any(is.na(library_size_vec))) {
    library_size_vec <- rowSums(dat)
  } else if(length(library_size_vec) == 1) {
    library_size_vec <- rep(library_size_vec, n)
  }

  library_size_vec / min(library_size_vec)
}
