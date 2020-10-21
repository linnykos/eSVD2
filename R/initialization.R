#' Initialize eSVD
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param k positive integer less than \code{min(c(nrow(dat), ncol(dat)))}
#' @param family character (\code{"gaussian"}, \code{"exponential"}, \code{"poisson"}, \code{"neg_binom"},
#' or \code{"curved_gaussian"})
#' @param nuisance_param_vec either \code{NA} or a single numeric or a length-\eqn{p}
#' vector of numerics representing nuisance parameters (for \code{family = "neg_binom"} or
#' \code{family = "curved_gausian"}).
#' @param library_size_vec either \code{NA} or a single numeric or
#' a length-\eqn{n} vector of numerics or \code{NULL}.
#' If \code{NA}, the library size will be estimated.
#' @param config additional parameters for the initialization, whose defaults can be
#' set with \code{eSVD2::initialization_default()}
#' @param verbose non-negative integer specifying level of printouts
#'
#' @return a list with elements \code{x_mat} and \code{y_mat}, representing the two
#' latent matrices
#' @export
initialize_esvd <- function(dat, k, family, nuisance_param_vec = NA, library_size_vec = NA,
                            config = initialization_options(), verbose = 0){
  if(family != "gaussian") stopifnot(all(dat[!is.na(dat)] >= 0))
  if(all(!is.na(nuisance_param_vec)) & length(nuisance_param_vec) == 1) {
    nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
  }

  # estimate library sizes if asked
  library_size_vec <- .estimate_library_size(dat, library_size_vec = library_size_vec,
                                        config = config)

  rescaled_dat <- t(sapply(1:nrow(dat), function(i){dat[i,]/library_size_vec[i]}))

  # determine initial matrix taking into account to missing values and library size
  rescaled_dat <- .matrix_completion(rescaled_dat, k = k)
  init_res <- .determine_initial_matrix(rescaled_dat, k = k, family = family,
                                        nuisance_param_vec = nuisance_param_vec,
                                        max_val = config$max_val,
                                        tol = config$tol)

  # project inital matrix into space of low-rank matrices
  nat_mat <- .initialize_nat_mat(init_res, k = k, config = config)

  # reparameterize
  res <- .factorize_matrix(nat_mat, k = k, equal_covariance = T)
  res <- .fix_rank_defficiency(res$x_mat, res$y_mat, domain = init_res$domain)

  structure(list(x_mat = res$x_mat, y_mat = res$y_mat,
                 library_size_vec = library_size_vec, nuisance_param_vec = nuisance_param_vec,
                 domain = init_res$domain), class = "eSVD")
}

#' Initialization defaults
#'
#' @param init_method character (\code{"kmean_rows"})
#' @param library_size_method string such as \code{"total_read"}, specificying how the library size of a cell will be calculated
#' @param max_val maximum magnitude of the inner product (positive numeric). This parameter
#' could be \code{NA}.
#' @param tol small positive value, which also controls (amongst other numerical things)
#' the smallest magnitude of the inner product
#'
#' @return a list of values, of class \code{initialization_param}
#' @export
initialization_options <- function(init_method = "kmean_rows",
                                 library_size_method = "total_read",
                                 max_val = NA, tol = 1e-3){
  stopifnot(init_method %in% c("kmean_rows"))
  stopifnot(library_size_method %in% c("total_read"))
  stopifnot(tol > 0, tol <= 1, (is.na(max_val) | max_val > 0))

  structure(list(init_method = init_method,
                 library_size_method = library_size_method,
                 max_val = max_val, tol = tol), class = "initialization_param")
}


################

.estimate_library_size <- function(dat, library_size_vec, config){
  if(all(!is.na(library_size_vec))){
    if(length(library_size_vec) == 1) library_size_vec <- rep(1, nrow(dat))
    stopifnot(length(library_size_vec) == nrow(dat))

  } else {
    if(config$library_size_method == "total_read"){
      library_size_vec <- rowSums(dat, na.rm = T)
      library_size_vec <- library_size_vec/min(library_size_vec)
    } else {
      stop("config library_size_method not found")
    }
  }

  library_size_vec
}

###################

.initialize_nat_mat <- function(init_res, k = k, config = config){
  if(config$init_method == "kmean_rows"){
    nat_mat <- .initialization_kmean(init_res$nat_mat, k = k, domain = init_res$domain,
                                     row = T)
  } else {
    stop("config init_method not found")
  }

  nat_mat
}

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
#' @param nuisance_param_vec length-\eqn{p}
#' vector of numerics representing nuisance parameters (for \code{family = "neg_binom"} or
#' \code{family = "curved_gausian"})
#' @param max_val maximum magnitude of the inner product
#' @param tol numeric
#' or \code{"curved gaussian"} for \code{family}
#'
#' @return \code{n} by \code{p} matrix
.determine_initial_matrix <- function(dat, k, family, nuisance_param_vec, max_val = NA, tol = 1e-3){
  stopifnot((is.na(max_val) || max_val >= 0), all(dat >= 0))

  domain <- .determine_domain(family, tol)
  if(!is.na(max_val)) domain <- .intersect_intervals(domain, c(-max_val, max_val))

  if(family != "bernoulli") dat[which(dat <= tol)] <- tol/2
  nat_mat <- .mean_transformation(dat, family, nuisance_param_vec = nuisance_param_vec)
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
#' Here, \code{domain} represents the domain where all the inner products in
#' \code{x_mat \%*\% t(y_mat)} are assumed to lie between
#'
#' @param x_mat a numeric matrix
#' @param y_mat a numeric matrix with the same number of columns as \code{x_mat}
#' @param domain a vector where \code{domain[1] < domain[2]}
#'
#' @return a list of \code{x_mat} and \code{y_mat}
.fix_rank_defficiency <- function(x_mat, y_mat, domain){
  k <- ncol(x_mat)
  nat_mat <- x_mat %*% t(y_mat)
  stopifnot(.check_domain(nat_mat, domain))
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

  # fix any remaining issue with lying within domain
  nat_mat <- x_mat %*% t(y_mat)
  mag <- max(abs(domain))
  if(max(abs(nat_mat)) > mag){
    nat_mat <- nat_mat * (mag/max(abs(nat_mat)))
    res <- .factorize_matrix(nat_mat, k = k, equal_covariance = T)
  } else {
    res <- .reparameterize(x_mat, y_mat, equal_covariance = T)
  }

  list(x_mat = res$x_mat, y_mat = res$y_mat)
}
