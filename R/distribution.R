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

.compute_mean_neg_binom <- function(nat_mat, nuisance_param_vec = NA,
                                    library_size_vec = rep(1, nrow(nat_mat))){
  stopifnot(length(nuisance_param_vec) == ncol(nat_mat))

  res <- exp(nat_mat)/(1-exp(nat_mat))
  res <- library_size_vec * (res * rep(nuisance_param_vec, each = nrow(nat_mat)))

  colnames(res) <- colnames(nat_mat)
  rownames(res) <- rownames(nat_mat)

  res
}

##############################
##############################

# natural parameter -> variance
.compute_variance <- function(nat_mat, family, nuisance_param_vec = NA,
                              library_size_vec = rep(1, nrow(nat_mat))){
  stopifnot(nrow(nat_mat) == length(library_size_vec))

  if(family == "gaussian") {
    stopifnot(ncol(nat_mat) == length(nuisance_param_vec))
    matrix(rep(library_size_vec, times = ncol(nat_mat)) * rep(nuisance_param_vec^2, each = nrow(nat_mat)),
           nrow = nrow(nat_mat), ncol = ncol(nat_mat))

  } else if(family == "poisson"){
    library_size_vec * exp(nat_mat)

  } else if(family == "neg_binom"){
    stopifnot(all(nat_mat < 0))
    .compute_variance_neg_binom(nat_mat, nuisance_param_vec = nuisance_param_vec,
                                library_size_vec = library_size_vec)

  } else if(family == "exponential"){
    stopifnot(all(nat_mat < 0))
    library_size_vec * (-1/nat_mat^2)

  } else if(family == "curved_gaussian"){
    stopifnot(all(nat_mat > 0))
    stopifnot(length(nuisance_param_vec) == ncol(nat_mat))

    library_size_vec * ((1/nat_mat)^2 * rep(1/nuisance_param_vec^2, each = nrow(nat_mat)))

  } else if(family == "bernoulli"){
    tmp <- stats::plogis(nat_mat)
    tmp * (1-tmp)

  } else {
    stop("family not found")
  }
}

.compute_variance_neg_binom <- function(nat_mat, nuisance_param_vec = NA,
                                        library_size_vec = rep(1, nrow(nat_mat))){
  stopifnot(length(nuisance_param_vec) == ncol(nat_mat))

  res <- sapply(1:ncol(nat_mat), function(i){
    exp(nat_mat[,i])/(1-exp(nat_mat[,i]))^2
  })
  res <- res * (rep(nuisance_param_vec, each = nrow(nat_mat)) * rep(library_size_vec, times = ncol(nat_mat)))

  colnames(res) <- colnames(nat_mat)
  rownames(res) <- rownames(nat_mat)

  res
}

##############################
##############################

# natural parameter -> canonical (for data generation convenience)

.convert_natural_to_canonical <- function(nat_mat, family) {
  if(family == "gaussian")        return(nat_mat)
  if(family == "curved_gaussian") return(1 / nat_mat)
  if(family == "exponential")     return(-nat_mat)
  if(family == "bernoulli")       return(stats::plogis(nat_mat))
  if(family %in% c("poisson", "neg_binom")) {
    if(any(nat_mat >= 10)) warning("Potential large values generated")
    return(exp(nat_mat))
  }
  stop("unknown distribution family")
}

.determine_domain <- function(family, tol = 1e-3){
 if(family %in% c("exponential", "neg_binom")) {
  domain <- c(-Inf, -tol)
 } else if(family %in% c("curved_gaussian")) {
  domain <-  c(tol, Inf)
 } else if(family %in% c("poisson", "gaussian", "bernoulli")) {
  domain <- c(-Inf, Inf)
 }

 domain
}

.check_domain <- function(nat_mat, domain){
  idx <- which.min(abs(domain))

  domain_modified <- domain
  domain_modified[idx] <- sign(domain_modified[idx])*(abs(domain_modified[idx])-.Machine$double.eps*1e2)
  if(abs(domain_modified[-idx]) != Inf){
    domain_modified[-idx] <- sign(domain_modified[-idx])*(abs(domain_modified[-idx])+.Machine$double.eps*1e2)
  }

  all(nat_mat >= domain_modified[1]) & all(nat_mat <= domain_modified[2])
}

##############################
##############################

# observed-value -> natural parameter (for initialization)

.mean_transformation <- function(dat, family, tol = 1e-3, ...){
  if(family == "exponential"){
    dat <- -1/(dat+1)
  } else if(family == "curved_gaussian"){
    dat <- 1/(dat+1)
  } else if(family == "poisson"){
    dat <- log(dat + tol)
  } else if(family == "neg_binom"){
    dat <- .mean_transformation_neg_binom(dat, tol, ...)
  } else if(family == "bernoulli"){
    stopifnot(all(dat[!is.na(dat)] %in% c(0,1)))
    dat[dat == 0] <- -1
  } else if(family != "gaussian") {
    stop("family not found")
  }

  dat
}

.mean_transformation_neg_binom <- function(dat, tol, nuisance_param_vec = NA, ...){
 if(any(is.na(nuisance_param_vec))) stop("No argument nuisance_param_vec provided for negative binomial")
 stopifnot(length(nuisance_param_vec) == ncol(dat))

 res <- sapply(1:ncol(dat), function(i){(dat[,i] + tol)/nuisance_param_vec[i]})
 res <- log(res/(1+res))

 colnames(res) <- colnames(dat)
 rownames(res) <- rownames(dat)

 res
}

.string_to_distr_funcs <- function(family){
  if(family == "exponential"){
    return(.exponential)
  } else if(family == "curved_gaussian"){
    return(.curved_gaussian)
  } else if(family == "poisson"){
    return(.poisson)
  } else if(family == "neg_binom"){
    return(.neg_binom)
  } else if(family == "bernoulli"){
    return(.bernoulli)
  } else if(family == "gaussian") {
    return(.gaussian)
  } else {
    stop("family not found")
  }
}
