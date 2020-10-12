#' Compute mean from natural parameter matrix
#'
#' @param nat_mat matrix
#' @param family character
#' @param ... additional parameters for distribution
#'
#' @return matrix
#' @export
compute_mean <- function(nat_mat, family, ...){
 if(family == "gaussian") {
  return(nat_mat)
 } else if(family == "poisson"){
  stopifnot(all(nat_mat > 0))
  exp(nat_mat)
 } else if(family == "neg_binom"){
  stopifnot(all(nat_mat < 0))
  .compute_mean_neg_binom(nat_mat, ...)
 } else if(family == "exponential"){
  stopifnot(all(nat_mat < 0))
  -1/nat_mat
 } else if(family == "curved_gaussian"){
  stopifnot(all(nat_mat > 0))
  1/nat_mat
 } else {
  stop("family not found")
 }
}

######

.convert_natural_to_canonical <- function(nat_mat, family) {
    if(family == "gaussian")        return(nat_mat)
    if(family == "curved_gaussian") return(1 / nat_mat)
    if(family == "exponential")     return(-nat_mat)
    if(family == "bernoulli")       return(plogis(nat_mat))
    if(family %in% c("poisson", "neg_binom")) {
        if(any(nat_mat >= 10)) warning("Potential large values generated")
        return(exp(nat_mat))
    }
    stop("unknown distribution family")
}

.compute_mean_neg_binom <- function(nat_mat, scalar, ...){
 if(is.na(scalar)) stop("No argument scalar provided for negative binomial")
 stopifnot(length(scalar) == 1)

 scalar*exp(nat_mat)/(1-exp(nat_mat))
}

.determine_domain <- function(family, tol = 1e-3){
 if(family %in% c("exponential", "neg_binom")) {
  domain <- c(-Inf, -tol)
 } else if(family %in% c("poisson", "curved_gaussian")) {
  domain <-  c(tol, Inf)
 } else if(family == "gaussian") {
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

.mean_transformation <- function(dat, family, tol = 1e-3, ...){
 if(family == "exponential"){
  dat <- -1/(dat+1)
 } else if(family == "curved_gaussian"){
  dat <- 1/(dat+1)
 } else if(family == "poisson"){
  dat <- log(dat + tol)
 } else if(family == "neg_binom"){
  dat <- .mean_transformation_neg_binom(dat, tol, ...)
 } else if(family != "gaussian") {
  stop("family not found")
 }

 dat
}

.mean_transformation_neg_binom <- function(dat, tol, scalar = NA, ...){
 if(is.na(scalar)) stop("No argument scalar provided for negative binomial")
 stopifnot(length(scalar) == 1)

 dat_new <- (dat + tol)/scalar
 log(dat_new / (1+dat_new))
}
