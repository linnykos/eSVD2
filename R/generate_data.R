#' Generate data
#'
#' @param nat_mat An \eqn{n\times p} matrix of natural parameters, where
#'                \eqn{n} rows represent cells and \eqn{p} columns represent genes.
#' @param family A character string, one of \code{"gaussian"}, \code{"exponential"},
#'               \code{"poisson"}, \code{"neg_binom"}, \code{"curved_gaussian"},
#'               and \code{"bernoulli"}.
#' @param nuisance_param_vec Either \code{NA} or a single numeric or a length-\eqn{p}
#'                           vector of numerics representing nuisance parameters
#'                           (for \code{family = "neg_binom"} and
#'                           \code{family = "curved_gausian"}).
#'                           It is only required if
#'                           \code{family \%in\% c("neg_binom", "curved_gaussian")}.
#' @param library_size_vec Either \code{NA} or a length-\eqn{n} vector of numerics
#' @param tol Small positive value to determine the smallest possible value in the output
#'            matrix, useful for only \code{family = "curved_gaussian"}.
#'
#' @return The generated data matrix
#' @export
generate_data <- function(
    nat_mat, family, nuisance_param_vec = NA, library_size_vec = 1, tol = 1e-3
) {
  family <- .string_to_distr_funcs(family)
  stopifnot(
    is.matrix(nat_mat),
    length(nuisance_param_vec) %in% c(1, ncol(nat_mat)),
    family$feasibility(nat_mat)
  )

  n <- nrow(nat_mat)
  if(length(library_size_vec) != 1) {
    stopifnot(
      length(library_size_vec) == n,
      all(!is.na(library_size_vec))
    )
  } else {
    library_size_vec <- rep(library_size_vec[1], n)
  }

  # library_size_vec is now a length-n vector
  dat <- .generate_values(nat_mat, family, nuisance_param_vec, library_size_vec)
  dim(dat) <- dim(nat_mat)

  if(family$name %in% c("curved_gaussian", "gaussian") && !is.na(tol))
    dat <- pmax(dat, tol)

  dat
}

.generate_values <- function(nat_mat, family, nuisance_param_vec, library_size_vec) {

  n <- nrow(nat_mat)
  p <- ncol(nat_mat)
  num_val <- length(nat_mat)

  stopifnot(length(nuisance_param_vec) %in% c(1, p))
  if(length(nuisance_param_vec) == 1)
    nuisance_param_vec <- rep(nuisance_param_vec[1], p)
  nuisance_param_na <- any(is.na(nuisance_param_vec))

  stopifnot(length(library_size_vec) == n)

  canon_mat <- family$nat_to_canon(nat_mat)

  if(family$name == "gaussian") {
    stopifnot(!nuisance_param_na)
    vec <- stats::rnorm(
      num_val, mean = canon_mat * library_size_vec,
      sd = rep(nuisance_param_vec, each = n) * sqrt(library_size_vec)
    )

  } else if(family$name == "curved_gaussian") {
    stopifnot(!nuisance_param_na)
    vec <-stats::rnorm(
      num_val, mean = canon_mat * library_size_vec,
      sd = sqrt(library_size_vec) * (canon_mat / rep(nuisance_param_vec, each = n))
    )

  } else if(family$name == "exponential") {
    vec <- stats::rgamma(num_val, shape = rep(library_size_vec, times = p), scale = canon_mat)

  } else if(family$name == "poisson") {
    vec <- stats::rpois(num_val, lambda = canon_mat * library_size_vec)

  } else if(family$name == "neg_binom") {
    stopifnot(!nuisance_param_na)
    vec <- stats::rnbinom(
      num_val, size = rep(nuisance_param_vec, each = n) * library_size_vec,
      prob = 1 - canon_mat
    )

  } else if(family$name == "neg_binom2") {
    stopifnot(!nuisance_param_na)
    vec <- stats::rnbinom(
      num_val, size = rep(nuisance_param_vec, each = n) * library_size_vec,
      mu = canon_mat
    )

  } else if(family$name == "bernoulli") {
    vec <- stats::rbinom(num_val, size = 1, prob = canon_mat)

  } else {
    stop("unknown distribution family")
  }

  vec
}
