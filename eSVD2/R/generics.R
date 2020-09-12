# Natural parameter: m_{ij} = u_i^Tv_j
# Data: a_{ij} ~ p(x|m_{ij})
#
# This function computes the negative log-likelihood value given the data
# `dat` may contain missing values
.evaluate_objective <- function(dat, u_mat, v_mat, ...) {
    stopifnot(
        is.matrix(dat),
        nrow(dat) == nrow(u_mat),
        ncol(dat) == nrow(v_mat),
        ncol(u_mat) == ncol(v_mat)
    )

    UseMethod(".evaluate_objective")
}

.evaluate_objective.default <- function(dat, u_mat, v_mat, ...) {
    .evaluate_objective.exponential(dat, u_mat, v_mat, ...)
}

##########

# This function computes the objective function on one row/column of data
# For example,
#   dat_vec = a_{i*} = (a_{i1}, ..., a_{ip})
#   current_vec = u_i
#   other_mat = V
# Or,
#   dat_vec = a_{*j} = (a_{1j}, ..., a_{nj})
#   current_vec = v_j
#   other_mat = U

.evaluate_objective_single <- function(dat_vec, current_vec, other_mat, ...) {
    stopifnot(
        !is.matrix(dat_vec),
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    UseMethod(".evaluate_objective_single")
}

.evaluate_objective_single.default <- function(dat_vec, current_vec, other_mat, ...) {
    .evaluate_objective_single.exponential(dat_vec, current_vec, other_mat, ...)
}

###########

# The gradient of .evaluate_objective_single with respect to current_vec
.gradient_vec <- function(dat_vec, current_vec, other_mat, ...) {
    stopifnot(
        !is.matrix(dat_vec),
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    UseMethod(".gradient_vec")
}

.gradient_vec.default <- function(dat_vec, current_vec, other_mat, ...) {
    .gradient_vec.exponential(dat_vec, current_vec, other_mat, ...)
}

###########

.evaluate_objective_mat <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)))

  UseMethod(".evaluate_objective_mat")
}

###########

#' Gradient of the objective function
#'
#' Computes the gradient for a particular model (based on the
#' class of \code{dat}) of the \code{.evaluate_objective_mat} function.
#'
#' Note, \code{dat} is NOT allowed to have any \code{NA} values for this
#' function.
#'
#' @param dat dataset where the \code{n} rows represent cells and \code{d} columns represent genes
#' @param nat_mat \code{n} by \code{d} matrix where each entry represents the natural parameter of the corresponding entry in \code{dat}
#' @param ... other parameters
#'
#' @return \code{n} by \code{d} matrix
.gradient_mat <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)))

  UseMethod(".gradient_mat")
}
