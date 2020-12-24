# Distribution: Gaussian

.evaluate_objective.gaussian <- function(
  dat, x_mat, y_mat, nuisance_param_vec, library_size_vec, ...
) {
  # Check dimensions
  n <- nrow(dat)
  p <- ncol(dat)
  stopifnot(
    length(nuisance_param_vec) == p,
    length(library_size_vec) == n,
    ncol(x_mat) == ncol(y_mat),
    nrow(x_mat) == n,
    nrow(y_mat) == p
  )

  # Compute natural parameters
  nat_mat <- tcrossprod(x_mat, y_mat)

  # Only compute likelihood on non-missing data
  idx <- which(!is.na(dat))
  stopifnot(length(idx) > 0)

  if(length(idx) == prod(dim(dat))){
    dat_vals <- (dat / library_size_vec)
    negloglik <- .mult_mat_vec((nat_mat - dat_vals)^2, 1/nuisance_param_vec^2)
    sum(negloglik) /  (n*p)
  } else {
    ## Check below lines ##
    # nat_vals <- nat_mat[idx]
    # dat_vals <- (dat / library_size_vec)[idx]
    # nuisance_vals <- rep(nuisance_param_vec, rep(nrow(dat), length(nuisance_param_vec)))[idx]
    # negloglik <- (nat_vals - dat_vals)^2/nuisance_vals^2
    # sum(negloglik) /  (n*p)
  }
}

# length(library_size_vec) == 1 if current vector is u
# length(library_size_vec) == n if current vector is v
.evaluate_objective_single.gaussian <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
      length(current_vec) == ncol(other_mat),
      length(dat_vec) == nrow(other_mat),
      ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
       | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
  )

  nat_vec <- c(other_mat %*% current_vec)
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  if(length(idx) == length(dat_vec)){
    dat_vals <- (dat_vec / library_size_vec)
    negloglik <- (nat_vec - dat_vals)^2/nuisance_param_vec^2
    sum(negloglik) / length(dat_vec)
  } else {
    ## Check below lines ##
    # nat_vals <- nat_vec[idx]
    # dat_vals <- (dat_vec / library_size_vec)[idx]
    # nuisance_vals <- nuisance_param_vec[idx]
    # negloglik <- (nat_vals - dat_vals)^2/nuisance_vals^2
    # sum(negloglik) / length(dat_vec)
  }
}

.gradient_vec.gaussian <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat),
    ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
     | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
  )

  nat_vec <- c(other_mat %*% current_vec)
  idx <- which(!is.na(dat_vec))
  stopifnot(length(idx) > 0)

  if(length(idx) == length(dat_vec)){
    dat_vals <- (dat_vec / library_size_vec)
    other_vals <- other_mat
    # Broadcast (nat_vals - dat_vals) to each column of other_vals
    grad <- 2 * other_vals * (nat_vec - dat_vals)/nuisance_param_vec^2
    colSums(grad) / length(dat_vec)
  } else {
    ## Check below lines ##
    # nat_vals <- nat_vec[idx]
    # dat_vals <- (dat_vec / library_size_vec)[idx]
    # other_vals <- other_mat[idx, , drop = FALSE]
    # # Broadcast (nat_vals - dat_vals) to each column of other_vals
    # grad <- 2 * other_vals * (nat_vals - dat_vals)
    # colSums(grad) / length(dat_vec)
  }
}

.hessian_vec.gaussian <- function(
  current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
  stopifnot(
    length(current_vec) == ncol(other_mat),
    length(dat_vec) == nrow(other_mat),
    ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
      | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
  )

  idx <- which(!is.na(dat_vec))

  if(length(idx) == length(dat_vec)){
    if(length(nuisance_param_vec) == 1){
      term1 <- 2 * crossprod(other_mat/nuisance_param_vec)
    } else {
      term1 <- 2 * crossprod(.mult_vec_mat(1/nuisance_param_vec, other_mat))
    }
    term1 / length(dat_vec)
  } else {
    ## Check below lines ##
    # other_vals <- other_mat[idx, , drop = FALSE]
    # term1 <- 2 * crossprod(other_vals)
    # term1 / length(dat_vec)
  }
}

.feasibility.gaussian <- function(current_vec, other_mat, ...) {
  TRUE
}

.gaussian <- structure(
  list(
    objfn_all = .evaluate_objective.gaussian,
    objfn     = .evaluate_objective_single.gaussian,
    grad      = .gradient_vec.gaussian,
    hessian   = .hessian_vec.gaussian,
    feas      = .feasibility.gaussian
  ),
  class = "esvd_family"
)
