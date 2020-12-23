# Distribution: exponential

.evaluate_objective.exponential <- function(
    dat, x_mat, y_mat, nuisance_param_vec = NA, library_size_vec, ...
) {
    # Check dimensions
    n <- nrow(dat)
    p <- ncol(dat)
    stopifnot(
        length(library_size_vec) == n,
        ncol(x_mat) == ncol(y_mat),
        nrow(x_mat) == n,
        nrow(y_mat) == p
    )

    # Compute natural parameters
    nat_mat <- tcrossprod(x_mat, y_mat)
    stopifnot(all(nat_mat < 0))

    # Only compute likelihood on non-missing data
    idx <- which(!is.na(dat))
    stopifnot(length(idx) > 0)

    if(length(idx) == prod(dim(dat))){
      negloglik <- -library_size_vec * log(-nat_mat) - dat * nat_mat
      sum(negloglik) / (n*p)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_mat[idx]
      # dat_vals <- (dat / library_size_vec)[idx]
      # negloglik <- -log(-nat_vals) - nat_vals * dat_vals
      # sum(negloglik) / (n*p)
    }
}

# length(library_size_vec) == 1 if current vector is u
# length(library_size_vec) == n if current vector is v
.evaluate_objective_single.exponential <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec = NA, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      negloglik <- -library_size_vec * log(-nat_vec) - dat_vec * nat_vec
      sum(negloglik) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # dat_vals <- (dat_vec / library_size_vec)[idx]
      # negloglik <- -log(-nat_vals) - nat_vals * dat_vals
      # sum(negloglik) / length(dat_vec)
    }
}

.gradient_vec.exponential <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec = NA, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      grad <- other_mat * (-library_size_vec/nat_vec - dat_vec)

      colSums(grad) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # dat_vals <- (dat_vec / library_size_vec)[idx]
      # other_vals <- other_mat[idx, , drop = FALSE]
      # grad <- other_vals * (-1 / nat_vals - dat_vals)
      #
      # colSums(grad) / length(dat_vec)
    }


}

.hessian_vec.exponential <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec = NA, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      term1 <-  crossprod(.mult_vec_mat(library_size_vec/nat_vec^2, other_mat), other_mat)
      term1 / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # dat_vals <- (dat_vec / library_size_vec)[idx]
      # other_vals <- other_mat[idx, , drop = FALSE]
      #
      # term1 <- t(other_vals) %*% diag(1 / nat_vals^2) %*% other_vals
      #
      # term1 / length(dat_vec)
    }
}

.feasibility.exponential <- function(current_vec, other_mat, ...) {
    nat_vec <- c(other_mat %*% current_vec)
    all(nat_vec < 0)
}

.exponential <- structure(
    list(
        objfn_all = .evaluate_objective.exponential,
        objfn     = .evaluate_objective_single.exponential,
        grad      = .gradient_vec.exponential,
        hessian   = .hessian_vec.exponential,
        feas      = .feasibility.exponential
    ),
    class = "esvd_family"
)
