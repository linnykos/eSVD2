# Distribution: one-parameter Gaussian where sd = mean/scalar

.evaluate_objective.curved_gaussian <- function(
    dat, x_mat, y_mat, nuisance_param_vec, library_size_vec, ...
) {
    # Check dimensions
    n <- nrow(dat)
    p <- ncol(dat)
    stopifnot(
        length(nuisance_param_vec) == p, length(library_size_vec) == n,
        ncol(x_mat) == ncol(y_mat),
        nrow(x_mat) == n,
        nrow(y_mat) == p
    )

    # Compute natural parameters
    nat_mat <- tcrossprod(x_mat, y_mat)
    stopifnot(all(nat_mat > 0))

    # Only compute likelihood on non-missing data
    idx <- which(!is.na(dat))
    stopifnot(length(idx) > 0)

    if(length(idx) == prod(dim(dat))){
      negloglik <- -log(nat_mat) +
        .mult_vec_mat(1/(2*library_size_vec), .mult_mat_vec(dat^2, nuisance_param_vec^2)) * nat_mat^2 -
        .mult_mat_vec(dat, nuisance_param_vec^2) * nat_mat
      sum(negloglik) / (n*p)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_mat[idx]
      # dat_vals <- (dat / library_size_vec)[idx]
      # negloglik <- -log(nat_vals) +
      #   scalar^2 * dat_vals^2 * nat_vals^2 / 2 -
      #   scalar^2 * dat_vals * nat_vals
      # sum(negloglik) / (n*p)
    }

}

# length(library_size_vec) == 1 if current vector is u
# length(library_size_vec) == n if current vector is v
.evaluate_objective_single.curved_gaussian <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat),
        ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
         | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec > 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      negloglik <- -log(nat_vec) +
        nuisance_param_vec^2 * dat_vec^2 * nat_vec^2 / (2*library_size_vec) -
        nuisance_param_vec^2 * dat_vec * nat_vec
      sum(negloglik) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # dat_vals <- (dat_vec / library_size_vec)[idx]
      # negloglik <- -log(nat_vals) +
      #   scalar^2 * dat_vals^2 * nat_vals^2 / 2 -
      #   scalar^2 * dat_vals * nat_vals
      # sum(negloglik) / length(dat_vec)
    }
}

.gradient_vec.curved_gaussian <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat),
        ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
         | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec > 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      grad <- other_mat * (-1 / nat_vec + nuisance_param_vec^2 * dat_vec^2 * nat_vec / library_size_vec -
                              nuisance_param_vec^2 * dat_vec)
      colSums(grad) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # dat_vals <- (dat_vec / library_size_vec)[idx]
      # other_vals <- other_mat[idx, , drop = FALSE]
      # grad <- other_vals * (-1 / nat_vals + scalar^2 * dat_vals^2 * nat_vals -
      #                         scalar^2 * dat_vals)
      # colSums(grad) / length(dat_vec)
    }

}

.hessian_vec.curved_gaussian <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat),
        ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
         | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec > 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      term1 <- crossprod(.mult_vec_mat(1/nat_vec^2, other_mat), other_mat)
      term2 <- crossprod(.mult_vec_mat(nuisance_param_vec^2*dat_vec^2/library_size_vec, other_mat), other_mat)

      (term1 + term2) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # dat_vals <- (dat_vec / library_size_vec)[idx]
      # other_vals <- other_mat[idx, , drop = FALSE]
      #
      # term1 <- t(other_vals) %*% diag(1 / nat_vals^2) %*% other_vals
      # term2 <- scalar^2 * t(other_vals) %*% diag(dat_vals^2) %*% other_vals
      #
      # (term1 + term2) / length(dat_vec)
    }

}

.feasibility.curved_gaussian <- function(current_vec, other_mat, ...) {
    nat_vec <- c(other_mat %*% current_vec)
    all(nat_vec > 0)
}

.curved_gaussian <- structure(
    list(
        objfn_all = .evaluate_objective.curved_gaussian,
        objfn     = .evaluate_objective_single.curved_gaussian,
        grad      = .gradient_vec.curved_gaussian,
        hessian   = .hessian_vec.curved_gaussian,
        feas      = .feasibility.curved_gaussian
    ),
    class = "esvd_family"
)

