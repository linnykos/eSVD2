# Distribution: Poisson

.evaluate_objective.poisson <- function(
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

    # Only compute likelihood on non-missing data
    idx <- which(!is.na(dat))
    stopifnot(length(idx) > 0)

    if(length(idx) == prod(dim(dat))){
      negloglik <- exp(nat_mat + log(library_size_vec)) - nat_mat * dat
      sum(negloglik) / (n*p)
    } else {
      ## Check below lines ##
      # log_library_size_vals <- rep(log(library_size_vec), p)[idx]
      # nat_vals <- nat_mat[idx] + log_library_size_vals
      # dat_vals <- dat[idx]
      # negloglik <- exp(nat_vals) - nat_vals * dat_vals
      # sum(negloglik) / (n*p)
    }
}

# length(library_size_vec) == 1 if current vector is u
# length(library_size_vec) == n if current vector is v
.evaluate_objective_single.poisson <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec = NA, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      negloglik <- exp(nat_vec + log(library_size_vec)) - nat_vec * dat_vec
      sum(negloglik) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- (nat_vec + log(library_size_vec))[idx]
      # dat_vals <- dat_vec[idx]
      # negloglik <- exp(nat_vals) - nat_vals * dat_vals
      # sum(negloglik) / length(dat_vec)
    }
}

.gradient_vec.poisson <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec = NA, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      grad <- other_mat * (exp(nat_vec + log(library_size_vec)) - dat_vec)
      colSums(grad) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- (nat_vec + log(library_size_vec))[idx]
      # dat_vals <- dat_vec[idx]
      # other_vals <- other_mat[idx, , drop = FALSE]
      # grad <- other_vals * (exp(nat_vals) - dat_vals)
      #
      # colSums(grad) / length(dat_vec)
    }
}

.hessian_vec.poisson <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec = NA, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat)
    )

    nat_vec <- c(other_mat %*% current_vec)
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      term1 <- crossprod(.mult_vec_mat(exp(nat_vec + log(library_size_vec)), other_mat), other_mat)

      term1 / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- (nat_vec + log(library_size_vec))[idx]
      # dat_vals <- dat_vec[idx]
      # other_vals <- other_mat[idx, , drop = FALSE]
      #
      # term1 <- t(other_vals) %*% diag(exp(nat_vals)) %*% other_vals
      #
      # term1 / length(dat_vec)
    }
}

.feasibility.poisson <- function(current_vec, other_mat, ...) {
    TRUE
}

.poisson <- structure(
    list(
        objfn_all = .evaluate_objective.poisson,
        objfn     = .evaluate_objective_single.poisson,
        grad      = .gradient_vec.poisson,
        hessian   = .hessian_vec.poisson,
        feas      = .feasibility.poisson
    ),
    class = "esvd_family"
)
