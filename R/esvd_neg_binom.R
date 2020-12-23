# Distribution: negative binomial

.evaluate_objective.neg_binom <- function(
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
    stopifnot(all(nat_mat < 0))

    # Only compute likelihood on non-missing data
    idx <- which(!is.na(dat))
    stopifnot(length(idx) > 0)

    if(length(idx) == prod(dim(dat))){
      negloglik <- .mult_mat_vec(.mult_vec_mat(-library_size_vec, log(1 - exp(nat_mat))), nuisance_param_vec) - nat_mat * dat
      sum(negloglik) / (n*p)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_mat[idx]
      # dat_vals <- dat[idx]
      # negloglik <- -scalar * log(1 - exp(nat_vals)) - nat_vals * dat_vals
      # sum(negloglik) / (n*p)
    }
}

.evaluate_objective_single.neg_binom <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat),
        ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
         | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      negloglik <- -library_size_vec * nuisance_param_vec * log(1 - exp(nat_vec)) - nat_vec * dat_vec
      sum(negloglik) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # dat_vals <- dat_vec[idx]
      # negloglik <- -scalar * log(1 - exp(nat_vals)) - nat_vals * dat_vals
      # sum(negloglik) / length(dat_vec)
    }
}


# f(x) = log(1 - exp(x))
# f'(x) = -exp(x) / (1 - exp(x)) = -1 / (exp(-x) - 1) = 1 / (1 - exp(-x))
# f''(x) = -exp(-x) / (1 - exp(-x))^2 = [1 / (1 - exp(-x))] * [-exp(-x) / (1 - exp(-x))]
#        = f'(x) * [1 - f'(x)]

.gradient_vec.neg_binom <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat),
        ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
         | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      grad <- other_mat * (library_size_vec * nuisance_param_vec / (exp(-nat_vec) - 1) - dat_vec)
      colSums(grad) / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # dat_vals <- dat_vec[idx]
      # other_vals <- other_mat[idx, , drop = FALSE]
      # grad <- other_vals * (scalar / (exp(-nat_vals) - 1) - dat_vals)
      # colSums(grad) / length(dat_vec)
    }

}

.hessian_vec.neg_binom <- function(
    current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
) {
    stopifnot(
        length(current_vec) == ncol(other_mat),
        length(dat_vec) == nrow(other_mat),
        ((length(library_size_vec) == 1 & length(nuisance_param_vec) == nrow(other_mat))
         | (length(library_size_vec) == nrow(other_mat) & length(nuisance_param_vec) == 1))
    )

    nat_vec <- c(other_mat %*% current_vec)
    stopifnot(all(nat_vec < 0))
    idx <- which(!is.na(dat_vec))
    stopifnot(length(idx) > 0)

    if(length(idx) == length(dat_vec)){
      tmp <- exp(nat_vec)
      diag_vec <- library_size_vec * nuisance_param_vec * tmp * (1/(1-tmp) + tmp/(1-tmp)^2)
      term1 <- crossprod(.mult_vec_mat(diag_vec, other_mat), other_mat)
      term1 / length(dat_vec)
    } else {
      ## Check below lines ##
      # nat_vals <- nat_vec[idx]
      # other_vals <- other_mat[idx, , drop = FALSE]
      #
      # fpx <- 1 / (1 - exp(nat_vals))
      # fppx <- fpx * (1 - fpx)
      # term1 <- t(other_vals) %*% diag(-scalar * fppx) %*% other_vals
      #
      # term1 / length(dat_vec)
    }
}

.feasibility.neg_binom <- function(current_vec, other_mat, ...) {
    nat_vec <- c(other_mat %*% current_vec)
    all(nat_vec < 0)
}

.neg_binom <- structure(
    list(
        objfn_all = .evaluate_objective.neg_binom,
        objfn     = .evaluate_objective_single.neg_binom,
        grad      = .gradient_vec.neg_binom,
        hessian   = .hessian_vec.neg_binom,
        feas      = .feasibility.neg_binom
    ),
    class = "esvd_family"
)
