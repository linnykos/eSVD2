# distribution: one-parameter Gaussian where var = mean/scalar
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: m_{ij} = 1/mu_{ij}
# optimization problem: -log(m_{ij}) - scalar^2*a_{ij}^2*(-m_{ij}^2)/n - scalar^2*a_{ij}*m_{ij}

.evaluate_objective.curved_gaussian <- function(dat, u_mat, v_mat, scalar = 2, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(u_mat) == nrow(dat), nrow(v_mat) == ncol(dat))

  n <- nrow(dat); p <- ncol(dat)
  nat_mat <- u_mat %*% t(v_mat)
  stopifnot(all(nat_mat > 0))

  idx <- which(!is.na(dat))

  1/(n*p) * sum(-log(nat_mat[idx]) -
        nat_mat[idx]*dat[idx]*scalar^2 +
        nat_mat[idx]^2*dat[idx]^2*scalar^2/2)
}

.evaluate_objective_single.curved_gaussian <- function(dat_vec, current_vec, other_mat, n, p,
                                                scalar = 2, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  stopifnot(all(pred_vec > 0))

  idx <- which(!is.na(dat_vec))

  1/(n*p) * sum(-log(pred_vec[idx]) -
        pred_vec[idx]*dat_vec[idx]*scalar^2 +
        pred_vec[idx]^2*dat_vec[idx]^2*scalar^2/2)
}

.gradient_vec.curved_gaussian <- function(dat_vec, current_vec, other_mat, n, p, scalar = 2, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  stopifnot(all(pred_vec > 0))

  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(- 1/pred_vec[j] - dat_vec[j]*scalar^2 +
                          pred_vec[j]*scalar^2*dat_vec[j]^2)
  })

  if(is.matrix(tmp)) 1/(n*p) * rowSums(tmp) else 1/(n*p) * sum(tmp)
}

.evaluate_objective_mat.curved_gaussian <- function(dat, nat_mat, scalar = 2, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat > 0))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum(-log(nat_mat[idx]) -
        nat_mat[idx]*dat[idx]*scalar^2 +
        nat_mat[idx]^2*dat[idx]^2*scalar^2/2)
}

.gradient_mat.curved_gaussian <- function(dat, nat_mat, scalar = 2, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat > 0))

  n <- nrow(dat); p <- ncol(dat)

  (-1/(nat_mat) - scalar^2*dat + scalar^2*dat^2*nat_mat)/(n*p)
}


