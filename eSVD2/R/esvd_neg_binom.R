# distribution: negative binomial
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: m_{ij} = log(p_{ij})
# optimization problem: (-r*log(1-exp(m_{ij}))) - a_{ij}*mu_{ij}

.evaluate_objective.neg_binom <- function(dat, u_mat, v_mat, scalar, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(u_mat) == nrow(dat), nrow(v_mat) == ncol(dat))

  n <- nrow(dat); p <- ncol(dat)
  nat_mat <- u_mat %*% t(v_mat)
  stopifnot(all(nat_mat < 0))

  idx <- which(!is.na(dat))

  1/(n*p) * sum(-scalar * log(1-exp(nat_mat[idx])) - nat_mat[idx]*dat[idx])
}

.evaluate_objective_single.neg_binom <- function(dat_vec, current_vec, other_mat, n, p, scalar, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  stopifnot(all(pred_vec < 0))

  idx <- which(!is.na(dat_vec))

  1/(n*p) * sum(-scalar * log(1-exp(pred_vec[idx])) - pred_vec[idx]*dat_vec[idx])
}

.gradient_vec.neg_binom <- function(dat_vec, current_vec, other_mat, n, p, scalar, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  stopifnot(all(pred_vec < 0))

  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    other_mat[j,,drop=F]*(scalar * exp(pred_vec[j])/(1-exp(pred_vec[j])) - dat_vec[j])
  })

  if(is.matrix(tmp)) 1/(n*p) * rowSums(tmp) else 1/(n*p) * sum(tmp)
}

.evaluate_objective_mat.neg_binom <- function(dat, nat_mat, scalar, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat < 0))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum(-scalar * log(1-exp(nat_mat[idx])) - nat_mat[idx]*dat[idx])
}

.gradient_mat.neg_binom <- function(dat, nat_mat, scalar, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)), all(nat_mat < 0))

  n <- nrow(dat); p <- ncol(dat)
  (scalar * exp(nat_mat)/(1-exp(nat_mat)) - dat)/(n*p)
}



