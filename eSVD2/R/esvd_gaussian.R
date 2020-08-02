# distribution: Gaussian
# natural parameter: m_{ij} = u_i^Tv_j
# relation to canonical parameters: m_{ij} = mu_{ij}
# optimization problem: (m_{ij} - a_{ij})^2

.evaluate_objective.gaussian <- function(dat, u_mat, v_mat, ...){
  stopifnot(ncol(u_mat) == ncol(v_mat), nrow(u_mat) == nrow(dat), nrow(v_mat) == ncol(dat))

  n <- nrow(dat); p <- ncol(dat)
  nat_mat <- u_mat %*% t(v_mat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum((nat_mat[idx] - dat[idx])^2)
}

.evaluate_objective_single.gaussian <- function(dat_vec, current_vec, other_mat, n, p, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  1/(n*p) * sum((pred_vec[idx] - dat_vec[idx])^2)
}

.gradient_vec.gaussian <- function(dat_vec, current_vec, other_mat, n, p, ...){
  stopifnot(length(current_vec) == ncol(other_mat), nrow(other_mat) == length(dat_vec))

  pred_vec <- other_mat %*% current_vec
  idx <- which(!is.na(dat_vec))

  tmp <- sapply(idx, function(j){
    2*other_mat[j,,drop=F]*(pred_vec[j] - dat_vec[j])
  })

  if(is.matrix(tmp)) 1/(n*p) * rowSums(tmp) else 1/(n*p) * sum(tmp)
}

.evaluate_objective_mat.gaussian <- function(dat, nat_mat, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)))

  n <- nrow(dat); p <- ncol(dat)
  idx <- which(!is.na(dat))

  1/(n*p) * sum((nat_mat[idx] - dat[idx])^2)
}

.gradient_mat.gaussian <- function(dat, nat_mat, scalar = 2, ...){
  stopifnot(all(dim(dat) == dim(nat_mat)))

  n <- nrow(dat); p <- ncol(dat)

  1/(n*p) * 2 * (nat_mat - dat)
}


