# Distribution: negative binomial

# .evaluate_objective.neg_binom <- function(
#   dat, x_mat, y_mat, nuisance_param_vec, library_size_vec, ...
# ) {
#   # Check dimensions
#   n <- nrow(dat)
#   p <- ncol(dat)
#   stopifnot(
#     length(nuisance_param_vec) == p,
#     length(library_size_vec) == n,
#     ncol(x_mat) == ncol(y_mat),
#     nrow(x_mat) == n,
#     nrow(y_mat) == p
#   )
#
#   # Compute natural parameters
#   nat_mat <- tcrossprod(x_mat, y_mat)
#   stopifnot(all(nat_mat < 0))
#
#   # Only compute likelihood on non-missing data
#   idx <- which(!is.na(dat))
#   stopifnot(length(idx) > 0)
#   # `negloglik` contains NAs with the same locations as those in `dat`
#   negloglik <- .mult_mat_vec(.mult_vec_mat(-library_size_vec, log(1 - exp(nat_mat))), nuisance_param_vec) - nat_mat * dat
#   mean(negloglik[idx])
# }
#
# # If current vector is x:
# #   length(library_size_vec) == 1, length(nuisance_param_vec) == p
# # If current vector is y:
# #   length(library_size_vec) == n, length(nuisance_param_vec) == 1
# .evaluate_objective_single.neg_binom <- function(
#   current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
# ) {
#   check_dim_x <- length(library_size_vec) == 1 &&
#                  length(nuisance_param_vec) == nrow(other_mat)
#   check_dim_y <- length(library_size_vec) == nrow(other_mat) &&
#                  length(nuisance_param_vec) == 1
#   stopifnot(
#     length(current_vec) == ncol(other_mat),
#     length(dat_vec) == nrow(other_mat),
#     check_dim_x || check_dim_y
#   )
#
#   nat_vec <- c(other_mat %*% current_vec)
#   stopifnot(all(nat_vec < 0))
#   idx <- which(!is.na(dat_vec))
#   stopifnot(length(idx) > 0)
#
#   negloglik <- -library_size_vec * nuisance_param_vec * log(1 - exp(nat_vec)) - nat_vec * dat_vec
#   mean(negloglik[idx])
# }
#
#
# # f(x) = log(1 - exp(x))
# # f'(x) = -exp(x) / (1 - exp(x)) = -1 / (exp(-x) - 1) = 1 / (1 - exp(-x))
# # f''(x) = -exp(-x) / (1 - exp(-x))^2 = [1 / (1 - exp(-x))] * [-exp(-x) / (1 - exp(-x))]
# #        = f'(x) * [1 - f'(x)]
#
# .gradient_vec.neg_binom <- function(
#   current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
# ) {
#   check_dim_x <- length(library_size_vec) == 1 &&
#                  length(nuisance_param_vec) == nrow(other_mat)
#   check_dim_y <- length(library_size_vec) == nrow(other_mat) &&
#                  length(nuisance_param_vec) == 1
#   stopifnot(
#     length(current_vec) == ncol(other_mat),
#     length(dat_vec) == nrow(other_mat),
#     check_dim_x || check_dim_y
#   )
#
#   nat_vec <- c(other_mat %*% current_vec)
#   stopifnot(all(nat_vec < 0))
#   idx <- which(!is.na(dat_vec))
#   stopifnot(length(idx) > 0)
#
#   grad <- other_mat * (library_size_vec * nuisance_param_vec / (exp(-nat_vec) - 1) - dat_vec)
#   # `grad` contains NA rows corresponding to the NAs in `dat_vec`
#   grad <- grad[idx, , drop = FALSE]
#   colMeans(grad)
# }
#
# .hessian_vec.neg_binom <- function(
#   current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
# ) {
#   check_dim_x <- length(library_size_vec) == 1 &&
#                  length(nuisance_param_vec) == nrow(other_mat)
#   check_dim_y <- length(library_size_vec) == nrow(other_mat) &&
#                  length(nuisance_param_vec) == 1
#   stopifnot(
#     length(current_vec) == ncol(other_mat),
#     length(dat_vec) == nrow(other_mat),
#     check_dim_x || check_dim_y
#   )
#
#   nat_vec <- c(other_mat %*% current_vec)
#   stopifnot(all(nat_vec < 0))
#   idx <- which(!is.na(dat_vec))
#   stopifnot(length(idx) > 0)
#
#   other_vals <- other_mat[idx, , drop = FALSE]
#   tmp <- exp(nat_vec)
#   diag_vec <- library_size_vec * nuisance_param_vec * tmp * (1 / (1 - tmp) + tmp / (1 - tmp)^2)
#   term1_vals <- .mult_vec_mat(diag_vec, other_mat)[idx, , drop = FALSE]
#   hess <- crossprod(term1_vals, other_vals) / length(idx)
#   hess
# }
#
# .feasibility.neg_binom <- function(current_vec, other_mat, ...) {
#   nat_vec <- c(other_mat %*% current_vec)
#   all(nat_vec < 0)
# }
#
# .neg_binom <- structure(
#   list(
#     objfn_all = .evaluate_objective.neg_binom,
#     objfn     = .evaluate_objective_single.neg_binom,
#     grad      = .gradient_vec.neg_binom,
#     hessian   = .hessian_vec.neg_binom,
#     feas      = .feasibility.neg_binom
#   ),
#   class = "esvd_family"
# )



# See eSVD2_writing/writeup/Writeup4
#
# Log-density for the whole data matrix [n x p]
.log_prob.neg_binom <- function(A, theta, s, gamma)
{
  slog1theta <- log(1 - exp(theta)) * s
  sweep(slog1theta, 2, gamma, "*") + A * theta
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.neg_binom <- function(Ai, thetai, si, gamma)
{
  Ai * thetai + si * gamma * log(1 - exp(thetai))
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.neg_binom <- function(Aj, thetaj, s, gammaj)
{
  Aj * thetaj + s * gammaj * log(1 - exp(thetaj))
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.neg_binom <- function(Ai, thetai, si, gamma)
{
  exptheta <- exp(thetai)
  Ai - si * gamma * exptheta / (1 - exptheta)
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.neg_binom <- function(Aj, thetaj, s, gammaj)
{
  exptheta <- exp(thetaj)
  Aj - s * gammaj * exptheta / (1 - exptheta)
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.neg_binom <- function(Ai, thetai, si, gamma)
{
  exptheta <- exp(thetai)
  -si * gamma * exptheta / (1 - exptheta)^2
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.neg_binom <- function(Aj, thetaj, s, gammaj)
{
  exptheta <- exp(thetaj)
  -s * gammaj * exptheta / (1 - exptheta)^2
}

# Feasibility of the natural parameter
.feasibility.neg_binom <- function(theta)
{
  all(theta < 0)
}

.dat_to_nat.neg_binom <- function(A, gamma, tol = 1e-3){
  res <- sapply(1:ncol(A), function(j){
    (A[,j] + tol)/gamma[j]
  })
  res <- log(res/(1+res))

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}

.nat_to_canon.neg_binom <- function(theta){
  exp(theta)
}

.esvd.neg_binom <- structure(
  list(
    name           = "neg_binom",
    log_prob       = .log_prob.neg_binom,
    log_prob_row   = .log_prob_row.neg_binom,
    log_prob_col   = .log_prob_col.neg_binom,
    dlog_prob_row  = .dlog_prob_row.neg_binom,
    dlog_prob_col  = .dlog_prob_col.neg_binom,
    d2log_prob_row = .d2log_prob_row.neg_binom,
    d2log_prob_col = .d2log_prob_col.neg_binom,
    feasibility    = .feasibility.neg_binom,
    feas_always    = FALSE,
    domain         = c(-Inf, 0),
    dat_to_nat     = .dat_to_nat.neg_binom,
    nat_to_canon   = .nat_to_canon.neg_binom
  ),
  class = "esvd_family"
)
.esvd.neg_binom <- list2env(.esvd.neg_binom)
