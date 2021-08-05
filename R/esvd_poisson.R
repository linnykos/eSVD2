# Distribution: Poisson

# .evaluate_objective.poisson <- function(
#   dat, x_mat, y_mat, nuisance_param_vec, library_size_vec, ...
# ) {
#   # Check dimensions
#   n <- nrow(dat)
#   p <- ncol(dat)
#   stopifnot(
#     length(library_size_vec) == n,
#     ncol(x_mat) == ncol(y_mat),
#     nrow(x_mat) == n,
#     nrow(y_mat) == p
#   )
#
#   # Compute natural parameters
#   nat_mat <- tcrossprod(x_mat, y_mat)
#
#   # Only compute likelihood on non-missing data
#   idx <- which(!is.na(dat))
#   stopifnot(length(idx) > 0)
#   # `negloglik` contains NAs with the same locations as those in `dat`
#   negloglik <- exp(nat_mat + log(library_size_vec)) - nat_mat * dat
#   mean(negloglik[idx])
# }
#
# # If current vector is x: length(library_size_vec) == 1
# # If current vector is y: length(library_size_vec) == n
# .evaluate_objective_single.poisson <- function(
#   current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
# ) {
#   stopifnot(
#     length(current_vec) == ncol(other_mat),
#     length(dat_vec) == nrow(other_mat)
#   )
#
#   nat_vec <- c(other_mat %*% current_vec)
#   idx <- which(!is.na(dat_vec))
#   stopifnot(length(idx) > 0)
#
#   negloglik <- exp(nat_vec + log(library_size_vec)) - nat_vec * dat_vec
#   mean(negloglik[idx])
# }
#
# .gradient_vec.poisson <- function(
#   current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
# ) {
#   stopifnot(
#     length(current_vec) == ncol(other_mat),
#     length(dat_vec) == nrow(other_mat)
#   )
#
#   nat_vec <- c(other_mat %*% current_vec)
#   idx <- which(!is.na(dat_vec))
#   stopifnot(length(idx) > 0)
#
#   grad <- other_mat * (exp(nat_vec + log(library_size_vec)) - dat_vec)
#   # `grad` contains NA rows corresponding to the NAs in `dat_vec`
#   grad <- grad[idx, , drop = FALSE]
#   colMeans(grad)
# }
#
# .hessian_vec.poisson <- function(
#   current_vec, other_mat, dat_vec, nuisance_param_vec, library_size_vec, ...
# ) {
#   stopifnot(
#     length(current_vec) == ncol(other_mat),
#     length(dat_vec) == nrow(other_mat)
#   )
#
#   nat_vec <- c(other_mat %*% current_vec)
#   idx <- which(!is.na(dat_vec))
#   stopifnot(length(idx) > 0)
#
#   exp_term <- exp(nat_vec + log(library_size_vec))
#   exp_other_term <- .mult_vec_mat(exp_term, other_mat)
#   exp_other_vals <- exp_other_term[idx, , drop = FALSE]
#   other_vals <- other_mat[idx, , drop = FALSE]
#   hess <- crossprod(exp_other_vals, other_vals) / length(idx)
#   hess
# }
#
# .feasibility.poisson <- function(current_vec, other_mat, ...) {
#   TRUE
# }
#
# .poisson <- structure(
#   list(
#     objfn_all = .evaluate_objective.poisson,
#     objfn     = .evaluate_objective_single.poisson,
#     grad      = .gradient_vec.poisson,
#     hessian   = .hessian_vec.poisson,
#     feas      = .feasibility.poisson
#   ),
#   class = "esvd_family"
# )



# See eSVD2_writing/writeup/Writeup4
#
# Log-density for the whole data matrix [n x p]
.log_prob.poisson <- function(A, theta, s, gamma)
{
  A * theta - exp(log(s) + theta)
}

# Log-density for the i-th row of the data matrix [p x 1]
.log_prob_row.poisson <- function(Ai, thetai, si, gamma)
{
  Ai * thetai - exp(log(si) + thetai)
}

# Log-density for the j-th column of the data matrix [n x 1]
.log_prob_col.poisson <- function(Aj, thetaj, s, gammaj)
{
  Aj * thetaj - exp(log(s) + thetaj)
}

# 1st derivative of log-density w.r.t. the i-th row of theta [p x 1]
.dlog_prob_row.poisson <- function(Ai, thetai, si, gamma)
{
  Ai - exp(log(si) + thetai)
}

# 1st derivative of log-density w.r.t. the j-th column of theta [n x 1]
.dlog_prob_col.poisson <- function(Aj, thetaj, s, gammaj)
{
  Aj - exp(log(s) + thetaj)
}

# 2nd derivative of log-density w.r.t. the i-th row of theta [p x 1]
.d2log_prob_row.poisson <- function(Ai, thetai, si, gamma)
{
  -exp(log(si) + thetai)
}

# 2nd derivative of log-density w.r.t. the j-th column of theta [n x 1]
.d2log_prob_col.poisson <- function(Aj, thetaj, s, gammaj)
{
  -exp(log(s) + thetaj)
}

# Feasibility of the natural parameter
.feasibility.poisson <- function(theta)
{
  TRUE
}

.dat_to_nat.poisson <- function(A, gamma, tol = 1e-3){
  res <- log(A + tol)

  if(length(rownames(A)) != 0) rownames(res) <- rownames(A)
  if(length(colnames(A)) != 0) colnames(res) <- colnames(A)

  res
}
.nat_to_canon.poisson <- function(theta){
  exp(theta)
}

.esvd.poisson <- structure(
  list(
    name           = "poisson",
    log_prob       = .log_prob.poisson,
    log_prob_row   = .log_prob_row.poisson,
    log_prob_col   = .log_prob_col.poisson,
    dlog_prob_row  = .dlog_prob_row.poisson,
    dlog_prob_col  = .dlog_prob_col.poisson,
    d2log_prob_row = .d2log_prob_row.poisson,
    d2log_prob_col = .d2log_prob_col.poisson,
    feasibility    = .feasibility.poisson,
    feas_always    = TRUE,
    domain         = c(-Inf, Inf),
    dat_to_nat     = .dat_to_nat.poisson,
    nat_to_canon   = .nat_to_canon.poisson
  ),
  class = "esvd_family"
)
.esvd.poisson <- list2env(.esvd.poisson)
