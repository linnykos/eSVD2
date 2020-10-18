initialize_nuisance_param <- function(dat, init_nat_mat, family,
                                      library_size_vec){
  stopifnot(all(dim(dat) == dim(init_nat_mat)), all(dat >= 0))
  stopifnot((family == "curved_gaussian" && all(init_nat_mat < 0)) ||
              (family == "neg_binom"))

  previous_family <- ifelse(family == "curved_gaussian", "exponential", "poisson")

  mean_mat <- compute_mean(init_nat_mat, family = previous_family,
                           library_size_vec = library_size_vec)

  param_vec <- sapply(1:ncol(dat), function(j){
    empirical_var_vec <- (dat[,j] - mean_mat[,j])^2

    if(family == "curved_gaussian"){
      tmp <- mean_mat[,j]/sqrt(empirical_var_vec)
      min_val <- max(min(tmp), 0); max_val <- max(max(tmp), 10)

      tryCatch(stats::uniroot(.root_curved_gaussian_closure(empirical_var_vec, mean_mat[,j]),
                              interval = c(min_val, max_val))$root,
               error = function(e){
                 1
               })

    } else {
      tmp <- (empirical_var_vec- mean_mat[,j])/(mean_mat[,j]^2)
      min_val <- max(min(tmp), 0); max_val <- max(max(tmp), 10)

      tryCatch(stats::uniroot(.root_neg_binom_closure(empirical_var_vec, mean_mat[,j]),
                              interval = c(min_val, max_val))$root,
               error = function(e){
                 1/(10*max(dat[,j]))
               })
    }
  })

  if(family == "neg_binom") param_vec <- 1/param_vec
  param_vec
}

.root_curved_gaussian_closure <- function(empirical_var_vec, theoretical_mean_vec){
  function(x){
    sum((empirical_var_vec - (theoretical_mean_vec/x)^2) / (theoretical_mean_vec/x^2))
  }
}

.root_neg_binom_closure <- function(empirical_var_vec, theoretical_mean_vec){
  function(x){
    sum((empirical_var_vec - theoretical_mean_vec * (1+ x * theoretical_mean_vec)) /
          (1+ x * theoretical_mean_vec))
  }
}
