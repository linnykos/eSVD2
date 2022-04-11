compute_posterior <- function(mat,
                              esvd_res,
                              nuisance_vec,
                              case_control_variable,
                              alpha_max = 50,
                              nuisance_lower_quantile = 0.01){
  offset_var <- setdiff(colnames(esvd_res$covariates), case_control_variable)

  nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
  nat_mat2 <- tcrossprod(esvd_res$covariates[,case_control_variable,drop = F],
                         esvd_res$b_mat[,case_control_variable,drop = F])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)
  library_mat <- exp(tcrossprod(
    esvd_res$covariates[,offset_var],
    esvd_res$b_mat[,offset_var]
  ))

  nuisance_vec <- pmax(nuisance_vec, quantile(nuisance_vec, probs = nuisance_lower_quantile))
  Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
                 STATS = nuisance_vec, FUN = "*")
  Alpha <- pmin(Alpha, alpha_max)
  AplusAlpha <- mat + Alpha
  SplusBeta <- sweep(library_mat, MARGIN = 2,
                     STATS = nuisance_vec, FUN = "+")
  posterior_mean_mat <- AplusAlpha/SplusBeta
  posterior_var_mat <- AplusAlpha/SplusBeta^2

  rownames(posterior_mean_mat) <- rownames(mat)
  rownames(posterior_var_mat) <- rownames(mat)
  colnames(posterior_mean_mat) <- colnames(mat)
  colnames(posterior_var_mat) <- colnames(mat)

  list(posterior_mean_mat = posterior_mean_mat,
       posterior_var_mat = posterior_var_mat)
}
