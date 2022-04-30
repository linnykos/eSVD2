initialize_esvd <- function(dat,
                            covariates = NULL,
                            case_control_variable = NULL,
                            k = 30,
                            lambda = 0.1,
                            mixed_effect_variables = NULL,
                            offset_variables = NULL,
                            p_val_thres = 0.05,
                            tol = 1e-3,
                            verbose = 0,
                            tmp_path = NULL){
  stopifnot(is.matrix(covariates),
            all(offset_variables %in% colnames(covariates)),
            length(case_control_variable) == 1,
            case_control_variable %in% colnames(covariates),
            k <= ncol(dat), k > 0, k %% 1 == 0,
            lambda <= 1e4)
  if(!all(is.null(offset_variables))){
    stopifnot(is.character(offset_variables),
              all(offset_variables %in% colnames(covariates)))
  }
  if(!all(is.null(case_control_variable))){
    stopifnot(is.character(case_control_variable),
              all(case_control_variable %in% colnames(covariates)))
  }
  if(!all(is.null(mixed_effect_variables))){
    stopifnot(is.character(mixed_effect_variables),
              all(mixed_effect_variables %in% colnames(covariates)))
  }

  n <- nrow(dat); p <- ncol(dat)
  dat[is.na(dat)] <- 0
  param <- .initalize_format_param(case_control_variable = case_control_variable,
                                   offset_variables = offset_variables,
                                   k = k,
                                   lambda = lambda,
                                   mixed_effect_variables = mixed_effect_variables,
                                   p_val_thres = p_val_thres)

  # [[note to self: make mixed_effect_variables easier to use]]

  if(verbose >= 1) print("Step 1: Performing GLMs")
  tmp <- .initialize_coefficient(case_control_variable = case_control_variable,
                                 covariates = covariates,
                                 dat = dat,
                                 lambda = lambda,
                                 mixed_effect_variables = mixed_effect_variables,
                                 offset_variables = offset_variables,
                                 p_val_thres = p_val_thres,
                                 verbose = verbose,
                                 tmp_path = tmp_path)
  z_mat <- tmp$z_mat
  log_pval_vec <- tmp$log_pval_vec
  if(!is.null(tmp_path)) save(z_mat, file = tmp_path)

  # do some cleanup
  if(verbose >= 1) print("Step 1b: Cleaning up coefficients")
  covariates <- cbind(rep(1, n), covariates)
  colnames(covariates)[1] <- "Intercept"
  col_idx <- sapply(colnames(covariates), function(i){which(colnames(z_mat) == i)})
  z_mat <- z_mat[,as.numeric(col_idx)]

  if(verbose >= 1) print("Step 2: Computing residuals")
  tmp <- .initialize_residuals(z_mat = z_mat,
                               covariates = covariates,
                               dat = dat,
                               k = k)

  structure(list(x_mat = tmp$x_mat, y_mat = tmp$y_mat,
                 z_mat = z_mat,
                 covariates = covariates,
                 log_pval_vec = log_pval_vec,
                 param = param),
            class = "eSVD")
}

#####################

.initialize_coefficient <- function(case_control_variable,
                                    covariates,
                                    dat,
                                    lambda,
                                    mixed_effect_variables,
                                    offset_variables,
                                    p_val_thres,
                                    verbose = 0,
                                    tmp_path = NULL){
  n <- nrow(dat); p <- ncol(dat)
  covariates_nooffset <- covariates[,which(!colnames(covariates) %in% offset_variables)]
  offset_vec <- Matrix::rowSums(covariates[,offset_variables,drop = F])

  log_pval_vec <- rep(NA, p)
  z_mat <- matrix(NA, nrow = p, ncol = ncol(covariates_nooffset)+1)
  rownames(z_mat) <- colnames(dat)
  colnames(z_mat) <- c("Intercept", colnames(covariates_nooffset))
  for(j in 1:p){
    if(verbose == 1 && p >= 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2){
      verbose_additional_msg <- paste0(" (", j , " of ", p, ")")
    } else {
      verbose_additional_msg <- ""
    }
    tmp <- .lrt_coefficient(case_control_variable = case_control_variable,
                            covariates = covariates_nooffset,
                            lambda = lambda,
                            mixed_effect_variables = mixed_effect_variables,
                            offset_vec = offset_vec,
                            p_val_thres = p_val_thres,
                            vec = as.numeric(dat[,j]),
                            verbose = verbose,
                            verbose_additional_msg = verbose_additional_msg,
                            verbose_gene_name = colnames(dat)[j])
    z_mat[j,] <- tmp$vec
    log_pval_vec[j] <- tmp$log_pval

    if(!is.null(tmp_path) && p >= 10 && floor(p/10) == 0){
      save(z_mat, file = tmp_path)
    }
  }

  z_mat <- .append_offset(z_mat = z_mat,
                          covariates = covariates)

  list(z_mat = z_mat,
       log_pval_vec = log_pval_vec)
}

# [[note to self: hard coded for Poisson at the moment]]
.lrt_coefficient <- function(case_control_variable,
                             covariates,
                             lambda,
                             mixed_effect_variables,
                             offset_vec,
                             p_val_thres,
                             vec,
                             verbose = 0,
                             verbose_additional_msg = "",
                             verbose_gene_name = ""){
  # see https://statisticaloddsandends.wordpress.com/2018/11/13/a-deep-dive-into-glmnet-penalty-factor/
  penalty_factor1 <- rep(0, ncol(covariates))
  penalty_factor1[colnames(covariates) %in% mixed_effect_variables] <- 1
  glm_fit1 <- glmnet::glmnet(x = covariates,
                             y = vec,
                             alpha = 0,
                             family = "poisson",
                             intercept = T,
                             lambda = exp(seq(log(1e4), log(lambda), length.out = 100)),
                             offset = offset_vec,
                             penalty.factor = penalty_factor1,
                             standardize = F)
  coef_vec1 <- c(glm_fit1$a0[length(glm_fit1$a0)], glm_fit1$beta[,ncol(glm_fit1$beta)])
  mean_vec1 <- exp(covariates %*% coef_vec1[-1] + offset_vec + coef_vec1[1])
  log_vec <- vec/mean_vec1; log_vec[vec != 0] <- log(log_vec[vec != 0])
  deviance1 <- 2*sum(vec*log_vec - (vec - mean_vec1))

  covariates2 <- covariates[,which(colnames(covariates) != case_control_variable), drop = F]
  penalty_factor2 <- rep(0, ncol(covariates2))
  penalty_factor2[colnames(covariates2) %in% mixed_effect_variables] <- 1
  glm_fit2 <- glmnet::glmnet(x = covariates2,
                             y = vec,
                             alpha = 0,
                             family = "poisson",
                             intercept = T,
                             lambda = exp(seq(log(1e4), log(lambda), length.out = 100)),
                             offset = offset_vec,
                             penalty.factor = penalty_factor2,
                             standardize = F)
  coef_vec2 <- c(glm_fit2$a0[length(glm_fit2$a0)], glm_fit2$beta[,ncol(glm_fit2$beta)])
  mean_vec2 <- exp(covariates2 %*% coef_vec2[-1] + offset_vec + coef_vec2[1])
  log_vec <- vec/mean_vec2; log_vec[vec != 0] <- log(log_vec[vec != 0])
  deviance2 <- 2*sum(vec*log_vec - (vec - mean_vec2))

  residual_deviance <- max(deviance2 - deviance1, 0)
  log_p_val <- stats::pchisq(residual_deviance, df = 1,
                             lower.tail = FALSE,
                             log.p = T)

  if(log_p_val <= log(p_val_thres)){
    names(coef_vec1) <- c("Intercept", colnames(covariates))
    if(verbose >= 2) print(paste0(verbose_gene_name, verbose_additional_msg,
                                  ": Significant (deviance=", round(residual_deviance,1), "), coefficent: ",
                                  round(coef_vec1[case_control_variable], 2)))
    return(list(vec = coef_vec1, log_pval = log_p_val))
  } else {
    names(coef_vec2) <- c("Intercept", colnames(covariates2))
    vec2 <- sapply(c("Intercept", colnames(covariates)), function(var){
      if(var %in% names(coef_vec2)) coef_vec2[var] else 0
    })
    names(vec2) <- c("Intercept", colnames(covariates))
    if(verbose >= 2) print(paste0(verbose_gene_name, verbose_additional_msg,
                                  ": Insignificant (deviance=", round(residual_deviance,1), ")"))
    return(list(vec = vec2, log_pval = log_p_val))
  }
}

.append_offset <- function(z_mat,
                           covariates){
  z_mat2 <- matrix(1, nrow = nrow(z_mat), ncol = ncol(covariates)+1)
  rownames(z_mat2) <- rownames(z_mat)
  colnames(z_mat2) <- c("Intercept", colnames(covariates))

  for(j in 1:ncol(z_mat2)){
    if(colnames(z_mat2)[j] %in% colnames(z_mat)){
      idx <- which(colnames(z_mat) == colnames(z_mat2)[j])
      z_mat2[,j] <- z_mat[,idx]
    } else if(colnames(z_mat2)[j] == "Intercept"){
      z_mat2[,j] <- z_mat[,"Intercept"]
    }
  }

  z_mat2
}

.initialize_residuals <- function(z_mat,
                                  covariates,
                                  dat,
                                  k){
  dat_transform <- log1p(dat)
  nat_mat <- tcrossprod(covariates, z_mat)
  residual_mat <- dat_transform - nat_mat

  svd_res <- .svd_safe(mat = residual_mat,
                       check_stability = T,
                       K = k,
                       mean_vec = NULL,
                       rescale = F,
                       scale_max = NULL,
                       sd_vec = NULL)
  x_mat <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_mat <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  rownames(x_mat) <- rownames(dat)
  rownames(y_mat) <- colnames(dat)

  list(x_mat = x_mat, y_mat = y_mat)
}

.initalize_format_param <- function(case_control_variable,
                                    offset_variables,
                                    k,
                                    lambda,
                                    mixed_effect_variables,
                                    p_val_thres) {
  list(case_control_variable = case_control_variable,
       offset_variables = offset_variables,
       k = k,
       lambda = lambda,
       mixed_effect_variables = mixed_effect_variables,
       p_val_thres = p_val_thres)
}
