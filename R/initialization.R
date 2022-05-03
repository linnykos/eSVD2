initialize_esvd <- function(dat,
                            bool_intercept = T,
                            covariates = NULL,
                            case_control_variable = NULL,
                            k = 30,
                            lambda = 0.1,
                            mixed_effect_variables = NULL,
                            offset_variables = NULL,
                            verbose = 0){
  stopifnot(inherits(dat, c("dgCMatrix", "matrix")),
            k <= ncol(dat), k > 0, k %% 1 == 0,
            lambda <= 1e4, lambda >= 1e-4)
  stopifnot(!all(is.null(covariates))) #[[note to self: remove this necessity]]

  if(!all(is.null(covariates))){
    stopifnot(is.matrix(covariates))

    if(!all(is.null(offset_variables))){
      stopifnot(is.character(offset_variables),
                all(offset_variables %in% colnames(covariates)))
    }
    if(!all(is.null(case_control_variable))){
      stopifnot(is.character(case_control_variable),
                all(case_control_variable %in% colnames(covariates)),
                length(unique(covariates[,case_control_variable])) == 2)
    }
    if(!all(is.null(mixed_effect_variables))){
      stopifnot(is.character(mixed_effect_variables),
                all(mixed_effect_variables %in% colnames(covariates)))
    }
  }

  n <- nrow(dat); p <- ncol(dat)
  dat[is.na(dat)] <- 0
  param <- .format_param_initialize(bool_intercept = bool_intercept,
                                    case_control_variable = case_control_variable,
                                    offset_variables = offset_variables,
                                    k = k,
                                    lambda = lambda,
                                    mixed_effect_variables = mixed_effect_variables)

  # [[note to self: make mixed_effect_variables easier to use]]

  if(verbose >= 1) print("Performing GLMs")
  initial_Reg <- .initialize_coefficient(bool_intercept = bool_intercept,
                                         case_control_variable = case_control_variable,
                                         covariates = covariates,
                                         dat = dat,
                                         lambda = lambda,
                                         mixed_effect_variables = mixed_effect_variables,
                                         offset_variables = offset_variables,
                                         verbose = verbose)

  if(bool_intercept){
    covariates <- cbind(rep(1, n), covariates)
    colnames(covariates)[1] <- "Intercept"
  }

  eSVD_obj <- list(dat = dat,
                   covariates = covariates,
                   initial_Reg = initial_Reg,
                   param = param)
  class(eSVD_obj) <- "eSVD"
  eSVD_obj
}

apply_initial_threshold <- function(eSVD_obj,
                                    pval_thres = 0.01,
                                    verbose = 0) {
  stopifnot(inherits(eSVD_obj, "eSVD"),
            c("initial_Reg" %in% names(eSVD_obj)),
            pval_thres >= 0, pval_thres <= 1)

  if(verbose >= 1) print("Assembling coefficents")
  z_mat1 <- .get_object(eSVD_obj = eSVD_obj,
                        what_obj = "z_mat1",
                        which_fit = "initial_Reg")
  z_mat2 <- .get_object(eSVD_obj = eSVD_obj,
                        what_obj = "z_mat2",
                        which_fit = "initial_Reg")
  log_pval <- .get_object(eSVD_obj = eSVD_obj,
                          what_obj = "log_pval",
                          which_fit = "initial_Reg")
  stopifnot(all(dim(z_mat1) == dim(z_mat2)), length(log_pval) == nrow(z_mat1))
  p <- nrow(z_mat1)
  z_mat <- matrix(NA, nrow = nrow(z_mat1), ncol = ncol(z_mat2))
  for(j in 1:p){
    if(log_pval[j] <= log(pval_thres)) z_mat[j,] <- z_mat1[j,] else z_mat[j,] <- z_mat2[j,]
  }

  dat <- .get_object(eSVD_obj = eSVD_obj, what_obj = "dat", which_fit = NULL)
  covariates <- .get_object(eSVD_obj = eSVD_obj, what_obj = "covariates", which_fit = NULL)
  k <- .get_object(eSVD_obj = eSVD_obj, what_obj = "init_k", which_fit = "param")

  if(verbose >= 1) print("Computing residuals")
  eSVD_obj[["fit_Init"]] <- .initialize_residuals(
    covariates = covariates,
    dat = dat,
    k = k,
    z_mat = z_mat
  )
  eSVD_obj[["latest_Fit"]] <- "fit_Init"
  eSVD_obj
}

#####################

.initialize_coefficient <- function(bool_intercept,
                                    case_control_variable,
                                    covariates,
                                    dat,
                                    lambda,
                                    mixed_effect_variables,
                                    offset_variables,
                                    verbose = 0){
  stopifnot(bool_intercept == T) # [[note to self: change this]]
  n <- nrow(dat); p <- ncol(dat)
  covariates_nooffset <- covariates[,which(!colnames(covariates) %in% offset_variables)]
  offset_vec <- Matrix::rowSums(covariates[,offset_variables,drop = F])

  log_pval <- rep(NA, p)
  z_mat1 <- matrix(NA, nrow = p, ncol = ncol(covariates_nooffset)+1)
  z_mat2 <- matrix(NA, nrow = p, ncol = ncol(covariates_nooffset)+1)

  rownames(z_mat1) <- colnames(dat)
  colnames(z_mat1) <- c("Intercept", colnames(covariates_nooffset))
  rownames(z_mat2) <- colnames(dat)
  colnames(z_mat2) <- c("Intercept", colnames(covariates_nooffset))

  for(j in 1:p){
    if(verbose == 1 && p >= 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(paste0("Finished variable ", j , " of ", p, ")"))
    tmp <- .lrt_coefficient(case_control_variable = case_control_variable,
                            covariates = covariates_nooffset,
                            lambda = lambda,
                            mixed_effect_variables = mixed_effect_variables,
                            offset_vec = offset_vec,
                            vec = as.numeric(dat[,j]),
                            verbose = verbose)
    z_mat1[j,] <- tmp$coef_vec1
    z_mat2[j,] <- tmp$coef_vec2
    log_pval[j] <- tmp$log_pval
  }

  z_mat1 <- .append_offset(bool_intercept = bool_intercept,
                           covariates = covariates,
                           z_mat = z_mat1)
  z_mat2 <- .append_offset(bool_intercept = bool_intercept,
                           covariates = covariates,
                           z_mat = z_mat2)

  # Do some final formatting
  col_idx <- sapply(c("Intercept", colnames(covariates)), function(i){which(colnames(z_mat1) == i)})
  z_mat1 <- z_mat1[,as.numeric(col_idx)]
  z_mat2 <- z_mat2[,as.numeric(col_idx)]

  structure(list(log_pval = log_pval,
                 z_mat1 = z_mat1,
                 z_mat2 = z_mat2),
            class = "initial_Reg")
}

# [[note to self: hard coded for Poisson at the moment]]
.lrt_coefficient <- function(case_control_variable,
                             covariates,
                             lambda,
                             mixed_effect_variables,
                             offset_vec,
                             vec,
                             verbose = 0){
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
  log_pval <- stats::pchisq(residual_deviance, df = 1,
                            lower.tail = FALSE,
                            log.p = T)

  # Tailor results and return
  names(coef_vec1) <- c("Intercept", colnames(covariates))
  names(coef_vec2) <- c("Intercept", colnames(covariates2))
  coef_vec2b <- sapply(c("Intercept", colnames(covariates)), function(var){
    if(var %in% names(coef_vec2)) coef_vec2[var] else 0
  })
  names(coef_vec2b) <- c("Intercept", colnames(covariates))

  list(coef_vec1 = coef_vec1,
       coef_vec2 = coef_vec2b,
       log_pval = log_pval)
}

.append_offset <- function(bool_intercept,
                           covariates,
                           z_mat){
  if(bool_intercept) p <- ncol(covariates)+1 else p <- ncol(covariates)
  z_mat2 <- matrix(1, nrow = nrow(z_mat), ncol = p)
  rownames(z_mat2) <- rownames(z_mat)
  if(bool_intercept){
    colnames(z_mat2) <- c("Intercept", colnames(covariates))
  } else {
    colnames(z_mat2) <- colnames(covariates)
  }

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

.initialize_residuals <- function(covariates,
                                  dat,
                                  k,
                                  z_mat){
  dat_transform <- log1p(as.matrix(dat))
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

  .form_esvd_fit(x_mat = x_mat, y_mat = y_mat, z_mat = z_mat)
}

.format_param_initialize <- function(bool_intercept,
                                     case_control_variable,
                                     offset_variables,
                                     k,
                                     lambda,
                                     mixed_effect_variables) {
  list(init_bool_intercept = bool_intercept,
       init_case_control_variable = case_control_variable,
       init_offset_variables = offset_variables,
       init_k = k,
       init_lambda = lambda,
       init_mixed_effect_variables = mixed_effect_variables)
}
