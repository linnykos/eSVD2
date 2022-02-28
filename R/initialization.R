#' Initialize eSVD
#'
#' @param dat                     dataset where the \eqn{n} rows represent cells and \eqn{p} columns represent genes
#' @param k                       positive integer less than \code{min(c(nrow(dat), ncol(dat)))}
#' @param family                  character (\code{"gaussian"}, \code{"exponential"}, \code{"poisson"}, \code{"neg_binom"},
#'                                or \code{"curved_gaussian"})
#' @param covariates              an \eqn{n \times d}{n × d} matrix representing the additional \eqn{d} covariates,
#'                                or \code{NULL} if no covariate is given
#' @param offset_vec              a vector of length-\eqn{n} that represents a constant amount added to each row of the
#'                                natural parameter matrix
#' @param column_set_to_one       vector of characters, all contained in \code{colnames(covariates)}
#' @param tol                     small positive number
#' @param verbose                 non-negative integer specifying level of printouts
#'
#' @return a list with elements \code{x_mat} and \code{y_mat} (and others), representing the two
#' latent matrices
#' @export
initialize_esvd <- function(dat,
                            covariates,
                            case_control_variable,
                            offset_variables,
                            family = "poisson",
                            k = 10,
                            lambda = 0.01,
                            p_val_thres = 0.05,
                            tol = 1e-3,
                            verbose = 0){
  stopifnot(is.character(family),
            family %in% c("gaussian", "poisson", "neg_binom2"),
            is.matrix(covariates),
            all(offset_variables %in% colnames(covariates)),
            length(case_control_variable) == 1,
            case_control_variable %in% colnames(covariates),
            k <= ncol(dat), k > 0, k %% 1 == 0,
            lambda <= 1e4)
  n <- nrow(dat); p <- ncol(dat)
  dat[is.na(dat)] <- 0

  offset_vec <- Matrix::rowSums(covariates[,offset_variables])

  b_mat <- .initialize_coefficient(case_control_variable = case_control_variable,
                                   covariates = covariates,
                                   dat = dat,
                                   lambda = lambda,
                                   offset_variables = offset_variables,
                                   offset_vec = offset_vec,
                                   p_val_thres = p_val_thres,
                                   verbose = verbose)
  if(include_intercept) {
    covariates <- cbind(rep(1, n), covariates)
    colnames(covariates) <- "Intercept"
  }

  tmp <- .initialize_residuals(b_mat = b_mat,
                               covariates = covariates,
                               dat = dat,
                               k = k)

  structure(list(x_mat = tmp$x_mat, y_mat = tmp$y_mat,
                 b_mat = b_mat,
                 covariates = covariates,
                 nuisance_param_vec = rep(0, ncol(dat))),
            class = "eSVD")
}

#####################

.initialize_coefficient <- function(case_control_variable,
                                    covariates,
                                    dat,
                                    lambda,
                                    offset_variables,
                                    offset_vec,
                                    p_val_thres,
                                    verbose = 0){
  n <- nrow(dat); p <- ncol(dat)
  covariates_nooffset <- covariates[,which(!colnames(covariates) %in% offset_variables)]

  b_mat <- sapply(1:p, function(j){
    .lrt_coefficient(case_control_variable = case_control_variable,
                     covariates = covariates_nooffset,
                     lambda = lambda,
                     offset_vec = offset_vec,
                     p_val_thres = p_val_thres,
                     vec = dat[,j],
                     verbose = verbose,
                     verbose_gene_name = colnames(dat)[j])
  })
  b_mat <- t(b_mat)
  rownames(b_mat) <- colnames(dat)
  colnames(b_mat) <- colnames(covariates_nooffset)

  b_mat <- .append_offset(b_mat = b_mat,
                          covariates = covariates,
                          offset_variables = offset_variables)

  b_mat
}

# [[note to self: hard coded for Poisson at the moment]]
.lrt_coefficient <- function(case_control_variable,
                             covariates,
                             lambda,
                             offset_vec,
                             p_val_thres,
                             vec,
                             verbose = 0,
                             verbose_gene_name = ""){
  glm_fit1 <- glmnet::glmnet(x = covariates,
                             y = vec,
                             family = "poisson",
                             offset = offset_vec,
                             alpha = 0,
                             standardize = F,
                             intercept = T,
                             lambda = exp(seq(log(1e4), log(lambda), length.out = 100)))
  coef_vec1 <- c(glm_fit1$a0[length(glm_fit1$a0)], glm_fit1$beta[,ncol(glm_fit1$beta)])
  mean_vec1 <- exp(covariates %*% coef_vec1[-1] + offset_vec + coef_vec1[1])
  log_vec <- vec/mean_vec1; log_vec[vec != 0] <- log(log_vec[vec != 0])
  deviance1 <- 2*sum(vec*log_vec - (vec - mean_vec1))

  covariates2 <- covariates[,which(colnames(covariates) != case_control_variable), drop = F]
  glm_fit2 <- glmnet::glmnet(x = covariates2,
                             y = vec,
                             family = "poisson",
                             offset = offset_vec,
                             alpha = 0,
                             standardize = F,
                             intercept = T,
                             lambda = exp(seq(log(1e4), log(lambda), length.out = 100)))
  coef_vec2 <- c(glm_fit2$a0[length(glm_fit2$a0)], glm_fit2$beta[,ncol(glm_fit2$beta)])
  mean_vec2 <- exp(covariates2 %*% coef_vec2[-1] + offset_vec + coef_vec2[1])
  log_vec <- vec/mean_vec2; log_vec[vec != 0] <- log(log_vec[vec != 0])
  deviance2 <- 2*sum(vec*log_vec - (vec - mean_vec2))

  residual_deviance <- max(deviance2 - deviance1, 0)
  p_val <- 1-stats::pchisq(residual_deviance, df = 1)

  if(p_val <= p_val_thres){
    names(coef_vec1) <- c("Intercept", colnames(covariates))
    if(verbose >= 2) print(paste0(verbose_gene_name, ": Significant (deviance=", round(residual_deviance,1), "), coefficent: ",
                             round(coef_vec1[case_control_variable], 2)))
    return(coef_vec1)
  } else {
    names(coef_vec2) <- c("Intercept", colnames(covariates2))
    vec2 <- sapply(c("Intercept", colnames(covariates)), function(var){
      if(var %in% names(coef_vec2)) coef_vec2[var] else 0
    })
    names(vec2) <- c("Intercept", colnames(covariates))
    if(verbose <= 2) print(paste0(verbose_gene_name, ": Insignificant (deviance=", round(residual_deviance,1), ")"))
    return(vec2)
  }
}

.append_offset <- function(b_mat,
                           covariates){
  b_mat2 <- matrix(1, nrow = nrow(b_mat), ncol = ncol(covariates))
  for(j in 1:ncol(covariates)){
    if(colnames(covariates)[j] %in% colnames(b_mat)){
      idx <- which(colnames(b_mat) == colnames(covariates)[j])
      b_mat2[,j] <- b_mat[,j]
    }
  }

  rownames(b_mat2) <- rownames(b_mat)
  colnames(b_mat2) <- colnames(covariates)

  b_mat2
}

.initialize_residuals <- function(b_mat,
                                  covariates,
                                  dat,
                                  k){
  dat_transform <- log1p(dat)
  nat_mat <- tcrossprod(covariates, b_mat)
  residual_mat <- dat_transform - nat_mat

  svd_res <- .svd_truncated(residual_mat,
                            K = k,
                            symmetric = F,
                            rescale = F,
                            mean_vec = NULL,
                            sd_vec = NULL,
                            K_full_rank = F)
  x_mat <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_mat <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  rownames(x_mat) <- rownames(dat)
  rownames(y_mat) <- colnames(dat)

  list(x_mat = x_mat, y_mat = y_mat)
}

# .initialize_nuisance <- function(){
#   init_theta <- tcrossprod(x_init, y_init) + nat_offset_mat
#   init_theta <- sweep(init_theta, 1, offset_vec, "+")
#
#   if(family$name %in% c("neg_binom", "neg_binom2")){
#     glmgampoi_init <- glmGamPoi::overdispersion_mle(y = Matrix::t(dat),
#                                                     mean = family$nat_to_canon(t(init_theta)),
#                                                     global_estimate = T)
#     nuisance_init <- rep(glmgampoi_init$estimate, p)
#   } else {
#     nuisance_init <- NULL
#   }
# }
