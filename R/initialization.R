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
#' @param verbose                 non-negative integer specifying level of printouts
#'
#' @return a list with elements \code{x_mat} and \code{y_mat} (and others), representing the two
#' latent matrices
#' @export
initialize_esvd <- function(dat,
                            k,
                            family,
                            covariates = NULL,
                            offset_vec = rep(0, nrow(dat)),
                            verbose = 0){
  stopifnot(is.character(family),
            family %in% c("gaussian", "poisson", "neg_binom2"),
            all(is.null(covariates)) || is.matrix(covariates))

  family <- .string_to_distr_funcs(family)
  if(family$name != "gaussian") stopifnot(all(dat[!is.na(dat)] >= 0))

  n <- nrow(dat); p <- ncol(dat)
  dat[is.na(dat)] <- 0

  if(!all(is.null(covariates))){
    b_init <- sapply(1:ncol(covariates), function(j){
      if(stats::sd(covariates[,j]) == 0) {
        log(matrixStats::colMeans2(dat))
      } else {
        rep(1, nrow(dat))
      }
    })
    colnames(b_init) <- colnames(covariates)
    nat_offset_mat <- tcrossprod(covariates, b_init)
  } else {
    b_init <- NULL
    nat_offset_mat <- 0
  }

  nat_mat <- family$dat_to_nat(dat)
  residual_mat <- nat_mat - nat_offset_mat
  residual_mat <- sweep(residual_mat, 1, offset_vec, "-")

  svd_res <- .svd_truncated(residual_mat,
                            K = k,
                            symmetric = F,
                            rescale = F,
                            mean_vec = NULL,
                            sd_vec = NULL,
                            K_full_rank = F)
  x_init <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_init <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  init_theta <- tcrossprod(x_init, y_init) + nat_offset_mat
  init_theta <- sweep(init_theta, 1, offset_vec, "+")

  if(family$name %in% c("neg_binom", "neg_binom2")){
    glmgampoi_init <- glmGamPoi::overdispersion_mle(y = Matrix::t(dat),
                                                    mean = family$nat_to_canon(t(init_theta)),
                                                    global_estimate = T)
    nuisance_init <- rep(glmgampoi_init$estimate, p)
  } else {
    nuisance_init <- NULL
  }

  rownames(x_init) <- rownames(dat)
  rownames(y_init) <- colnames(dat)

  structure(list(x_mat = x_init, y_mat = y_init, b_mat = b_init,
                 covariates = covariates,
                 nuisance_param_vec = nuisance_init,
                 offset_vec = offset_vec),
            class = "eSVD")
}
