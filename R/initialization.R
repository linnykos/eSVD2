#' Initialize eSVD
#'
#' For each gene, this function estimates two ridge-regression penalized GLMs (using the
#' Poisson model) -- one using the \code{case_control_variable} and one without, and
#' both sets of coefficients as well as the p-value (according to a deviance test) is returned.
#' This p-value is on the log10-scale.
#'
#' @param dat                      Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
#'                                 and \eqn{p} columns represent genes.
#'                                 The rows and columns of the matrix should be named.
#' @param covariates               \code{matrix} object with \eqn{n} rows with the same rownames as \code{dat} where the columns
#'                                 represent the different covariates.
#'                                 Notably, this should contain only numerical columns (i.e., all categorical
#'                                 variables should have already been split into numerous indicator variables), and all the columns
#'                                 in \code{covariates} will (strictly speaking) be included in the eSVD matrix factorization model.
#' @param metadata_individual      \code{factor} vector of length \eqn{n} that denotes which cell originates from which individual.
#' @param bool_intercept           Boolean on whether or not an intercept will be included as a covariate.
#' @param case_control_variable    A string of the column name of \code{covariates} which depicts the case-control
#'                                 status of each cell. Notably, this should be a binary variable where a \code{1}
#'                                 is hard-coded to describe case, and a \code{0} to describe control.
#' @param k                        Number of latent dimensions.
#' @param lambda                   Penalty of the \code{mixed_effect_variables} when using \code{glmnet::glmnet} to
#'                                 initialize the coefficients.
#' @param library_size_variable    A string of the variable name (which must be in \code{covariates}) of which variable denotes the sequenced (i.e., observed) library size.
#' @param metadata_case_control    (Optional) vector of length \eqn{n} with values strictly 0 or 1 that denotes if a cell is from cases or controls.
#'                                 By default, this is set to \code{NULL} since the code will extract this information from \code{covariates}.
#' @param offset_variables         A vector of strings depicting which column names in \code{covariate} will
#'                                 be set to have a coefficient of \code{1} automatically (i.e., there will be no estimation
#'                                 of their coefficient).
#' @param verbose                  Integer
#'
#' @return \code{eSVD} object with elements \code{dat}, \code{covariates},
#' \code{initial_Reg} and \code{param}
#' @export
initialize_esvd <- function(dat,
                            covariates,
                            metadata_individual,
                            bool_intercept = F,
                            case_control_variable = NULL,
                            k = 30,
                            lambda = 0.01,
                            library_size_variable = "Log_UMI",
                            offset_variables = "Log_UMI",
                            metadata_case_control = NULL,
                            verbose = 0){
  stopifnot(inherits(dat, c("dgCMatrix", "matrix")),
            nrow(dat) == nrow(covariates),
            is.matrix(covariates),
            k <= ncol(dat), k > 0, k %% 1 == 0,
            lambda <= 1e4, lambda >= 1e-4,
            library_size_variable %in% colnames(covariates),
            is.null(case_control_variable) || case_control_variable %in% colnames(covariates),
            "Intercept" %in% colnames(covariates),
            all(is.null(metadata_case_control)) || (all(is.numeric(metadata_case_control)) && all(metadata_case_control %in% c(0,1))),
            all(is.factor(metadata_individual)))
  stopifnot(all(is.null(offset_variables)) ||
              (all(offset_variables %in% colnames(covariates)) && !"Intercept" %in% offset_variables))

  n <- nrow(dat); p <- ncol(dat)
  if(is.matrix(dat)) dat[is.na(dat)] <- 0
  param <- .format_param_initialize(bool_intercept = bool_intercept,
                                    case_control_variable = case_control_variable,
                                    k = k,
                                    lambda = lambda,
                                    library_size_variable = library_size_variable,
                                    offset_variables = offset_variables)

  if(verbose >= 1) print("Performing GLMs")
  z_mat <- .initialize_coefficient(bool_intercept = bool_intercept,
                                   covariates = covariates,
                                   dat = dat,
                                   lambda = lambda,
                                   offset_variables = offset_variables,
                                   verbose = verbose)

  eSVD_obj <- structure(list(dat = dat,
                             covariates = covariates,
                             param = param),
                        class = "eSVD")

  if(verbose >= 1) print("Computing residuals")
  eSVD_obj[["fit_Init"]] <- .initialize_residuals(
    covariates = covariates,
    dat = dat,
    k = k,
    z_mat = z_mat
  )

  eSVD_obj[["latest_Fit"]] <- "fit_Init"

  if(all(is.null(metadata_case_control)) & !is.null(case_control_variable)){
    metadata_case_control <- covariates[,case_control_variable]
  }
  eSVD_obj[["case_control"]] <- metadata_case_control
  eSVD_obj[["individual"]] <- metadata_individual

  eSVD_obj
}

#####################

.initialize_coefficient <- function(bool_intercept,
                                    covariates,
                                    dat,
                                    lambda,
                                    offset_variables,
                                    verbose = 0){
  n <- nrow(dat); p <- ncol(dat)
  covariates_tmp <- covariates[,which(colnames(covariates) != "Intercept"), drop = F]
  if(!is.null(offset_variables)){
    covariates_tmp <- covariates_tmp[,which(!colnames(covariates_tmp) %in% offset_variables), drop=F]
    offset_vec <- Matrix::rowSums(covariates[,offset_variables,drop = F])
  } else {
    offset_vec <- NULL
  }

  z_mat <- matrix(1, nrow = p, ncol = ncol(covariates))
  colnames(z_mat) <- colnames(covariates)
  rownames(z_mat) <- colnames(dat)

  for(j in 1:p){
    if(verbose == 1 && p >= 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(paste0("Working on variable ", j , " of ", p))

    if(ncol(covariates_tmp) > 1){
      glm_fit <- glmnet::glmnet(x = covariates_tmp,
                                y = as.numeric(dat[,j]),
                                family = "poisson",
                                offset = offset_vec,
                                alpha = 0,
                                standardize = F,
                                intercept = bool_intercept,
                                lambda = exp(seq(log(1e4), log(lambda), length.out = 100)))

      if(bool_intercept){
        z_mat[j, c("Intercept", colnames(covariates_tmp))] <- c(glm_fit$a0[length(glm_fit$a0)], glm_fit$beta[,ncol(glm_fit$beta)])
      } else {
        z_mat[j, c("Intercept", colnames(covariates_tmp))] <- c(0, glm_fit$beta[,ncol(glm_fit$beta)])
      }
    } else {
      # Handle corner case when there are no covariates to adjust for

      df <- as.data.frame(cbind(
        as.numeric(dat[,j]), covariates_tmp
      ))
      colnames(df) <- "y"

      if(bool_intercept){
        glm_fit <- stats::glm(y ~ ., data = df, family = stats::poisson)
        z_mat[j, c("Intercept", colnames(covariates_tmp))] <- stats::coef(glm_fit)
      } else {
        glm_fit <- stats::glm(y ~ . - 1, data = df, family = stats::poisson)
        z_mat[j, c("Intercept", colnames(covariates_tmp))] <- c(0, stats::coef(glm_fit))
      }
    }
  }

  z_mat
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
  rownames(z_mat) <- colnames(dat)
  colnames(z_mat) <- colnames(covariates)

  .form_esvd_fit(x_mat = x_mat, y_mat = y_mat, z_mat = z_mat)
}

.format_param_initialize <- function(bool_intercept,
                                     case_control_variable,
                                     k,
                                     lambda,
                                     library_size_variable,
                                     offset_variables) {
  list(init_bool_intercept = bool_intercept,
       init_case_control_variable = case_control_variable,
       init_k = k,
       init_lambda = lambda,
       init_library_size_variable = library_size_variable,
       init_offset_variables = offset_variables)
}
