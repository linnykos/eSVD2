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
#' \code{param}, \code{fit_Init} (an \code{eSVD_Fit} with \code{x_mat},
#' \code{y_mat}, \code{z_mat}), \code{latest_Fit} (\code{"fit_Init"}),
#' \code{case_control} and \code{individual}
#' @examples
#' set.seed(10)
#' sim <- generate_null(cell_per_person = 15, num_genes = 40,
#'                      num_individuals = 8)
#' esvd_obj <- initialize_esvd(dat = sim$obs_mat,
#'                             covariates = sim$covariates,
#'                             metadata_individual = sim$metadata_individual,
#'                             case_control_variable = "CC",
#'                             bool_intercept = TRUE,
#'                             k = 2,
#'                             lambda = 0.1)
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_Init",
#'                                                omitted_variables = "Log_UMI")
#' names(esvd_obj)
#' dim(esvd_obj$fit_Init$x_mat)
#' @export
initialize_esvd <- function(dat,
                            covariates,
                            metadata_individual,
                            bool_intercept = FALSE,
                            case_control_variable = NULL,
                            k = 30,
                            lambda = 0.01,
                            library_size_variable = "Log_UMI",
                            offset_variables = "Log_UMI",
                            metadata_case_control = NULL,
                            verbose = 0){
  stopifnot(inherits(dat, c("dgCMatrix", "matrix")),
            nrow(dat) == nrow(covariates),
            is.matrix(covariates))
  if(length(k) != 1 || k %% 1 != 0 || k <= 0 || k > ncol(dat)){
    stop("`k` = ", paste0(k, collapse = ", "), " must be a positive integer ",
         "no larger than the number of genes, ncol(dat) = ", ncol(dat))
  }
  stopifnot(lambda <= 1e4, lambda >= 1e-4,
            library_size_variable %in% colnames(covariates),
            is.null(case_control_variable) || case_control_variable %in% colnames(covariates),
            "Intercept" %in% colnames(covariates),
            all(is.null(metadata_case_control)) || (all(is.numeric(metadata_case_control)) && all(metadata_case_control %in% c(0,1))),
            all(is.factor(metadata_individual)))
  stopifnot(all(is.null(offset_variables)) ||
              (all(offset_variables %in% colnames(covariates)) && !"Intercept" %in% offset_variables))

  n <- nrow(dat); p <- ncol(dat)
  if(anyNA(covariates) || any(!is.finite(covariates))){
    stop("`covariates` contains NA or non-finite entries")
  }
  qr_res <- qr(covariates)
  if(qr_res$rank < ncol(covariates)){
    # A rank-deficient design is fit by glmnet's ridge without complaint,
    # but the reparameterization step cannot regress on it, and the paper
    # warns that including every individual's indicator makes the design
    # collinear with the intercept. Refuse it here with the column names.
    aliased_vec <- colnames(covariates)[qr_res$pivot[-seq_len(qr_res$rank)]]
    stop("`covariates` is rank deficient (rank ", qr_res$rank, " for ",
         ncol(covariates), " columns): `", paste0(aliased_vec, collapse = "`, `"),
         "` can be written as a combination of the other columns. Remove ",
         "one variable from each collinear set. Including an indicator for ",
         "every individual, or for every level of a factor, causes this")
  }
  # NAs are zeroed on both storage types. The sparse branch used to be
  # skipped, and a single NA in a dgCMatrix then errored inside glmnet.
  if(is.matrix(dat)){
    dat[is.na(dat)] <- 0
  } else if(anyNA(dat@x)){
    dat@x[is.na(dat@x)] <- 0
    dat <- Matrix::drop0(dat)
  }

  # An all-zero gene makes glmnet warn and return the wrong lambda's fit,
  # sends gamma_rate to its clamp, and gives a 0/0 test statistic. The
  # helper `eSVD_helper` removes such genes before the pipeline; reaching
  # here with one is refused rather than silently mis-fitted (Q-STATUS-1).
  all_zero_idx <- .which_all_zero(dat)
  if(length(all_zero_idx) > 0){
    stop(length(all_zero_idx), " gene(s) are all zero (",
         paste0(utils::head(colnames(dat)[all_zero_idx], 5), collapse = ", "),
         if(length(all_zero_idx) > 5) ", ..." else "",
         "); remove them before initialization, or run `eSVD_helper`, ",
         "which removes them and records the removal in `gene_status`")
  }
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

  if(all(is.null(metadata_case_control)) && !is.null(case_control_variable)){
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
  covariates_tmp <- covariates[,which(colnames(covariates) != "Intercept"), drop = FALSE]
  if(!is.null(offset_variables)){
    covariates_tmp <- covariates_tmp[,which(!colnames(covariates_tmp) %in% offset_variables), drop=F]
    offset_vec <- Matrix::rowSums(covariates[,offset_variables,drop = FALSE])
  } else {
    offset_vec <- NULL
  }

  z_mat <- matrix(1, nrow = p, ncol = ncol(covariates))
  colnames(z_mat) <- colnames(covariates)
  rownames(z_mat) <- colnames(dat)

  for(j in seq_len(p)){
    if(verbose == 1 && p >= 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(paste0("Working on variable ", j , " of ", p))

    if(ncol(covariates_tmp) > 1){
      glm_fit <- glmnet::glmnet(x = covariates_tmp,
                                y = as.numeric(dat[,j]),
                                family = "poisson",
                                offset = offset_vec,
                                alpha = 0,
                                standardize = FALSE,
                                intercept = bool_intercept,
                                lambda = exp(seq(log(1e4), log(lambda), length.out = 100)))

      if(bool_intercept){
        z_mat[j, c("Intercept", colnames(covariates_tmp))] <- c(glm_fit$a0[length(glm_fit$a0)], glm_fit$beta[,ncol(glm_fit$beta)])
      } else {
        z_mat[j, c("Intercept", colnames(covariates_tmp))] <- c(0, glm_fit$beta[,ncol(glm_fit$beta)])
      }
    } else {
      # `glmnet` needs at least two predictor columns, so with zero or one
      # covariate left to estimate the fit is an unpenalized `stats::glm`.
      # The offset (the library size) has to be carried into this branch
      # too: it used to be dropped here, so a design with only the
      # case-control indicator initialized every gene's intercept without
      # adjusting for sequencing depth.
      y_vec <- as.numeric(dat[,j])
      glm_offset_vec <- if(is.null(offset_vec)) rep(0, n) else offset_vec

      if(ncol(covariates_tmp) == 1){
        x_vec <- covariates_tmp[,1]
        if(bool_intercept){
          glm_fit <- stats::glm(y_vec ~ x_vec + offset(glm_offset_vec),
                                family = stats::poisson)
          coef_vec <- stats::coef(glm_fit)
        } else {
          glm_fit <- stats::glm(y_vec ~ 0 + x_vec + offset(glm_offset_vec),
                                family = stats::poisson)
          coef_vec <- c(0, stats::coef(glm_fit))
        }
      } else {
        # Nothing to estimate beyond the intercept.
        if(bool_intercept){
          glm_fit <- stats::glm(y_vec ~ 1 + offset(glm_offset_vec),
                                family = stats::poisson)
          coef_vec <- stats::coef(glm_fit)
        } else {
          coef_vec <- 0
        }
      }
      z_mat[j, c("Intercept", colnames(covariates_tmp))] <- unname(coef_vec)
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
                       check_stability = TRUE,
                       K = k,
                       mean_vec = NULL,
                       rescale = FALSE,
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
