.identification <- function(cov_x, cov_y, check = FALSE, tol = 1e-6){
  stopifnot(all(dim(cov_x) == dim(cov_y)), nrow(cov_x) == ncol(cov_x))
  if(nrow(cov_x) == 1){
    return(matrix((as.numeric(cov_y)/as.numeric(cov_x))^(1/4), 1, 1))
  }

  eigen_x <- eigen(cov_x)
  eigen_y <- eigen(cov_y)

  Vx <- eigen_x$vectors
  Vy <- eigen_y$vectors

  if(any(eigen_x$values <= tol) || any(eigen_y$values <= tol))
    warning("Detecting rank deficiency in reparameterization step")

  # Proceed after the warning (question Q-REP-1). An eigenvalue at or below
  # `tol` is floored at `tol` rather than inverted: `Dx^(-1/2)` of an exact
  # zero is Inf and of a rounding-error negative is NaN, and either used to
  # kill the `eigen()` below with "infinite or missing values". The floored
  # direction carries no signal in `x_mat` (that is what the zero eigenvalue
  # means), so the product `x_mat %*% t(y_mat)` is still preserved.
  Dx <- pmax(eigen_x$values, tol)
  Dy <- pmax(eigen_y$values, tol)

  # form R
  tmp <- crossprod(.mult_mat_vec(Vy, sqrt(Dy)), .mult_mat_vec(Vx, sqrt(Dx)))
  svd_tmp <- svd(tmp)
  Q <- tcrossprod(svd_tmp$u, svd_tmp$v)

  # run a check
  if(check){
    sym_mat <- crossprod(Q, tmp)
    stopifnot(sum(abs(sym_mat - t(sym_mat))) <= 1e-6)
  }

  # now form the symmetric matrix to later factorize
  sym_prod <- tcrossprod(tcrossprod(.mult_mat_vec(Vx, Dx^(-1/2)), Q), .mult_mat_vec(Vy, sqrt(Dy)))
  sym_prod[which(abs(sym_prod) <= tol)] <- 0

  if(check){
    stopifnot(sum(abs(sym_prod - t(sym_prod))) <= 1e-6)
  }

  eigen_sym <- eigen(sym_prod)
  # `sym_prod` is positive semi-definite in exact arithmetic; after the
  # flooring above a rounding-error negative eigenvalue would give NaN.
  W_mat <- .mult_vec_mat(sqrt(pmax(eigen_sym$values, 0)), t(eigen_sym$vectors))

  if(check){
    mat1 <- tcrossprod(W_mat %*% cov_x, W_mat)

    W_mat_inv <- solve(W_mat)
    mat2 <- crossprod(W_mat_inv, cov_y) %*% W_mat_inv
    stopifnot(sum(abs(mat1 - mat2)) <= 1e-6)
  }

  # adjust the transformation so it yields a diagonal matrix
  eig_res <- eigen(tcrossprod(W_mat %*% cov_x, W_mat))

  crossprod(eig_res$vectors, W_mat)
}


#' Function to reparameterize two matrices
#'
#' Designed to output matrices of the same dimension as \code{x_mat}
#' and \code{y_mat}, but linearly transformed so \code{x_mat \%*\% t(y_mat)}
#' is preserved but either the \eqn{k \times k} matrix \code{t(x_mat) \%*\% x_mat} is diagonal and equal to
#' \code{t(y_mat) \%*\% y_mat} (if \code{equal_covariance} is \code{FALSE})
#' or \code{t(x_mat) \%*\% x_mat/nrow(x_mat)} is diagonal and equal to
#' \code{t(y_mat) \%*\% y_mat/nrow(y_mat)} (if \code{equal_covariance} is \code{TRUE})
#'
#' @param x_mat matrix of dimension \code{n} by \code{k}
#' @param y_mat matrix of dimension \code{p} by \code{k}
#' @param equal_covariance boolean
#'
#' @return list of two matrices, \code{x_mat} and \code{y_mat}
#' @noRd
.reparameterize <- function(x_mat, y_mat, equal_covariance){
  stopifnot(ncol(x_mat) == ncol(y_mat))
  n <- nrow(x_mat); p <- nrow(y_mat)

  res <- .identification(crossprod(x_mat), crossprod(y_mat))

  if(equal_covariance){
    list(x_mat = (n/p)^(1/4)*tcrossprod(x_mat, res), y_mat = (p/n)^(1/4)*y_mat %*% solve(res))
  } else {
    list(x_mat = tcrossprod(x_mat, res), y_mat = y_mat %*% solve(res))
  }
}

.factorize_matrix <- function(mat, k, equal_covariance){
  stopifnot(k <= min(dim(mat)))

  svd_res <- .svd_safe(mat = mat,
                       check_stability = TRUE,
                       K = k,
                       mean_vec = NULL,
                       rescale = FALSE,
                       scale_max = NULL,
                       sd_vec = NULL)
  x_mat <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_mat <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  .reparameterize(x_mat, y_mat, equal_covariance = equal_covariance)
}

#' Reparameterize eSVD object
#'
#' @param input_obj           \code{eSVD} object, after either \code{initialize_esvd()} or
#'                            \code{opt_esvd()}
#' @param fit_name            The name of the fit in \code{input_obj} that you wish to reparameterize. This should be in \code{names(input_obj)}
#' @param omitted_variables   Either \code{NULL} (the default) or variables in \code{input_obj$covariates} that should not be reparameterized
#' @param verbose             Integer.
#'
#' @return \code{eSVD} object after adjusting the fit in \code{fit_name}.
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
#' # x_mat is now orthogonal to the retained covariates
#' max(abs(crossprod(esvd_obj$fit_Init$x_mat,
#'                   esvd_obj$covariates[, "Sex"])))
#' @export
reparameterization_esvd_covariates <- function(input_obj,
                                               fit_name,
                                               omitted_variables = NULL,
                                               verbose = 0){
  if(length(fit_name) != 1 || !fit_name %in% names(input_obj)){
    available_fits <- names(input_obj)[sapply(input_obj, inherits, "eSVD_Fit")]
    stop("`fit_name` = \"", paste0(fit_name, collapse = "\", \""),
         "\" is not a fit in `input_obj`; the available fits are \"",
         paste0(available_fits, collapse = "\", \""), "\"")
  }
  stopifnot(is.null(omitted_variables) || all(omitted_variables %in% colnames(input_obj$covariates)))

  x_mat <- input_obj[[fit_name]]$x_mat
  y_mat <- input_obj[[fit_name]]$y_mat
  z_mat <- input_obj[[fit_name]]$z_mat
  k <- ncol(x_mat)
  p <- nrow(y_mat)

  covariate_mat <- input_obj$covariates
  stopifnot("Intercept" %in% colnames(covariate_mat))
  covariate_mat <- covariate_mat[,which(!colnames(covariate_mat) %in% omitted_variables),drop = FALSE]

  # Each latent factor is regressed on the retained covariates (intercept
  # included) and replaced by its residual; the fitted part is folded into
  # the covariate coefficients so that x_mat %*% t(y_mat) + C %*% t(z_mat)
  # is unchanged. This is plain least squares on the design matrix. It used
  # to go through `stats::lm(x ~ ., data = as.data.frame(...))`, which
  # (a) rewrote non-syntactic covariate names such as "Diagnosis_ASD (severe)"
  # via make.names(), so `z_mat[, names(coef)]` was a subscript error, and
  # (b) returned NA coefficients for aliased columns, which then propagated
  # silently into z_mat and every posterior (additional_context/CRAN_READINESS.md
  # section 1.3, not shipped with the package).
  qr_res <- qr(covariate_mat)
  if(qr_res$rank < ncol(covariate_mat)){
    aliased_vec <- colnames(covariate_mat)[qr_res$pivot[-seq_len(qr_res$rank)]]
    stop("the covariates retained for reparameterization are collinear: `",
         paste0(aliased_vec, collapse = "`, `"),
         "` can be written as a combination of the others. Remove one ",
         "variable from each collinear set (or name it in ",
         "`omitted_variables`)")
  }
  coef_mat <- qr.coef(qr_res, x_mat)         # ncol(covariate_mat) x k
  rownames(coef_mat) <- colnames(covariate_mat)
  stopifnot(all(is.finite(coef_mat)))
  fitted_mat <- qr.fitted(qr_res, x_mat)
  if(verbose > 0){
    r2_vec <- 1 - colSums((x_mat - fitted_mat)^2) /
      colSums(sweep(x_mat, 2, colMeans(x_mat))^2)
    for(ell in seq_len(k)){
      print(paste0(ell, ": R2 of ", round(r2_vec[ell], 2)))
    }
  }

  # z_mat[, retained] += y_mat %*% t(coef_mat), matched by name.
  z_mat[, rownames(coef_mat)] <- z_mat[, rownames(coef_mat), drop = FALSE] +
    tcrossprod(y_mat, coef_mat)
  x_mat <- x_mat - fitted_mat

  res <- .reparameterize(x_mat, y_mat, equal_covariance = TRUE)
  x_mat <- res$x_mat; y_mat <- res$y_mat

  input_obj[[fit_name]]$x_mat <- x_mat
  input_obj[[fit_name]]$y_mat <- y_mat
  input_obj[[fit_name]]$z_mat <- z_mat

  input_obj
}


#########

.svd_safe <- function(mat,
                      check_stability, # boolean
                      K, # positive integer
                      mean_vec, # boolean, NULL or vector
                      rescale, # boolean
                      scale_max, # NULL or positive integer
                      sd_vec){ # boolean, NULL or vector
  if(is.na(K)) K <- min(dim(mat))
  stopifnot(min(dim(mat)) >= K)

  mean_vec <- .compute_matrix_mean(mat, mean_vec)
  sd_vec <- .compute_matrix_sd(mat, sd_vec)

  res <- .svd_in_sequence(check_stability = check_stability,
                          K = K,
                          mat = mat,
                          mean_vec = mean_vec,
                          scale_max = scale_max,
                          sd_vec = sd_vec)
  res <- list(d = res$d, u = res$u, v = res$v, method = res$method)
  class(res) <- "svd"

  # pass row-names and column-names
  if(length(rownames(mat)) > 0) rownames(res$u) <- rownames(mat)
  if(length(colnames(mat)) > 0) rownames(res$v) <- colnames(mat)

  # useful only if your application requires only the singular vectors
  # if the number of rows or columns is too large, the singular vectors themselves
  # are often a bit too small numerically
  if(rescale){
    n <- nrow(mat); p <- ncol(mat)
    res$u <- res$u * sqrt(n)
    res$v <- res$v * sqrt(p)
    res$d <- res$d / (sqrt(n)*sqrt(p))
  }

  res
}

.svd_in_sequence <- function(check_stability,
                             K,
                             mat,
                             mean_vec,
                             scale_max,
                             sd_vec){
  res <- tryCatch({
    .irlba_custom(check_stability = check_stability,
                  K = K,
                  mat = mat,
                  mean_vec = mean_vec,
                  scale_max = scale_max,
                  sd_vec = sd_vec)
  },
  warning = function(e){NULL},
  error = function(e){NULL})
  if(!all(is.null(res))) {res$method <- "irlba"; return(res)}

  ##

  res <- tryCatch({
    .rpsectra_custom(check_stability = check_stability,
                     K = K,
                     mat = mat,
                     mean_vec = mean_vec,
                     scale_max = scale_max,
                     sd_vec = sd_vec)
  },
  warning = function(e){NULL},
  error = function(e){NULL})
  if(!all(is.null(res))) {res$method <- "RSpectra"; return(res)}

  ##

  if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
  if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
  if(!is.null(scale_max)){
    mat[mat > abs(scale_max)] <- abs(scale_max)
    mat[mat < -abs(scale_max)] <- -abs(scale_max)
  }
  res <- svd(mat)
  res$method <- "base"
  res
}

#' Deterministic starting vector for the iterative SVD solvers
#'
#' \code{irlba::irlba} and \code{RSpectra::svds} both start from a random
#' vector drawn from the user's RNG stream. Two runs of the pipeline on
#' identical input therefore differed by about 1e-9 after initialization,
#' and the alternating optimization amplified that to O(1) differences in
#' the test statistics of the most strongly DE genes. A fixed start removes
#' the run-to-run variation (though not the sensitivity it exposed).
#'
#' A golden-ratio (Weyl) sequence pushed through \code{qnorm} gives a
#' pseudo-random-looking vector that is deterministic and touches no RNG
#' state, so a user's \code{set.seed()} is not consumed.
#'
#' @param n Length.
#'
#' @returns Unit-norm numeric vector of length \code{n}.
#' @noRd
.svd_start_vector <- function(n){
  uniform_vec <- (seq_len(n) * (sqrt(5) - 1) / 2 + 0.5) %% 1
  start_vec <- stats::qnorm(uniform_vec)

  start_vec / sqrt(sum(start_vec^2))
}

.irlba_custom <- function(check_stability,
                          K,
                          mat,
                          mean_vec,
                          scale_max,
                          sd_vec){
  start_vec <- .svd_start_vector(ncol(mat))

  if(inherits(mat, "dgCMatrix")){
    if(!all(is.null(scale_max))) warning("scale_max does not work with sparse matrices when using irlba")
    tmp <- irlba::irlba(A = mat,
                        nv = K,
                        work = min(c(K + 10, dim(mat))),
                        scale = sd_vec,
                        center = mean_vec,
                        v = start_vec)

    if(check_stability && K > 5) {
      tmp2 <- irlba::irlba(A = mat,
                           nv = 5,
                           scale = sd_vec,
                           center = mean_vec,
                           v = start_vec)
      ratio_vec <- tmp2$d/tmp$d[1:5]
      if(any(ratio_vec > 2) || any(ratio_vec < 1/2)) warning("irlba is potentially unstable")
    }

    return(tmp)

  } else {
    if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
    if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
    if(!is.null(scale_max)){
      mat[mat > abs(scale_max)] <- abs(scale_max)
      mat[mat < -abs(scale_max)] <- -abs(scale_max)
    }

    tmp <- irlba::irlba(A = mat, nv = K, v = start_vec)

    if(check_stability && K > 5) {
      tmp2 <- irlba::irlba(A = mat, nv = 5, v = start_vec)
      ratio_vec <- tmp2$d/tmp$d[1:5]
      if(any(ratio_vec > 2) || any(ratio_vec < 1/2)) warning("irlba is potentially unstable")
    }

    return(tmp)
  }
}

.rpsectra_custom <- function(check_stability,
                             K,
                             mat,
                             mean_vec,
                             scale_max,
                             sd_vec){

  if(inherits(mat, "dgCMatrix")){
    if(!all(is.null(mean_vec))) warning("mean_vec does not work with sparse matrices when using RSpectra")
    if(!all(is.null(sd_vec))) warning("sd_vec does not work with sparse matrices when using RSpectra")
    if(!all(is.null(scale_max))) warning("scale_max does not work with sparse matrices when using RSpectra")
  } else {
    if(!all(is.null(mean_vec))) mat <- sweep(mat, MARGIN = 2, STATS = mean_vec, FUN = "-")
    if(!all(is.null(sd_vec))) mat <- sweep(mat, MARGIN = 2, STATS = sd_vec, FUN = "/")
    if(!is.null(scale_max)){
      mat[mat > abs(scale_max)] <- abs(scale_max)
      mat[mat < -abs(scale_max)] <- -abs(scale_max)
    }
  }

  # RSpectra works on the smaller Gram matrix, so the start has that length.
  svds_opts <- list(initvec = .svd_start_vector(min(dim(mat))))
  tmp <- RSpectra::svds(A = mat, k = K, opts = svds_opts)

  if(check_stability && K > 5) {
    tmp2 <- RSpectra::svds(A = mat, k = 5, opts = svds_opts)
    ratio_vec <- tmp2$d/tmp$d[1:5]
    if(any(ratio_vec > 2) || any(ratio_vec < 1/2)) warning("RSpectra is potentially unstable")
  }

  tmp
}

# `mean_vec` is either NULL, a single logical (TRUE: compute the column means;
# FALSE: do not center) or a full numeric vector of length ncol(mat). The
# guard tests `is.logical()` rather than `length() == 1`, because `0.5` and
# `TRUE` are both length 1 and `if(0.5)` would silently coerce (Q-SVD-2).
.compute_matrix_mean <- function(mat, mean_vec){
  if(is.null(mean_vec)) return(NULL)

  if(is.logical(mean_vec) && length(mean_vec) == 1){
    if(mean_vec){
      mean_vec <- Matrix::colMeans(mat)
    } else{
      mean_vec <- NULL
    }
  } else if(!is.numeric(mean_vec) || length(mean_vec) != ncol(mat)){
    stop("`mean_vec` must be NULL, a single TRUE/FALSE, or a numeric vector ",
         "of length ncol(mat) = ", ncol(mat), "; got a ", class(mean_vec)[1],
         " of length ", length(mean_vec))
  }

  mean_vec
}

.compute_matrix_sd <- function(mat, sd_vec){
  if(is.null(sd_vec)) return(NULL)

  if(is.logical(sd_vec) && length(sd_vec) == 1){
    if(sd_vec){
      if(inherits(x = mat, what = 'dgCMatrix')){
        sd_vec <- .sparse_col_sds(mat)
      } else {
        sd_vec <- matrixStats::colSds(as.matrix(mat))
      }
    } else{
      sd_vec <- NULL
    }
  } else if(!is.numeric(sd_vec) || length(sd_vec) != ncol(mat)){
    stop("`sd_vec` must be NULL, a single TRUE/FALSE, or a numeric vector ",
         "of length ncol(mat) = ", ncol(mat), "; got a ", class(sd_vec)[1],
         " of length ", length(sd_vec))
  }

  sd_vec
}

#' Column standard deviations of a sparse matrix, without densifying
#'
#' Replaces \code{sparseMatrixStats::colSds}, which is Bioconductor-only.
#' The naive two-pass form
#' \code{sqrt((colSums(x^2) - n*colMeans(x)^2)/(n-1))} suffers catastrophic
#' cancellation on a column with a large mean and a small variance and
#' returns \code{NaN}; this sums squared deviations over the stored non-zeros
#' and adds the contribution of the implied zeros exactly:
#' \code{sum_i (x_i - m)^2 = sum_stored (x_i - m)^2 + (n - nnz) * m^2}.
#'
#' @param mat  \code{dgCMatrix}.
#'
#' @returns Numeric vector of length \code{ncol(mat)}.
#' @noRd
.sparse_col_sds <- function(mat){
  stopifnot(inherits(mat, "dgCMatrix"))

  n <- nrow(mat)
  col_mean_vec <- Matrix::colMeans(mat)
  num_stored_vec <- diff(mat@p)
  col_idx_vec <- rep(seq_len(ncol(mat)), num_stored_vec)
  deviation_vec <- (mat@x - col_mean_vec[col_idx_vec])^2
  stored_ss_vec <- as.numeric(tapply(deviation_vec,
                                     factor(col_idx_vec,
                                            levels = seq_len(ncol(mat))),
                                     sum))
  stored_ss_vec[is.na(stored_ss_vec)] <- 0

  sqrt((stored_ss_vec + (n - num_stored_vec) * col_mean_vec^2) / (n - 1))
}


