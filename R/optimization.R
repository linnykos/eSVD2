# Optimization for eSVD

#' Main optimization function for eSVD
#'
#' @param x_init             initial estimate of a \eqn{n \times k}{n × k} matrix for the \eqn{k}-dimensional
#'                           embedding for the \eqn{n} cells
#' @param y_init             initial estimate of a \eqn{p \times k}{p × k} matrix for the \eqn{k}-dimensional
#'                           embedding for the \eqn{p} genes
#' @param dat                dataset where the \eqn{n} rows represent cells and \eqn{p} columns represent genes
#' @param family             a character string, one of \code{"gaussian"}, \code{"exponential"},
#'                           \code{"poisson"}, \code{"neg_binom"}, \code{"curved_gaussian"},
#'                           and \code{"bernoulli"}
#' @param method             a character string indicating the optimization method,
#'                           either \code{"newton"} or \code{"lbfgs"}
#' @param b_init             initial estimate of a \eqn{n \times d}{n × d} matrix for the \eqn{d}-dimensional
#'                           embedding for the \eqn{n} cells
#' @param covariates         an \eqn{n \times d}{n × d} matrix representing the additional \eqn{d} covariates,
#'                           or \code{NULL} if no covariate is given
#' @param nuisance_param_vec either \code{NA} or a single numeric or a length-\eqn{p}
#'                           vector of numerics representing nuisance parameters (for \code{family = "neg_binom"}
#'                           or \code{family = "curved_gausian"})
#' @param library_size_vec   either \code{NA} or a single numeric (default is \code{1}) or
#'                           a length-\eqn{n} vector of numerics. If \code{NA}, the library size will be estimated
#' @param offset_vec         a vector of length-\eqn{n} that represents a constant amount added to each row of the
#'                           natural parameter matrix
#' @param l2pen              the ridge penalty parameter for (X, Y, B)
#' @param reparameterize     reparameterize \code{x_mat} and \code{y_mat} after every iteration
#' @param reestimate_nuisance a boolean for whether the nuisance parameter is reestimate after every iteration
#' @param global_estimate    a boolean for whether or not the same nuisance parameter is used for all genes
#' @param min_nuisance       a small positive number for the minimum nuisance parameter value
#' @param max_nuisance       a large positive number for the maximum nuisance parameter value
#' @param max_iter           a positive integer giving the maximum number of iterations
#' @param tol                a small positive number for the tolerance of optimization error
#' @param verbose            a non-negative integer to indicate the verbosity of messages
#' @param ...                additional parameters, currently not used
#'
#' @return A list with elements \code{x_mat} and \code{y_mat} (and others), representing the two
#'         latent matrices.
#' @export
opt_esvd <- function(x_init,
                     y_init,
                     dat,
                     family = "gaussian",
                     b_init = NULL,
                     covariates = NULL,
                     library_size_vec = NA,
                     nuisance_param_vec = NA,
                     offset_vec = rep(0, nrow(x_init)),
                     bool_run_cpp = T,
                     gene_group_factor = factor(rep("1", ncol(dat))),
                     l2pen = 0,
                     max_cell_subsample = 10*nrow(dat),
                     max_iter = 100,
                     method = c("newton", "lbfgs"),
                     nuisance_value_lower = NA,
                     nuisance_value_upper = NA,
                     reparameterize = F,
                     reestimate_nuisance = F,
                     reestimate_nuisance_per_iteration = 1,
                     tol = 1e-6,
                     verbose = 0,
                     ...)
{
  n <- nrow(dat)
  p <- ncol(dat)
  k <- ncol(x_init)
  stopifnot(
    nrow(x_init) == n, nrow(y_init) == p, ncol(y_init) == k,
    is.character(family), sum(!is.na(dat)) > 0
  )

  # Whether to use C++ code
  if(bool_run_cpp)
  {
    load_cpp_code()
    # C++ code needs double type
    storage.mode(x_init) <- "double"
    storage.mode(y_init) <- "double"
    storage.mode(dat)    <- "double"
  }

  param <- .opt_esvd_format_param(family = family,
                                  gene_group_factor = gene_group_factor,
                                  l2pen = l2pen,
                                  max_cell_subsample = max_cell_subsample,
                                  max_iter = max_iter,
                                  method = method,
                                  nuisance_value_lower = nuisance_value_lower,
                                  nuisance_value_upper = nuisance_value_upper,
                                  reparameterize = reparameterize,
                                  reestimate_nuisance = reestimate_nuisance,
                                  reestimate_nuisance_per_iteration = reestimate_nuisance_per_iteration,
                                  tol = tol,
                                  verbose = verbose)

  # Convert family string to internal family object, e.g. `.esvd.poisson` and `.esvd.neg_binom2`
  family <- .string_to_distr_funcs(family)
  # Parse library size and nuisance parameter
  library_size_vec <- .parse_library_size(dat, library_size_vec)
  if(all(!is.na(nuisance_param_vec)) && length(nuisance_param_vec) == 1) {
    nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
  }
  library_size_vec <- as.numeric(library_size_vec)
  nuisance_param_vec <- as.numeric(nuisance_param_vec)
  # Parse optimization method
  method <- match.arg(method)
  opt_fun <- if(method == "newton") constr_newton else constr_lbfgs

  # Initialize embedding matrices
  b_mat <- .opt_esvd_setup_b_mat(b_init, covariates)
  x_mat <- x_init
  y_mat <- y_init

  losses <- c()
  for(i in seq_len(param$max_iter))
  {
    if(verbose >= 1)
      cat("========== eSVD Iter ", i, " ==========\n\n", sep = "")
    # Optimize X given Y and B
    x_mat <- opt_x(X0 = x_mat,
                   Y = y_mat,
                   B = b_mat,
                   Z = covariates,
                   A = dat,
                   family = family,
                   s = library_size_vec,
                   gamma = nuisance_param_vec,
                   offset_vec = offset_vec,
                   l2pen = param$l2pen,
                   opt_fun = opt_fun,
                   gene_group_factor = gene_group_factor,
                   verbose = verbose, ...)

    # Optimize Y and B given X
    yb_mat <- cbind(y_mat, b_mat)
    xz_mat <- cbind(x_mat, covariates)
    yb_mat <- opt_yb(yb_mat,
                     XZ = xz_mat,
                     A = dat,
                     family = family,
                     s = library_size_vec,
                     gamma = nuisance_param_vec,
                     offset_vec = offset_vec,
                     l2pen = param$l2pen,
                     opt_fun = opt_fun,
                     verbose = verbose, ...)

    # Split Y and B
    if(is.null(covariates))
    {
      y_mat <- yb_mat
      b_mat <- NULL
    } else {
      y_mat <- yb_mat[, 1:k, drop = FALSE]
      b_mat <- yb_mat[, -(1:k), drop = FALSE]
    }

    if(param$reestimate_nuisance &
       family$name == "neg_binom2" &
       i %% param$reestimate_nuisance_per_iteration == 0){
      nuisance_param_vec <- .opt_nuisance(covariates = covariates,
                                          dat = dat,
                                          gene_group_factor = param$gene_group_factor,
                                          max_cell_subsample = param$max_cell_subsample,
                                          offset_vec = offset_vec,
                                          x_mat = x_mat,
                                          yb_mat = yb_mat,
                                          value_lower = param$nuisance_value_lower,
                                          value_upper = param$nuisance_value_upper,
                                          verbose = verbose)
    }

    if(param$reparameterize){
      tmp <- tryCatch(.reparameterize(x_mat, y_mat, equal_covariance = T),
                      error = function(e){list(x_mat = x_mat, y_mat = y_mat)})
      x_mat <- tmp$x_mat
      y_mat <- tmp$y_mat
    }

    # Loss function
    loss <- objfn_all(X = x_mat,
                      Y = y_mat,
                      B = b_mat,
                      Z = covariates,
                      A = dat,
                      family = family,
                      s = library_size_vec,
                      gamma = nuisance_param_vec,
                      offset = offset_vec,
                      l2pen = param$l2pen)
    losses <- c(losses, loss)
    if(verbose >= 1)
      cat("========== eSVD Iter ", i, ", loss = ", loss, " ==========\n\n", sep = "")

    # Convergence test
    resid <- abs(losses[i] - losses[i - 1])
    thresh <- param$tol * max(1, abs(losses[i - 1]))
    if(i >=2 && resid <= thresh)
      break
  }

  # update nuisances one last time
  if(param$reestimate_nuisance & family$name == "neg_binom2"){
    nuisance_param_vec <- .opt_nuisance(covariates = covariates,
                                        dat = dat,
                                        gene_group_factor = param$gene_group_factor,
                                        max_cell_subsample = param$max_cell_subsample,
                                        offset_vec = offset_vec,
                                        x_mat = x_mat,
                                        yb_mat = yb_mat,
                                        value_lower = param$nuisance_value_lower,
                                        value_upper = param$nuisance_value_upper,
                                        verbose = verbose)
  }

  tmp <- tryCatch(.reparameterize(x_mat, y_mat, equal_covariance = T),
                  error = function(e){list(x_mat = x_mat, y_mat = y_mat)})
  x_mat <- tmp$x_mat
  y_mat <- tmp$y_mat

  tmp <- .opt_esvd_format_matrices(b_mat, covariates, dat, x_mat, y_mat)

  list(x_mat = tmp$x_mat,
       y_mat = tmp$y_mat,
       covariates = covariates,
       b_mat = tmp$b_mat,
       loss = losses,
       nuisance_param_vec = nuisance_param_vec,
       library_size_vec = library_size_vec,
       offset_vec = offset_vec,
       param = param)
}
