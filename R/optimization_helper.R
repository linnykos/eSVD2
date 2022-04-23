#' Optimize X given C, Y and Z
#'
#' @param XC_init Initial value for the `[X C]` matrix, X [n x k], C [n x r]
#' @param YZ      The `[Y Z]` matrix, Y [p x k], Z [p x r]
#' @param k       Number of columns in X and Y
#' @param loader  The data loader, typically returned by data_loader()
#' @param family  A family object, typically returned by esvd_family()
#' @param s       The library size vector, [n x 1]
#' @param gamma   The nuisance parameter vector, [p x 1]
#' @param l2penx  The l2 penalty parameter for each row of X, [n x 1]
#' @param verbose Verbosity parameter
#' @param inplace Whether the input XC_init will be modified and returned
#'
opt_x <- function(XC_init, YZ, k, loader, family, s, gamma, l2penx,
                  verbose = 0, inplace = FALSE, ...)
{
  storage.mode(YZ) <- "double"
  .opt_x(XC_init, YZ, k, loader, family, s, gamma, l2penx, verbose, inplace)
}

#' Optimize Y and Z given X and C
#'
#' @param YZ_init    Initial value for the `[Y Z]` matrix, Y [p x k], Z [p x r]
#' @param XC         The `[X C]` matrix, X [n x k], C [n x r]
#' @param k          Number of columns in X and Y
#' @param fixed_cols Which columns in YZ need to be fixed
#' @param loader     The data loader, typically returned by data_loader()
#' @param family     A family object, typically returned by esvd_family()
#' @param s          The library size vector, [n x 1]
#' @param gamma      The nuisance parameter vector, [p x 1]
#' @param l2peny     The l2 penalty parameter for each row of Y, [p x 1]
#' @param l2penz     The l2 penalty parameter for each row of Z, [p x 1]
#' @param verbose    Verbosity parameter
#' @param inplace    Whether the input XC_init will be modified and returned
#'
opt_yz <- function(YZ_init, XC, k, fixed_cols, loader, family, s, gamma,
                   l2peny, l2penz, verbose = 0, inplace = FALSE, ...)
{
  # YZind will be passed to C++, should be zero-based
  YZind <- setdiff(1:ncol(YZ_init), fixed_cols) - 1
  storage.mode(XC) <- "double"
  .opt_yz(YZ_init, XC, k, YZind, loader, family, s, gamma, l2peny, l2penz,
          verbose, inplace)
}

#############################

# .opt_nuisance <- function(covariates,
#                           dat,
#                           gene_group_factor,
#                           max_cell_subsample,
#                           offset_vec,
#                           x_mat,
#                           yb_mat,
#                           value_lower,
#                           value_upper,
#                           verbose) {
#   stopifnot(is.factor(gene_group_factor), length(gene_group_factor) == ncol(dat),
#             nrow(x_mat) == nrow(dat), nrow(yb_mat) == ncol(dat),
#             max_cell_subsample > 0,
#             length(value_lower) == length(levels(gene_group_factor)),
#             length(value_upper) == length(levels(gene_group_factor)))
#   if(all(!is.null(covariates))){
#     stopifnot(nrow(covariates) == nrow(dat),
#     ncol(x_mat) + ncol(covariates) == ncol(yb_mat))
#   }
#   p <- ncol(dat); n <- nrow(dat)
#   nuisance_param_vec <- rep(NA, p)
#
#   theta_mat <- tcrossprod(cbind(x_mat, covariates), yb_mat)
#   if(all(!is.null(offset_vec))) {
#     stopifnot(length(offset_vec) == nrow(x_mat))
#     theta_mat <- sweep(theta_mat, 1, offset_vec, "+")
#   }
#
#   gene_groups <- levels(gene_group_factor)
#   for(i in 1:length(gene_groups)){
#     if(verbose == 1 & length(gene_groups) > 10 & i %% floor(length(gene_groups)/10) == 0){
#       cat('*')
#     } else if(verbose >= 2) {
#       print(paste0("Updating nuisance parameter for group ", gene_groups[i]))
#     }
#     gene_idx <- which(gene_group_factor == gene_groups[i])
#
#     tmp_dat <- dat[,gene_idx]
#     tmp_theta <- theta_mat[,gene_idx]
#     if(length(tmp_dat) > max_cell_subsample){
#       cell_idx <- sample(1:length(tmp_dat), size = max_cell_subsample, replace = F)
#       y_vec <- as.numeric(tmp_dat[cell_idx])
#       mean_vec <- as.numeric(exp(tmp_theta[cell_idx]))
#     } else {
#       y_vec <- as.numeric(tmp_dat)
#       mean_vec <- as.numeric(exp(tmp_theta))
#     }
#
#     # seems more robust than glmGamPoi::overdispersion_mle
#     val <- MASS::theta.mm(y = y_vec,
#                           mu = mean_vec,
#                           dfr = length(y_vec)-1)
#
#     # apply the constraints
#     if(!is.na(value_lower[i])){
#       val <- max(val, value_lower[i])
#     }
#     if(!is.na(value_upper[i])){
#       val <- min(val, value_upper[i])
#     }
#
#     nuisance_param_vec[gene_idx] <- val
#   }
#
#   nuisance_param_vec
# }

##########################

.opt_esvd_format_param <- function(family,
                                   l2pen,
                                   max_cell_subsample,
                                   max_iter,
                                   method,
                                   nuisance_value_lower,
                                   nuisance_value_upper,
                                   reparameterize,
                                   reestimate_nuisance,
                                   reestimate_nuisance_per_iteration,
                                   tol,
                                   verbose) {
  list(family = family,
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
}

# Initialize the B matrix according to covariates
.opt_esvd_setup_b_mat <- function(b_init, covariates, p) {
  if(is.null(covariates))
  {
    b_mat <- NULL
  } else {
    if(is.null(b_init))
    {
      b_mat <- matrix(0, nrow = p, ncol = ncol(covariates))
    } else {
      b_mat <- b_init
    }
  }
}

# Set row and column names of output matrices
.opt_esvd_format_matrices <- function(b_mat, covariates, dat, x_mat, y_mat) {
  rownames(x_mat) <- rownames(dat)
  rownames(y_mat) <- colnames(dat)
  colnames(x_mat) <- paste0("latent_", 1:ncol(x_mat))
  colnames(y_mat) <- paste0("latent_", 1:ncol(y_mat))
  rownames(b_mat) <- colnames(dat)
  colnames(b_mat) <- colnames(covariates)

  list(b_mat = b_mat, x_mat = x_mat, y_mat = y_mat)
}
