
# Optimize X given Y and B
opt_x <- function(X0, Y, B, Z, A,
                  family, s, gamma, offset_vec, l2pen,
                  opt_fun,
                  bool_run_cpp,
                  gene_group_factor,
                  gene_ignore_excessive_zero,
                  verbose = 0, ...)
{
  if(bool_run_cpp && !is.null(family$cpp_functions) && identical(opt_fun, constr_newton))
  {
    return(opt_x_cpp(X0, Y, B, Z, A, family, s, gamma, offset_vec, l2pen, verbose))
  }

  n <- nrow(A)
  X <- X0
  # Optimize each row of X
  for(i in 1:n)
  {
    if(verbose >= 2)
      cat("===== Optimizing Row ", i, " of X =====\n", sep = "")
    Zi <- if(is.null(Z)) NULL else Z[i, ]
    if(gene_ignore_excessive_zero){
      gene_idx <- sort(unique(c(which(A[i,] > 0), grep("normal", gene_group_factor))))
    } else {
      gene_idx <- 1:ncol(A)
    }

    opt <- opt_fun(
      x0 = X0[i, ],
      f = objfn_Xi,
      gr = grad_Xi,
      hn = hessian_Xi,
      direc = direction_Xi,
      feas = feas_Xi,
      eps_rel = 1e-3,
      verbose = (verbose >= 3),
      Y = Y[gene_idx,,drop=F], B = B[gene_idx,,drop=F], Zi = Zi,
      Ai = A[i, gene_idx],
      family = family, si = s[i], gamma = gamma[gene_idx],
      offseti = offset_vec[i], l2pen = l2pen, ...
    )

    X[i, ] <- opt$x
    if(verbose >= 4) {
      print(paste0("Iteration ", i))
      print(X[i, ])
    }
    if(verbose >= 3) cat("==========\n\n")
  }
  X
}

# Optimize Y and B given X
opt_yb <- function(YB0, XZ, A,
                   family, s, gamma, offset_vec, l2pen,
                   opt_fun,
                   bool_run_cpp,
                   verbose = 0, ...)
{
  if(bool_run_cpp && !is.null(family$cpp_functions) && identical(opt_fun, constr_newton))
  {
    return(opt_yb_cpp(YB0, XZ, A, family, s, gamma, offset_vec, l2pen, verbose))
  }

  p <- ncol(A)
  YB <- YB0
  # Optimize each row of Y and B
  for(j in 1:p)
  {
    if(verbose >= 2)
      cat("===== Optimizing Row ", j, " of Y =====\n", sep = "")

    opt <- opt_fun(
      x0 = YB0[j, ],
      f = objfn_Yj,
      gr = grad_Yj,
      hn = hessian_Yj,
      direc = direction_Yj,
      feas = feas_Yj,
      eps_rel = 1e-3,
      verbose = (verbose >= 3),
      X = XZ, Bj = NULL, Z = NULL, Aj = A[, j],
      family = family, s = s, gammaj = gamma[j], offset = offset_vec, l2pen = l2pen, ...
    )

    YB[j, ] <- opt$x
    if(verbose >= 4) {
      print(paste0("Iteration ", i))
      print(X[i, ])
    }
    if(verbose >= 3)  cat("==========\n\n")
  }
  YB
}

#############################

.opt_nuisance <- function(covariates,
                          dat,
                          gene_group_factor,
                          max_cell_subsample,
                          offset_vec,
                          x_mat,
                          yb_mat,
                          value_lower,
                          value_upper,
                          verbose){
  print(gene_group_factor)
  print(class(gene_group_factor))
  stopifnot(is.factor(gene_group_factor), length(gene_group_factor) == ncol(dat),
            nrow(x_mat) == nrow(dat), nrow(yb_mat) == ncol(dat),
            max_cell_subsample > 0,
            length(value_lower) == length(levels(gene_group_factor)),
            length(value_upper) == length(levels(gene_group_factor)))
  if(all(!is.null(covariates))){
    stopifnot(nrow(covariates) == nrow(dat),
    ncol(x_mat) + ncol(covariates) == ncol(yb_mat))
  }
  p <- ncol(dat); n <- nrow(dat)
  nuisance_param_vec <- rep(NA, p)

  theta_mat <- tcrossprod(cbind(x_mat, covariates), yb_mat)
  if(all(!is.null(offset_vec))) {
    stopifnot(length(offset_vec) == nrow(x_mat))
    theta_mat <- sweep(theta_mat, 1, offset_vec, "+")
  }

  gene_groups <- levels(gene_group_factor)
  print(gene_groups)
  for(i in 1:length(gene_groups)){
    if(verbose >= 1) print(paste0("Updating nuisance parameter for group ", gene_groups[i]))
    gene_idx <- which(gene_group_factor == gene_groups[i])

    tmp_dat <- dat[,gene_idx]
    tmp_theta <- theta_mat[,gene_idx]
    if(length(tmp_dat) > max_cell_subsample){
      cell_idx <- sample(1:length(tmp_dat), size = max_cell_subsample, replace = F)
      y_vec <- as.numeric(tmp_dat[cell_idx])
      mean_vec <- as.numeric(exp(tmp_theta[cell_idx]))
    } else {
      y_vec <- as.numeric(tmp_dat)
      mean_vec <- as.numeric(exp(tmp_theta))
    }

    # seems more robust than glmGamPoi::overdispersion_mle
    val <- MASS::theta.mm(y = y_vec,
                          mu = mean_vec,
                          dfr = length(y_vec)-1)

    # apply the constraints
    if(!is.na(value_lower[i])){
      val <- max(val, value_lower[i])
    }
    if(!is.na(value_upper[i])){
      val <- min(val, value_upper[i])
    }

    nuisance_param_vec[gene_idx] <- val
  }

  nuisance_param_vec
}

##########################

.opt_esvd_format_param <- function(bool_run_cpp,
                                   family,
                                   gene_group_factor,
                                   gene_ignore_excessive_zero,
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
                                   verbose){
  list(bool_run_cpp = bool_run_cpp,
       family = family,
       gene_group_factor = gene_group_factor,
       gene_ignore_excessive_zero = gene_ignore_excessive_zero,
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

.opt_esvd_setup_b_mat <- function(b_init, covariates){
  if(is.null(covariates))
  {
    b_mat <- NULL
  } else {
    if(is.null(b_init))
    {
      r <- ncol(covariates)
      b_mat <- matrix(0, p, r)
    } else {
      b_mat <- b_init
    }
  }
}

.opt_esvd_format_matrices <- function(b_mat, covariates,
                                      dat, x_mat, y_mat){
  rownames(x_mat) <- rownames(dat)
  rownames(y_mat) <- colnames(dat)
  colnames(x_mat) <- paste0("latent_", 1:ncol(x_mat))
  colnames(y_mat) <- paste0("latent_", 1:ncol(y_mat))
  rownames(b_mat) <- colnames(dat)
  colnames(b_mat) <- colnames(covariates)

  list(b_mat = b_mat, x_mat = x_mat, y_mat = y_mat)
}
