#' @export
estimate_nuisance <- function(input_obj, ...) UseMethod("estimate_nuisance")

#' @export
estimate_nuisance.eSVD <- function(input_obj,
                                   verbose = 0, ...){
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"))

  dat <- .get_object(eSVD_obj = input_obj, what_obj = "dat", which_fit = NULL)
  covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
  latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
  case_control_variable <- .get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
  x_mat <-.get_object(eSVD_obj = input_obj, what_obj = "x_mat", which_fit = latest_Fit)
  y_mat <-.get_object(eSVD_obj = input_obj, what_obj = "y_mat", which_fit = latest_Fit)
  z_mat <-.get_object(eSVD_obj = input_obj, what_obj = "z_mat", which_fit = latest_Fit)

  nat_mat1 <- tcrossprod(x_mat, y_mat)
  nat_mat2 <- tcrossprod(covariates[,case_control_variable,drop = F],
                         z_mat[,case_control_variable,drop = F])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)

  offset_variable <- setdiff(colnames(covariates), case_control_variable)
  library_mat <- exp(tcrossprod(
    covariates[,offset_variable], z_mat[,offset_variable]
  ))

  nuisance_vec <- estimate_nuisance.default(
    input_obj = dat,
    mean_mat = mean_mat_nolib,
    library_mat = library_mat,
    verbose = verbose
  )

  input_obj[[latest_Fit]]$nuisance_vec <- nuisance_vec
  input_obj
}

#' @export
estimate_nuisance.default <- function(input_obj,
                                      mean_mat,
                                      library_mat,
                                      verbose = 0, ...){
  stopifnot(inherits(input_obj, c("matrix", "dgCMatrix")),
            is.matrix(mean_mat), is.matrix(library_mat),
            all(dim(mean_mat) == dim(input_obj)),
            all(dim(library_mat) == dim(input_obj)))

  p <- ncol(input_obj)
  nuisance_vec <- sapply(1:p, function(j){
    if(verbose ==1 && p > 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(j)

    val <- tryCatch(
      exp(log_gamma_rate(x = as.numeric(input_obj[,j]),
                     mu = mean_mat[,j],
                     s = library_mat[,j])),
      error = function(c) {
        if(verbose > 0) print(paste0("Nuisance estimation failed at variable ", j))
        0
      })

    val
  })
  if(length(colnames(input_obj)) > 0) names(nuisance_vec) <- colnames(input_obj)

  nuisance_vec
}

