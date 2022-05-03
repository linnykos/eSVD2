.get_object <- function(eSVD_obj,
                        what_obj,
                        which_fit){
  if(is.null(what_obj)){
    stopifnot(inherits(eSVD_obj[[which_fit]], "eSVD_Fit"))
    eSVD_obj[[which_fit]]

  } else {
    if(what_obj == "dat"){
      return(eSVD_obj$dat)

    } else if(what_obj == "covariates"){
      return(eSVD_obj$covariates)

    } else if(what_obj == "z_mat1"){
      stopifnot(which_fit == "initial_Reg")
      return(eSVD_obj$initial_Reg$z_mat1)

    } else if(what_obj == "z_mat2"){
      stopifnot(which_fit == "initial_Reg")
      return(eSVD_obj$initial_Reg$z_mat2)

    } else if(what_obj == "log_pval"){
      stopifnot(which_fit == "initial_Reg")
      return(eSVD_obj$initial_Reg$log_pval)

    } else if(what_obj == "x_mat"){
      stopifnot(which_fit %in% names(eSVD_obj))
      return(eSVD_obj[[which_fit]]$x_mat)

    } else if(what_obj == "y_mat"){
      stopifnot(which_fit %in% names(eSVD_obj))
      return(eSVD_obj[[which_fit]]$y_mat)

    } else if(what_obj == "z_mat"){
      stopifnot(which_fit %in% names(eSVD_obj))
      return(eSVD_obj[[which_fit]]$z_mat)

    } else if(what_obj == "posterior_mean_mat"){
      return(eSVD_obj[[which_fit]]$posterior_mean_mat)

    } else if(what_obj == "posterior_var_mat"){
      return(eSVD_obj[[which_fit]]$posterior_var_mat)

    } else if(what_obj == "nuisance"){
      stopifnot(which_fit %in% names(eSVD_obj))
      return(eSVD_obj[[which_fit]]$nuisance_vec)

    } else if(what_obj == "latest_Fit"){
      return(eSVD_obj[["latest_Fit"]])

    } else if(which_fit == "param"){
      param <- eSVD_obj$param
      stopifnot(what_obj %in% names(param))
      return(param[[what_obj]])

    } else {
      stopifnot("what_obj is not found")
    }
  }
}

.form_esvd_fit <- function(x_mat,
                           y_mat,
                           z_mat, ...){
  structure(list(x_mat = x_mat,
                 y_mat = y_mat,
                 z_mat = z_mat, ...),
            class = "eSVD_Fit")
}
