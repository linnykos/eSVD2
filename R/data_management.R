.get_object <- function(esvd_obj,
                        what_obj,
                        which_fit){
  if(what_obj == "data"){
    esvd_obj$data

  } else if(what_obj == "covariates"){
    esvd_obj$covariates

  } else if(what_obj == "z_mat1"){
    stopifnot(which_fit == "initial_Reg")
    esvd_obj$initial_Reg$z_mat1

  } else if(what_obj == "z_mat2"){
    stopifnot(which_fit == "initial_Reg")
    esvd_obj$initial_Reg$z_mat2

  } else if(what_obj == "log_pval"){
    stopifnot(which_fit == "initial_Reg")
    esvd_obj$initial_Reg$log_pval

  } else if(what_obj == "x_mat"){
    stopifnot(which_fit %in% names(esvd_obj))
    esvd_obj[[which_fit]]$x_mat

  } else if(what_obj == "y_mat"){
    stopifnot(which_fit %in% names(esvd_obj))
    esvd_obj[[which_fit]]$y_mat

  } else if(what_obj == "z_mat"){
    stopifnot(which_fit %in% names(esvd_obj))
    esvd_obj[[which_fit]]$z_mat

  } else if(which_fit == "param"){
    param <- esvd_obj$param
    stopifnot(what_obj %in% names(param))
    param$what_obj

  } else {
    stopifnot("what_obj is not found")
  }
}

.form_esvd_fit <- function(x_mat,
                           y_mat,
                           z_mat, ...){
  structure(list(x_mat = x_mat,
                 y_mat = y_mat,
                 z_mat = z_mat, ...),
            class = "eSVD_fit")
}
