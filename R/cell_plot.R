cell_plot <- function(input_obj,
                      variable,
                      what_1,
                      what_2 = NULL,
                      bool_include_diagonal = T,
                      bool_jitter_x = F,
                      bool_jitter_y = F,
                      col_case = 2,
                      col_control = 3,
                      indiv_cases = NULL,
                      indiv_controls = NULL,
                      indiv_vec = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      ...){
  what_vec <- c("dat", "relative_expression",
                "fit", "fit_nolibrary",
                "posterior_mean_nonadjusted", "posterior_mean",
                "posterior_variance")
  stopifnot(class(input_obj) == "eSVD",
            what_1 %in% what_vec,
            is.null(what_2) || what_2 %in% what_vec,
            variable %in% colnames(input_obj$dat))

  # first compute the relevant vectors
  vec_1 <- .compute_what_cell_vector(input_obj = input_obj,
                                     variable = variable,
                                     what = what_1)
  vec_2 <- .compute_what_cell_vector(input_obj = input_obj,
                                     variable = variable,
                                     what = what_2)
  if(bool_jitter_x) vec_1 <- jitter(vec_1)
  if(all(!is.null(vec_2))) {
    if(bool_jitter_y) vec_2 <- jitter(vec_2)
    vec_all <- cbind(vec_1, vec_2)
  } else {
    vec_all <- matrix(vec_1, nrow = length(vec_1), ncol = 1)
  }

  if(all(!is.null(indiv_cases)) & all(!is.null(indiv_controls)) & all(!is.null(indiv_vec))){
    stopifnot(is.factor(indiv_vec), length(indiv_vec) == nrow(input_obj$dat),
              length(intersect(indiv_cases, indiv_controls)) == 0,
              all(indiv_cases %in% levels(indiv_vec)),
              all(indiv_controls %in% levels(indiv_vec)))

    case_indiv_idx <- lapply(indiv_cases, function(indiv){
      which(indiv_vec == indiv)
    })
    control_indiv_idx <- lapply(indiv_controls, function(indiv){
      which(indiv_vec == indiv)
    })

    case_means <- sapply(1:length(case_indiv_idx), function(i){
      Matrix::colMeans(vec_all[case_indiv_idx[[i]],,drop = F])
    })
    control_means <- sapply(1:length(control_indiv_idx), function(i){
      Matrix::colMeans(vec_all[control_indiv_idx[[i]],,drop = F])
    })

    # do some formatting
    if(is.matrix(case_means) && nrow(case_means) == 2){
      case_means <- t(case_means); control_means <- t(control_means)
    } else {
      case_means <- matrix(case_means, nrow = length(case_means), ncol = 1)
      control_means <- matrix(control_means, nrow = length(control_means), ncol = 1)
    }

    means_all <- cbind(
      rbind(case_means, control_means),
      c(rep(1, length(case_indiv_idx)), rep(2, length(control_indiv_idx)))
    )

    shuff_idx <- sample(1:nrow(means_all))
    means_all <- means_all[shuff_idx,,drop = F]
    color_vec <- sapply(means_all[,ncol(means_all)], function(x){
      ifelse(x == 1, col_case, col_control)
    })
  } else {
    means_all <- NULL
  }

  # now plot
  if(all(!is.null(vec_2))){
    if(is.null(xlab)) xlab <- what_1
    if(is.null(ylab)) ylab <- what_2

    graphics::plot(x = vec_1, y = vec_2,
                   cex = 1,
                   pch = 16,
                   col = rgb(0.5, 0.5, 0.5, 0.2),
                   xlab = xlab, ylab = ylab,
                   ...)

    if(bool_include_diagonal){
      graphics::lines(x = c(min(c(vec_1, vec_2)), max(c(vec_1, vec_2))),
                      y = c(min(c(vec_1, vec_2)), max(c(vec_1, vec_2))),
                      lwd = 2, lty = 2,
                      col = 2)
    }

    if(all(!is.null(means_all))){
      graphics::points(x = means_all[,1], y = means_all[,2],
                       col = "white", pch = 16, cex = 3)
      graphics::points(x = means_all[,1], y = means_all[,2],
                       col = color_vec, pch = 16, cex = 2.5)

      for(idx in shuff_idx){
        graphics::rug(means_all[idx,1], col = color_vec[idx], side = 3,
                      lwd = 2)
        graphics::rug(means_all[idx,2], col = color_vec[idx], side = 4,
                      lwd = 2)
      }
    }

  } else {
    if(is.null(xlab)) xlab <- what_1
    graphics::hist(vec_1,
                   xlab = xlab,
                   ...)

    if(all(!is.null(means_all))){
      for(idx in shuff_idx){
        graphics::rug(means_all[idx,1], col = color_vec[idx], side = 1,
                      lwd = 2)
      }
    }
  }

  invisible()
}

.compute_what_cell_vector <- function(input_obj,
                                      variable,
                                      what){
  if(is.null(what)) return(NULL)
  latest_Fit <- .get_object(eSVD_obj = input_obj,
                            what_obj = "latest_Fit",
                            which_fit = NULL)

  if(what == "dat") {
    return(as.numeric(input_obj$dat[,variable]))

  } else if (what == "relative_expression"){
    library_size <- Matrix::rowSums(input_obj$dat)
    return(as.numeric(input_obj$dat[,variable])/library_size)

  } else if (what == "posterior_mean"){
    return(input_obj[[latest_Fit]]$posterior_mean_mat[,variable])

  } else if (what == "posterior_variance"){
    return(input_obj[[latest_Fit]]$posterior_var_mat[,variable])

  } else {
    covariates <- .get_object(eSVD_obj = input_obj,
                              what_obj = "covariates",
                              which_fit = NULL)
    x_mat <- .get_object(eSVD_obj = input_obj,
                         what_obj = "x_mat",
                         which_fit = latest_Fit)
    y_mat <- .get_object(eSVD_obj = input_obj,
                         what_obj = "y_mat",
                         which_fit = latest_Fit)
    z_mat <- .get_object(eSVD_obj = input_obj,
                         what_obj = "z_mat",
                         which_fit = latest_Fit)
    nat_mat1 <- tcrossprod(x_mat, y_mat)

    if(what == "fit"){
      nat_mat2 <- tcrossprod(covariates, z_mat)
      return(exp(nat_mat1[,variable] + nat_mat2[,variable]))

    } else if (what == "fit_nolibrary"){
      library_size_variable <- input_obj$param$init_library_size_variable
      library_idx <- which(colnames(covariates) == library_size_variable)

      nat_mat2 <- tcrossprod(covariates[,-library_idx,drop = F],
                             z_mat[,-library_idx,drop = F])
      nat_mat_nolib <- nat_mat1 + nat_mat2
      mean_mat_nolib <- exp(nat_mat_nolib)
      return(mean_mat_nolib[,variable])

    } else if(what == "posterior_mean_nonadjusted"){
      library_size_variable <- input_obj$param$init_library_size_variable
      library_idx <- which(colnames(covariates) == library_size_variable)
      nat_mat2 <- tcrossprod(covariates[,-library_idx,drop = F],
                             z_mat[,-library_idx,drop = F])
      nat_mat_nolib <- nat_mat1 + nat_mat2
      mean_mat_nolib <- exp(nat_mat_nolib)
      library_mat <- exp(tcrossprod(
        covariates[,library_idx], z_mat[,library_idx]
      ))
      colnames(library_mat) <- rownames(z_mat)

      nuisance_vec <-.get_object(eSVD_obj = input_obj,
                                 what_obj = "nuisance",
                                 which_fit = latest_Fit)

      alpha <- as.numeric(input_obj$dat[,variable]) + mean_mat_nolib[,variable]*nuisance_vec[variable]
      beta <- library_mat[,variable]+nuisance_vec[variable]
      return(alpha/beta)
    } else {
      stop("what not found")
    }
  }

}
