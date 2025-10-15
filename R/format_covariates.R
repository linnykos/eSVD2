#' Format covariates
#'
#' Mainly, this method splits the categorical variables (which should be `factor` variables)
#' into indicator variables (i.e.,
#' one-hot encoding), dropping the last level, and then rescales
#' all the numerical variables (but does not center them),
#' and computes the \code{"Log_UMI"} (i.e., log total counts) for each cell.
#' \code{"Log_UMI"} is added as its own column.
#'
#' @param dat                         Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
#'                                    and \eqn{p} columns represent genes.
#'                                    The rows and columns of the matrix should be named.
#' @param covariate_df                \code{data.frame} where each row represents a cell, and the
#'                                    columns are the different categorical or numerical variables that you wish to adjust for
#' @param bool_center                 Boolean if the numerical variables should be centered around zero, default is \code{FALSE}
#' @param rescale_numeric_variables   A vector of strings denoting the column names in \code{covariate_df} that are numerical and you wish to rescale
#' @param variables_enumerate_all     If not \code{NULL}, this allows you to control specifically which \code{factor} variables
#'                                    in \code{covariate_df} you would like to split into indicators. By default, this is \code{NULL}, meaning all the \code{factor} variables are split into indicators
#'
#' @return a \code{matrix} with the same number of rows as \code{dat}
#' @export
format_covariates <- function(dat,
                              covariate_df,
                              bool_center = FALSE,
                              rescale_numeric_variables = NULL,
                              variables_enumerate_all = NULL){
  stopifnot(nrow(dat) == nrow(covariate_df), is.data.frame(covariate_df),
            all(is.null(rescale_numeric_variables)) || all(rescale_numeric_variables %in% colnames(covariate_df)),
            all(is.null(variables_enumerate_all)) || all(variables_enumerate_all %in% colnames(covariate_df)))
  n <- nrow(covariate_df)

  factor_vec <- colnames(covariate_df)[sapply(covariate_df, is.factor)]

  numeric_vec <- setdiff(colnames(covariate_df), factor_vec)
  if(length(numeric_vec) > 0){
    covariate_df2 <- covariate_df[,numeric_vec,drop = F]
    colnames(covariate_df2) <- numeric_vec

    if(!all(is.null(rescale_numeric_variables))){
      stopifnot(all(rescale_numeric_variables %in% numeric_vec))

      for(var in rescale_numeric_variables){
        covariate_df2[,var] <- scale(covariate_df[,var], center = bool_center, scale = T)
      }
    }
  } else {
    covariate_df2 <- matrix(0, nrow = n, ncol = 0)
  }
  rownames(covariate_df2) <- rownames(covariate_df)

  logumi_vec <- log(Matrix::rowSums(dat))
  covariate_df2 <- cbind(logumi_vec, covariate_df2)
  colnames(covariate_df2)[1] <- "Log_UMI"

  for(var in factor_vec){
    vec <- covariate_df[,var]
    vec <- droplevels(vec)
    if(!all(is.null(variables_enumerate_all)) && var %in% variables_enumerate_all){
      uniq_level <- levels(vec)
    } else {
      stopifnot(length(levels(vec)) > 1)
      uniq_level <- levels(vec)[-1]
    }

    for(lvl in uniq_level){
      tmp <- rep(0, n)
      tmp[which(vec == lvl)] <- 1

      var_name <- paste0(var, "_", lvl)
      covariate_df2 <- cbind(covariate_df2, tmp)
      colnames(covariate_df2)[ncol(covariate_df2)] <- var_name
    }
  }

  covariate_df2 <- cbind(1, covariate_df2)
  colnames(covariate_df2)[1] <- "Intercept"

  as.matrix(covariate_df2)
}
