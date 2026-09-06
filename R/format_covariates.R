#' Format covariates
#'
#' Mainly, this method splits the categorical variables (which should be `factor` variables)
#' into indicator variables (i.e.,
#' one-hot encoding), dropping the first level of each factor as the
#' reference, then rescales the numerical variables named in
#' \code{rescale_numeric_variables} (only those; the others are left as they
#' are), and computes the \code{"Log_UMI"} (i.e., log total counts) for each cell.
#' \code{"Log_UMI"} is added as its own column.
#'
#' The rescaling divides each named column by its root-mean-square
#' (\code{scale(x, center = FALSE, scale = TRUE)}), not by its standard
#' deviation, so the result does not have unit variance unless the column is
#' also centered (\code{bool_center = TRUE}). What it does deliver is
#' invariance to the unit the covariate was recorded in: age in years and age
#' in months give identical columns.
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
#' @return a \code{matrix} with the same number of rows as \code{dat}, whose
#' first two columns are \code{"Intercept"} and \code{"Log_UMI"}, followed by
#' the numerical variables and then one indicator column
#' \code{"<variable>_<level>"} per retained factor level.
#' @examples
#' set.seed(10)
#' dat <- matrix(stats::rpois(60 * 5, lambda = 3), nrow = 60, ncol = 5)
#' dimnames(dat) <- list(paste0("cell", 1:60), paste0("gene", 1:5))
#' covariate_df <- data.frame(
#'   CC = factor(rep(c("control", "case"), each = 30),
#'               levels = c("control", "case")),
#'   Sex = factor(rep(c("F", "M"), times = 30)),
#'   Age = stats::rnorm(60, mean = 40, sd = 10)
#' )
#' covariates <- format_covariates(dat = dat,
#'                                 covariate_df = covariate_df,
#'                                 rescale_numeric_variables = "Age")
#' colnames(covariates)
#' utils::head(covariates, 3)
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
    covariate_df2 <- covariate_df[,numeric_vec,drop = FALSE]
    colnames(covariate_df2) <- numeric_vec

    if(!all(is.null(rescale_numeric_variables))){
      stopifnot(all(rescale_numeric_variables %in% numeric_vec))

      for(var in rescale_numeric_variables){
        covariate_df2[,var] <- scale(covariate_df[,var], center = bool_center, scale = TRUE)
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
      if(length(levels(vec)) < 2){
        stop("factor variable `", var, "` has only one level (\"",
             paste0(levels(vec), collapse = "\", \""),
             "\") after dropping unused levels, so it cannot be split into ",
             "indicators; remove it from `covariate_df`")
      }
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
