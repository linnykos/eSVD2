format_covariates <- function(dat,
                              covariate_df,
                              bool_center = F,
                              bool_library_size = T,
                              subject_variable = NULL){
  stopifnot(nrow(dat) == nrow(covariate_df), is.data.frame(covariate_df),
            is.null(subject_variable) || subject_variable %in% colnames(covariate_df))
  n <- nrow(covariate_df)

  factor_vec <- colnames(covariate_df)[sapply(covariate_df, is.factor)]
  numeric_vec <- setdiff(colnames(covariate_df), factor_vec)
  if(length(numeric_vec) > 0){
    covariate_df2 <- sapply(numeric_vec, function(var){
      scale(covariate_df[,var], center = bool_center, scale = T)
    })
  } else {
    covariate_df2 <- matrix(0, nrow = n, ncol = 0)
  }
  rownames(covariate_df2) <- rownames(covariate_df)
  colnames(covariate_df2) <- numeric_vec

  if(bool_library_size){
    logumi_vec <- log1p(Matrix::rowSums(dat))
    covariate_df2 <- cbind(logumi_vec, covariate_df2)
    colnames(covariate_df2)[1] <- "Log_UMI"
  }

  for(var in factor_vec){
    vec <- covariate_df[,var]
    if(!is.null(subject_variable) && var == subject_variable){
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

  # move all the subject variables to the back
  if(!is.null(subject_variable)){
    col_idx <- grep(subject_variable, colnames(covariate_df2))
    covariate_df2 <- cbind(covariate_df2[,-col_idx,drop = F],
                           covariate_df2[,col_idx,drop = F])
  }

  as.matrix(covariate_df2)
}
