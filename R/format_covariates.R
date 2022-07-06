format_covariates <- function(dat,
                              covariate_df,
                              bool_center = F,
                              rescale_numeric_variables = NULL){
  stopifnot(nrow(dat) == nrow(covariate_df), is.data.frame(covariate_df),
            all(is.null(rescale_numeric_variables)) || all(rescale_numeric_variables %in% colnames(covariate_df)))
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
    stopifnot(length(levels(vec)) > 1)
    uniq_level <- levels(vec)[-1]

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
