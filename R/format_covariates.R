format_covariates <- function(dat, covariate_df,
                              mixed_effect_variables = c()){
  stopifnot(nrow(dat) == nrow(covariate_df),
            is.data.frame(covariate_df),
            all(mixed_effect_variables %in% colnames(covariate_df)),
            all(sapply(mixed_effect_variables, function(var){
              is.factor(covariate_df[,var])
            })))

  n <- nrow(covariate_df)

  factor_vec <- colnames(covariate_df)[sapply(covariate_df, is.factor)]
  numeric_vec <- setdiff(colnames(covariate_df), factor_vec)
  covariate_df2 <- sapply(numeric_vec, function(var){
    scale(covariate_df[,var], center = T, scale = T)
  })
  logumi_vec <- log(Matrix::rowSums(dat))
  covariate_df2 <- cbind(logumi_vec, covariate_df2[,numeric_vec])
  colnames(covariate_df2)[1] <- "Log_UMI"

  for(var in factor_vec){
    vec <- covariate_df[,var]
    if(var %in% mixed_effect_variables){
      uniq_level <- levels(vec)
    } else {
      uniq_level <- levels(vec)[2]
    }

    for(lvl in uniq_level){
      tmp <- rep(0, n)
      tmp[which(vec == lvl)] <- 1

      var_name <- paste0(var, "_", lvl)
      covariate_df2 <- cbind(covariate_df2, tmp)
      colnames(covariate_df2)[ncol(covariate_df2)] <- var_name
    }
  }

  as.matrix(covariate_df2)
}
