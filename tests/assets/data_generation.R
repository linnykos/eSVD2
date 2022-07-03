rm(list=ls())

set.seed(123)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

num_indiv <- 20
n_per_indiv <- 100
p <- 150
k <- 5
n <- num_indiv * n_per_indiv
x_mat <- matrix(abs(rnorm(n * k))*.5, nrow = n, ncol = k)
y_mat <- matrix(abs(rnorm(p * k))*.5, nrow = p, ncol = k)
covariate_df <- cbind(
  c(rep(0, n/2), rep(1, n/2)),
  c(rep(0, n/4), rep(1, n/4), rep(0, n/4), rep(1, n/4)),
  sapply(1:n, function(i){rnorm(1, mean = i/n*3, sd = 0.5)})
)
colnames(covariate_df) <- c("case_control", "gender", "Log_UMI")
covariate_df <- as.data.frame(covariate_df)
indiv_vec <- rep(paste0("individual_", 1:num_indiv), each = n_per_indiv)
covariate_df <- cbind(covariate_df, indiv_vec)
colnames(covariate_df)[ncol(covariate_df)] <- "individual"
covariate_df[,"case_control"] <- as.factor(covariate_df[,"case_control"])
covariate_df[,"gender"] <- as.factor(covariate_df[,"gender"])
covariate_df[,"individual"] <- as.factor(covariate_df[,"individual"])
covariates <- format_covariates(dat = abs(matrix(rnorm(n*p), nrow = n, ncol = p)),
                               covariate_df = covariate_df[,which(colnames(covariate_df) != "Log_UMI")],
                               subject_variable = "individual")
covariates[,"Log_UMI"] <- covariate_df[,"Log_UMI"]

z_mat <- cbind(0.5,
               rep(1,p),
               c(rep(0, .9*p), rep(1, 2/3*(.1*p)), rep(2, 1/3*(.1*p))),
               rnorm(p)
)
for(i in 1:num_indiv){
  tmp <- rnorm(p, mean = 0, sd = .2)
  z_mat <- cbind(z_mat, tmp)
}
colnames(z_mat) <-  colnames(covariates)
true_cc_status <- ifelse(z_mat[,"case_control_1"] > 1e-6, 2, 1)
case_control_variable <- "case_control_1"
case_control_idx <- which(colnames(z_mat) == case_control_variable)
library_idx <- which(colnames(z_mat) == "Log_UMI")

nat_mat_nolib <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates[,case_control_idx], z_mat[,case_control_idx])
library_mat <- exp(tcrossprod(covariates[,-library_idx], z_mat[,-library_idx]))
nuisance_vec <- rep(c(5, 1, 1/5), times = p/3)

# Simulate data
gamma_mat <- matrix(NA, nrow = n, ncol = p)
dat <- matrix(NA, nrow = n, ncol = p)
for(i in 1:n){
  for(j in 1:p){
    gamma_mat[i,j] <- stats::rgamma(n = 1,
                                    shape = nuisance_vec[j]*exp(nat_mat_nolib[i,j]),
                                    rate = nuisance_vec[j])
    dat[i,j] <- stats::rpois(n = 1, lambda = library_mat[i,j] * gamma_mat[i,j])
  }
}
dat <- pmin(dat, 100)
# quantile(dat)
# length(which(dat == 0))/prod(dim(dat))
# image(dat)
dat <- Matrix::Matrix(dat, sparse = T)
rownames(dat) <- paste0("c", 1:n)
colnames(dat) <- paste0("g", 1:p)
metadata <- data.frame(individual = factor(rep(1:num_indiv, each = n_per_indiv)))
rownames(metadata) <- rownames(dat)
covariates[,"Log_UMI"] <- log(Matrix::rowSums(dat))

# fit eSVD
eSVD_obj <- eSVD2::initialize_esvd(dat = dat,
                                   covariates = covariates,
                                   case_control_variable = case_control_variable,
                                   k = 5,
                                   subject_variables = colnames(covariates)[grep("individual", colnames(covariates))],
                                   verbose = 1)
eSVD_obj <- eSVD2::opt_esvd(input_obj = eSVD_obj,
                            max_iter = 50,
                            verbose = 1)
# pred_mat <- exp(tcrossprod(eSVD_obj$fit_First$x_mat, eSVD_obj$fit_First$y_mat) + tcrossprod(eSVD_obj$covariates, eSVD_obj$fit_First$z_mat)); image(pred_mat)
# plot(eSVD_obj$fit_First$z_mat[,"case_control"], col = true_cc_status)
eSVD_obj <- eSVD2::estimate_nuisance(input_obj = eSVD_obj,
                                     bool_library_includes_interept = T,
                                     verbose = 1)
# plot(eSVD_obj$fit_First$nuisance_vec, jitter(nuisance_vec), asp = T)
eSVD_obj <- eSVD2::compute_posterior(input_obj = eSVD_obj,
                                     bool_adjust_covariates = F,
                                     bool_covariates_as_library = T)
eSVD_obj <- eSVD2::compute_test_statistic(input_obj = eSVD_obj,
                                          covariate_individual = "individual",
                                          metadata = metadata)

save(dat, covariates, metadata, nuisance_vec,
     nat_mat_nolib, gamma_mat, library_mat,
     true_cc_status,
     eSVD_obj,
     date_of_run, session_info,
     file = "tests/assets/synthetic_data.RData")
