# Provenance script for tests/assets/synthetic_data.RData, the fixture used by
# the legacy test files (test_initialization.R, test_optimization.R,
# test_nuisance.R, test_posterior.R, test_compute_test_statistic.R,
# test_compute_test_per_gene.R, test_gamma_rate.R, test_report_results.R).
#
# Run from the package root with the package loaded (devtools::load_all()).
# It is NOT run at check time.
#
# The fixture used to be 2000 cells x 150 genes and carried a fully fitted
# 8.7 MB eSVD object, three 2000 x 150 generation by-products and a
# session_info; the tarball was 13 MB against CRAN's 5 MB limit. It is now
# 20 individuals x 20 cells = 400 cells and 60 genes, saved with xz, and
# keeps only what a test reads: `dat`, `covariates`, `metadata`,
# `nuisance_vec`, `library_mat`, `true_cc_status`, `eSVD_obj` and
# `date_of_run`.

rm(list = ls())

set.seed(123)
date_of_run <- Sys.time()

num_indiv <- 20
n_per_indiv <- 20
p <- 60
k <- 5
n <- num_indiv * n_per_indiv
x_mat <- matrix(abs(stats::rnorm(n * k)) * 0.5, nrow = n, ncol = k)
y_mat <- matrix(abs(stats::rnorm(p * k)) * 0.5, nrow = p, ncol = k)
covariate_df <- cbind(
  c(rep(0, n / 2), rep(1, n / 2)),
  c(rep(0, n / 4), rep(1, n / 4), rep(0, n / 4), rep(1, n / 4)),
  sapply(seq_len(n), function(i){stats::rnorm(1, mean = i / n * 3, sd = 0.5)})
)
colnames(covariate_df) <- c("case_control", "gender", "Log_UMI")
covariate_df <- as.data.frame(covariate_df)
indiv_vec <- rep(paste0("individual_", seq_len(num_indiv)), each = n_per_indiv)
covariate_df <- cbind(covariate_df, indiv_vec)
colnames(covariate_df)[ncol(covariate_df)] <- "individual"
covariate_df[, "case_control"] <- as.factor(covariate_df[, "case_control"])
covariate_df[, "gender"] <- as.factor(covariate_df[, "gender"])
covariate_df[, "individual"] <- as.factor(covariate_df[, "individual"])
covariates <- format_covariates(
  dat = abs(matrix(stats::rnorm(n * p), nrow = n, ncol = p)),
  covariate_df = covariate_df[, which(colnames(covariate_df) != "Log_UMI")]
)
covariates[, "Log_UMI"] <- covariate_df[, "Log_UMI"]

# 10% of the genes are truly DE (the last two blocks of the case-control
# coefficient), which needs 0.9 * p, 2/3 * 0.1 * p and 1/3 * 0.1 * p to be
# integers: p = 60 gives 54, 4 and 2.
z_mat <- cbind(0.5,
               rep(1, p),
               c(rep(0, 0.9 * p), rep(1, 2 / 3 * (0.1 * p)), rep(2, 1 / 3 * (0.1 * p))),
               stats::rnorm(p))
for(i in 2:num_indiv){
  tmp <- stats::rnorm(p, mean = 0, sd = 0.2)
  z_mat <- cbind(z_mat, tmp)
}
colnames(z_mat) <- colnames(covariates)
true_cc_status <- ifelse(z_mat[, "case_control_1"] > 1e-6, 2, 1)
case_control_variable <- "case_control_1"
case_control_idx <- which(colnames(z_mat) == case_control_variable)
library_idx <- which(colnames(z_mat) == "Log_UMI")

nat_mat_nolib <- tcrossprod(x_mat, y_mat) +
  tcrossprod(covariates[, case_control_idx], z_mat[, case_control_idx])
library_mat <- exp(tcrossprod(covariates[, -library_idx], z_mat[, -library_idx]))
nuisance_vec <- rep(c(5, 1, 1 / 5), times = p / 3)

# Simulate data from the Gamma-Poisson hierarchy; `nuisance_vec` is the rate.
gamma_mat <- matrix(stats::rgamma(n = n * p,
                                  shape = as.numeric(exp(nat_mat_nolib)) *
                                    rep(nuisance_vec, each = n),
                                  rate = rep(nuisance_vec, each = n)),
                    nrow = n, ncol = p)
dat <- matrix(stats::rpois(n = n * p, lambda = as.numeric(library_mat * gamma_mat)),
              nrow = n, ncol = p)
dat <- pmin(dat, 100)
dat <- Matrix::Matrix(dat, sparse = TRUE)
rownames(dat) <- paste0("c", seq_len(n))
colnames(dat) <- paste0("g", seq_len(p))
metadata <- data.frame(individual = factor(rep(seq_len(num_indiv), each = n_per_indiv)))
rownames(metadata) <- rownames(dat)
covariates[, "Log_UMI"] <- log(Matrix::rowSums(dat))
stopifnot(length(.which_all_zero(dat)) == 0)

########################

id_var <- "individual"

# fit eSVD, through the posterior and the test statistic
eSVD_obj <- initialize_esvd(dat = dat,
                            covariates = covariates[, -grep(id_var, colnames(covariates))],
                            case_control_variable = case_control_variable,
                            bool_intercept = TRUE,
                            k = 5,
                            lambda = 0.1,
                            metadata_case_control = covariates[, case_control_variable],
                            metadata_individual = covariate_df[, id_var],
                            verbose = 0)

eSVD_obj <- opt_esvd(input_obj = eSVD_obj,
                     max_iter = 50,
                     verbose = 0)

eSVD_obj <- estimate_nuisance(input_obj = eSVD_obj,
                              bool_covariates_as_library = TRUE,
                              bool_library_includes_interept = TRUE,
                              bool_use_log = FALSE,
                              verbose = 0)

eSVD_obj <- compute_posterior(input_obj = eSVD_obj,
                              bool_adjust_covariates = FALSE,
                              alpha_max = 2 * max(dat@x),
                              bool_covariates_as_library = TRUE,
                              bool_stabilize_underdispersion = TRUE,
                              library_min = 0.1,
                              pseudocount = 0)

eSVD_obj <- compute_test_statistic(input_obj = eSVD_obj,
                                   verbose = 0)

save(dat, covariates, metadata, nuisance_vec, library_mat,
     true_cc_status, eSVD_obj, date_of_run,
     file = "tests/assets/synthetic_data.RData",
     compress = "xz")
