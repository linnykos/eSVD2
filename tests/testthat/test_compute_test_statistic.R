context("Test compute_test_statistic")

## compute_test_statistic is correct

test_that("compute_test_statistic works when using the true mean and variance", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))*.5, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))*.5, nrow = p, ncol = k)
  covariates <- cbind(c(rep(0, n/2), rep(1, n/2)),
                      matrix(abs(rnorm(n * 3, mean = 1, sd = 0.1)), nrow = n, ncol = 3))
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(2, p/2)), rep(1,p), rep(1,p), rep(1,p))
  colnames(z_mat) <-  colnames(covariates)
  case_control_variable <- "covariate_1"
  case_control_idx <- which(colnames(z_mat) == case_control_variable)
  nat_mat_nolib <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates[,case_control_idx], z_mat[,case_control_idx])
  library_mat <- exp(tcrossprod(covariates[,-case_control_idx], z_mat[,-case_control_idx]))
  nuisance_vec <- rep(c(5, 1, 1/5), times = 50)

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
  dat <- pmin(dat, 200)
  dat <- Matrix::Matrix(dat, sparse = T)
  rownames(dat) <- paste0("c", 1:n)
  colnames(dat) <- paste0("g", 1:p)
  metadata <- data.frame(individual = factor(rep(1:4, each = n/4)))
  rownames(metadata) <- rownames(dat)

  posterior_mean_mat <- gamma_mat
  posterior_var_mat <- sweep(gamma_mat, MARGIN = 2, STATS = nuisance_vec, FUN = "/")
  res <- compute_test_statistic(input_obj = posterior_mean_mat,
                                posterior_var_mat = posterior_var_mat,
                                case_individuals = c("3", "4"),
                                control_individuals = c("1", "2"),
                                covariate_individual = "individual",
                                metadata = metadata)

  expect_true(2*mean(abs(res[1:p/2])) < mean(abs(res[(p/2+1):p])))
})

test_that("compute_test_statistic works", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))*.5, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))*.5, nrow = p, ncol = k)
  covariates <- cbind(c(rep(0, n/2), rep(2, n/2)),
                      matrix(abs(rnorm(n * 3, mean = 1, sd = 0.1)), nrow = n, ncol = 3))
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(2, p/2)), rep(1,p), rep(1,p), rep(1,p))
  colnames(z_mat) <-  colnames(covariates)
  case_control_variable <- "covariate_1"
  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)
  expected_lib <- rowSums(nat_mat)
  nat_mat <- .mult_vec_mat(100/expected_lib, nat_mat)
  library_mat <- matrix(25, nrow = n, ncol = p)

  # mean_mat <- exp(nat_mat)
  # expected_lib <- rowSums(mean_mat)
  # tmp <- median(expected_lib)/expected_lib
  # library_mat <- matrix(10*tmp, nrow = n, ncol = p)

  nuisance_vec <- rep(c(5, 1, 1/5), times = 50)

  # Simulate data
  gamma_mat <- matrix(NA, nrow = n, ncol = p)
  dat <- matrix(NA, nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      gamma_mat[i,j] <- stats::rgamma(n = 1,
                                      shape = nuisance_vec[j]*exp(nat_mat[i,j]),
                                      rate = nuisance_vec[j])
      dat[i,j] <- stats::rpois(n = 1, lambda = library_mat[i,j] * gamma_mat[i,j])
    }
  }
  dat <- pmin(dat, 200)
  dat <- Matrix::Matrix(dat, sparse = T)
  rownames(dat) <- paste0("c", 1:n)
  colnames(dat) <- paste0("g", 1:p)
  covariates <- cbind(covariates, log(Matrix::rowSums(dat)))
  colnames(covariates)[ncol(covariates)] <- "Log_UMI"
  metadata <- data.frame(individual = factor(rep(1:4, each = n/4)))
  rownames(metadata) <- rownames(dat)

  # fit eSVD
  eSVD_obj <- initialize_esvd(dat = dat,
                              bool_intercept = T,
                              covariates = covariates,
                              case_control_variable = case_control_variable,
                              k = 5,
                              lambda = 0.01,
                              mixed_effect_variables = NULL,
                              offset_variables = NULL)
  logp_quant <- quantile(eSVD_obj$initial_Reg$log_pval, probs = 0.25)
  eSVD_obj <- apply_initial_threshold(eSVD_obj = eSVD_obj,
                                      pval_thres = exp(logp_quant))
  eSVD_obj <- opt_esvd(input_obj = eSVD_obj,
                       max_iter = 10)
  eSVD_obj <- estimate_nuisance(input_obj = eSVD_obj,
                                verbose = 0)
  eSVD_obj <- compute_posterior(input_obj = eSVD_obj)

  res <- compute_test_statistic(input_obj = eSVD_obj,
                                covariate_individual = "individual",
                                metadata = metadata)
  # plot(res[["teststat_vec"]])

  expect_true(is.numeric(res[["teststat_vec"]]))
  expect_true(length(res[["teststat_vec"]]) == ncol(dat))
  expect_true(2*mean(abs(res[["teststat_vec"]][1:p/2])) < mean(abs(res[["teststat_vec"]][(p/2+1):p])))
})

######################

## .determine_individual_indices is correct

test_that(".determine_individual_indices works", {
  n <- 100
  metadata <- data.frame(individual = factor(rep(1:4, each = n/4)))
  rownames(metadata) <- paste0("c", 1:n)
  res <- .determine_individual_indices(case_individuals = c("1", "2"),
                                       control_individuals = c("3", "4"),
                                       covariate_individual = "individual",
                                       metadata = metadata)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("case_indiv_idx", "control_indiv_idx"))))

  tmp <- c(res$case_indiv_idx, res$control_indiv_idx)
  expect_true(is.list(tmp))
  expect_true(length(tmp) == 4)
  expect_true(length(unique(unlist(tmp))) == n)
})

#######################

## .construct_averaging_matrix is correct

test_that(".construct_averaging_matrix works", {
  n <- 100
  metadata <- data.frame(individual = factor(rep(1:4, each = n/4)))
  rownames(metadata) <- paste0("c", 1:n)
  tmp <- .determine_individual_indices(case_individuals = c("1", "2"),
                                       control_individuals = c("3", "4"),
                                       covariate_individual = "individual",
                                       metadata = metadata)
  all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
  res <- .construct_averaging_matrix(idx_list = all_indiv_idx,
                                     n = n)

  expect_true(inherits(res, "dgCMatrix"))
  expect_true(all(dim(res) == c(4,n)))
  res2 <- as.matrix(res)
  for(i in 1:4){
    idx <- which(metadata[,"individual"] == as.character(i))
    expect_true(all(res2[i,idx] == 1))
    expect_true(all(res2[i,-idx] == 0))
  }
})

###############################

## .compute_mixture_gaussian_variance is correct

test_that(".compute_mixture_gaussian_variance works", {
  set.seed(10)
  n <- 5
  p <- 10
  avg_posterior_mean_mat <- matrix(runif(n*p), nrow = n, ncol = p)
  avg_posterior_var_mat <- matrix(runif(n*p), nrow = n, ncol = p)

  res <- .compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat,
    avg_posterior_var_mat = avg_posterior_var_mat
  )
  res2 <- sapply(1:p, function(j){
    mean(avg_posterior_var_mat[,j]) + mean(avg_posterior_mean_mat[,j]^2) - mean(avg_posterior_mean_mat[,j])^2
  })

  expect_true(sum(abs(res - res2)) <= 1e-6)
})
