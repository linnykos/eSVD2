context("Test compute_test_statistic")

## compute_test_statistic is correct

test_that("compute_test_statistic works", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  covariates <- .get_object(eSVD_obj = eSVD_obj, what_obj = "covariates", which_fit = NULL)
  cc_vec <- covariates[,"case_control"]
  cc_levels <- sort(unique(cc_vec), decreasing = F)
  control_idx <- which(cc_vec == cc_levels[1])
  case_idx <- which(cc_vec == cc_levels[2])

  individual_vec <- metadata[,"individual"]
  control_individuals <- as.character(unique(individual_vec[control_idx]))
  case_individuals <- as.character(unique(individual_vec[case_idx]))

  res <- compute_test_statistic(input_obj = eSVD_obj$fit_First$posterior_mean_mat,
                                posterior_var_mat = eSVD_obj$fit_First$posterior_var_mat,
                                case_individuals = case_individuals,
                                control_individuals = control_individuals,
                                covariate_individual = "individual",
                                metadata = metadata)

  expect_true(2*mean(abs(res[true_cc_status == 1])) < mean(abs(res[true_cc_status == 2])))
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
  set.seed(10)
  vec <- numeric(0)
  for(i in 1:4){
    vec <- c(vec, rep(i, round(50*runif(1))))
  }
  metadata <- data.frame(individual = factor(vec))
  n <- length(vec)
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
    expect_true(all(res2[i,-idx] == 0))
  }

  response_vec <- runif(n)
  avg1 <- as.numeric(res %*% response_vec)
  avg2 <- sapply(1:4, function(i){
    idx <- which(metadata[,"individual"] == as.character(i))
    mean(response_vec[idx])
  })
  expect_true(sum(abs(avg1 - avg2)) <= 1e-5)
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
