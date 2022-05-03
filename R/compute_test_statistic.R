#' @export
compute_test_statistic <- function(input_obj, ...) UseMethod("compute_test_statistic")

#' @export
compute_test_statistic.eSVD <- function(input_obj,
                                        covariate_individual,
                                        metadata,
                                        verbose = 0){
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"),
            all(rownames(metadata) == rownames(input_obj$dat)),
            covariate_individual %in% colnames(metadata),
            is.factor(metadata[,covariate_individual]))

  case_control_variable <- .get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
  covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
  cc_vec <- covariates[,case_control_variable]
  cc_levels <- sort(unique(cc_vec), decreasing = F)
  stopifnot(length(cc_levels) == 2)
  control_idx <- which(cc_vec == cc_levels[1])
  case_idx <- which(cc_vec == cc_levels[2])

  individual_vec <- metadata[,covariate_individual]
  control_individuals <- unique(individual_vec[control_idx])
  case_individuals <- unique(individual_vec[case_idx])
  stopifnot(length(intersect(control_individuals, case_individuals)) == 0)

  param <- .format_param_test_statistic(case_individuals = case_individuals,
                                        control_individuals = control_individuals)
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)

  latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
  posterior_mean_mat <- .get_object(eSVD_obj = input_obj, what_obj = "posterior_mean_mat", which_fit = latest_Fit)
  posterior_var_mat <- .get_object(eSVD_obj = input_obj, what_obj = "posterior_var_mat", which_fit = latest_Fit)

  teststat_vec <- compute_test_statistic.default(
    input_obj = posterior_mean_mat,
    posterior_var_mat = posterior_var_mat,
    case_individuals = case_individuals,
    control_individuals = control_individuals,
    covariate_individual = covariate_individual,
    metadata = metadata,
    verbose = verbose
  )

  input_obj[["teststat_vec"]] <- teststat_vec
  input_obj
}

#' @export
compute_test_statistic.default <- function(input_obj,
                                           posterior_var_mat,
                                           case_individuals,
                                           control_individuals,
                                           covariate_individual,
                                           metadata,
                                           verbose = 0) {
  posterior_mean_mat <- input_obj
  stopifnot(all(dim(posterior_mean_mat) == dim(posterior_var_mat)),
            covariate_individual %in% colnames(metadata),
            all(rownames(metadata) == rownames(posterior_mean_mat)))

  p <- ncol(posterior_mean_mat)
  individual_stats <- lapply(1:p, function(j){
    if(verbose > 0 && p > 10 && j %% floor(p/10) == 0) cat('*')

    # next find the cells, then compute one gaussian per individual
    case_gaussians <- sapply(case_individuals, function(indiv){
      cell_names <- rownames(metadata)[which(metadata[,covariate_individual] == indiv)]
      cell_idx <- which(rownames(posterior_mean_mat) %in% cell_names)

      mean_val <- mean(posterior_mean_mat[cell_idx,j])
      var_val <- mean(posterior_var_mat[cell_idx,j])
      c(mean_val = mean_val, var_val = var_val)
    })

    control_gaussians <- sapply(control_individuals, function(indiv){
      cell_names <- rownames(metadata)[which(metadata[,covariate_individual] == indiv)]
      cell_idx <- which(rownames(posterior_mean_mat) %in% cell_names)

      mean_val <- mean(posterior_mean_mat[cell_idx,j])
      var_val <- mean(posterior_var_mat[cell_idx,j])
      c(mean_val = mean_val, var_val = var_val)
    })

    list(case_gaussians = case_gaussians,
         control_gaussians = control_gaussians)
  })

  # see https://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
  group_stats <- lapply(1:p, function(j){
    case_gaussians <- individual_stats[[j]]$case_gaussians
    control_gaussians <- individual_stats[[j]]$control_gaussians

    case_gaussian <- list(mean_val = mean(case_gaussians[1,]),
                          var_val = mean(case_gaussians[2,]) + mean(case_gaussians[1,]^2) - (mean(case_gaussians[1,]))^2,
                          n = ncol(case_gaussians))
    control_gaussian <- list(mean_val = mean(control_gaussians[1,]),
                             var_val = mean(control_gaussians[2,]) + mean(control_gaussians[1,]^2) - (mean(control_gaussians[1,]))^2,
                             n = ncol(control_gaussians))

    list(case_gaussian = case_gaussian,
         control_gaussian = control_gaussian)
  })

  teststat_vec <- sapply(1:p, function(j){
    case_gaussian <- group_stats[[j]]$case_gaussian
    control_gaussian <- group_stats[[j]]$control_gaussian

    n1 <- control_gaussian$n; n2 <- case_gaussian$n
    mean1 <- control_gaussian$mean_val; mean2 <- case_gaussian$mean_val
    cov1 <- control_gaussian$var_val; cov2 <- control_gaussian$var_val

    combined_cov <- cov1/n1 + cov2/n2
    (mean2 - mean1)/sqrt(combined_cov)
  })
  names(teststat_vec) <- colnames(posterior_mean_mat)

  teststat_vec
}

.format_param_test_statistic <- function(case_individuals,
                                         control_individuals){
  list(test_case_individuals = case_individuals,
       test_control_individuals = control_individuals)
}
