#' Compute test statistics
#'
#' Generic function interface
#'
#' @param input_obj Main object
#' @param ...       Additional parameters
#'
#' @return Output dependent on class of \code{input_obj}
#' @export
compute_test_statistic <- function(input_obj, ...) {UseMethod("compute_test_statistic")}

#' Compute test statistics for eSVD object
#'
#' @param input_obj             \code{eSVD} object outputed from \code{compute_posterior.eSVD}.
#' @param covariate_individual  A string of the column name of \code{metadata} which depicts the individual
#'                              each cell originates from. Notably, this column in \code{metadata} should be a factor vector.
#' @param metadata              \code{data frame} object with \eqn{n} rows with the same row names as \code{input_obj$dat}
#'                              where the columns represent the different covariates.
#'                              Notably, this should can contain categorical variables.
#' @param verbose               Integer.
#' @param ...                   Additional parameters.
#'
#' @return \code{eSVD} object with added element \code{"teststat_vec"}
#' @export
compute_test_statistic.eSVD <- function(input_obj,
                                        covariate_individual,
                                        metadata,
                                        verbose = 0,
                                        ...){
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

  res <- compute_test_statistic.default(
    input_obj = posterior_mean_mat,
    posterior_var_mat = posterior_var_mat,
    case_individuals = case_individuals,
    control_individuals = control_individuals,
    covariate_individual = covariate_individual,
    metadata = metadata,
    verbose = verbose
  )

  input_obj[["teststat_vec"]] <- res$teststat_vec
  input_obj[["diffmean_vec"]] <- res$diffmean_vec
  input_obj
}

#' Compute test statistics for matrices
#'
#' @param input_obj            Posterior mean matrix (a \code{matrix}) where the \eqn{n} rows represent cells
#'                             and \eqn{p} columns represent genes.
#'                             The rows and columns of the matrix should be named.
#' @param posterior_var_mat    Posterior variance matrix (a \code{matrix}) where the \eqn{n} rows represent cells
#'                             and \eqn{p} columns represent genes.
#'                             The rows and columns of the matrix should be the same as those in \code{input_obj}.
#' @param case_individuals     Vector of strings representing the individuals in \code{metadata[,covariate_individual]}
#'                             that are the case individuals.
#' @param control_individuals  Vector of strings representing the individuals in \code{metadata[,covariate_individual]}
#'                             that are the control individuals.
#' @param covariate_individual A string of the column name of \code{metadata} which depicts the individual
#'                             each cell originates from. Notably, this column in \code{metadata} should be a factor vector.
#' @param metadata             \code{data frame} object with \eqn{n} rows with the same row names as \code{input_obj$dat}
#'                             where the columns represent the different covariates.
#'                             Notably, this should can contain categorical variables.
#' @param verbose              Integer.
#' @param ...                  Additional parameters.
#'
#' @return A vector of test statistics of length \code{ncol(input_obj)}
#' @export
compute_test_statistic.default <- function(input_obj,
                                           posterior_var_mat,
                                           case_individuals,
                                           control_individuals,
                                           covariate_individual,
                                           metadata,
                                           verbose = 0,
                                           ...) {
  stopifnot(inherits(input_obj, "matrix"))

  posterior_mean_mat <- input_obj
  stopifnot(all(dim(posterior_mean_mat) == dim(posterior_var_mat)),
            covariate_individual %in% colnames(metadata),
            all(rownames(metadata) == rownames(posterior_mean_mat)))

  p <- ncol(posterior_mean_mat)

  if(verbose >= 1) print("Computing individual-level statistics")
  tmp <- .determine_individual_indices(case_individuals = case_individuals,
                                       control_individuals = control_individuals,
                                       covariate_individual = covariate_individual,
                                       metadata = metadata)
  all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
  avg_mat <- .construct_averaging_matrix(idx_list = all_indiv_idx,
                                         n = nrow(posterior_mean_mat))
  avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
  avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

  # see https://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
  if(verbose >= 1) print("Computing group-level statistics")
  case_row_idx <- 1:length(case_individuals)
  control_row_idx <- (length(case_individuals)+1):nrow(avg_posterior_mean_mat)
  case_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[case_row_idx,,drop = F])
  control_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[control_row_idx,,drop = F])
  case_gaussian_var <- .compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[case_row_idx,,drop = F],
    avg_posterior_var_mat = avg_posterior_var_mat[case_row_idx,,drop = F]
  )
  control_gaussian_var <- .compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[control_row_idx,,drop = F],
    avg_posterior_var_mat = avg_posterior_var_mat[control_row_idx,,drop = F]
  )

  if(verbose >= 1) print("Computing test statistics")
  n1 <- length(case_individuals)
  n2 <- length(control_individuals)
  teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
    (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))
  names(teststat_vec) <- colnames(posterior_mean_mat)
  diffmean_vec <- case_gaussian_mean - control_gaussian_mean
  names(diffmean_vec) <- colnames(posterior_mean_mat)

  list(teststat_vec = teststat_vec,
       diffmean_vec = diffmean_vec)
}

.determine_individual_indices <- function(case_individuals,
                                          control_individuals,
                                          covariate_individual,
                                          metadata){
  case_indiv_idx <- lapply(case_individuals, function(indiv){
    which(metadata[,covariate_individual] == indiv)
  })
  control_indiv_idx <- lapply(control_individuals, function(indiv){
    which(metadata[,covariate_individual] == indiv)
  })

  list(case_indiv_idx = case_indiv_idx,
       control_indiv_idx = control_indiv_idx)
}

.construct_averaging_matrix <- function(idx_list,
                                        n){
  tmp <- unlist(idx_list)
  stopifnot(max(table(tmp)) == 1, max(tmp) <= n, min(tmp) >= 1, all(tmp %% 1 == 0))

  averaging_indices <- do.call(rbind, lapply(1:length(idx_list), function(i){
    cbind(rep(i, length(idx_list[[i]])), idx_list[[i]], rep(1/length(idx_list[[i]]), length(idx_list[[i]])))
  }))
  Matrix::sparseMatrix(i = averaging_indices[,1],
                       j = averaging_indices[,2],
                       x = averaging_indices[,3],
                       dims = c(length(idx_list), n))
}

.compute_mixture_gaussian_variance <- function(avg_posterior_mean_mat,
                                               avg_posterior_var_mat){
  Matrix::colMeans(avg_posterior_var_mat) + Matrix::colMeans(avg_posterior_mean_mat^2) -
    Matrix::colMeans(avg_posterior_mean_mat)^2
}

.format_param_test_statistic <- function(case_individuals,
                                         control_individuals){
  list(test_case_individuals = case_individuals,
       test_control_individuals = control_individuals)
}
