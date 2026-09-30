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
#' @param input_obj             \code{eSVD} object output from \code{compute_posterior.eSVD}.
#' @param min_cells_per_individual  Minimum number of cells an individual must
#'                              contribute; see \code{compute_test_statistic.default}.
#' @param verbose               Integer.
#' @param ...                   Additional parameters.
#'
#' @return \code{eSVD} object with added elements, each a named numeric
#' vector with one entry per gene: \code{teststat_vec} (the Welch statistic),
#' \code{case_mean} and \code{control_mean} (the mean over the individuals of
#' each arm of their mean posterior expression), \code{case_var} and
#' \code{control_var} (the variance of each arm, the two quantities the
#' statistic divides by), and \code{log2fc_vec} and \code{log2fc_se_vec} (the
#' log2 fold change and its standard error; see
#' \code{compute_log_fold_change}). The individuals of each arm and
#' \code{min_cells_per_individual} are recorded in \code{param}.
#' @examples
#' set.seed(10)
#' sim <- generate_null(cell_per_person = 15, num_genes = 40,
#'                      num_individuals = 8)
#' esvd_obj <- initialize_esvd(dat = sim$obs_mat,
#'                             covariates = sim$covariates,
#'                             metadata_individual = sim$metadata_individual,
#'                             case_control_variable = "CC",
#'                             bool_intercept = TRUE,
#'                             k = 2,
#'                             lambda = 0.1)
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_Init",
#'                                                omitted_variables = "Log_UMI")
#' esvd_obj <- opt_esvd(input_obj = esvd_obj,
#'                      max_iter = 5,
#'                      offset_variables = setdiff(colnames(esvd_obj$covariates), "CC"),
#'                      fit_name = "fit_First",
#'                      fit_previous = "fit_Init")
#' esvd_obj <- reparameterization_esvd_covariates(input_obj = esvd_obj,
#'                                                fit_name = "fit_First",
#'                                                omitted_variables = "Log_UMI")
#' esvd_obj <- estimate_nuisance(input_obj = esvd_obj)
#' esvd_obj <- compute_posterior(input_obj = esvd_obj,
#'                               alpha_max = 2 * max(sim$obs_mat))
#' esvd_obj <- compute_test_statistic(input_obj = esvd_obj)
#' utils::head(esvd_obj$teststat_vec)
#' @export
compute_test_statistic.eSVD <- function(input_obj,
                                        min_cells_per_individual = 3,
                                        verbose = 0,
                                        ...){
  if(is.null(input_obj[["dat"]])){
    stop("`input_obj$dat` is missing. An object built with ",
         "`eSVD(bool_diet = TRUE)` has no count matrix or posterior matrices ",
         "and already carries its test statistics; the matrix path needs ",
         "`bool_diet = FALSE`")
  }
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"),
            all(!is.null(input_obj[["case_control"]])) && all(input_obj[["case_control"]] %in% c(0,1)) && length(input_obj[["case_control"]]) == nrow(input_obj[["dat"]]),
            all(!is.null(input_obj[["individual"]])) && all(is.factor(input_obj[["individual"]])) && length(input_obj[["individual"]]) == nrow(input_obj[["dat"]]))

  cc_vec <- input_obj[["case_control"]]
  cc_levels <- sort(unique(cc_vec), decreasing = FALSE)
  stopifnot(length(cc_levels) == 2)
  control_idx <- which(cc_vec == cc_levels[1])
  case_idx <- which(cc_vec == cc_levels[2])

  individual_vec <- input_obj[["individual"]]
  control_individuals <- unique(individual_vec[control_idx])
  case_individuals <- unique(individual_vec[case_idx])
  stopifnot(length(intersect(control_individuals, case_individuals)) == 0)

  param <- .format_param_test_statistic(case_individuals = case_individuals,
                                        control_individuals = control_individuals,
                                        min_cells_per_individual = min_cells_per_individual)
  input_obj$param <- .combine_two_named_lists(input_obj$param, param)
  # `.combine_two_named_lists` keeps an entry that is already there, so a
  # rerun on a changed cohort would leave the previous individuals in place
  # and `compute_log_fold_change` would divide by the wrong number of them.
  input_obj$param[names(param)] <- param

  latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
  posterior_mean_mat <- .get_object(eSVD_obj = input_obj, what_obj = "posterior_mean_mat", which_fit = latest_Fit)
  posterior_var_mat <- .get_object(eSVD_obj = input_obj, what_obj = "posterior_var_mat", which_fit = latest_Fit)

  res <- compute_test_statistic.default(
    input_obj = posterior_mean_mat,
    posterior_var_mat = posterior_var_mat,
    case_individuals = case_individuals,
    control_individuals = control_individuals,
    individual_vec = individual_vec,
    min_cells_per_individual = min_cells_per_individual,
    verbose = verbose
  )

  input_obj[["teststat_vec"]] <- res$teststat_vec
  input_obj[["case_mean"]] <- res$case_mean
  input_obj[["control_mean"]] <- res$control_mean
  input_obj[["case_var"]] <- res$case_var
  input_obj[["control_var"]] <- res$control_var
  input_obj[["log2fc_vec"]] <- res$log2fc_vec
  input_obj[["log2fc_se_vec"]] <- res$log2fc_se_vec
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
#' @param case_individuals     Vector of strings: the values of \code{individual_vec}
#'                             that are the case individuals.
#' @param control_individuals  Vector of strings: the values of \code{individual_vec}
#'                             that are the control individuals.
#' @param individual_vec       Vector of strings of length \eqn{n} (i.e., the number of cells) that denote which cell originates from which individual.
#' @param min_cells_per_individual  Minimum number of cells an individual must
#'                             contribute for the test to be computed. Individuals
#'                             with fewer are expected to have been dropped
#'                             upstream by \code{eSVD_helper}; reaching here with
#'                             one is an error rather than a silently noisy
#'                             statistic. Set to \code{0} to disable the check.
#' @param verbose              Integer.
#' @param ...                  Additional parameters.
#'
#' @return A list of named numeric vectors, each of length
#' \code{ncol(input_obj)} and carrying \code{colnames(input_obj)}:
#' \code{teststat_vec}, \code{case_mean}, \code{control_mean},
#' \code{case_var}, \code{control_var}, \code{log2fc_se_vec} and
#' \code{log2fc_vec}. The two variances are those of the mixture over
#' individuals that summarizes each arm, \eqn{\mathrm{mean}_i(\bar V_i) +
#' \mathrm{var}_i(m_i)}, where \eqn{m_i} and \eqn{\bar V_i} are individual
#' \eqn{i}'s average posterior mean and variance and the variance over
#' individuals divides by their number, with no Bessel correction. The
#' statistic is \code{(case_mean - control_mean) / sqrt(case_var / n1 +
#' control_var / n0)} for \code{n1} case and \code{n0} control individuals.
#' See \code{compute_log_fold_change} for the last two.
#' @export
compute_test_statistic.default <- function(input_obj,
                                           posterior_var_mat,
                                           case_individuals,
                                           control_individuals,
                                           individual_vec,
                                           min_cells_per_individual = 3,
                                           verbose = 0,
                                           ...) {
  stopifnot(inherits(input_obj, "matrix"))

  posterior_mean_mat <- input_obj
  stopifnot(all(dim(posterior_mean_mat) == dim(posterior_var_mat)),
            length(individual_vec) == nrow(posterior_mean_mat))

  .check_cohort_is_testable(case_individuals = case_individuals,
                            control_individuals = control_individuals,
                            individual_vec = individual_vec,
                            min_cells_per_individual = min_cells_per_individual)

  p <- ncol(posterior_mean_mat)

  if(verbose >= 1) print("Computing individual-level statistics")
  tmp <- .determine_individual_indices(case_individuals = case_individuals,
                                       control_individuals = control_individuals,
                                       individual_vec = individual_vec)
  all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
  avg_mat <- .construct_averaging_matrix(idx_list = all_indiv_idx,
                                         n = nrow(posterior_mean_mat))
  avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
  avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

  # see https://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
  if(verbose >= 1) print("Computing group-level statistics")
  case_row_idx <- 1:length(case_individuals)
  control_row_idx <- (length(case_individuals)+1):nrow(avg_posterior_mean_mat)
  case_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[case_row_idx,,drop = FALSE])
  control_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[control_row_idx,,drop = FALSE])
  case_gaussian_var <- .compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[case_row_idx,,drop = FALSE],
    avg_posterior_var_mat = avg_posterior_var_mat[case_row_idx,,drop = FALSE]
  )
  control_gaussian_var <- .compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[control_row_idx,,drop = FALSE],
    avg_posterior_var_mat = avg_posterior_var_mat[control_row_idx,,drop = FALSE]
  )

  if(verbose >= 1) print("Computing test statistics")
  n1 <- length(case_individuals)
  n2 <- length(control_individuals)
  teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
    (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))
  names(teststat_vec) <- colnames(posterior_mean_mat)

  lfc_res <- .compute_log2_fold_change(case_mean = case_gaussian_mean,
                                       control_mean = control_gaussian_mean,
                                       case_var = case_gaussian_var,
                                       control_var = control_gaussian_var,
                                       num_case = n1,
                                       num_control = n2)

  list(teststat_vec = teststat_vec,
       case_mean = case_gaussian_mean,
       control_mean = control_gaussian_mean,
       case_var = case_gaussian_var,
       control_var = control_gaussian_var,
       log2fc_se_vec = lfc_res$log2fc_se_vec,
       log2fc_vec = lfc_res$log2fc_vec)
}

.determine_individual_indices <- function(case_individuals,
                                          control_individuals,
                                          individual_vec){
  case_indiv_idx <- lapply(case_individuals, function(indiv){
    which(individual_vec == indiv)
  })
  control_indiv_idx <- lapply(control_individuals, function(indiv){
    which(individual_vec == indiv)
  })

  list(case_indiv_idx = case_indiv_idx,
       control_indiv_idx = control_indiv_idx)
}

.construct_averaging_matrix <- function(idx_list,
                                        n){
  # An empty index set is reachable through an unused factor level. Its row
  # of the averaging matrix would be all zeros, which "averages" to 0 and
  # silently biases the group mean.
  empty_idx <- which(sapply(idx_list, length) == 0)
  if(length(empty_idx) > 0){
    stop("`idx_list` has ", length(empty_idx), " empty element(s) at position(s) ",
         paste0(empty_idx, collapse = ", "),
         "; every individual must contribute at least one cell. ",
         "An unused factor level in the individual vector causes this")
  }
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

#' Check that a cohort can support a two-sample test
#'
#' Two conditions, both of which otherwise fail silently or obscurely much
#' further downstream.
#'
#' Individuals with very few cells are meant to have been dropped upstream by
#' the cohort filter; reaching here with one means the filter was skipped, and
#' the honest response is to stop rather than to compute a statistic from a
#' single cell's posterior.
#'
#' An arm with fewer than two individuals is a harder failure: the
#' Welch-Satterthwaite denominator contains \code{(v/n)^2/(n-1)}, so
#' \code{n = 1} makes the degrees of freedom \code{0} and every downstream
#' \code{stats::pt()} call returns \code{NaN}.
#'
#' @param case_individuals          Vector of individuals in the case arm.
#' @param control_individuals       Vector of individuals in the control arm.
#' @param individual_vec            \code{factor} of length \eqn{n} naming each
#'                                  cell's individual.
#' @param min_cells_per_individual  Minimum number of cells an individual must
#'                                  contribute. Set to \code{0} to disable.
#'
#' @return \code{invisible(TRUE)}; called for the error.
#' @noRd
.check_cohort_is_testable <- function(case_individuals,
                                      control_individuals,
                                      individual_vec,
                                      min_cells_per_individual = 3){
  if(length(case_individuals) < 2 || length(control_individuals) < 2){
    stop("each arm needs at least 2 individuals to form a two-sample test; ",
         "found ", length(case_individuals), " case and ",
         length(control_individuals), " control. ",
         "With one individual in an arm the Welch degrees of freedom are 0 ",
         "and every p-value is NaN")
  }

  if(min_cells_per_individual > 0){
    tested_individuals <- c(as.character(case_individuals),
                            as.character(control_individuals))
    cell_count_vec <- table(as.character(individual_vec))
    cell_count_vec <- cell_count_vec[names(cell_count_vec) %in%
                                       tested_individuals]

    sparse_idx <- which(cell_count_vec < min_cells_per_individual)
    if(length(sparse_idx) > 0){
      stop("individual(s) ",
           paste0(names(cell_count_vec)[sparse_idx], " (",
                  as.integer(cell_count_vec)[sparse_idx], " cells)",
                  collapse = ", "),
           " have fewer than ", min_cells_per_individual,
           " cells. Remove them upstream (see `eSVD_helper`), or pass ",
           "`min_cells_per_individual = 0` to proceed anyway")
    }
  }

  invisible(TRUE)
}

.format_param_test_statistic <- function(case_individuals,
                                         control_individuals,
                                         min_cells_per_individual){
  list(test_case_individuals = case_individuals,
       test_control_individuals = control_individuals,
       test_min_cells_per_individual = min_cells_per_individual)
}
