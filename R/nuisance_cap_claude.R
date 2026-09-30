#' Check the multiplier of the cap on the nuisance rate
#'
#' @param cap_multiplier  The value to check.
#'
#' @returns \code{invisible()}; called for its error.
#' @noRd
.check_cap_multiplier <- function(cap_multiplier){
  if(!is.numeric(cap_multiplier) || length(cap_multiplier) != 1 ||
     is.na(cap_multiplier) || cap_multiplier <= 0){
    stop("`cap_multiplier` must be one positive number (`Inf` for no cap); ",
         "received ",
         if(is.null(cap_multiplier)) "NULL" else paste0(cap_multiplier, collapse = ", "))
  }

  invisible()
}

#' Check the floor of the nuisance rate against the cap
#'
#' Both are multiples of the median library size of the gene, so a floor at
#' or above the cap would leave every gene above its cap.
#'
#' @param min_val         The value to check.
#' @param cap_multiplier  The cap, already checked by
#'                        \code{.check_cap_multiplier}.
#'
#' @returns \code{invisible()}; called for its error.
#' @noRd
.check_min_val <- function(min_val, cap_multiplier){
  if(!is.numeric(min_val) || length(min_val) != 1 || is.na(min_val) ||
     !is.finite(min_val) || min_val <= 0){
    stop("`min_val` must be one positive finite number; received ",
         if(is.null(min_val)) "NULL" else paste0(min_val, collapse = ", "))
  }
  if(min_val >= cap_multiplier){
    stop("`min_val` (", min_val, ") must be below `cap_multiplier` (",
         cap_multiplier, "); both are multiples of the median library size")
  }

  invisible()
}

# Around the Poisson limit (beta -> Inf) the marginal log-likelihood of one
# gene is
#   l(beta) = l_Poisson + D / (2 * beta) + O(beta^-2),
#   D = sum_i [ (A_i - m_i)^2 - A_i ] / mu_i,   m_i = mu_i * s_i,
# from the expansion of the negative binomial with size r_i = mu_i * beta in
# 1 / r_i. D is the score statistic for over-dispersion, see
# Dean and Lawless (1989), https://doi.org/10.1080/01621459.1989.10478792

#' Compute the score statistic that decides whether a finite rate exists
#'
#' For one gene. A gene with \code{D <= 0} has counts that are no more
#' variable around the fit than a Poisson, so its likelihood rises all the
#' way to an infinite rate and \code{gamma_rate} returns wherever its
#' iterations stopped. Called inside the per-gene loop of
#' \code{.estimate_nuisance_matrix}, beside the estimate, so that each
#' column of the counts is extracted once.
#'
#' @param x_vec   Numeric vector, the gene's counts over cells.
#' @param mu_vec  Numeric vector, its fitted mean without the library size.
#' @param s_vec   Numeric vector, its library size.
#'
#' @returns One number. It is not finite when \code{mu_vec} holds a zero or
#' a non-finite value.
#' @noRd
.compute_boundary_statistic <- function(x_vec,
                                        mu_vec,
                                        s_vec){
  sum(((x_vec - mu_vec * s_vec)^2 - x_vec) / mu_vec)
}

#' Apply the cap to the maximum-likelihood rates and label every gene
#'
#' The one place the rule lives, shared by \code{estimate_nuisance} and
#' \code{recompute_pvalue}: \code{max(min(MLE, cap_multiplier * m), min_val *
#' m)}, with \code{m} the median library size of the gene and the MLE of a
#' gene at the boundary being infinite.
#'
#' @param nuisance_mle_vec    The maximum-likelihood rates, one per gene,
#'                            already floored at \code{min_val} times
#'                            \code{library_median_vec}. Its names are
#'                            carried to the output.
#' @param library_median_vec  The median over cells of each gene's library
#'                            size.
#' @param bool_boundary_vec   Logical, one per gene: is the gene at the
#'                            boundary (see \code{.compute_boundary_statistic})?
#' @param bool_failed_vec     Logical, one per gene: did both estimation
#'                            routes fail?
#' @param cap_multiplier      One positive number, \code{Inf} for no cap.
#' @param min_val             Minimum value of the rate, as a multiple of
#'                            \code{library_median_vec}; below
#'                            \code{cap_multiplier}.
#'
#' @returns List with \code{num_boundary}, \code{num_capped} (every gene whose
#' rate the cap replaced, the boundary genes among them included),
#' \code{nuisance_status} (factor with levels \code{estimated}, \code{capped},
#' \code{boundary}, \code{failed}) and \code{nuisance_vec}.
#' @noRd
.apply_nuisance_cap <- function(nuisance_mle_vec,
                                library_median_vec,
                                bool_boundary_vec,
                                bool_failed_vec,
                                cap_multiplier,
                                min_val){
  p <- length(nuisance_mle_vec)
  stopifnot(length(library_median_vec) == p,
            length(bool_boundary_vec) == p,
            length(bool_failed_vec) == p)

  # `Inf * 0` is NaN, and a library size can underflow to 0.
  if(is.infinite(cap_multiplier)){
    cap_vec <- rep(Inf, p)
  } else {
    cap_vec <- cap_multiplier * as.numeric(library_median_vec)
  }

  # A boundary gene has no finite maximum-likelihood rate, so min(MLE, cap)
  # is the cap itself. What the optimizer returned for it is where its
  # iterations stopped (about 1e7 on one route, exp(10) on the other), which
  # is below the cap when the library size is in the thousands.
  bool_at_cap_vec <- bool_boundary_vec & !bool_failed_vec & is.finite(cap_vec)
  # A failed gene sits at the floor by assignment; the cap did not act on it.
  bool_capped_vec <- !bool_failed_vec &
    (nuisance_mle_vec > cap_vec | bool_at_cap_vec)
  nuisance_vec <- pmin(nuisance_mle_vec, cap_vec)
  nuisance_vec[bool_at_cap_vec] <- cap_vec[bool_at_cap_vec]
  # The floor is in the units of the cap and below it (`.check_min_val`), so
  # it can only move a failed gene, whose value is 0, and no gene ends above
  # its cap.
  nuisance_vec <- pmax(nuisance_vec, min_val * as.numeric(library_median_vec))
  names(nuisance_vec) <- names(nuisance_mle_vec)

  status_vec <- rep("estimated", p)
  status_vec[bool_capped_vec] <- "capped"
  status_vec[bool_boundary_vec] <- "boundary"
  status_vec[bool_failed_vec] <- "failed"
  nuisance_status <- factor(status_vec, levels = .nuisance_status_levels())
  names(nuisance_status) <- names(nuisance_mle_vec)

  list(num_boundary = sum(nuisance_status == "boundary"),
       num_capped = sum(bool_capped_vec),
       nuisance_status = nuisance_status,
       nuisance_vec = nuisance_vec)
}

.nuisance_status_levels <- function(){
  c("estimated", "capped", "boundary", "failed")
}

#' Which covariates make up the library size
#'
#' The one place the rule lives. Shared by \code{estimate_nuisance.eSVD},
#' \code{compute_posterior.default}, \code{compute_test_per_gene} and
#' \code{plot_fitted_vs_observed}, which must all use the same library.
#'
#' @param covariates             Covariate matrix, columns named.
#' @param case_control_variable  Name of the case-control column, or
#'                               \code{NULL} or an empty vector for none.
#' @param library_size_variable  Name of the column of the observed library
#'                               size.
#' @param bool_covariates_as_library      See \code{estimate_nuisance.eSVD}.
#' @param bool_library_includes_interept  See \code{estimate_nuisance.eSVD}.
#'
#' @returns Integer vector of column indices of \code{covariates}.
#' @noRd
.nuisance_library_idx <- function(covariates,
                                  case_control_variable,
                                  library_size_variable,
                                  bool_covariates_as_library,
                                  bool_library_includes_interept){
  if(is.null(case_control_variable)) case_control_variable <- character(0)

  library_size_variables <- library_size_variable
  if(bool_covariates_as_library){
    library_size_variables <- c(library_size_variables,
                                setdiff(colnames(covariates),
                                        c("Intercept", case_control_variable)))
  }
  if(bool_library_includes_interept){
    library_size_variables <- c("Intercept", library_size_variables)
  }

  which(colnames(covariates) %in% library_size_variables)
}
