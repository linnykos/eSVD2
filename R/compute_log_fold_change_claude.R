#' Compute the log2 fold change and its standard error
#'
#' Reports, for every gene, the log2 fold change between the case and the
#' control individuals and a standard error for it, on the scale and with the
#' unit of replication that 'DESeq2', 'dreamlet' and 'NEBULA' report theirs.
#' \code{compute_test_statistic} and \code{compute_test_per_gene} call this
#' for you, so the two vectors are part of the ordinary output of \code{eSVD};
#' it is exported for objects on which they need to be recomputed.
#'
#' \strong{The estimate.} Write \eqn{\bar m_{1j}} and \eqn{\bar m_{0j}} for
#' \code{case_mean} and \code{control_mean}: gene \eqn{j}'s posterior mean
#' expression, averaged over the cells of each individual and then over the
#' individuals of the arm. The fold change is
#' \deqn{\widehat{\beta}_j = \log_2(\bar m_{1j} / \bar m_{0j}).}
#'
#' \strong{The standard error.} Write \eqn{v_{1j}} and \eqn{v_{0j}} for
#' \code{case_var} and \code{control_var}, the variances the Welch statistic
#' divides by, and \eqn{n_1}, \eqn{n_0} for the number of case and control
#' \emph{individuals}. Treating the two arms as independent, the delta method
#' gives
#' \deqn{\widehat{\mathrm{SE}}(\widehat{\beta}_j) = \frac{1}{\ln 2}
#' \sqrt{\frac{v_{1j}}{n_1 \bar m_{1j}^2} + \frac{v_{0j}}{n_0 \bar m_{0j}^2}}.}
#' Multiply by \code{log(2)} for the standard error of the natural-log fold
#' change, which is the scale 'NEBULA' reports on. The standard error of the
#' difference \eqn{\bar m_{1j} - \bar m_{0j}} itself, the denominator of
#' \code{teststat_vec}, is \code{sqrt(case_var / n1 + control_var / n0)}.
#'
#' \strong{What it is a standard error of.} The target is the log2 ratio of
#' the arm-level mean expression, each arm mean being an average over
#' individuals of their covariate- and depth-adjusted expression. The
#' variance is divided by the number of individuals and never by the number
#' of cells, so sequencing more cells from the same individuals does not
#' shrink it. That is the property it shares with a pseudobulk analysis in
#' 'DESeq2' or 'dreamlet' and with the standard error 'NEBULA' reports for a
#' subject-level predictor. Things to know when setting the numbers side by
#' side:
#' \itemize{
#'   \item It is a ratio of arithmetic means, as in 'DESeq2' and 'NEBULA'.
#'         'dreamlet' reports a difference of mean log expression, which is a
#'         different quantity whenever the two arms differ in spread.
#'   \item Each arm's variance is the variance of a mixture over individuals,
#'         \eqn{v_g = \mathrm{mean}_i(\bar V_i) + \mathrm{var}_i(m_i)}, where
#'         \eqn{m_i} and \eqn{\bar V_i} are individual \eqn{i}'s average
#'         posterior mean and variance. The first term is an average per-cell
#'         posterior variance, is not divided by the number of cells, and is
#'         usually most of \eqn{v_g}. The standard error is therefore larger,
#'         often several times larger, than the standard deviation of the
#'         fold change when individuals are resampled and the fit is held
#'         fixed.
#'   \item It treats the fitted factorization and the nuisance parameters as
#'         known, as 'DESeq2' treats its dispersions. This matters most for a
#'         gene whose estimated nuisance rate (\code{nuisance_vec}) is very
#'         large. The posterior of such a gene concentrates on the fitted
#'         prediction, the first term of \eqn{v_g} tends to zero, and the
#'         standard error reflects only how much the fitted predictions
#'         differ between individuals. It can then be much smaller than the
#'         error in the fold change. Inspect \code{nuisance_vec} before
#'         relying on an unusually small standard error.
#'   \item Expression is adjusted for each cell's observed sequencing depth,
#'         so the fold change is relative to total expression. If the genes
#'         that differ between the arms make up a large share of the counts,
#'         every gene's fold change is shifted by the same amount in the
#'         opposite direction.
#'   \item \code{log2fc_vec / log2fc_se_vec} does not reproduce the p-value in
#'         \code{pvalue_list}. That p-value comes from the statistic on the
#'         linear scale, its Welch-Satterthwaite degrees of freedom, and an
#'         empirical null fitted across genes.
#' }
#' The fold change is a ratio of posterior means, which are shrunk toward the
#' fitted low-rank prediction; it is therefore a shrunken estimate, closer in
#' kind to the output of \code{lfcShrink} in 'DESeq2' than to an unshrunken
#' maximum likelihood estimate.
#'
#' @param input_obj  \code{eSVD} object output from
#'                   \code{compute_test_statistic} or
#'                   \code{compute_test_per_gene}. It must carry
#'                   \code{case_mean}, \code{control_mean}, \code{case_var}
#'                   and \code{control_var} (named numeric vectors, one entry
#'                   per gene) and, in \code{param}, the individuals of each
#'                   arm. The count matrix and the posterior matrices are not
#'                   needed, so an object built with
#'                   \code{eSVD(bool_diet = TRUE)} is accepted.
#' @param verbose    Integer; \code{0} is silent.
#'
#' @returns The \code{eSVD} object with two named numeric vectors added or
#' replaced, one entry per gene in the order of \code{case_mean}:
#' \code{log2fc_vec}, the log2 fold change (positive when expression is
#' higher in the cases), and \code{log2fc_se_vec}, its standard error. A gene
#' whose \code{case_mean} or \code{control_mean} is \code{NA}, as for the
#' all-zero genes reinserted by \code{eSVD_helper}, is \code{NA} in both. A
#' gene whose mean is not positive is also \code{NA} in both, with a warning.
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
#' # compute_test_statistic has already stored both vectors
#' utils::head(cbind(log2fc = esvd_obj$log2fc_vec,
#'                   se = esvd_obj$log2fc_se_vec))
#' # the standard error on the natural-log scale
#' utils::head(esvd_obj$log2fc_se_vec * log(2))
#' esvd_obj <- compute_log_fold_change(input_obj = esvd_obj)
#' @export
compute_log_fold_change <- function(input_obj,
                                    verbose = 0){
  stopifnot(inherits(input_obj, "eSVD"))

  # `[[` throughout: `$` on a list matches names partially.
  required_vec <- c("case_mean", "control_mean", "case_var", "control_var")
  missing_vec <- required_vec[!required_vec %in% names(input_obj)]
  if(length(missing_vec) > 0){
    stop("`input_obj` has no `", paste0(missing_vec, collapse = "`, `"),
         "`. Run `compute_test_statistic` or `compute_test_per_gene` first. ",
         "An object built by eSVD2 before version 1.1.0 stores the arm means ",
         "but not the arm variances, so it has to be rebuilt: the variances ",
         "cannot be recovered from an object made with ",
         "`eSVD(bool_diet = TRUE)`")
  }

  case_individuals <- input_obj[["param"]][["test_case_individuals"]]
  control_individuals <- input_obj[["param"]][["test_control_individuals"]]
  if(length(case_individuals) == 0 || length(control_individuals) == 0){
    stop("`input_obj$param` does not record `test_case_individuals` and ",
         "`test_control_individuals`, so the number of individuals in each ",
         "arm is unknown. Run `compute_test_statistic` or ",
         "`compute_test_per_gene` first")
  }

  if(verbose > 0) print("Computing the log2 fold change and its standard error")
  res <- .compute_log2_fold_change(
    case_mean = input_obj[["case_mean"]],
    control_mean = input_obj[["control_mean"]],
    case_var = input_obj[["case_var"]],
    control_var = input_obj[["control_var"]],
    num_case = length(case_individuals),
    num_control = length(control_individuals)
  )

  input_obj[["log2fc_vec"]] <- res$log2fc_vec
  input_obj[["log2fc_se_vec"]] <- res$log2fc_se_vec
  input_obj
}

#' Log2 fold change and its delta-method standard error
#'
#' The one place the formula lives; \code{compute_test_statistic.default},
#' \code{compute_test_per_gene} and \code{compute_log_fold_change} all call
#' it, so the matrix path and the per-gene path cannot drift.
#'
#' With \eqn{g(a, b) = \log_2(a / b)}, the gradient is
#' \eqn{(1 / (a \ln 2), -1 / (b \ln 2))}, so for independent arms
#' \eqn{\mathrm{Var}(g) \approx \mathrm{Var}(a) / (a \ln 2)^2 +
#' \mathrm{Var}(b) / (b \ln 2)^2}, with \eqn{\mathrm{Var}(a) =}
#' \code{case_var / num_case} and \eqn{\mathrm{Var}(b) =}
#' \code{control_var / num_control}.
#'
#' A gene whose four inputs are all \code{NA} is \code{NA} in both outputs,
#' silently: that is how the all-zero genes padded by \code{.reinsert_genes}
#' arrive. Any other gene with an input that is not finite, a mean that is
#' not positive, or a variance that is negative is set to \code{NA} in both
#' with one warning, because \code{log2()} and \code{sqrt()} would otherwise
#' return \code{-Inf} or \code{NaN} without a word. A posterior mean from
#' \code{compute_posterior} is positive by construction, so this is reachable
#' only through the matrix method with a caller's own matrices.
#'
#' @param case_mean     Named numeric vector, one entry per gene: the mean
#'                      over case individuals of their mean posterior
#'                      expression.
#' @param control_mean  The same for the control individuals.
#' @param case_var      Named numeric vector, one entry per gene: the mixture
#'                      variance of the case arm, as returned by
#'                      \code{.compute_mixture_gaussian_variance}.
#' @param control_var   The same for the control individuals.
#' @param num_case      Number of case individuals (not cells).
#' @param num_control   Number of control individuals (not cells).
#'
#' @returns A list with \code{log2fc_se_vec} and \code{log2fc_vec}, both
#' carrying \code{names(case_mean)}.
#' @noRd
.compute_log2_fold_change <- function(case_mean,
                                      control_mean,
                                      case_var,
                                      control_var,
                                      num_case,
                                      num_control){
  stopifnot(length(control_mean) == length(case_mean),
            length(case_var) == length(case_mean),
            length(control_var) == length(case_mean),
            length(num_case) == 1, num_case >= 1,
            length(num_control) == 1, num_control >= 1)

  # A gene whose four inputs are ALL NA is a gene that was deliberately
  # padded, and it passes through as NA without a warning. One NA among valid
  # inputs is not that, and neither is NaN (`is.na(NaN)` is TRUE, so it is
  # separated out here).
  is_padded <- function(vec){is.na(vec) & !is.nan(vec)}
  padded_idx <- which(is_padded(case_mean) & is_padded(control_mean) &
                        is_padded(case_var) & is_padded(control_var))
  invalid_idx <- which(!is.finite(case_mean) | !is.finite(control_mean) |
                         !is.finite(case_var) | !is.finite(control_var) |
                         case_mean <= 0 | control_mean <= 0 |
                         case_var < 0 | control_var < 0)
  invalid_idx <- setdiff(invalid_idx, padded_idx)
  valid_idx <- setdiff(seq_along(case_mean), c(padded_idx, invalid_idx))
  if(length(invalid_idx) > 0){
    warning(length(invalid_idx), " gene(s) have a case or control mean that ",
            "is not positive, or a mean or variance that is negative or not ",
            "finite; their log2 fold change and its standard error are set ",
            "to NA")
  }

  # Computed on the valid genes only, so that neither `log2()` nor `sqrt()`
  # ever sees an input it would turn into NaN, and so that every other gene
  # is NA_real_ and not a platform-dependent mix of NA and NaN.
  log2fc_vec <- rep(NA_real_, length(case_mean))
  log2fc_se_vec <- rep(NA_real_, length(case_mean))
  log2fc_vec[valid_idx] <- log2(case_mean[valid_idx] /
                                  control_mean[valid_idx])
  log2fc_se_vec[valid_idx] <- sqrt(
    case_var[valid_idx] / (num_case * case_mean[valid_idx]^2) +
      control_var[valid_idx] / (num_control * control_mean[valid_idx]^2)
  ) / log(2)
  names(log2fc_vec) <- names(case_mean)
  names(log2fc_se_vec) <- names(case_mean)

  list(log2fc_se_vec = log2fc_se_vec,
       log2fc_vec = log2fc_vec)
}
