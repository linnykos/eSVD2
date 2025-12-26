#' Gene-wise posterior + test statistic + p-value without storing posterior matrices
#'
#' This function reproduces the behavior of
#' \code{compute_posterior()}, \code{compute_test_statistic()}, and
#' \code{compute_pvalue()}, but it computes everything one gene at a time
#' without ever allocating an n x p posterior mean/variance matrix.
#'
#' @param input_obj  eSVD object after nuisance estimation.
#' @param alpha_max  Maximum value of the prior mean (same role as in
#'                   \code{compute_posterior.eSVD}); default 1e3.
#' @param bool_adjust_covariates Boolean; if TRUE, adjust the prior by
#'                   confounding covariates (same as in
#'                   \code{compute_posterior.default}).
#' @param bool_covariates_as_library Boolean; if TRUE, include non–case-control
#'                   covariates in the covariate-adjusted library size.
#' @param bool_stabilize_underdispersion Boolean; if TRUE, mean-center
#'                   log10(nuisance_vec) when it suggests under-dispersion.
#' @param library_min Minimum value for the covariate-adjusted library size.
#' @param nuisance_lower_quantile Lower quantile at which to floor nuisance_vec.
#' @param pseudocount Numeric; additional count added to each entry in the
#'                   count matrix when forming the posterior.
#' @param verbose    Integer; controls printed messages.
#'
#' @return The input \code{eSVD} object with:
#'   \itemize{
#'     \item \code{teststat_vec}  — Welch t-statistics (length p).
#'     \item \code{case_mean}    — case Gaussian means per gene.
#'     \item \code{control_mean} — control Gaussian means per gene.
#'     \item \code{pvalue_list}  — list with \code{df_vec}, \code{fdr_vec},
#'           \code{gaussian_teststat}, \code{log10pvalue},
#'           \code{null_mean}, \code{null_sd}.
#'   }
#' @export
compute_test_per_gene <- function(input_obj,
                                  alpha_max = 1e3,
                                  bool_adjust_covariates = FALSE,
                                  bool_covariates_as_library = TRUE,
                                  bool_stabilize_underdispersion = TRUE,
                                  library_min = 1e-2,
                                  nuisance_lower_quantile = 0.01,
                                  pseudocount = 0,
                                  verbose = 0) {

  ## ------------------------------------------------------------
  ## 0. Basic checks and pull eSVD pieces
  ## ------------------------------------------------------------
  if(verbose > 0) print("Basic checks")
  stopifnot(
    inherits(input_obj, "eSVD"),
    "latest_Fit" %in% names(input_obj)
  )
  latest_Fit_name <- input_obj[["latest_Fit"]]
  stopifnot(
    latest_Fit_name %in% names(input_obj),
    inherits(input_obj[[latest_Fit_name]], "eSVD_Fit")
  )

  dat        <- .get_object(eSVD_obj = input_obj,
                            what_obj = "dat",
                            which_fit = NULL)
  covariates <- .get_object(eSVD_obj = input_obj,
                            what_obj = "covariates",
                            which_fit = NULL)
  esvd_res   <- .get_object(eSVD_obj = input_obj,
                            what_obj = NULL,
                            which_fit = latest_Fit_name)
  nuisance_vec <- .get_object(eSVD_obj = input_obj,
                              what_obj = "nuisance",
                              which_fit = latest_Fit_name)
  stopifnot(all(nuisance_vec >= 0))

  # case/control & individual info (same logic as compute_test_statistic.eSVD / .compute_df)
  stopifnot(
    !is.null(input_obj[["case_control"]]),
    length(input_obj[["case_control"]]) == nrow(dat),
    !is.null(input_obj[["individual"]]),
    length(input_obj[["individual"]]) == nrow(dat)
  )
  cc_vec         <- input_obj[["case_control"]]
  individual_vec <- input_obj[["individual"]]

  cc_levels <- sort(unique(cc_vec), decreasing = FALSE)
  stopifnot(length(cc_levels) == 2)
  control_idx <- which(cc_vec == cc_levels[1])
  case_idx    <- which(cc_vec == cc_levels[2])

  control_individuals <- unique(individual_vec[control_idx])
  case_individuals    <- unique(individual_vec[case_idx])
  stopifnot(length(intersect(control_individuals, case_individuals)) == 0)

  tmp_idx   <- .determine_individual_indices(
    case_individuals    = case_individuals,
    control_individuals = control_individuals,
    individual_vec      = individual_vec
  )
  all_indiv_idx <- c(tmp_idx$case_indiv_idx, tmp_idx$control_indiv_idx)
  avg_mat <- .construct_averaging_matrix(
    idx_list = all_indiv_idx,
    n        = nrow(dat)
  )

  n1 <- length(case_individuals)
  n2 <- length(control_individuals)

  case_row_idx    <- seq_len(n1)
  control_row_idx <- (n1 + 1):nrow(avg_mat)

  ## ------------------------------------------------------------
  ## 1. Setup covariate indices for posterior construction
  ##    (mirrors compute_posterior.default)
  ## ------------------------------------------------------------
  if(verbose > 0) print("Setting up posterior calculations")

  # These are stored in param by initialize_esvd / nuisance estimation
  case_control_variable <- .get_object(
    eSVD_obj  = input_obj,
    which_fit = "param",
    what_obj  = "init_case_control_variable"
  )
  library_size_variable <- .get_object(
    eSVD_obj  = input_obj,
    which_fit = "param",
    what_obj  = "init_library_size_variable"
  )
  bool_library_includes_interept <- .get_object(
    eSVD_obj  = input_obj,
    which_fit = "param",
    what_obj  = "nuisance_bool_library_includes_interept"
  )

  if (!is.null(pseudocount) && pseudocount > 0) {
    # we'll add this gene-wise, but record here for clarity
    # (no need to actually allocate dat + pseudocount)
  }

  if (is.null(case_control_variable)) {
    case_control_variable <- numeric(0)
  }
  case_control_idx <- which(colnames(covariates) == case_control_variable)

  library_size_variables <- library_size_variable
  if (bool_covariates_as_library) {
    library_size_variables <- unique(c(
      library_size_variables,
      setdiff(colnames(covariates),
              c("Intercept", case_control_variable))
    ))
  }
  if (bool_library_includes_interept) {
    library_size_variables <- unique(c("Intercept", library_size_variables))
  }
  library_idx <- which(colnames(covariates) %in% library_size_variables)
  idx_vec     <- c(case_control_idx, library_idx)

  # Pre-split covariates to avoid repeatedly indexing inside the loop
  cov_nolib      <- covariates[, -library_idx, drop = FALSE]
  cov_lib        <- covariates[,  library_idx, drop = FALSE]
  confounder_idx <- setdiff(seq_len(ncol(covariates)), idx_vec)
  cov_confounder <- if (length(confounder_idx) > 0) {
    covariates[, confounder_idx, drop = FALSE]
  } else {
    NULL
  }

  # Process nuisance_vec once (same as compute_posterior.default)
  if (!is.null(nuisance_lower_quantile)) {
    lower_bound <- stats::quantile(nuisance_vec, probs = nuisance_lower_quantile)
    nuisance_vec <- pmax(nuisance_vec, lower_bound)
  }
  if (bool_stabilize_underdispersion &&
      mean(log10(nuisance_vec)) > 0) {
    # center log10 nuisance if it suggests under-dispersion
    nuisance_vec <- 10^(scale(log10(nuisance_vec), center = TRUE, scale = FALSE))
  }

  ## ------------------------------------------------------------
  ## 2. Prepare containers for per-gene outputs
  ## ------------------------------------------------------------
  if(verbose > 0) print("Preparing per-gene outputs")

  n  <- nrow(dat)
  p  <- ncol(dat)
  gn <- colnames(dat)

  teststat_vec   <- numeric(p)
  case_mean_vec  <- numeric(p)
  control_mean_vec <- numeric(p)
  df_vec         <- numeric(p)

  names(teststat_vec)    <- gn
  names(case_mean_vec)   <- gn
  names(control_mean_vec) <- gn
  names(df_vec)          <- gn

  x_mat <- esvd_res$x_mat
  y_mat <- esvd_res$y_mat
  z_mat <- esvd_res$z_mat

  ## ------------------------------------------------------------
  ## 3. Gene-wise loop: posterior -> individual-level Gaussian -> t + df
  ## ------------------------------------------------------------

  if (verbose >= 1) message("Looping over ", p, " genes in compute_test_per_gene()")

  for (j in seq(p)) {
    if (verbose == 2 && p > 100 && (j %% 100 == 0)) {
      message("  gene ", j, "/", p, " (", gn[j], ")")
    }
    if (verbose >= 3) {
      message("  gene ", j, "/", p, " (", gn[j], ")")
    }

    # Counts for this gene
    if(verbose >= 4) print("Extracting gene")

    y_j <- as.numeric(dat[, j])
    if (!is.null(pseudocount) && pseudocount > 0) {
      y_j <- y_j + pseudocount
    }

    # 3a. Gamma-Poisson posterior *for this gene only*
    if(verbose >= 4) print("Computing natural parameters")
    nat1_j <- as.vector(tcrossprod(x_mat, y_mat[j, , drop = FALSE]))
    nat2_j <- as.vector(tcrossprod(cov_nolib, z_mat[j, -library_idx, drop = FALSE]))

    nat_nolib_j   <- nat1_j + nat2_j
    mean_nolib_j  <- exp(nat_nolib_j)
    if (!is.null(alpha_max)) {
      mean_nolib_j <- pmin(mean_nolib_j, alpha_max)
    }

    # gene-specific nuisance (column sweep in original code)
    nuisance_j <- nuisance_vec[j]

    # Alpha = mean_mat_nolib * nuisance_vec (column-wise)
    Alpha_j      <- mean_nolib_j * nuisance_j
    AplusAlpha_j <- y_j + Alpha_j

    # Confounder adjustment (if requested)
    if (bool_adjust_covariates && !is.null(cov_confounder)) {
      nat_confounder_j <- as.vector(
        tcrossprod(cov_confounder, z_mat[j, confounder_idx, drop = FALSE])
      )
      tmp              <- log(AplusAlpha_j)
      AplusAlpha_j     <- exp(tmp - nat_confounder_j)
    }

    # library_mat[, j] = exp( covariates[,library_idx] %*% t(z_mat[j, library_idx]) )
    if(verbose >= 4) print("Computing library size")
    library_j <- as.vector(
      tcrossprod(cov_lib, z_mat[j, library_idx, drop = FALSE])
    )
    library_j <- exp(library_j)
    if (!is.null(library_min)) {
      library_j <- pmax(library_j, library_min)
    }

    # SplusBeta = library_mat + nuisance_vec (column-wise)
    if(verbose >= 4) print("Computing posteriors")
    SplusBeta_j <- library_j + nuisance_j

    posterior_mean_j <- AplusAlpha_j / SplusBeta_j
    posterior_var_j  <- AplusAlpha_j / (SplusBeta_j^2)

    # 3b. Average posterior mean/var over individuals (one gene)
    if(verbose >= 4) print("Averaging over people")
    avg_mean_j <- as.vector(avg_mat %*% posterior_mean_j)
    avg_var_j  <- as.vector(avg_mat %*% posterior_var_j)

    # Split into case vs control individuals
    case_mean_indiv    <- avg_mean_j[case_row_idx]
    case_var_indiv     <- avg_var_j[case_row_idx]
    control_mean_indiv <- avg_mean_j[control_row_idx]
    control_var_indiv  <- avg_var_j[control_row_idx]

    # Mixture Gaussian means (same as Matrix::colMeans in the original 1-col case)
    if(verbose >= 4) print("Computing Gaussian distributions")
    case_gaussian_mean    <- mean(case_mean_indiv)
    control_gaussian_mean <- mean(control_mean_indiv)

    # Mixture Gaussian variances: E[var] + E[mean^2] - (E[mean])^2
    case_gaussian_var <- mean(case_var_indiv) +
      mean(case_mean_indiv^2) - case_gaussian_mean^2

    control_gaussian_var <- mean(control_var_indiv) +
      mean(control_mean_indiv^2) - control_gaussian_mean^2

    # 3c. Welch t-statistic for this gene (same formula as compute_test_statistic.default)
    if(verbose >= 4) print("Computing test statistic")
    t_j <- (case_gaussian_mean - control_gaussian_mean) /
      sqrt(case_gaussian_var / n1 + control_gaussian_var / n2)

    # 3d. Gene-specific df_j (same as .compute_df(), but per gene)
    if(verbose >= 4) print("Computing degree of freedom")
    numerator_j <- (case_gaussian_var / n1 + control_gaussian_var / n2)^2
    denominator_j <- (case_gaussian_var / n1)^2 / (n1 - 1) +
      (control_gaussian_var / n2)^2 / (n2 - 1)
    df_j <- numerator_j / denominator_j

    # 3e. Store per-gene summaries
    teststat_vec[j]      <- t_j
    df_vec[j]            <- df_j
    case_mean_vec[j]     <- case_gaussian_mean
    control_mean_vec[j]  <- control_gaussian_mean
  }

  ## ------------------------------------------------------------
  ## 4. Convert t + df to Gaussian test stats and empirical-null p-values
  ##    (mirrors compute_pvalue)
  ## ------------------------------------------------------------
  if(verbose > 0) print("Computing Gaussianized test statistics")
  # t + df -> Gaussian statistic (vectorized version of the original sapply)
  gaussian_teststat <- stats::qnorm(stats::pt(teststat_vec, df = df_vec))
  names(gaussian_teststat) <- names(teststat_vec)

  # empirical null & FDR (unchanged)
  if(verbose > 0) print("Computing multiple-testing adjusted via empirical null")
  fdr_res <- multtest(gaussian_teststat)
  fdr_vec <- fdr_res$fdr_vec
  names(fdr_vec) <- names(gaussian_teststat)

  null_mean <- fdr_res$null_mean
  null_sd   <- fdr_res$null_sd

  # two-sided log p-values using symmetric tail around null_mean
  if(verbose > 0) print("Computing p-values")
  logp_vec <- sapply(gaussian_teststat, function(x) {
    if (x < null_mean) {
      Rmpfr::pnorm(x,
                   mean  = null_mean,
                   sd    = null_sd,
                   log.p = TRUE)
    } else {
      Rmpfr::pnorm(null_mean - (x - null_mean),
                   mean  = null_mean,
                   sd    = null_sd,
                   log.p = TRUE)
    }
  })
  # convert log p to -log10(p) and account for two-sided test
  log10pvalue_vec <- -(logp_vec / log(10) + log10(2))
  names(log10pvalue_vec) <- names(teststat_vec)

  pvalue_list <- list(
    df_vec            = df_vec,
    fdr_vec           = fdr_vec,
    gaussian_teststat = gaussian_teststat,
    log10pvalue       = log10pvalue_vec,
    null_mean         = null_mean,
    null_sd           = null_sd
  )

  ## ------------------------------------------------------------
  ## 5. Attach per-gene outputs to object and return
  ## ------------------------------------------------------------
  if(verbose > 0) print("Returning")
  input_obj[["teststat_vec"]]  <- teststat_vec
  input_obj[["case_mean"]]     <- case_mean_vec
  input_obj[["control_mean"]]  <- control_mean_vec
  input_obj[["pvalue_list"]]   <- pvalue_list

  # NOTE: we *deliberately* do NOT create / store posterior_mean_mat or posterior_var_mat
  # in input_obj[[latest_Fit_name]]

  if (verbose >= 1) message("compute_test_per_gene(): finished.")

  input_obj
}
