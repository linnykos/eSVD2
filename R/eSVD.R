#' Run the full eSVD-DE pipeline on a Seurat object
#'
#' The end-to-end wrapper: \code{format_covariates}, \code{initialize_esvd},
#' two rounds of \code{opt_esvd} (the case-control coefficient held fixed,
#' then freed) each followed by \code{reparameterization_esvd_covariates},
#' \code{estimate_nuisance}, and then either the fused
#' \code{compute_test_per_gene} (\code{bool_diet = TRUE}) or
#' \code{compute_posterior}, \code{compute_test_statistic} and
#' \code{compute_pvalue} in sequence.
#'
#' \strong{This function refuses inputs the model cannot handle}, rather than
#' filtering them: a gene whose counts are all zero, an individual with fewer
#' than \code{min_cells_per_individual} cells, or an arm with fewer than two
#' individuals each produce an error naming the offending genes or
#' individuals. The recommended entry point, \code{eSVD_helper}, applies the
#' cohort filters and removes all-zero genes before calling this function,
#' and reinserts the removed genes afterwards with a \code{gene_status}
#' record. It does \emph{not} refuse an underpowered cohort (few cells or few
#' individuals overall); those minima are \code{eSVD_helper}'s policy.
#'
#' The model is a Gamma-Poisson hierarchy: \eqn{A_{ji} \mid \lambda_{ji} \sim
#' \mathrm{Poisson}(\ell_{ji} \lambda_{ji})} with \eqn{\lambda_{ji} \sim
#' \mathrm{Gamma}(\mathrm{mean} = \mu_{ji}, \mathrm{var} = \gamma_j \mu_{ji})},
#' where \eqn{\mu_{ji}} is a low-rank function of the cell and gene latent
#' vectors plus the case-control covariate, and \eqn{\ell_{ji}} is the
#' covariate-adjusted sequencing depth. Note that the package stores the Gamma
#' \emph{rate} \eqn{\beta_j = 1/\gamma_j} as \code{nuisance_vec}, so larger
#' values mean \emph{less} over-dispersion. Differential expression is a
#' difference in mean expression between the case and control individuals,
#' tested with a Welch statistic on the per-individual posterior means. See
#' Lin, Qiu and Roeder (2024) \doi{10.1186/s12859-024-05724-7}.
#'
#' @param batch_var_prefix  \code{NULL}, or a regular expression matched against
#'                          the names of \code{categorical_vars}; the matching
#'                          covariates are treated as batch variables and their
#'                          coefficients are not reparameterized.
#' @param case_control_levels  Character vector of length 2 giving the two
#'                          levels of \code{case_control_var}, \strong{control
#'                          first, then case}. Reversing them flips the sign of
#'                          every test statistic.
#' @param case_control_var  Name of the column in \code{seurat_obj@meta.data}
#'                          holding the case-control status.
#' @param categorical_vars  Character vector of column names in
#'                          \code{seurat_obj@meta.data} to adjust for as factors
#'                          (one-hot encoded, first level dropped), or \code{NULL}.
#'                          A variable with a single level is dropped.
#' @param id_var            Name of the column in \code{seurat_obj@meta.data}
#'                          naming each cell's individual (donor).
#' @param numerical_vars    Character vector of column names in
#'                          \code{seurat_obj@meta.data} to adjust for as
#'                          numerical covariates, or \code{NULL}. They are rescaled
#'                          by \code{format_covariates}.
#' @param seurat_obj        A \code{Seurat} object whose \code{"RNA"} assay has a
#'                          \code{"counts"} layer (genes by cells) and whose
#'                          \code{meta.data} holds the variables above.
#' @param alpha_max         Maximum value of the prior mean in the posterior;
#'                          \code{NULL} (the default) uses twice the largest count.
#' @param bool_adjust_covariates  See \code{compute_posterior}.
#' @param bool_covariates_as_library  See \code{estimate_nuisance} and
#'                          \code{compute_posterior}.
#' @param bool_diet         If \code{TRUE} (the default), compute the posterior,
#'                          test statistics and p-values one gene at a time
#'                          without allocating cell-by-gene posterior matrices,
#'                          and drop \code{dat}, \code{covariates} and the
#'                          intermediate fits \code{fit_Init} and
#'                          \code{fit_First} from the returned object. The final
#'                          fit (\code{fit_Second}: \code{x_mat}, \code{y_mat},
#'                          \code{z_mat}, \code{nuisance_vec}) is kept either way.
#' @param bool_intercept    Whether the per-gene GLMs include an intercept.
#' @param bool_library_includes_interept  See \code{estimate_nuisance}.
#' @param bool_stabilize_underdispersion  See \code{compute_posterior}.
#' @param bool_use_log      See \code{estimate_nuisance}.
#' @param intermediate_save \code{NULL}, or a file path at which the object is
#'                          saved after each major stage.
#' @param k                 Number of latent dimensions; must not exceed the
#'                          number of genes.
#' @param l2pen             L2 penalty for \code{opt_esvd}.
#' @param lambda            Ridge penalty for \code{initialize_esvd}.
#' @param library_min       See \code{compute_posterior}.
#' @param max_iter          Maximum iterations for each \code{opt_esvd} call.
#' @param min_cells_per_individual  Minimum number of cells every individual
#'                          must contribute; an individual with fewer is an
#'                          error. \code{0} disables the check. See
#'                          \code{compute_test_statistic}.
#' @param pseudocount       See \code{compute_posterior}.
#' @param tol               Convergence tolerance for \code{opt_esvd}.
#' @param verbose           Integer; \code{0} is silent.
#'
#' @returns An \code{eSVD} object with elements \code{teststat_vec},
#' \code{case_mean}, \code{control_mean}, \code{pvalue_list} (see
#' \code{compute_pvalue}), \code{param}, \code{case_control},
#' \code{individual}, \code{latest_Fit} and the fit it names (an
#' \code{eSVD_Fit} with \code{x_mat}, \code{y_mat}, \code{z_mat},
#' \code{nuisance_vec}), plus \code{dat}, \code{covariates}, the earlier fits
#' and the posterior matrices when \code{bool_diet = FALSE}.
#' \code{report_results} turns it into a data frame.
#' @export
eSVD <- function(batch_var_prefix, # a variable inside categorical_vars. Can be NULL
                 case_control_levels, # Control and then Case
                 case_control_var,
                 categorical_vars,
                 id_var,
                 numerical_vars,
                 seurat_obj,
                 alpha_max = NULL,
                 bool_adjust_covariates = FALSE,
                 bool_covariates_as_library = TRUE,
                 bool_diet = TRUE,
                 bool_intercept = TRUE,
                 bool_library_includes_interept = TRUE,
                 bool_stabilize_underdispersion = TRUE,
                 bool_use_log = FALSE,
                 intermediate_save = NULL, # NULL or filepath to save intermediary results
                 k = 30,
                 l2pen = 0.1,
                 lambda = 0.1,
                 library_min = 0.1,
                 max_iter = 100,
                 min_cells_per_individual = 3,
                 pseudocount = 0,
                 tol = 1e-6,
                 verbose = 0){
  # SeuratObject is in Suggests: CRAN requires conditional use.
  if(!requireNamespace("SeuratObject", quietly = TRUE)){
    stop("package `SeuratObject` is required by `eSVD()` to read the count ",
         "matrix from `seurat_obj`; install it with ",
         "install.packages(\"SeuratObject\")")
  }
  stopifnot(inherits(seurat_obj, "Seurat"))
  stopifnot(all(is.null(categorical_vars)) || (length(unique(categorical_vars)) == length(categorical_vars) && all(is.character(categorical_vars))))
  stopifnot(all(is.null(numerical_vars)) || (length(unique(numerical_vars)) == length(numerical_vars) && all(is.character(numerical_vars))))
  stopifnot(length(case_control_var) == 1,
            is.character(case_control_var),
            length(id_var) == 1,
            is.character(id_var),
            length(case_control_levels) == 2,
            all(is.character(case_control_levels)))
  stopifnot(length(k) == 1, k > 0, k %% 1 == 0,
            length(min_cells_per_individual) == 1, min_cells_per_individual >= 0)

  # make sure there's an appropriate batch variable
  if(!is.null(batch_var_prefix) &&
     length(grep(batch_var_prefix, categorical_vars)) == 0){
    stop("`batch_var_prefix` = \"", batch_var_prefix,
         "\" matches none of `categorical_vars`")
  }

  missing_vars <- setdiff(c(case_control_var, id_var, categorical_vars, numerical_vars),
                          colnames(seurat_obj@meta.data))
  if(length(missing_vars) > 0){
    stop("variable(s) `", paste0(missing_vars, collapse = "`, `"),
         "` are not columns of `seurat_obj@meta.data`")
  }
  cc_raw_vec <- as.character(seurat_obj@meta.data[,case_control_var])
  if(!all(cc_raw_vec %in% case_control_levels)){
    stop("`case_control_var` = \"", case_control_var, "\" takes value(s) \"",
         paste0(setdiff(unique(cc_raw_vec), case_control_levels), collapse = "\", \""),
         "\" outside `case_control_levels`; subset `seurat_obj` to the two ",
         "levels first")
  }

  # extract count matrix
  mat <- Matrix::t(SeuratObject::LayerData(seurat_obj,
                                           assay = "RNA",
                                           layer = "counts"))

  # The three refusals (Q-COH-7): the wrapper `eSVD_helper` filters, this
  # function refuses, so a helper bug fails loudly here rather than as an
  # obscure numerical failure further down.
  all_zero_idx <- .which_all_zero(mat)
  if(length(all_zero_idx) > 0){
    stop(length(all_zero_idx), " gene(s) are all zero (",
         paste0(utils::head(colnames(mat)[all_zero_idx], 5), collapse = ", "),
         if(length(all_zero_idx) > 5) ", ..." else "",
         "); `eSVD()` does not filter genes. Use `eSVD_helper()`, which removes ",
         "them and records the removal in `gene_status`, or remove them yourself")
  }
  if(k > ncol(mat)){
    stop("`k` = ", k, " exceeds the number of genes to analyze, ", ncol(mat),
         "; a `k` that was legal before filtering all-zero genes can become ",
         "illegal after it")
  }
  individual_vec <- seurat_obj@meta.data[,id_var]
  case_individuals <- unique(individual_vec[cc_raw_vec == case_control_levels[2]])
  control_individuals <- unique(individual_vec[cc_raw_vec == case_control_levels[1]])
  if(length(intersect(as.character(case_individuals),
                      as.character(control_individuals))) > 0){
    stop("individual(s) `",
         paste0(intersect(as.character(case_individuals),
                          as.character(control_individuals)), collapse = "`, `"),
         "` appear in both the case and the control arm")
  }
  .check_cohort_is_testable(case_individuals = case_individuals,
                            control_individuals = control_individuals,
                            individual_vec = individual_vec,
                            min_cells_per_individual = min_cells_per_individual)

  if(verbose > 0) print("Processing the covariates")
  if(is.factor(seurat_obj@meta.data[,id_var])) seurat_obj@meta.data[,id_var] <- droplevels(seurat_obj@meta.data[,id_var])
  for(variable in categorical_vars){
    if(is.factor(seurat_obj@meta.data[,variable])) seurat_obj@meta.data[,variable] <- droplevels(seurat_obj@meta.data[,variable])
  }

  if(length(categorical_vars) >= 1){
    categorical_vars_subset <- categorical_vars[sapply(categorical_vars, function(x){
      length(levels(droplevels(seurat_obj@meta.data[,x]))) > 1
    })]
  } else {
    categorical_vars_subset <- NULL
  }
  covariate_dat <- seurat_obj@meta.data[,c(id_var, categorical_vars_subset, numerical_vars), drop = FALSE]
  covariate_df <- data.frame(covariate_dat)

  covariate_df[,case_control_var] <- factor(seurat_obj@meta.data[,case_control_var],
                                            levels = case_control_levels)
  for(variable in setdiff(c(id_var, categorical_vars_subset), case_control_var)){
    covariate_df[,variable] <- factor(covariate_df[,variable],
                                      levels = names(sort(table(covariate_df[,variable]),
                                                          decreasing = TRUE)))
  }

  covariates <- format_covariates(dat = mat,
                                  covariate_df = covariate_df,
                                  rescale_numeric_variables = numerical_vars)

  case_control_variable <- paste0(case_control_var, "_", case_control_levels[2])

  ############

  # The individual indicators are removed from the design by their exact
  # names, `<id_var>_<level>`. A regex `grep(id_var, ...)` here used to drop
  # every column containing `id_var` as a substring, e.g. `donor_age` for
  # `id_var = "donor"`, silently un-adjusting the model for it.
  individual_columns <- paste0(id_var, "_", levels(covariate_df[,id_var]))
  keep_columns <- which(!colnames(covariates) %in% individual_columns)

  if(verbose > 0) print("Initialization")
  eSVD_obj <- initialize_esvd(dat = mat,
                              covariates = covariates[, keep_columns, drop = FALSE],
                              case_control_variable = case_control_variable,
                              bool_intercept = bool_intercept,
                              k = k,
                              lambda = lambda,
                              metadata_case_control = covariates[,case_control_variable],
                              metadata_individual = covariate_df[,id_var],
                              verbose = verbose - 1)

  if(!is.null(batch_var_prefix)){
    omitted_variables <- colnames(eSVD_obj$covariates)[grep(batch_var_prefix, colnames(eSVD_obj$covariates))]
  } else {
    omitted_variables <- NULL
  }

  if(!is.null(intermediate_save)){
    save(eSVD_obj,
         file = intermediate_save)
  }

  if(verbose > 0) print("Doing the first reparameterization")
  eSVD_obj <- reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_Init",
    omitted_variables = c("Log_UMI", omitted_variables)
  )

  if(verbose > 0)  print("First fit")
  eSVD_obj <- opt_esvd(input_obj = eSVD_obj,
                       l2pen = l2pen,
                       max_iter = max_iter,
                       offset_variables = setdiff(colnames(eSVD_obj$covariates), case_control_variable),
                       tol = tol,
                       verbose = verbose - 1,
                       fit_name = "fit_First",
                       fit_previous = "fit_Init")

  if(!is.null(intermediate_save)){
    save(eSVD_obj,
         file = intermediate_save)
  }

  eSVD_obj <- reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_First",
    omitted_variables = c("Log_UMI", omitted_variables)
  )

  if(bool_diet){
    eSVD_obj[["fit_Init"]] <- NULL
  }

  if(verbose > 0) print("Second fit")
  eSVD_obj <- opt_esvd(input_obj = eSVD_obj,
                       l2pen = l2pen,
                       max_iter = max_iter,
                       offset_variables = NULL,
                       tol = tol,
                       verbose = verbose - 1,
                       fit_name = "fit_Second",
                       fit_previous = "fit_First")

  if(!is.null(intermediate_save)){
    save(eSVD_obj,
         file = intermediate_save)
  }

  eSVD_obj <- reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_Second",
    omitted_variables = omitted_variables
  )

  if(bool_diet){
    eSVD_obj[["fit_First"]] <- NULL
  }

  if(verbose > 0) print("Nuisance estimation")
  eSVD_obj <- estimate_nuisance(input_obj = eSVD_obj,
                                bool_covariates_as_library = bool_covariates_as_library,
                                bool_library_includes_interept = bool_library_includes_interept,
                                bool_use_log = bool_use_log,
                                verbose = verbose - 1)

  if(!is.null(intermediate_save)){
    save(eSVD_obj,
         file = intermediate_save)
  }

  if(is.null(alpha_max)){
    alpha_max <- 2*max(eSVD_obj$dat@x)
    stopifnot(alpha_max > 0)
  }

  if(bool_diet){
    if(verbose > 0) print("Computing posterior and p-values, one gene at a time")
    eSVD_obj <- compute_test_per_gene(input_obj = eSVD_obj,
                                      alpha_max = alpha_max,
                                      bool_adjust_covariates = bool_adjust_covariates,
                                      bool_covariates_as_library = bool_covariates_as_library,
                                      bool_stabilize_underdispersion = bool_stabilize_underdispersion,
                                      library_min = library_min,
                                      min_cells_per_individual = min_cells_per_individual,
                                      pseudocount = pseudocount,
                                      verbose = verbose - 1)

  } else {
    if(verbose > 0) print("Computing posterior")
    eSVD_obj <- compute_posterior(input_obj = eSVD_obj,
                                  bool_adjust_covariates = bool_adjust_covariates,
                                  alpha_max = alpha_max,
                                  bool_covariates_as_library = bool_covariates_as_library,
                                  bool_stabilize_underdispersion = bool_stabilize_underdispersion,
                                  library_min = library_min,
                                  pseudocount = pseudocount)

    if(!is.null(intermediate_save)){
      save(eSVD_obj,
           file = intermediate_save)
    }

    if(verbose > 0) print("Computing p-values")
    eSVD_obj <- compute_test_statistic(input_obj = eSVD_obj,
                                       min_cells_per_individual = min_cells_per_individual,
                                       verbose = verbose - 1)
    eSVD_obj <- compute_pvalue(input_obj = eSVD_obj,
                               min_cells_per_individual = min_cells_per_individual)
  }

  if(verbose > 0) print("Finalizing")
  if(bool_diet){
    # The final fit is kept: it is small (n x k and p x (k + r)) and it is the
    # model. Only the data copies and the superseded fits go.
    eSVD_obj$covariates <- NULL
    eSVD_obj$fit_Init <- NULL
    eSVD_obj$fit_First <- NULL
    eSVD_obj$dat <- NULL
  }

  eSVD_obj
}
