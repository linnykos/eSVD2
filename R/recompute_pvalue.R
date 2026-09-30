#' Redo the test at another cap on the nuisance rate, without fitting again
#'
#' The factorization (\code{x_mat}, \code{y_mat}, \code{z_mat}) is the
#' expensive part of an analysis and does not depend on the cap on the
#' nuisance rate; see \code{estimate_nuisance}. This function takes a fitted
#' and tested \code{eSVD} object, applies another \code{cap_multiplier} to the
#' maximum-likelihood rates stored on it, and repeats what comes after the
#' rates: the posterior, the test statistic, the log2 fold change and the
#' p-values. Nothing is estimated again, and every other setting is the one
#' recorded in \code{input_obj$param} when the object was built.
#'
#' The posterior needs the counts and the covariates. An object built with
#' \code{eSVD(bool_diet = FALSE)}, or stage by stage, carries them. An object
#' built with \code{bool_diet = TRUE} (the default of \code{eSVD} and
#' \code{eSVD_helper}) does not, so pass \code{seurat_obj}, the Seurat object
#' the analysis was run on. The cells and genes of the fit are looked up in it
#' \strong{by name}, so it may hold more cells and genes than were analyzed
#' (the individuals and the all-zero genes \code{eSVD_helper} removed), and
#' the covariates are rebuilt from its \code{meta.data} with the variable
#' names \code{eSVD} recorded. Summaries of the rebuilt counts and
#' covariates are compared with those recorded at the fit, and the function
#' stops if they differ, for example because the \code{meta.data} was edited
#' since.
#'
#' The result is what \code{eSVD} would have returned had it been called with
#' this \code{cap_multiplier}.
#'
#' @param input_obj       \code{eSVD} object that has been tested, by
#'                        \code{eSVD}, \code{eSVD_helper},
#'                        \code{compute_test_per_gene} or \code{compute_pvalue},
#'                        with version 1.2.0 or later.
#' @param cap_multiplier  One positive number: no gene's nuisance rate may
#'                        exceed this many times the median of its library
#'                        size. \code{Inf} means no cap. See
#'                        \code{estimate_nuisance}.
#' @param seurat_obj      \code{NULL}, or the \code{Seurat} object the analysis
#'                        was run on. Required when \code{input_obj} has no
#'                        counts (\code{input_obj$dat}) or covariates, and
#'                        ignored otherwise.
#' @param verbose         Integer; \code{0} is silent.
#'
#' @returns The \code{eSVD} object with \code{teststat_vec}, \code{case_mean},
#' \code{control_mean}, \code{case_var}, \code{control_var},
#' \code{log2fc_vec}, \code{log2fc_se_vec} and \code{pvalue_list} replaced,
#' and, in the fit named by \code{latest_Fit}, \code{nuisance_vec},
#' \code{nuisance_status} and the posterior matrices (if the object had them)
#' replaced. \code{param$nuisance_cap_multiplier},
#' \code{param$nuisance_num_capped} and \code{param$nuisance_num_boundary}
#' describe the new cap. The factorization and \code{nuisance_mle_vec} are
#' unchanged, so the function can be called again at any other cap. An object
#' that had no counts is returned without them.
#' @export
recompute_pvalue <- function(input_obj,
                             cap_multiplier,
                             seurat_obj = NULL,
                             verbose = 0){
  stopifnot(inherits(input_obj, "eSVD"), "latest_Fit" %in% names(input_obj),
            input_obj[["latest_Fit"]] %in% names(input_obj),
            inherits(input_obj[[input_obj[["latest_Fit"]]]], "eSVD_Fit"))
  .check_cap_multiplier(cap_multiplier)
  latest_Fit <- input_obj[["latest_Fit"]]

  missing_elements <- setdiff(c("nuisance_mle_vec",
                                "nuisance_library_median_vec",
                                "nuisance_status"),
                              names(input_obj[[latest_Fit]]))
  if(length(missing_elements) > 0){
    stop("`input_obj$", latest_Fit, "` has no `",
         paste0(missing_elements, collapse = "`, `"),
         "`, which hold the nuisance rates before the cap. The object was ",
         "built before eSVD2 1.2.0; run `estimate_nuisance` on it again ",
         "(which needs its counts), or `eSVD` on the data")
  }
  .check_fit_is_named(input_obj)
  if(is.null(input_obj[["teststat_vec"]]) || is.null(input_obj[["pvalue_list"]])){
    stop("`input_obj` has no `teststat_vec` or no `pvalue_list`, so there is ",
         "no test to redo; run `compute_test_per_gene`, or ",
         "`compute_posterior`, `compute_test_statistic` and `compute_pvalue`")
  }

  if(verbose > 0) print("Gathering the counts and covariates")
  input_list <- .get_counts_and_covariates(input_obj = input_obj,
                                           seurat_obj = seurat_obj)

  # What the object looked like, to return it the same way.
  bool_has_dat <- !is.null(input_obj[["dat"]])
  bool_has_covariates <- !is.null(input_obj[["covariates"]])
  bool_has_posterior <- !is.null(input_obj[[latest_Fit]]$posterior_mean_mat)
  gene_status <- input_obj[["gene_status"]]
  dat_original <- input_obj[["dat"]]

  # Work on the analyzed genes only: the genes `eSVD_helper` reinserted are
  # NA in the fit, and the empirical null must not see them.
  eSVD_obj <- .remove_reinserted_genes(eSVD_obj = input_obj,
                                       gene_vec = colnames(input_list$dat))
  eSVD_obj[["dat"]] <- input_list$dat
  eSVD_obj[["covariates"]] <- input_list$covariates

  if(verbose > 0) print("Applying the cap")
  fit <- eSVD_obj[[latest_Fit]]
  min_val <- .get_param_or_default(eSVD_obj, "nuisance_min_val", 1e-4)
  .check_min_val(min_val, cap_multiplier)
  cap_res <- .apply_nuisance_cap(
    nuisance_mle_vec = fit$nuisance_mle_vec,
    library_median_vec = fit$nuisance_library_median_vec,
    bool_boundary_vec = as.character(fit$nuisance_status) == "boundary",
    bool_failed_vec = as.character(fit$nuisance_status) == "failed",
    cap_multiplier = cap_multiplier,
    min_val = min_val
  )
  eSVD_obj[[latest_Fit]]$nuisance_vec <- cap_res$nuisance_vec
  eSVD_obj[[latest_Fit]]$nuisance_status <- cap_res$nuisance_status
  eSVD_obj$param[c("nuisance_cap_multiplier",
                   "nuisance_num_boundary",
                   "nuisance_num_capped")] <- list(cap_multiplier,
                                                   cap_res$num_boundary,
                                                   cap_res$num_capped)

  posterior_args <- .get_posterior_args(eSVD_obj)
  min_cells_per_individual <- .get_param_or_default(
    eSVD_obj, "test_min_cells_per_individual", 3
  )

  if(bool_has_posterior){
    if(verbose > 0) print("Computing posterior")
    eSVD_obj <- do.call(compute_posterior,
                        c(list(input_obj = eSVD_obj), posterior_args))

    if(verbose > 0) print("Computing p-values")
    eSVD_obj <- compute_test_statistic(
      input_obj = eSVD_obj,
      min_cells_per_individual = min_cells_per_individual,
      verbose = verbose - 1
    )
    eSVD_obj <- compute_pvalue(
      input_obj = eSVD_obj,
      min_cells_per_individual = min_cells_per_individual
    )
  } else {
    if(verbose > 0) print("Computing posterior and p-values, one gene at a time")
    posterior_args$bool_return_components <- NULL
    eSVD_obj <- do.call(compute_test_per_gene,
                        c(list(input_obj = eSVD_obj,
                               min_cells_per_individual = min_cells_per_individual,
                               verbose = verbose - 1),
                          posterior_args))
  }

  if(verbose > 0) print("Finalizing")
  if(!is.null(gene_status)){
    eSVD_obj <- .reinsert_genes(eSVD_obj = eSVD_obj,
                                gene_status = gene_status,
                                count_mat = dat_original)
    eSVD_obj[["gene_status"]] <- gene_status
  }
  if(bool_has_dat){
    eSVD_obj[["dat"]] <- dat_original
  } else {
    eSVD_obj[["dat"]] <- NULL
  }
  if(!bool_has_covariates) eSVD_obj[["covariates"]] <- NULL

  # Removing an element and adding it back moves it to the end; return the
  # object in the layout it came in.
  name_vec <- c(intersect(names(input_obj), names(eSVD_obj)),
                setdiff(names(eSVD_obj), names(input_obj)))
  structure(unclass(eSVD_obj)[name_vec], class = class(input_obj))
}

#' The counts and covariates an eSVD object was fitted with
#'
#' Taken from the object when it carries them, and rebuilt from the Seurat
#' object otherwise. Either way they are restricted to the genes that were
#' analyzed, which for an object of \code{eSVD_helper} excludes the all-zero
#' genes it reinserted.
#'
#' A rebuilt pair is checked against records of the fit before it is
#' returned: each gene's mean count (\code{gene_mean_count_vec}, stored by
#' \code{estimate_nuisance}), each covariate's column sum, and the weighted
#' column sums of the counts and of the covariates (the three stored in
#' \code{param} by \code{eSVD}; see \code{.compute_weighted_colsum} for
#' why the plain sums are not enough). The tolerance is relative,
#' \code{1e-8}. These are summaries: they notice counts or metadata that
#' were edited, reordered or exchanged between cells, and they are not a
#' proof that the matrices are equal.
#'
#' @param input_obj   \code{eSVD} object.
#' @param seurat_obj  \code{NULL} or a \code{Seurat} object.
#'
#' @returns List with \code{covariates} (numeric matrix, rows are cells) and
#' \code{dat} (counts, rows are cells in the order of \code{covariates},
#' columns are the analyzed genes in the order of the fit).
#' @noRd
.get_counts_and_covariates <- function(input_obj, seurat_obj){
  latest_Fit <- input_obj[["latest_Fit"]]
  fit <- input_obj[[latest_Fit]]
  cell_vec <- rownames(fit$x_mat)
  gene_vec <- .get_analyzed_genes(input_obj)

  if(!is.null(input_obj[["dat"]]) && !is.null(input_obj[["covariates"]])){
    return(list(covariates = input_obj[["covariates"]],
                dat = input_obj[["dat"]][, gene_vec, drop = FALSE]))
  }

  if(is.null(seurat_obj)){
    stop("`input_obj` has no counts or no covariates (it was built with ",
         "`bool_diet = TRUE`), so `seurat_obj` is needed: pass the Seurat ",
         "object the analysis was run on")
  }
  if(!inherits(seurat_obj, "Seurat")){
    stop("`seurat_obj` must be a Seurat object; received an object of class `",
         paste0(class(seurat_obj), collapse = "`, `"), "`")
  }
  if(!requireNamespace("SeuratObject", quietly = TRUE)){
    stop("package `SeuratObject` is required to read the count matrix from ",
         "`seurat_obj`; install it with install.packages(\"SeuratObject\")")
  }
  esvd_names <- c("esvd_case_control_levels", "esvd_case_control_var",
                  "esvd_categorical_vars", "esvd_count_weighted_colsum_vec",
                  "esvd_covariate_colsum_vec",
                  "esvd_covariate_weighted_colsum_vec",
                  "esvd_id_var", "esvd_numerical_vars")
  if(!all(esvd_names %in% names(input_obj$param))){
    stop("`input_obj$param` does not record the variables `eSVD` was called ",
         "with (`", paste0(setdiff(esvd_names, names(input_obj$param)),
                           collapse = "`, `"),
         "`), so its covariates cannot be rebuilt from `seurat_obj`. The ",
         "object was not built by `eSVD` or `eSVD_helper` of eSVD2 1.2.0 or ",
         "later")
  }
  # Subset BEFORE building the covariates: `Log_UMI` and the rescaled
  # numerical variables depend on which cells and genes are present. The
  # function stops if a cell or a gene of the fit is not in `seurat_obj`.
  mat <- .extract_count_matrix(seurat_obj,
                               cell_vec = cell_vec,
                               gene_vec = gene_vec)
  metadata_df <- seurat_obj@meta.data[cell_vec, , drop = FALSE]

  param <- input_obj$param
  covariate_list <- .prepare_esvd_covariates(
    mat = mat,
    metadata_df = metadata_df,
    case_control_levels = param[["esvd_case_control_levels"]],
    case_control_var = param[["esvd_case_control_var"]],
    categorical_vars = param[["esvd_categorical_vars"]],
    id_var = param[["esvd_id_var"]],
    numerical_vars = param[["esvd_numerical_vars"]]
  )
  covariates <- covariate_list$covariates

  # The mean notices a count that changed, and the weighted sum one that
  # moved to another cell.
  bool_equal_vec <- .is_relatively_equal(Matrix::colMeans(mat),
                                         fit$gene_mean_count_vec[gene_vec],
                                         bool_all = FALSE) &
    .is_relatively_equal(.compute_weighted_colsum(mat),
                         param[["esvd_count_weighted_colsum_vec"]][gene_vec],
                         bool_all = FALSE)
  if(!all(bool_equal_vec)){
    stop("the counts in `seurat_obj` are not the counts `input_obj` was ",
         "fitted with: the counts of ", sum(!bool_equal_vec),
         " gene(s) differ from what was recorded at the fit")
  }
  if(!identical(colnames(covariates), colnames(fit$z_mat))){
    stop("the covariates rebuilt from `seurat_obj` (",
         paste0(colnames(covariates), collapse = ", "),
         ") are not the covariates of the fit (",
         paste0(colnames(fit$z_mat), collapse = ", "), ")")
  }
  covariate_vec <- colnames(covariates)
  bool_equal_vec <- .is_relatively_equal(
    Matrix::colSums(covariates),
    param[["esvd_covariate_colsum_vec"]][covariate_vec],
    bool_all = FALSE
  ) & .is_relatively_equal(
    .compute_weighted_colsum(covariates),
    param[["esvd_covariate_weighted_colsum_vec"]][covariate_vec],
    bool_all = FALSE
  )
  if(!all(bool_equal_vec)){
    stop("the covariates rebuilt from `seurat_obj` are not the covariates ",
         "`input_obj` was fitted with: `",
         paste0(colnames(covariates)[!bool_equal_vec], collapse = "`, `"),
         "` differ from what was recorded at the fit. Has the `meta.data` ",
         "changed?")
  }

  list(covariates = covariates,
       dat = mat)
}

#' Stop if the fit of an eSVD object does not name its cells and genes
#'
#' \code{initialize_esvd} accepts a count matrix without names, and
#' everything that works on a fitted object matches cells and genes by name.
#' \code{dat[, NULL]} selects no column and raises nothing, so the absence is
#' refused here.
#'
#' @param input_obj   \code{eSVD} object.
#' @param bool_cells  Whether the cell names are needed as well.
#'
#' @returns \code{invisible()}; called for its error.
#' @noRd
.check_fit_is_named <- function(input_obj, bool_cells = TRUE){
  latest_Fit <- input_obj[["latest_Fit"]]
  fit <- input_obj[[latest_Fit]]

  if(is.null(rownames(fit$y_mat))){
    stop("`input_obj$", latest_Fit, "$y_mat` has no gene names (row names), ",
         "and the genes are matched by name; the count matrix given to ",
         "`initialize_esvd` had no column names")
  }
  if(bool_cells && is.null(rownames(fit$x_mat))){
    stop("`input_obj$", latest_Fit, "$x_mat` has no cell names (row names), ",
         "and the cells are matched by name; the count matrix given to ",
         "`initialize_esvd` had no row names")
  }

  invisible()
}

# A sum over cells does not change when values are exchanged between cells:
# the `Age` of two individuals with as many cells each, or the counts of two
# cells. With a weight that differs from cell to cell it does. The weights
# are fixed and depend on the position of the cell alone.
#
#' Weighted column sums of a matrix whose rows are cells
#'
#' @param mat  \code{matrix} or \code{dgCMatrix}, rows are cells.
#'
#' @returns Numeric vector with one entry per column, named by
#' \code{colnames(mat)}: \code{sum_i cos(sqrt(2) * i) * mat[i, j]}.
#' @noRd
.compute_weighted_colsum <- function(mat){
  weight_vec <- cos(sqrt(2) * seq_len(nrow(mat)))

  res <- as.numeric(Matrix::crossprod(mat, weight_vec))
  names(res) <- colnames(mat)
  res
}

# Relative comparison of two numeric vectors, element by element.
.is_relatively_equal <- function(vec1, vec2, bool_all = TRUE, tol = 1e-8){
  stopifnot(length(vec1) == length(vec2))
  vec1 <- as.numeric(vec1)
  vec2 <- as.numeric(vec2)

  bool_vec <- abs(vec1 - vec2) <= tol * pmax(abs(vec1), abs(vec2), 1)
  bool_vec[is.na(bool_vec)] <- FALSE

  if(bool_all) all(bool_vec) else bool_vec
}

#' The genes of an eSVD object that were analyzed
#'
#' @param eSVD_obj  \code{eSVD} object.
#'
#' @returns Character vector: the genes of the fit, without those
#' \code{eSVD_helper} reinserted (whose \code{gene_status} is not
#' \code{"analyzed"}).
#' @noRd
.get_analyzed_genes <- function(eSVD_obj){
  gene_vec <- rownames(eSVD_obj[[eSVD_obj[["latest_Fit"]]]]$y_mat)
  gene_status <- eSVD_obj[["gene_status"]]
  if(is.null(gene_status)) return(gene_vec)

  stopifnot(identical(names(gene_status), gene_vec))
  gene_vec[as.character(gene_status) == "analyzed"]
}

#' Restrict an eSVD object to the genes that were analyzed
#'
#' The inverse of \code{.reinsert_genes}: every per-gene element it pads is
#' subset here, \strong{by name}.
#'
#' @param eSVD_obj  \code{eSVD} object.
#' @param gene_vec  Character vector, the genes to keep.
#'
#' @returns The \code{eSVD} object without \code{gene_status}, every per-gene
#' element holding the genes of \code{gene_vec} in that order.
#' @noRd
.remove_reinserted_genes <- function(eSVD_obj, gene_vec){
  subset_vector <- function(vec){
    if(is.null(vec)) return(NULL)
    stopifnot(all(gene_vec %in% names(vec)))
    vec[gene_vec]
  }

  for(element_name in c("teststat_vec", "case_mean", "control_mean",
                        "case_var", "control_var", "log2fc_vec",
                        "log2fc_se_vec")){
    eSVD_obj[[element_name]] <- subset_vector(eSVD_obj[[element_name]])
  }

  pvalue_list <- eSVD_obj[["pvalue_list"]]
  for(element_name in c("df_vec", "gaussian_teststat", "log10pvalue",
                        "fdr_vec")){
    pvalue_list[[element_name]] <- subset_vector(pvalue_list[[element_name]])
  }
  eSVD_obj[["pvalue_list"]] <- pvalue_list

  fit_names <- names(eSVD_obj)[sapply(eSVD_obj, inherits, "eSVD_Fit")]
  for(fit_name in fit_names){
    fit <- eSVD_obj[[fit_name]]
    fit$y_mat <- fit$y_mat[gene_vec, , drop = FALSE]
    fit$z_mat <- fit$z_mat[gene_vec, , drop = FALSE]
    for(element_name in c(.nuisance_numeric_elements(), "nuisance_status")){
      fit[[element_name]] <- subset_vector(fit[[element_name]])
    }
    for(element_name in c("posterior_mean_mat", "posterior_var_mat")){
      if(!is.null(fit[[element_name]])){
        fit[[element_name]] <- fit[[element_name]][, gene_vec, drop = FALSE]
      }
    }
    eSVD_obj[[fit_name]] <- fit
  }

  if(!is.null(eSVD_obj[["dat"]])){
    eSVD_obj[["dat"]] <- eSVD_obj[["dat"]][, gene_vec, drop = FALSE]
  }
  eSVD_obj[["gene_status"]] <- NULL

  eSVD_obj
}

# An entry of `param`, or `default_val` when the object does not record it.
.get_param_or_default <- function(eSVD_obj, param_name, default_val){
  if(!param_name %in% names(eSVD_obj$param)) return(default_val)

  eSVD_obj$param[[param_name]]
}

#' The settings of the posterior recorded on an eSVD object
#'
#' @param eSVD_obj  \code{eSVD} object.
#'
#' @returns Named list of the arguments of \code{compute_posterior.eSVD},
#' each read from \code{param$posterior_<argument>}.
#' @noRd
.get_posterior_args <- function(eSVD_obj){
  arg_names <- c("alpha_max", "bool_adjust_covariates",
                 "bool_covariates_as_library", "bool_return_components",
                 "bool_stabilize_underdispersion", "library_min",
                 "nuisance_lower_quantile", "pseudocount")
  param_names <- paste0("posterior_", arg_names)
  if(!all(param_names %in% names(eSVD_obj$param))){
    stop("`input_obj$param` does not record the settings of the posterior (`",
         paste0(setdiff(param_names, names(eSVD_obj$param)), collapse = "`, `"),
         "`), so the test cannot be repeated as it was run. The object was ",
         "tested with a version before eSVD2 1.2.0")
  }

  # Single brackets: a setting that is `NULL` stays in the list.
  arg_list <- eSVD_obj$param[param_names]
  names(arg_list) <- arg_names
  arg_list
}
