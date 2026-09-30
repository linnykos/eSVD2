#' Filter a cohort before running eSVD-DE
#'
#' Applies the cohort-level filters that \code{eSVD_helper} runs before
#' \code{eSVD}, in this order:
#' \enumerate{
#'   \item \strong{Drop} every individual with fewer than
#'         \code{min_cells_per_id} cells, with a warning naming them.
#'   \item \strong{Check} the filtered cohort against \code{min_cells},
#'         \code{min_ids}, \code{min_ids_per_arm} and
#'         \code{min_cells_casecontrol}; on the first failure, warn and
#'         return \code{NA}.
#' }
#' Dropping first means every threshold describes the cohort that will
#' actually be analyzed.
#'
#' \code{min_cells}, \code{min_ids} and \code{min_cells_casecontrol} reject
#' at \code{<=}, and \code{min_cells_per_id} and \code{min_ids_per_arm} at
#' \code{<}, so setting one to \code{0} disables it. \code{min_cells_per_id} is the one exception: a value
#' of \code{1} or \code{2} would not disable the drop but weaken it, and a
#' surviving one-cell individual has no within-individual variance, so it
#' must be \code{0} (off) or at least \code{3}.
#'
#' The four rejections return the same bare \code{NA}, so the warning text
#' is what tells a caller which filter fired; each names its own threshold.
#'
#' @param seurat_obj        A \code{Seurat} object.
#' @param id_var            Column of \code{seurat_obj@meta.data} naming each
#'                          cell's individual.
#' @param case_control_var  Column of \code{seurat_obj@meta.data} holding the
#'                          case-control status (two levels).
#' @param min_cells         Reject a cohort with this many cells or fewer,
#'                          after the drop.
#' @param min_cells_casecontrol  Reject a cohort in which either arm has this
#'                          many cells or fewer.
#' @param min_cells_per_id  Drop every individual with fewer than this many
#'                          cells; \code{0} or at least \code{3}.
#' @param min_ids           Reject a cohort with this many individuals or fewer,
#'                          pooled across arms.
#' @param min_ids_per_arm   Reject a cohort in which either arm has fewer than
#'                          this many individuals. The default \code{2} is the
#'                          smallest number for which the Welch degrees of
#'                          freedom are defined.
#' @param verbose           Integer; \code{0} is silent.
#'
#' @returns Either \code{NA} (the cohort was rejected; see the warning), or a
#' list with \code{dropped_individuals} (character vector, possibly empty) and
#' \code{seurat_obj} (the object with those individuals' cells removed).
#' @examples
#' if(requireNamespace("SeuratObject", quietly = TRUE)){
#'   set.seed(10)
#'   sim <- generate_null(cell_per_person = 15, num_genes = 40,
#'                        num_individuals = 8)
#'   meta_df <- data.frame(sim$covariates[, c("CC", "Sex", "Age")])
#'   meta_df$Sex <- factor(meta_df$Sex)
#'   meta_df$Individual <- sim$metadata_individual
#'   rownames(meta_df) <- rownames(sim$obs_mat)
#'   # Seurat rewrites "gene_1" as "gene-1" (with a warning); rename first
#'   count_mat <- Matrix::t(sim$obs_mat)
#'   rownames(count_mat) <- gsub("_", "", rownames(count_mat))
#'   seurat_obj <- SeuratObject::CreateSeuratObject(counts = count_mat,
#'                                                  meta.data = meta_df)
#'
#'   # drop two cells from one individual: with min_cells_per_id = 15 that
#'   # individual is removed, with a warning naming it
#'   seurat_small <- seurat_obj[, -(1:2)]
#'   res <- filter_cohort(seurat_obj = seurat_small,
#'                        id_var = "Individual",
#'                        case_control_var = "CC",
#'                        min_cells_per_id = 15)
#'   res$dropped_individuals
#'   ncol(res$seurat_obj)
#' }
#' @export
filter_cohort <- function(seurat_obj,
                          id_var,
                          case_control_var,
                          min_cells = 20,
                          min_cells_casecontrol = 20,
                          min_cells_per_id = 3,
                          min_ids = 4,
                          min_ids_per_arm = 2,
                          verbose = 0){
  if(!requireNamespace("SeuratObject", quietly = TRUE)){
    stop("package `SeuratObject` is required by `filter_cohort()`; ",
         "install it with install.packages(\"SeuratObject\")")
  }
  stopifnot(inherits(seurat_obj, "Seurat"),
            length(id_var) == 1, is.character(id_var),
            length(case_control_var) == 1, is.character(case_control_var),
            all(c(id_var, case_control_var) %in% colnames(seurat_obj@meta.data)),
            min_cells >= 0, min_cells_casecontrol >= 0, min_ids >= 0,
            min_ids_per_arm >= 0)
  if(!(min_cells_per_id == 0 || min_cells_per_id >= 3)){
    stop("`min_cells_per_id` = ", min_cells_per_id, " must be 0 (disable the ",
         "drop) or at least 3. A value of 1 or 2 weakens the drop instead of ",
         "disabling it, and a surviving individual with one cell has no ",
         "within-individual variance")
  }

  # Step 1: drop individuals with too few cells.
  id_vec <- as.character(seurat_obj@meta.data[,id_var])
  indiv_count <- table(id_vec)
  dropped_individuals <- names(indiv_count)[indiv_count < min_cells_per_id]
  if(length(dropped_individuals) > 0){
    warning("dropping ", length(dropped_individuals), " individual(s) with ",
            "fewer than `min_cells_per_id` = ", min_cells_per_id, " cells: ",
            paste0(dropped_individuals, " (",
                   as.integer(indiv_count[dropped_individuals]), " cells)",
                   collapse = ", "))

    keep_vec <- !id_vec %in% dropped_individuals
    num_cells <- sum(keep_vec)
    # Seurat refuses to subset to zero cells; that cohort is rejected below.
    if(num_cells > 0){
      seurat_obj <- seurat_obj[, SeuratObject::Cells(seurat_obj)[keep_vec]]
      if(is.factor(seurat_obj@meta.data[,id_var])){
        seurat_obj@meta.data[,id_var] <- droplevels(seurat_obj@meta.data[,id_var])
      }
      if(is.factor(seurat_obj@meta.data[,case_control_var])){
        seurat_obj@meta.data[,case_control_var] <- droplevels(seurat_obj@meta.data[,case_control_var])
      }
    }
  } else {
    keep_vec <- rep(TRUE, length(id_vec))
    num_cells <- length(id_vec)
  }
  if(verbose > 0) print(paste0("Cohort has ", num_cells, " cells after the drop"))

  # Step 2: check the filtered cohort.
  if(num_cells <= min_cells){
    warning("Not enough cells: ", num_cells, " cells remain, which is at or ",
            "below `min_cells` = ", min_cells, "; returning NA")
    return(NA)
  }

  id_kept_vec <- id_vec[keep_vec]
  cc_kept_vec <- as.character(seurat_obj@meta.data[,case_control_var])
  num_ids <- length(unique(id_kept_vec))
  if(num_ids <= min_ids){
    warning("Not enough individuals: ", num_ids, " individuals remain, which ",
            "is at or below `min_ids` = ", min_ids, "; returning NA")
    return(NA)
  }

  ids_per_arm <- tapply(id_kept_vec, cc_kept_vec, function(x){length(unique(x))})
  if(length(ids_per_arm) < 2 || any(ids_per_arm < min_ids_per_arm)){
    warning("Not enough individuals in an arm: ",
            paste0(names(ids_per_arm), " = ", as.integer(ids_per_arm),
                   collapse = ", "),
            " individuals, and each arm needs at least `min_ids_per_arm` = ",
            min_ids_per_arm, "; returning NA")
    return(NA)
  }

  cells_per_arm <- table(cc_kept_vec)
  if(any(cells_per_arm <= min_cells_casecontrol)){
    warning("Not enough cells in the case or control arm: ",
            paste0(names(cells_per_arm), " = ", as.integer(cells_per_arm),
                   collapse = ", "),
            " cells, and each arm needs more than `min_cells_casecontrol` = ",
            min_cells_casecontrol, "; returning NA")
    return(NA)
  }

  list(dropped_individuals = dropped_individuals,
       seurat_obj = seurat_obj)
}

#' Run eSVD-DE on a Seurat object, with cohort filtering and gene status
#'
#' The recommended entry point. It wraps \code{eSVD} in the five steps a
#' real dataset needs and \code{eSVD} itself refuses to perform:
#' \enumerate{
#'   \item Drop individuals with fewer than \code{min_cells_per_id} cells,
#'         with a warning naming them (\code{filter_cohort}).
#'   \item Check \code{min_cells}, \code{min_ids}, \code{min_ids_per_arm} and
#'         \code{min_cells_casecontrol} on the filtered cohort; on failure,
#'         warn and return \code{NA} (\code{filter_cohort}).
#'   \item Label every gene's \code{gene_status} from the filtered counts and
#'         remove the genes whose counts are all zero. This runs
#'         \emph{after} the drop, because a gene expressed only in a dropped
#'         individual becomes all-zero.
#'   \item Call \code{eSVD} on what remains; every further argument goes
#'         through \code{...}.
#'   \item Reinsert the removed genes at their original positions.
#' }
#'
#' \strong{\code{gene_status}} is a named factor over the original genes, in
#' their original order, with levels \code{"analyzed"} (1) and
#' \code{"all_zero"} (2). An \code{all_zero} gene comes back with \code{NA}
#' in every per-gene estimate (\code{teststat_vec}, \code{case_mean},
#' \code{control_mean}, \code{case_var}, \code{control_var},
#' \code{log2fc_vec}, \code{log2fc_se_vec}, \code{pvalue_list$df_vec},
#' \code{pvalue_list$gaussian_teststat}, and the rows of \code{y_mat},
#' \code{z_mat}, \code{nuisance_vec} and the other per-gene vectors of
#' \code{estimate_nuisance} in the final fit), a
#' \code{pvalue_list$log10pvalue} of \code{0} (that is, a p-value of 1) and a
#' \code{pvalue_list$fdr_vec} of \code{1}. The empirical null and the
#' Benjamini-Hochberg adjustment are computed on the analyzed genes only,
#' which is guaranteed by construction: \code{eSVD} never sees an all-zero
#' gene.
#'
#' Removing all-zero genes changes the factorization of the retained genes,
#' because an all-zero column is not a zero column of the residual matrix
#' the initialization decomposes. Results for the retained genes therefore
#' differ from a run that kept the all-zero genes; they equal a run of
#' \code{eSVD} on the filtered object exactly.
#'
#' @inheritParams eSVD
#' @inheritParams filter_cohort
#' @param ...  Further arguments to \code{eSVD}, such as \code{k},
#'             \code{bool_diet}, \code{cap_multiplier} or
#'             \code{intermediate_save}.
#'
#' @returns Either \code{NA} (the cohort was rejected; see the warning), or
#' the \code{eSVD} object returned by \code{eSVD} with the removed genes
#' reinserted as described above and an added element \code{gene_status}.
#' @examples
#' \donttest{
#' if(requireNamespace("SeuratObject", quietly = TRUE)){
#'   set.seed(10)
#'   sim <- generate_null(cell_per_person = 15, num_genes = 40,
#'                        num_individuals = 8)
#'   meta_df <- data.frame(sim$covariates[, c("CC", "Sex", "Age")])
#'   meta_df$Sex <- factor(meta_df$Sex)
#'   meta_df$Individual <- sim$metadata_individual
#'   rownames(meta_df) <- rownames(sim$obs_mat)
#'   # Seurat rewrites "gene_1" as "gene-1" (with a warning); rename first
#'   count_mat <- Matrix::t(sim$obs_mat)
#'   rownames(count_mat) <- gsub("_", "", rownames(count_mat))
#'   seurat_obj <- SeuratObject::CreateSeuratObject(counts = count_mat,
#'                                                  meta.data = meta_df)
#'
#'   esvd_obj <- eSVD_helper(batch_var_prefix = NULL,
#'                           case_control_levels = c("0", "1"),
#'                           case_control_var = "CC",
#'                           categorical_vars = "Sex",
#'                           id_var = "Individual",
#'                           numerical_vars = "Age",
#'                           seurat_obj = seurat_obj,
#'                           k = 2,
#'                           max_iter = 5)
#'   table(esvd_obj$gene_status)
#'   utils::head(report_results(esvd_obj))
#' }
#' }
#' @export
eSVD_helper <- function(batch_var_prefix, # a variable inside categorical_vars. Can be NULL
                        case_control_levels, # Control and then Case
                        case_control_var,
                        categorical_vars,
                        id_var,
                        numerical_vars,
                        seurat_obj,
                        min_cells = 20,
                        min_cells_casecontrol = 20,
                        min_cells_per_id = 3,
                        min_ids = 4,
                        min_ids_per_arm = 2,
                        verbose = 0,
                        ...){
  # Steps 1 and 2.
  if(verbose > 0) print("Filtering the cohort")
  filter_res <- filter_cohort(seurat_obj = seurat_obj,
                              id_var = id_var,
                              case_control_var = case_control_var,
                              min_cells = min_cells,
                              min_cells_casecontrol = min_cells_casecontrol,
                              min_cells_per_id = min_cells_per_id,
                              min_ids = min_ids,
                              min_ids_per_arm = min_ids_per_arm,
                              verbose = verbose - 1)
  if(!is.list(filter_res)) return(NA)
  seurat_obj <- filter_res$seurat_obj

  # Step 3: gene status, from the DONOR-FILTERED counts.
  if(verbose > 0) print("Labeling gene status")
  count_mat <- .extract_count_matrix(seurat_obj)
  gene_vec <- colnames(count_mat)
  all_zero_idx <- .which_all_zero(count_mat)
  gene_status <- factor(rep("analyzed", length(gene_vec)),
                        levels = c("analyzed", "all_zero"))
  gene_status[all_zero_idx] <- "all_zero"
  names(gene_status) <- gene_vec

  num_analyzed <- sum(gene_status == "analyzed")
  if(num_analyzed == 0){
    stop("every one of the ", length(gene_vec), " genes is all zero after the ",
         "cohort filter, so 0 genes remain to analyze")
  }
  if(length(all_zero_idx) > 0){
    if(verbose > 0){
      print(paste0("Temporarily removing ", length(all_zero_idx),
                   " all-zero gene(s); ", num_analyzed, " remain"))
    }
    # A Seurat FEATURE subset. `seurat_obj[genes, ]` is not supported for
    # this; `subset(features = )` is.
    seurat_obj <- subset(seurat_obj, features = gene_vec[-all_zero_idx])
  }

  # Step 4. `eSVD()` re-checks the per-individual minimum; forward this
  # function's threshold so `min_cells_per_id = 0` really does disable it
  # rather than tripping `eSVD()`'s default of 3 downstream.
  esvd_args <- list(batch_var_prefix = batch_var_prefix,
                    case_control_levels = case_control_levels,
                    case_control_var = case_control_var,
                    categorical_vars = categorical_vars,
                    id_var = id_var,
                    numerical_vars = numerical_vars,
                    seurat_obj = seurat_obj,
                    verbose = verbose)
  dot_args <- list(...)
  if(!"min_cells_per_individual" %in% names(dot_args)){
    dot_args$min_cells_per_individual <- min_cells_per_id
  }
  eSVD_obj <- do.call(eSVD, c(esvd_args, dot_args))

  # Step 5.
  if(verbose > 0) print("Reinserting the removed genes")
  eSVD_obj <- .reinsert_genes(eSVD_obj = eSVD_obj,
                              gene_status = gene_status,
                              count_mat = count_mat)
  eSVD_obj[["gene_status"]] <- gene_status

  eSVD_obj
}

#' Reinsert the removed genes into an eSVD object
#'
#' Every per-gene element is expanded from the analyzed genes to all genes
#' in \code{names(gene_status)}, matched \strong{by name} rather than by
#' position, so it does not depend on the feature order a Seurat subset
#' happens to return. Elements are only padded if present, so both
#' \code{bool_diet} paths are handled by the same code.
#'
#' @param eSVD_obj     Output of \code{eSVD} on the analyzed genes.
#' @param gene_status  Named factor over all genes; see \code{eSVD_helper}.
#' @param count_mat    The full cells-by-genes count matrix (all genes), used
#'                     to restore \code{dat} when it is present.
#'
#' @returns The padded \code{eSVD} object.
#' @noRd
.reinsert_genes <- function(eSVD_obj, gene_status, count_mat){
  gene_vec <- names(gene_status)
  analyzed_vec <- gene_vec[gene_status == "analyzed"]
  stopifnot(setequal(names(eSVD_obj$teststat_vec), analyzed_vec))

  pad_vector <- function(vec, fill){
    if(is.null(vec)) return(NULL)
    stopifnot(!is.null(names(vec)), all(names(vec) %in% analyzed_vec))
    out <- rep(fill, length(gene_vec))
    names(out) <- gene_vec
    out[names(vec)] <- vec
    out
  }
  # A factor assigned into a numeric vector leaves its integer codes behind.
  pad_factor <- function(vec){
    if(is.null(vec)) return(NULL)
    stopifnot(is.factor(vec), !is.null(names(vec)),
              all(names(vec) %in% analyzed_vec))
    out <- rep(NA_character_, length(gene_vec))
    names(out) <- gene_vec
    out[names(vec)] <- as.character(vec)
    out_factor <- factor(out, levels = levels(vec))
    names(out_factor) <- gene_vec
    out_factor
  }
  pad_rows <- function(mat){
    if(is.null(mat)) return(NULL)
    stopifnot(!is.null(rownames(mat)), all(rownames(mat) %in% analyzed_vec))
    out <- matrix(NA_real_, nrow = length(gene_vec), ncol = ncol(mat),
                  dimnames = list(gene_vec, colnames(mat)))
    out[rownames(mat), ] <- mat
    out
  }
  pad_columns <- function(mat){
    if(is.null(mat)) return(NULL)
    stopifnot(!is.null(colnames(mat)), all(colnames(mat) %in% analyzed_vec))
    out <- matrix(NA_real_, nrow = nrow(mat), ncol = length(gene_vec),
                  dimnames = list(rownames(mat), gene_vec))
    out[, colnames(mat)] <- mat
    out
  }

  eSVD_obj$teststat_vec <- pad_vector(eSVD_obj$teststat_vec, fill = NA_real_)
  eSVD_obj$case_mean <- pad_vector(eSVD_obj$case_mean, fill = NA_real_)
  eSVD_obj$control_mean <- pad_vector(eSVD_obj$control_mean, fill = NA_real_)
  for(element_name in c("case_var", "control_var", "log2fc_vec",
                        "log2fc_se_vec")){
    eSVD_obj[[element_name]] <- pad_vector(eSVD_obj[[element_name]],
                                           fill = NA_real_)
  }

  pvalue_list <- eSVD_obj$pvalue_list
  pvalue_list$df_vec <- pad_vector(pvalue_list$df_vec, fill = NA_real_)
  pvalue_list$gaussian_teststat <- pad_vector(pvalue_list$gaussian_teststat,
                                              fill = NA_real_)
  # Deliberately not NA (Q-STATUS-3): a downstream `which(fdr < 0.05)` must
  # not have to think about missingness. p = 1 is -log10(p) = 0.
  pvalue_list$log10pvalue <- pad_vector(pvalue_list$log10pvalue, fill = 0)
  pvalue_list$fdr_vec <- pad_vector(pvalue_list$fdr_vec, fill = 1)
  eSVD_obj$pvalue_list <- pvalue_list

  # The factorization: the gene dimension of every fit that survived.
  # `x_mat` is cells by k and is left alone.
  fit_names <- names(eSVD_obj)[sapply(eSVD_obj, inherits, "eSVD_Fit")]
  for(fit_name in fit_names){
    fit <- eSVD_obj[[fit_name]]
    fit$y_mat <- pad_rows(fit$y_mat)
    fit$z_mat <- pad_rows(fit$z_mat)
    for(element_name in .nuisance_numeric_elements()){
      fit[[element_name]] <- pad_vector(fit[[element_name]], fill = NA_real_)
    }
    fit$nuisance_status <- pad_factor(fit$nuisance_status)
    fit$posterior_mean_mat <- pad_columns(fit$posterior_mean_mat)
    fit$posterior_var_mat <- pad_columns(fit$posterior_var_mat)
    eSVD_obj[[fit_name]] <- fit
  }

  # With `bool_diet = FALSE` the counts are still on the object; restore the
  # full matrix so its columns line up with the padded `y_mat` rows.
  if(!is.null(eSVD_obj$dat)){
    eSVD_obj$dat <- count_mat
  }

  eSVD_obj
}

# The per-gene numeric vectors `estimate_nuisance` stores on a fit, in one
# place so that padding them and stripping the padding cannot drift.
.nuisance_numeric_elements <- function(){
  c("nuisance_vec", "nuisance_mle_vec", "nuisance_library_median_vec",
    "gene_mean_count_vec", "gene_sparsity_vec")
}
