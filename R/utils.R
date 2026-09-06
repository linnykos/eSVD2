.l2norm <- function(x){sqrt(sum(x^2))}

.rotate <- function(a) { t(a[nrow(a):1,]) }

# for diag(vec) %*% mat
.mult_vec_mat <- function(vec, mat){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == nrow(mat))
  vec * mat
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == ncol(mat))
  mat * rep(vec, rep(nrow(mat), length(vec)))
}

## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))

  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]

  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

# see https://stackoverflow.com/questions/7944809/assigning-null-to-a-list-element-in-r
.combine_two_named_lists <- function(list1, list2){
  idx <- which(!names(list2) %in% names(list1))
  for(i in idx){
    if(all(is.null(list2[[i]]))){
      list1 <- c(list1, list(TEMP_NAME = NULL))
      names(list1)[which(names(list1) == "TEMP_NAME")] <- names(list2)[i]
    } else {
      list1[[names(list2)[i]]] <- list2[[i]]
    }
  }

  list1
}

#' Which genes are all zero?
#'
#' The one definition of "all zero" shared by \code{eSVD_helper} (which
#' labels and removes such genes), \code{eSVD} and \code{initialize_esvd}
#' (which both refuse them). Three copies of this predicate would drift, and
#' the helper would then start passing genes that \code{eSVD} rejects.
#'
#' \code{NA} entries count as zero, consistent with \code{initialize_esvd}
#' zeroing them on entry: a gene that is all \code{NA} is all zero.
#'
#' @param dat  \code{matrix} or \code{dgCMatrix}, cells by genes.
#'
#' @returns Integer vector of column indices whose column sum is zero.
#' @noRd
.which_all_zero <- function(dat){
  stopifnot(inherits(dat, c("matrix", "dgCMatrix")))

  which(Matrix::colSums(dat, na.rm = TRUE) == 0)
}

# An integer attribute set by the C++ optimizers, or 0 when absent (a matrix
# that never went through `opt_x`/`opt_yz`).
.attr_or_zero <- function(mat, attr_name){
  val <- attr(mat, attr_name, exact = TRUE)
  if(is.null(val)) return(0)
  as.numeric(val)
}
