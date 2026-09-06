#' Internal data loader
#'
#' Wraps a dense \code{matrix} or a \code{dgCMatrix} in a C++ iterator used
#' by the optimizer. The pointer is only valid in the R session that created
#' it; a restored copy (e.g. from \code{readRDS}) errors when used.
#'
#' @param mat \code{matrix} (integer or double) or \code{dgCMatrix}, cells by
#'            genes.
#'
#' @returns An external pointer of class \code{esvd_data_loader}.
#' @keywords internal
data_loader <- function(mat)
{
  ptr <- .data_loader(mat)
  structure(ptr, class = "esvd_data_loader")
}

#' Print a data loader
#'
#' @param x   An \code{esvd_data_loader}, from \code{data_loader}.
#' @param ... Ignored.
#'
#' @returns \code{x}, invisibly; called for its printed one-line description.
#' @keywords internal
#' @export
print.esvd_data_loader <- function(x, ...)
{
  data_loader_description(x)
  invisible(x)
}
