#' Internal data loader function
#'
#' @param mat matrix
#'
#' @return pointer
data_loader <- function(mat)
{
  ptr <- .data_loader(mat)
  structure(ptr, class = "esvd_data_loader")
}

#' Internal data loader function
#'
#' @param x object
#' @param ... additional arguments
#'
#' @returns description
#' @export
print.esvd_data_loader <- function(x, ...)
{
  data_loader_description(x)
}
