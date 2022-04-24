data_loader <- function(mat)
{
  ptr <- .data_loader(mat)
  structure(ptr, class = "esvd_data_loader")
}

#' @export
print.esvd_data_loader <- function(x, ...)
{
  data_loader_description(x)
}
