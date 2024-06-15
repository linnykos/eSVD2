#' @useDynLib eSVD2
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL
# see https://stackoverflow.com/questions/74172476/devtoolsdocument-wont-include-usedylib-in-namespace

.onUnload <- function(libpath) {
  library.dynam.unload("eSVD2", libpath)
}
