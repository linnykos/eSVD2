#' @useDynLib eSVD2
#' @importFrom Rcpp evalCpp
NULL
# see https://stackoverflow.com/questions/74172476/devtoolsdocument-wont-include-usedylib-in-namespace
#
# There is deliberately no `@exportPattern` here. It used to export every
# non-dot-prefixed object, including the raw Rcpp bindings that take external
# pointers; the public API is now whatever carries an explicit `@export`.

.onUnload <- function(libpath) {
  library.dynam.unload("eSVD2", libpath)
}
