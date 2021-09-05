.onLoad <- function(libname, pkgname) {
  # library.dynam("eSVD2", pkgname, libname)

  .esvd.neg_binom$cpp_functions <- distribution_neg_binom()
  .esvd.poisson$cpp_functions <- distribution_poisson()
}

.onUnload <- function(libpath) {
  .esvd.neg_binom$cpp_functions <- NULL
  .esvd.poisson$cpp_functions <- NULL
  gc()

  library.dynam.unload("eSVD2", libpath)
}
