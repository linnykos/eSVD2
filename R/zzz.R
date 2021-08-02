.onLoad <- function(libname, pkgname) {
  library.dynam("eSVD2", pkgname, libname)

  .esvd.poisson$cpp_functions <- distribution_poisson()
}

.onUnload <- function(libpath) {
  .esvd.poisson$cpp_functions <- NULL
  gc()

  library.dynam.unload("eSVD2", libpath)
}
