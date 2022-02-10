# Environment to store package-wide variables
.pkg.env <- new.env()

# Flag to indicate whether the C++ code has been loaded
.pkg.env$cpp_loaded <- FALSE

# Function to load C++ code. Will be called in the optimization function
load_cpp_code <- function()
{
  # This makes sure that the code is loaded only once
  if(.pkg.env$cpp_loaded)
    return(invisible(NULL))

  .esvd.neg_binom$cpp_functions <- distribution_neg_binom()
  .esvd.neg_binom2$cpp_functions <- distribution_neg_binom2()
  .esvd.poisson$cpp_functions <- distribution_poisson()

  .pkg.env$cpp_loaded <- TRUE
  invisible(NULL)
}

# Unload C++ code. Will be called when the package is unloaded
unload_cpp_code <- function()
{
  .esvd.neg_binom$cpp_functions <- NULL
  .esvd.neg_binom2$cpp_functions <- NULL
  .esvd.poisson$cpp_functions <- NULL
  gc()
}
