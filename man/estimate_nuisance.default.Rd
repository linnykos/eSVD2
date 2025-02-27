% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nuisance.R
\name{estimate_nuisance.default}
\alias{estimate_nuisance.default}
\title{Estimate nuisance values for matrix or sparse matrices.}
\usage{
\method{estimate_nuisance}{default}(
  input_obj,
  mean_mat,
  library_mat,
  bool_use_log = F,
  min_val = 1e-04,
  verbose = 0,
  ...
)
}
\arguments{
\item{input_obj}{Dataset (either \code{matrix} or \code{dgCMatrix}) where the \eqn{n} rows represent cells
and \eqn{p} columns represent genes.
The rows and columns of the matrix should be named.}

\item{mean_mat}{A \code{matrix} of \eqn{n} rows and \eqn{p} columns that represents the
expected value of each entry.}

\item{library_mat}{A \code{matrix} of \eqn{n} rows and \eqn{p} columns that represents the
library size of each entry.}

\item{bool_use_log}{Boolean if the nuisance (i.e., over-dispersion) parameter should be estimated on the log scale, default is \code{FALSE}.}

\item{min_val}{Minimum value of the nuisance parameter.}

\item{verbose}{Integer.}

\item{...}{Additional parameters.}
}
\value{
Vector of length \eqn{p}
}
\description{
Estimate nuisance values for matrix or sparse matrices.
}
