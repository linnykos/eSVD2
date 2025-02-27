% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_pvalue.R
\name{compute_pvalue}
\alias{compute_pvalue}
\title{Compute p-values}
\usage{
compute_pvalue(input_obj, verbose = 0, ...)
}
\arguments{
\item{input_obj}{\code{eSVD} object outputed from \code{compute_test_statistic}.}

\item{verbose}{Integer.}

\item{...}{Additional parameters.}
}
\value{
\code{eSVD} object with added element \code{"pvalue_list"}
}
\description{
Compute p-values
}
