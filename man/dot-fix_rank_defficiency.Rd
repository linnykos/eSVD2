% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialization.R
\name{.fix_rank_defficiency}
\alias{.fix_rank_defficiency}
\title{Fix rank defficiency among two matrices}
\usage{
.fix_rank_defficiency(x_mat, y_mat, domain)
}
\arguments{
\item{x_mat}{a numeric matrix}

\item{y_mat}{a numeric matrix with the same number of columns as \code{x_mat}}

\item{domain}{a vector where \code{domain[1] < domain[2]}}
}
\value{
a list of \code{x_mat} and \code{y_mat}
}
\description{
Given two matrices, \code{x_mat} and \code{y_mat} (both with the same number of columns),
adjust these two matrices so the \code{x_mat \%*\% t(y_mat)} is actually of the desired
rank. This function is needed since sometimes upstream, the matrices \code{x_mat} and \code{y_mat}
do not actually have the rank equal to the number of columns (i.e., empirically we have
observed that \code{x_mat} might have a column that is all constant).
}
\details{
Here, \code{domain} represents the domain where all the inner products in
\code{x_mat \%*\% t(y_mat)} are assumed to lie between
}
