#ifndef ESVD2_FAMILY_H
#define ESVD2_FAMILY_H

#include <RcppEigen.h>
#include "distribution.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Environment;

double objfn_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_, double offseti
);

double objfn_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_, NumericVector offset
);

NumericVector grad_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_, double offseti
);

NumericVector grad_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_, NumericVector offset
);

NumericMatrix hessian_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_, double offseti
);

NumericMatrix hessian_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_, NumericVector offset
);

List direction_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_, double offseti
);

List direction_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_, NumericVector offset
);

bool feas_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_, Environment family_, double offseti
);

bool feas_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_, Environment family_, NumericVector offset
);

#endif  // ESVD2_FAMILY_H
