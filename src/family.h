#ifndef ESVD2_FAMILY_H
#define ESVD2_FAMILY_H

#include <RcppEigen.h>
#include "distribution.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Environment;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;

double objfn_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    MapVec Ai, Environment family,
    double si, MapVec gamma, double offseti
);

double objfn_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    MapVec Aj, Environment family,
    MapVec s, double gammaj, MapVec offset
);

NumericVector grad_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    MapVec Ai, Environment family,
    double si, MapVec gamma, double offseti
);

NumericVector grad_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    MapVec Aj, Environment family,
    MapVec s, double gammaj, MapVec offset
);

NumericMatrix hessian_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    MapVec Ai, Environment family,
    double si, MapVec gamma, double offseti
);

NumericMatrix hessian_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    MapVec Aj, Environment family,
    MapVec s, double gammaj, MapVec offset
);

List direction_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    MapVec Ai, Environment family,
    double si, MapVec gamma, double offseti
);

List direction_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    MapVec Aj, Environment family,
    MapVec s, double gammaj, MapVec offset
);

bool feas_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi, Environment family, double offseti
);

bool feas_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z, Environment family, MapVec offset
);

#endif  // ESVD2_FAMILY_H
