#ifndef ESVD2_FAMILY_H
#define ESVD2_FAMILY_H

#include <RcppEigen.h>
#include "data_loader.h"
#include "distribution.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Environment;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using MapMat = Eigen::Map<MatrixXd>;
using MapVec = Eigen::Map<VectorXd>;

double objfn_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    SEXP loader, int row_ind, Environment family,
    double si, MapVec gamma, double offseti, double l2pen
);

double objfn_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    SEXP loader, int col_ind, Environment family,
    MapVec s, double gammaj, MapVec offset, double l2pen
);

NumericVector grad_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    SEXP loader, int row_ind, Environment family,
    double si, MapVec gamma, double offseti, double l2pen
);

NumericVector grad_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    SEXP loader, int col_ind, Environment family,
    MapVec s, double gammaj, MapVec offset, double l2pen
);

NumericMatrix hessian_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    SEXP loader, int row_ind, Environment family,
    double si, MapVec gamma, double offseti, double l2pen
);

NumericMatrix hessian_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    SEXP loader, int col_ind, Environment family,
    MapVec s, double gammaj, MapVec offset, double l2pen
);

List direction_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    SEXP loader, int row_ind, Environment family,
    double si, MapVec gamma, double offseti, double l2pen
);

List direction_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    SEXP loader, int col_ind, Environment family,
    MapVec s, double gammaj, MapVec offset, double l2pen
);

bool feas_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi, Environment family, double offseti
);

bool feas_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z, Environment family, MapVec offset
);

#endif  // ESVD2_FAMILY_H
