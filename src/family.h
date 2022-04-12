#ifndef ESVD2_FAMILY_H
#define ESVD2_FAMILY_H

#include <RcppEigen.h>
#include "data_loader.h"
#include "distribution.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::List;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using MapMat = Eigen::Map<MatrixXd>;
using MapVec = Eigen::Map<VectorXd>;

double objfn_Xi(
    const MapVec& XCi, const MapMat& YZ, int k,
    DataLoader* loader, int row_ind, const Distribution* distr,
    double si, const MapVec& gamma, double l2penx
);

NumericVector grad_Xi(
    const MapVec& XCi, const MapMat& YZ, int k,
    DataLoader* loader, int row_ind, const Distribution* distr,
    double si, const MapVec& gamma, double l2penx
);

NumericMatrix hessian_Xi(
    const MapVec& XCi, const MapMat& YZ, int k,
    DataLoader* loader, int row_ind, const Distribution* distr,
    double si, const MapVec& gamma, double l2penx
);

List direction_Xi(
    const MapVec& XCi, const MapMat& YZ, int k,
    DataLoader* loader, int row_ind, const Distribution* distr,
    double si, const MapVec& gamma, double l2penx
);

bool feas_Xi(const MapVec& XCi, const MapMat& YZ, const Distribution* distr);



double objfn_YZj(
    const MapMat& XC, const MapVec& YZj, int k, IntegerVector YZind,
    DataLoader* loader, int col_ind, const Distribution* distr,
    const MapVec& s, double gammaj, double l2peny, double l2penz
);

NumericVector grad_YZj(
    const MapMat& XC, const MapVec& YZj, int k, IntegerVector YZind,
    DataLoader* loader, int col_ind, const Distribution* distr,
    const MapVec& s, double gammaj, double l2peny, double l2penz
);

NumericMatrix hessian_YZj(
    const MapMat& XC, const MapVec& YZj, int k, IntegerVector YZind,
    DataLoader* loader, int col_ind, const Distribution* distr,
    const MapVec& s, double gammaj, double l2peny, double l2penz
);

List direction_YZj(
    const MapMat& XC, const MapVec& YZj, int k, IntegerVector YZind,
    DataLoader* loader, int col_ind, const Distribution* distr,
    const MapVec& s, double gammaj, double l2peny, double l2penz
);

bool feas_YZj(const MapMat& XC, const MapVec& YZj, const Distribution* distr);


#endif  // ESVD2_FAMILY_H
