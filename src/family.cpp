#include "family.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Environment;
using Rcpp::Function;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using MapMat = Eigen::Map<MatrixXd>;
using MapVec = Eigen::Map<VectorXd>;
using DistrXPtr = Rcpp::XPtr<Distribution>;
using DataLoaderXPtr = Rcpp::XPtr<DataLoader>;

// theta = U * Vi + P * Qi
// U and P are matrices
// Vi and Qi are vectors
// P and Qi may be NULL
inline VectorXd compute_theta(MapMat U, MapVec Vi, SEXP P_, SEXP Qi_)
{
    VectorXd theta = U * Vi;
    if(P_ == R_NilValue || Qi_ == R_NilValue)
        return theta;

    MapMat P = Rcpp::as<MapMat>(P_);
    MapVec Qi = Rcpp::as<MapVec>(Qi_);
    if(P.cols() < 1 || Qi.size() < 1)
        return theta;

    theta.noalias() += P * Qi;
    return theta;
}

// offset is a scalar
template <typename T>
inline VectorXd compute_theta_with_offset(MapMat U, MapVec Vi, SEXP P, SEXP Qi, T offset)
{
  VectorXd theta = compute_theta(U, Vi, P, Qi);
  theta.array() += offset;
  return theta;
}

// offset is a vector
template <>
inline VectorXd compute_theta_with_offset(MapMat U, MapVec Vi, SEXP P, SEXP Qi, MapVec offset)
{
  VectorXd theta = compute_theta(U, Vi, P, Qi);
  theta.noalias() += offset;
  return theta;
}



// [[Rcpp::export]]
double objfn_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    SEXP loader, int row_ind, Environment family,
    double si, MapVec gamma, double offseti, double l2pen
)
{
    const int p = Y.rows();
    VectorXd thetai = compute_theta_with_offset(Y, Xi, B, Zi, offseti);

    DistrXPtr distr = family["cpp_functions"];
    DataLoaderXPtr data_loader(loader);
    VectorXd log_prob(p);
    int non_na = distr->log_prob_row(data_loader, row_ind, thetai.data(), si, gamma.data(), log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    return (-log_prob.sum() + l2pen * Xi.squaredNorm()) / double(non_na);
}

// [[Rcpp::export]]
double objfn_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    SEXP loader, int col_ind, Environment family,
    MapVec s, double gammaj, MapVec offset, double l2pen
)
{
    const int n = X.rows();
    VectorXd thetaj = compute_theta_with_offset(X, Yj, Z, Bj, offset);

    DistrXPtr distr = family["cpp_functions"];
    DataLoaderXPtr data_loader(loader);
    VectorXd log_prob(n);
    int non_na = distr->log_prob_col(data_loader, col_ind, thetaj.data(), s.data(), gammaj, log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    return (-log_prob.sum() + l2pen * Yj.squaredNorm()) / double(non_na);
}

// [[Rcpp::export]]
NumericVector grad_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    SEXP loader, int row_ind, Environment family,
    double si, MapVec gamma, double offseti, double l2pen
)
{
    const int k = Xi.size();
    const int p = Y.rows();
    NumericVector res_ = NumericVector(Rcpp::no_init_vector(k));
    MapVec res = Rcpp::as<MapVec>(res_);
    VectorXd thetai = compute_theta_with_offset(Y, Xi, B, Zi, offseti);

    DistrXPtr distr = family["cpp_functions"];
    DataLoaderXPtr data_loader(loader);
    VectorXd dlog_prob(p);
    int non_na = distr->d12log_prob_row(data_loader, row_ind, thetai.data(), si, gamma.data(), dlog_prob.data(), NULL);
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    res.noalias() = -Y.transpose() * dlog_prob + 2.0 * l2pen * Xi;
    res /= double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericVector grad_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    SEXP loader, int col_ind, Environment family,
    MapVec s, double gammaj, MapVec offset, double l2pen
)
{
    const int k = Yj.size();
    const int n = X.rows();
    NumericVector res_ = NumericVector(Rcpp::no_init_vector(k));
    MapVec res = Rcpp::as<MapVec>(res_);
    VectorXd thetaj = compute_theta_with_offset(X, Yj, Z, Bj, offset);

    DistrXPtr distr = family["cpp_functions"];
    DataLoaderXPtr data_loader(loader);
    VectorXd dlog_prob(n);
    int non_na = distr->d12log_prob_col(data_loader, col_ind, thetaj.data(), s.data(), gammaj, dlog_prob.data(), NULL);
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    res.noalias() = -X.transpose() * dlog_prob + 2.0 * l2pen * Yj;
    res /= double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericMatrix hessian_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    SEXP loader, int row_ind, Environment family,
    double si, MapVec gamma, double offseti, double l2pen
)
{
    const int k = Xi.size();
    const int p = Y.rows();
    NumericMatrix res_ = NumericMatrix(Rcpp::no_init_matrix(k, k));
    MapMat res = Rcpp::as<MapMat>(res_);
    VectorXd thetai = compute_theta_with_offset(Y, Xi, B, Zi, offseti);

    DistrXPtr distr = family["cpp_functions"];
    DataLoaderXPtr data_loader(loader);
    VectorXd d2log_prob(p);
    int non_na = distr->d12log_prob_row(data_loader, row_ind, thetai.data(), si, gamma.data(), NULL, d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    res.noalias() = -Y.transpose() * d2log_prob.asDiagonal() * Y;
    res.diagonal().array() += 2.0 * l2pen;
    res /= double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericMatrix hessian_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    SEXP loader, int col_ind, Environment family,
    MapVec s, double gammaj, MapVec offset, double l2pen
)
{
    const int k = Yj.size();
    const int n = X.rows();
    NumericMatrix res_ = NumericMatrix(Rcpp::no_init_matrix(k, k));
    MapMat res = Rcpp::as<MapMat>(res_);
    VectorXd thetaj = compute_theta_with_offset(X, Yj, Z, Bj, offset);

    DistrXPtr distr = family["cpp_functions"];
    DataLoaderXPtr data_loader(loader);
    VectorXd d2log_prob(n);
    int non_na = distr->d12log_prob_col(data_loader, col_ind, thetaj.data(), s.data(), gammaj, NULL, d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    res.noalias() = -X.transpose() * d2log_prob.asDiagonal() * X;
    res.diagonal().array() += 2.0 * l2pen;
    res /= double(non_na);
    return res_;
}

// [[Rcpp::export]]
List direction_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi,
    SEXP loader, int row_ind, Environment family,
    double si, MapVec gamma, double offseti, double l2pen
)
{
    const int k = Xi.size();
    const int p = Y.rows();
    VectorXd thetai = compute_theta_with_offset(Y, Xi, B, Zi, offseti);

    DistrXPtr distr = family["cpp_functions"];
    DataLoaderXPtr data_loader(loader);
    VectorXd dlog_prob(p), d2log_prob(p);
    int non_na = distr->d12log_prob_row(data_loader, row_ind, thetai.data(), si, gamma.data(),
                                        dlog_prob.data(), d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");

    VectorXd g = Y.transpose() * (-dlog_prob) + 2.0 * l2pen * Xi;
    MatrixXd H = Y.transpose() * (-d2log_prob).asDiagonal() * Y;
    H.diagonal().array() += 2.0 * l2pen;

    VectorXd direc(k);
    Eigen::LLT<MatrixXd> solver(H);
    if(solver.info() != Eigen::Success)
    {
        // Fall back to gradient direction if Hessian is singular
        direc.noalias() = -g / double(non_na);
    } else {
        direc.noalias() = -solver.solve(g);
    }
    g /= double(non_na);
    return List::create(
        Rcpp::Named("grad") = g,
        Rcpp::Named("direction") = direc
    );
}

// [[Rcpp::export]]
List direction_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z,
    SEXP loader, int col_ind, Environment family,
    MapVec s, double gammaj, MapVec offset, double l2pen
)
{
    const int k = Yj.size();
    const int n = X.rows();
    VectorXd thetaj = compute_theta_with_offset(X, Yj, Z, Bj, offset);

    DistrXPtr distr = family["cpp_functions"];
    DataLoaderXPtr data_loader(loader);
    VectorXd dlog_prob(n), d2log_prob(n);
    int non_na = distr->d12log_prob_col(data_loader, col_ind, thetaj.data(), s.data(), gammaj,
                                        dlog_prob.data(), d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");

    VectorXd g = X.transpose() * (-dlog_prob) + 2.0 * l2pen * Yj;
    MatrixXd H = X.transpose() * (-d2log_prob).asDiagonal() * X;
    H.diagonal().array() += 2.0 * l2pen;

    VectorXd direc(k);
    Eigen::LLT<MatrixXd> solver(H);
    if(solver.info() != Eigen::Success)
    {
        // Fall back to gradient direction if Hessian is singular
        direc.noalias() = -g / double(non_na);
    } else {
        direc.noalias() = -solver.solve(g);
    }
    g /= double(non_na);
    return List::create(
        Rcpp::Named("grad") = g,
        Rcpp::Named("direction") = direc
    );
}

bool feas_Xi_impl(
    MapVec Xi, MapMat Y, SEXP B, SEXP Zi, Environment family, double offseti
)
{
    DistrXPtr distr = family["cpp_functions"];
    if(distr->feas_always())
        return true;

    VectorXd thetai = compute_theta_with_offset(Y, Xi, B, Zi, offseti);
    return distr->feasibility(thetai.size(), thetai.data());
}

bool feas_Yj_impl(
    MapVec Yj, MapMat X, SEXP Bj, SEXP Z, Environment family, MapVec offset
)
{
    DistrXPtr distr = family["cpp_functions"];
    if(distr->feas_always())
        return true;

    VectorXd thetaj = compute_theta_with_offset(X, Yj, Z, Bj, offset);
    return distr->feasibility(thetaj.size(), thetaj.data());
}
