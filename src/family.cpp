#include <RcppEigen.h>
#include "distribution.h"
#include "family.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Environment;
using Rcpp::Function;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;

// theta = U * Vi + P * Qi
// U and P are matrices
// Vi and Qi are vectors
// P and Qi may be NULL
inline VectorXd compute_theta(NumericMatrix U_, NumericVector Vi_, SEXP P_, SEXP Qi_)
{
    MapMat U = Rcpp::as<MapMat>(U_);
    MapVec Vi = Rcpp::as<MapVec>(Vi_);
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

// [[Rcpp::export]]
double objfn_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_
)
{
    const int p = Y_.nrow();
    VectorXd thetai = compute_theta(Y_, Xi_, B_, Zi_);

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd log_prob(p);
    int non_na = distr->log_prob_row(p, Ai_.begin(), thetai.data(), si_, gamma_.begin(), log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    return -log_prob.sum() / double(non_na);
}

// [[Rcpp::export]]
double objfn_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_
)
{
    const int n = X_.nrow();
    VectorXd thetaj = compute_theta(X_, Yj_, Z_, Bj_);

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd log_prob(n);
    int non_na = distr->log_prob_col(n, Aj_.begin(), thetaj.data(), s_.begin(), gammaj_, log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    return -log_prob.sum() / double(non_na);
}

// [[Rcpp::export]]
NumericVector grad_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_
)
{
    const int k = Xi_.length();
    const int p = Y_.nrow();
    NumericVector res_ = NumericVector(Rcpp::no_init_vector(k));
    MapVec res = Rcpp::as<MapVec>(res_);
    MapMat Y = Rcpp::as<MapMat>(Y_);
    VectorXd thetai = compute_theta(Y_, Xi_, B_, Zi_);

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd dlog_prob(p);
    int non_na = distr->d12log_prob_row(p, Ai_.begin(), thetai.data(), si_, gamma_.begin(), dlog_prob.data(), NULL);
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    res.noalias() = Y.transpose() * dlog_prob;
    res /= -double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericVector grad_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_
)
{
    const int k = Yj_.length();
    const int n = X_.nrow();
    NumericVector res_ = NumericVector(Rcpp::no_init_vector(k));
    MapVec res = Rcpp::as<MapVec>(res_);
    MapMat X = Rcpp::as<MapMat>(X_);
    VectorXd thetaj = compute_theta(X_, Yj_, Z_, Bj_);

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd dlog_prob(n);
    int non_na = distr->d12log_prob_col(n, Aj_.begin(), thetaj.data(), s_.begin(), gammaj_, dlog_prob.data(), NULL);
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    res.noalias() = X.transpose() * dlog_prob;
    res /= -double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericMatrix hessian_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_
)
{
    const int k = Xi_.length();
    const int p = Y_.nrow();
    NumericMatrix res_ = NumericMatrix(Rcpp::no_init_matrix(k, k));
    MapMat res = Rcpp::as<MapMat>(res_);
    MapMat Y = Rcpp::as<MapMat>(Y_);
    VectorXd thetai = compute_theta(Y_, Xi_, B_, Zi_);

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd d2log_prob(p);
    int non_na = distr->d12log_prob_row(p, Ai_.begin(), thetai.data(), si_, gamma_.begin(), NULL, d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    res.noalias() = Y.transpose() * d2log_prob.asDiagonal() * Y;
    res /= -double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericMatrix hessian_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_
)
{
    const int k = Yj_.length();
    const int n = X_.nrow();
    NumericMatrix res_ = NumericMatrix(Rcpp::no_init_matrix(k, k));
    MapMat res = Rcpp::as<MapMat>(res_);
    MapMat X = Rcpp::as<MapMat>(X_);
    VectorXd thetaj = compute_theta(X_, Yj_, Z_, Bj_);

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd d2log_prob(n);
    int non_na = distr->d12log_prob_col(n, Aj_.begin(), thetaj.data(), s_.begin(), gammaj_, NULL, d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    res.noalias() = X.transpose() * d2log_prob.asDiagonal() * X;
    res /= -double(non_na);
    return res_;
}

// [[Rcpp::export]]
List direction_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_
)
{
    const int p = Y_.nrow();
    MapMat Y = Rcpp::as<MapMat>(Y_);
    VectorXd thetai = compute_theta(Y_, Xi_, B_, Zi_);

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd dlog_prob(p), d2log_prob(p);
    int non_na = distr->d12log_prob_row(p, Ai_.begin(), thetai.data(), si_, gamma_.begin(),
                                        dlog_prob.data(), d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");

    VectorXd g = Y.transpose() * (-dlog_prob);
    MatrixXd H = Y.transpose() * (-d2log_prob).asDiagonal() * Y;
    Eigen::LLT<MatrixXd> solver(H);
    if(solver.info() != Eigen::Success)
        Rcpp::stop("the Hessian matrix of Xi is singular");
    VectorXd direc = -solver.solve(g);
    g /= double(non_na);
    return List::create(
        Rcpp::Named("grad") = g,
        Rcpp::Named("direction") = direc
    );
}

// [[Rcpp::export]]
List direction_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_
)
{
    const int n = X_.nrow();
    MapMat X = Rcpp::as<MapMat>(X_);
    VectorXd thetaj = compute_theta(X_, Yj_, Z_, Bj_);

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd dlog_prob(n), d2log_prob(n);
    int non_na = distr->d12log_prob_col(n, Aj_.begin(), thetaj.data(), s_.begin(), gammaj_,
                                        dlog_prob.data(), d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");

    VectorXd g = X.transpose() * (-dlog_prob);
    MatrixXd H = X.transpose() * (-d2log_prob).asDiagonal() * X;
    Eigen::LLT<MatrixXd> solver(H);
    if(solver.info() != Eigen::Success)
        Rcpp::stop("the Hessian matrix of Yj is singular");
    VectorXd direc = -solver.solve(g);
    g /= double(non_na);
    return List::create(
        Rcpp::Named("grad") = g,
        Rcpp::Named("direction") = direc
    );
}

bool feas_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, SEXP B_, SEXP Zi_, Environment family_
)
{
    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    if(distr->feas_always())
        return true;

    VectorXd thetai = compute_theta(Y_, Xi_, B_, Zi_);
    return distr->feasibility(thetai.size(), thetai.data());
}

bool feas_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, SEXP Bj_, SEXP Z_, Environment family_
)
{
    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    if(distr->feas_always())
        return true;

    VectorXd thetaj = compute_theta(X_, Yj_, Z_, Bj_);
    return distr->feasibility(thetaj.size(), thetaj.data());
}
