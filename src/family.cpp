#include <RcppEigen.h>
#include "distribution.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Environment;
using Rcpp::Function;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;

// [[Rcpp::export]]
double objfn_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, NumericMatrix B_, NumericVector Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_
)
{
    const int p = Y_.nrow();
    const int r = B_.ncol();
    NumericVector thetai_ = NumericVector(Rcpp::no_init_vector(p));

    MapVec thetai = Rcpp::as<MapVec>(thetai_);
    MapVec Xi = Rcpp::as<MapVec>(Xi_);
    MapMat Y = Rcpp::as<MapMat>(Y_);
    MapMat B = Rcpp::as<MapMat>(B_);
    MapVec Zi = Rcpp::as<MapVec>(Zi_);

    thetai.noalias() = Y * Xi;
    if(r > 0)
        thetai.noalias() += B * Zi;

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd log_prob(p);
    int non_na = distr->log_prob_row(p, Ai_.begin(), thetai_.begin(), si_, gamma_.begin(), log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    return -log_prob.sum() / double(non_na);
}

// [[Rcpp::export]]
double objfn_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, NumericVector Bj_, NumericMatrix Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_
)
{
    const int n = X_.nrow();
    const int r = Z_.ncol();
    NumericVector thetaj_ = NumericVector(Rcpp::no_init_vector(n));

    MapVec thetaj = Rcpp::as<MapVec>(thetaj_);
    MapMat X = Rcpp::as<MapMat>(X_);
    MapVec Yj = Rcpp::as<MapVec>(Yj_);
    MapVec Bj = Rcpp::as<MapVec>(Bj_);
    MapMat Z = Rcpp::as<MapMat>(Z_);

    thetaj.noalias() = X * Yj;
    if(r > 0)
        thetaj.noalias() += Z * Bj;

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd log_prob(n);
    int non_na = distr->log_prob_col(n, Aj_.begin(), thetaj_.begin(), s_.begin(), gammaj_, log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    return -log_prob.sum() / double(non_na);
}

// [[Rcpp::export]]
NumericVector grad_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, NumericMatrix B_, NumericVector Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_
)
{
    const int k = Xi_.length();
    const int p = Y_.nrow();
    const int r = B_.ncol();
    NumericVector res_ = NumericVector(Rcpp::no_init_vector(k));
    NumericVector thetai_ = NumericVector(Rcpp::no_init_vector(p));

    MapVec res = Rcpp::as<MapVec>(res_);
    MapVec thetai = Rcpp::as<MapVec>(thetai_);
    MapVec Xi = Rcpp::as<MapVec>(Xi_);
    MapMat Y = Rcpp::as<MapMat>(Y_);
    MapMat B = Rcpp::as<MapMat>(B_);
    MapVec Zi = Rcpp::as<MapVec>(Zi_);

    thetai.noalias() = Y * Xi;
    if(r > 0)
        thetai.noalias() += B * Zi;

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd dlog_prob(p);
    int non_na = distr->dlog_prob_row(p, Ai_.begin(), thetai_.begin(), si_, gamma_.begin(), dlog_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    res.noalias() = Y.transpose() * dlog_prob;
    res /= -double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericVector grad_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, NumericVector Bj_, NumericMatrix Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_
)
{
    const int k = Yj_.length();
    const int n = X_.nrow();
    const int r = Z_.ncol();
    NumericVector res_ = NumericVector(Rcpp::no_init_vector(k));
    NumericVector thetaj_ = NumericVector(Rcpp::no_init_vector(n));

    MapVec res = Rcpp::as<MapVec>(res_);
    MapVec thetaj = Rcpp::as<MapVec>(thetaj_);
    MapMat X = Rcpp::as<MapMat>(X_);
    MapVec Yj = Rcpp::as<MapVec>(Yj_);
    MapVec Bj = Rcpp::as<MapVec>(Bj_);
    MapMat Z = Rcpp::as<MapMat>(Z_);

    thetaj.noalias() = X * Yj;
    if(r > 0)
        thetaj.noalias() += Z * Bj;

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd dlog_prob(n);
    int non_na = distr->dlog_prob_col(n, Aj_.begin(), thetaj_.begin(), s_.begin(), gammaj_, dlog_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    res.noalias() = X.transpose() * dlog_prob;
    res /= -double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericMatrix hessian_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, NumericMatrix B_, NumericVector Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_
)
{
    const int k = Xi_.length();
    const int p = Y_.nrow();
    const int r = B_.ncol();
    NumericMatrix res_ = NumericMatrix(Rcpp::no_init_matrix(k, k));
    NumericVector thetai_ = NumericVector(Rcpp::no_init_vector(p));

    MapMat res = Rcpp::as<MapMat>(res_);
    MapVec thetai = Rcpp::as<MapVec>(thetai_);
    MapVec Xi = Rcpp::as<MapVec>(Xi_);
    MapMat Y = Rcpp::as<MapMat>(Y_);
    MapMat B = Rcpp::as<MapMat>(B_);
    MapVec Zi = Rcpp::as<MapVec>(Zi_);

    thetai.noalias() = Y * Xi;
    if(r > 0)
        thetai.noalias() += B * Zi;

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd d2log_prob(p);
    int non_na = distr->d2log_prob_row(p, Ai_.begin(), thetai_.begin(), si_, gamma_.begin(), d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    res.noalias() = Y.transpose() * d2log_prob.asDiagonal() * Y;
    res /= -double(non_na);
    return res_;
}

// [[Rcpp::export]]
NumericMatrix hessian_Yj_impl(
    NumericVector Yj_, NumericMatrix X_, NumericVector Bj_, NumericMatrix Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_
)
{
    const int k = Yj_.length();
    const int n = X_.nrow();
    const int r = Z_.ncol();
    NumericMatrix res_ = NumericMatrix(Rcpp::no_init_matrix(k, k));
    NumericVector thetaj_ = NumericVector(Rcpp::no_init_vector(n));

    MapMat res = Rcpp::as<MapMat>(res_);
    MapVec thetaj = Rcpp::as<MapVec>(thetaj_);
    MapMat X = Rcpp::as<MapMat>(X_);
    MapVec Yj = Rcpp::as<MapVec>(Yj_);
    MapVec Bj = Rcpp::as<MapVec>(Bj_);
    MapMat Z = Rcpp::as<MapMat>(Z_);

    thetaj.noalias() = X * Yj;
    if(r > 0)
        thetaj.noalias() += Z * Bj;

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd d2log_prob(n);
    int non_na = distr->d2log_prob_col(n, Aj_.begin(), thetaj_.begin(), s_.begin(), gammaj_, d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    res.noalias() = X.transpose() * d2log_prob.asDiagonal() * X;
    res /= -double(non_na);
    return res_;
}

// [[Rcpp::export]]
List direction_Xi_impl(
    NumericVector Xi_, NumericMatrix Y_, NumericMatrix B_, NumericVector Zi_,
    NumericVector Ai_, Environment family_,
    double si_, NumericVector gamma_
)
{
    const int p = Y_.nrow();
    const int r = B_.ncol();
    NumericVector thetai_ = NumericVector(Rcpp::no_init_vector(p));

    MapVec thetai = Rcpp::as<MapVec>(thetai_);
    MapVec Xi = Rcpp::as<MapVec>(Xi_);
    MapMat Y = Rcpp::as<MapMat>(Y_);
    MapMat B = Rcpp::as<MapMat>(B_);
    MapVec Zi = Rcpp::as<MapVec>(Zi_);

    thetai.noalias() = Y * Xi;
    if(r > 0)
        thetai.noalias() += B * Zi;

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd dlog_prob(p), d2log_prob(p);
    int non_na = distr->d12log_prob_row(p, Ai_.begin(), thetai_.begin(), si_, gamma_.begin(),
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
    NumericVector Yj_, NumericMatrix X_, NumericVector Bj_, NumericMatrix Z_,
    NumericVector Aj_, Environment family_,
    NumericVector s_, double gammaj_
)
{
    const int n = X_.nrow();
    const int r = Z_.ncol();
    NumericVector thetaj_ = NumericVector(Rcpp::no_init_vector(n));

    MapVec thetaj = Rcpp::as<MapVec>(thetaj_);
    MapMat X = Rcpp::as<MapMat>(X_);
    MapVec Yj = Rcpp::as<MapVec>(Yj_);
    MapVec Bj = Rcpp::as<MapVec>(Bj_);
    MapMat Z = Rcpp::as<MapMat>(Z_);

    thetaj.noalias() = X * Yj;
    if(r > 0)
        thetaj.noalias() += Z * Bj;

    Rcpp::XPtr<Distribution> distr = family_["cpp_functions"];
    VectorXd dlog_prob(n), d2log_prob(n);
    int non_na = distr->d12log_prob_col(n, Aj_.begin(), thetaj_.begin(), s_.begin(), gammaj_,
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
