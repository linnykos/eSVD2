#include "family2.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::List;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using MapMat = Eigen::Map<MatrixXd>;
using MapVec = Eigen::Map<VectorXd>;

// Penalty for X
inline double l2_penalty(const MapVec& XCi, int k, double l2penx)
{
    if(l2penx == double(0))
        return 0.0;

    return l2penx * XCi.head(k).squaredNorm();
}

// Xi [k], Ci [r], Y [p x k], Z [p x r]
double objfn_Xi(
    MapVec XCi, MapMat YZ, int k,
    DataLoader* loader, int row_ind, const Distribution* distr,
    double si, MapVec gamma, double l2penx
)
{
    const int p = YZ.rows();
    VectorXd thetai = YZ * XCi;  // [p x 1]

    VectorXd log_prob(p);
    int non_na = distr->log_prob_row(loader, row_ind, thetai.data(), si, gamma.data(), log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");
    return (-log_prob.sum() + l2_penalty(XCi, k, l2penx)) / double(non_na);
}

NumericVector grad_Xi(
    MapVec XCi, MapMat YZ, int k,
    DataLoader* loader, int row_ind, const Distribution* distr,
    double si, MapVec gamma, double l2penx
)
{
    const int p = YZ.rows();
    NumericVector res_ = NumericVector(Rcpp::no_init_vector(k));
    MapVec res = Rcpp::as<MapVec>(res_);
    VectorXd thetai = YZ * XCi;  // [p x 1]

    VectorXd dlog_prob(p);
    int non_na = distr->d12log_prob_row(loader, row_ind, thetai.data(), si, gamma.data(), dlog_prob.data(), nullptr);
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");

    res.noalias() = -YZ.leftCols(k).transpose() * dlog_prob;
    if(l2penx != double(0))
        res.noalias() += 2.0 * l2penx * XCi.head(k);
    res /= double(non_na);
    return res_;
}

NumericMatrix hessian_Xi(
    MapVec XCi, MapMat YZ, int k,
    DataLoader* loader, int row_ind, const Distribution* distr,
    double si, MapVec gamma, double l2penx
)
{
    const int p = YZ.rows();
    NumericMatrix res_ = NumericMatrix(Rcpp::no_init_matrix(k, k));
    MapMat res = Rcpp::as<MapMat>(res_);
    VectorXd thetai = YZ * XCi;  // [p x 1]

    VectorXd d2log_prob(p);
    int non_na = distr->d12log_prob_row(loader, row_ind, thetai.data(), si, gamma.data(), nullptr, d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");

    res.noalias() = -YZ.leftCols(k).transpose() * d2log_prob.asDiagonal() * YZ.leftCols(k);
    if(l2penx != double(0))
        res.diagonal().array() += 2.0 * l2penx;
    res /= double(non_na);
    return res_;
}

List direction_Xi(
    MapVec XCi, MapMat YZ, int k,
    DataLoader* loader, int row_ind, const Distribution* distr,
    double si, MapVec gamma, double l2penx
)
{
    const int p = YZ.rows();
    VectorXd thetai = YZ * XCi;  // [p x 1]

    VectorXd dlog_prob(p), d2log_prob(p);
    int non_na = distr->d12log_prob_row(loader, row_ind, thetai.data(), si, gamma.data(),
                                        dlog_prob.data(), d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Ai are NA");

    VectorXd g = YZ.leftCols(k).transpose() * (-dlog_prob);
    MatrixXd H = YZ.leftCols(k).transpose() * (-d2log_prob).asDiagonal() * YZ.leftCols(k);
    if(l2penx != double(0))
    {
        g.noalias() += 2.0 * l2penx * XCi.head(k);
        H.diagonal().array() += 2.0 * l2penx;
    }

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

bool feas_Xi(MapVec XCi, MapMat YZ, const Distribution* distr)
{
    if(distr->feas_always())
        return true;

    VectorXd thetai = YZ * XCi;  // [p x 1]
    return distr->feasibility(thetai.size(), thetai.data());
}



// Penalty for Y and Z
inline double l2_penalty(const MapVec& YZj, int k, double l2peny, double l2penz)
{
    const int r = YZj.size() - k;
    double res = 0.0;
    if(l2peny != double(0))
        res += l2peny * YZj.head(k).squaredNorm();
    if((l2penz != double(0)) && (r > 0))
        res += l2penz * YZj.tail(r).squaredNorm();
    return res;
}

// Penalty term in the gradient, for Y and Z
inline void add_l2_penalty_grad(const MapVec& YZj, int k, double l2peny, double l2penz, VectorXd& res)
{
    const int r = YZj.size() - k;
    if(l2peny != double(0))
        res.head(k).noalias() += 2.0 * l2peny * YZj.head(k);
    if((l2penz != double(0)) && (r > 0))
        res.tail(r).noalias() += 2.0 * l2penz * YZj.tail(r);
}

// Penalty term in the Hessian matrix, for Y and Z
inline void add_l2_penalty_hessian(const MapVec& YZj, int k, double l2peny, double l2penz, MatrixXd& res)
{
    const int r = YZj.size() - k;
    if(l2peny != double(0))
        res.diagonal().array().head(k) += 2.0 * l2peny;
    if((l2penz != double(0)) && (r > 0))
        res.diagonal().array().tail(r) += 2.0 * l2penz;
}

// vec[ind], ind is 0-based
inline NumericVector subset_vector(const VectorXd& vec, IntegerVector ind)
{
    const int m = ind.length();
    NumericVector res(m);
    const double* src = vec.data();
    double* dest = res.begin();
    for(int i = 0; i < m; i++)
        dest[i] = src[ind[i]];
    return res;
}

// mat[ind, ind], ind is 0-based
inline NumericMatrix subset_matrix(const MatrixXd& mat, IntegerVector ind)
{
    const int n = mat.rows();
    const int m = ind.length();
    NumericMatrix res(m, m);
    for(int j = 0; j < m; j++)
    {
        const double* src_head = mat.data() + ind[j] * n;
        double* dest_head = res.begin() + j * m;
        for(int i = 0; i < m; i++)
        {
            dest_head[i] = src_head[ind[i]];
        }
    }
    return res;
}

// X [n x k], C [n x r], Yj [k], Zj [r]
// YZind is the index vector indicating which variables in YZj to update, 0-based
double objfn_YZj(
    MapMat XC, MapVec YZj, int k, IntegerVector YZind,
    DataLoader* loader, int col_ind, const Distribution* distr,
    MapVec s, double gammaj, double l2peny, double l2penz
)
{
    const int n = XC.rows();
    VectorXd thetaj = XC * YZj;  // [n x 1]

    VectorXd log_prob(n);
    int non_na = distr->log_prob_col(loader, col_ind, thetaj.data(), s.data(), gammaj, log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");
    return (-log_prob.sum() + l2_penalty(YZj, k, l2peny, l2penz)) / double(non_na);
}

NumericVector grad_YZj(
    MapMat XC, MapVec YZj, int k, IntegerVector YZind,
    DataLoader* loader, int col_ind, const Distribution* distr,
    MapVec s, double gammaj, double l2peny, double l2penz
)
{
    const int n = XC.rows();
    VectorXd thetaj = XC * YZj;  // [n x 1]

    VectorXd dlog_prob(n);
    int non_na = distr->d12log_prob_col(loader, col_ind, thetaj.data(), s.data(), gammaj, dlog_prob.data(), nullptr);
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");

    VectorXd res = -XC.transpose() * dlog_prob;
    add_l2_penalty_grad(YZj, k, l2peny, l2penz, res);
    res /= double(non_na);
    return subset_vector(res, YZind);
}

NumericMatrix hessian_YZj(
    MapMat XC, MapVec YZj, int k, IntegerVector YZind,
    DataLoader* loader, int col_ind, const Distribution* distr,
    MapVec s, double gammaj, double l2peny, double l2penz
)
{
    const int n = XC.rows();
    VectorXd thetaj = XC * YZj;  // [n x 1]

    VectorXd d2log_prob(n);
    int non_na = distr->d12log_prob_col(loader, col_ind, thetaj.data(), s.data(), gammaj, nullptr, d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");

    MatrixXd res = -XC.transpose() * d2log_prob.asDiagonal() * XC;
    add_l2_penalty_hessian(YZj, k, l2peny, l2penz, res);
    res /= double(non_na);
    return subset_matrix(res, YZind);
}

List direction_YZj(
    MapMat XC, MapVec YZj, int k, IntegerVector YZind,
    DataLoader* loader, int col_ind, const Distribution* distr,
    MapVec s, double gammaj, double l2peny, double l2penz
)
{
    const int n = XC.rows();
    VectorXd thetaj = XC * YZj;  // [n x 1]

    VectorXd dlog_prob(n), d2log_prob(n);
    int non_na = distr->d12log_prob_col(loader, col_ind, thetaj.data(), s.data(), gammaj,
                                        dlog_prob.data(), d2log_prob.data());
    if(non_na < 1)
        Rcpp::stop("all elements in Aj are NA");

    VectorXd g = XC.transpose() * (-dlog_prob);
    add_l2_penalty_grad(YZj, k, l2peny, l2penz, g);
    MatrixXd H = XC.transpose() * (-d2log_prob).asDiagonal() * XC;
    add_l2_penalty_hessian(YZj, k, l2peny, l2penz, H);

    NumericVector gsub = subset_vector(g, YZind);
    NumericMatrix Hsub = subset_matrix(H, YZind);

    VectorXd direc(YZind.length());
    Eigen::LLT<MatrixXd> solver(Rcpp::as<MapMat>(Hsub));
    if(solver.info() != Eigen::Success)
    {
        // Fall back to gradient direction if Hessian is singular
        direc.noalias() = -Rcpp::as<MapVec>(gsub) / double(non_na);
    } else {
        direc.noalias() = -solver.solve(Rcpp::as<MapVec>(gsub));
    }
    return List::create(
        Rcpp::Named("grad") = gsub / double(non_na),
        Rcpp::Named("direction") = direc
    );
}

bool feas_YZj(MapMat XC, MapVec YZj, const Distribution* distr)
{
    if(distr->feas_always())
        return true;

    VectorXd thetaj = XC * YZj;  // [n x 1]
    return distr->feasibility(thetaj.size(), thetaj.data());
}
