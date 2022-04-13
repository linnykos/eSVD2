#include <RcppEigen.h>
#include "objective.h"
#include "family.h"

using Rcpp::NumericVector;
using Rcpp::List;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using MapMat = Eigen::Map<MatrixXd>;
using MapVec = Eigen::Map<VectorXd>;

// Defined in constrained_newton.cpp
List constr_newton(
    const VectorXd& x0, Objective& objective, int max_iter = 100,
    int max_linesearch = 30, double eps_rel = 1e-5, bool verbose = false
);

// vec[ind], ind is 0-based
inline VectorXd subset_vector(const MapVec& vec, IntegerVector ind)
{
    const int m = ind.length();
    VectorXd res(m);
    const double* src = vec.data();
    double* dest = res.data();
    for(int i = 0; i < m; i++)
        dest[i] = src[ind[i]];
    return res;
}

// res[ind] <- vec, ind is 0-based
inline void subset_vector_assign(MapVec res, NumericVector vec, IntegerVector ind)
{
    const int m = ind.length();
    const double* src = vec.begin();
    double* dest = res.data();
    for(int i = 0; i < m; i++)
        dest[ind[i]] = src[i];
}

// Objective function, gradient, and Hessian for Xi
class ObjectiveX: public Objective
{
private:
    MapVec              m_XCi;
    const MapMat&       m_YZ;
    const int           m_k;
    DataLoader*         m_loader;
    const int           m_row_ind;
    const Distribution* m_distr;
    const double        m_si;
    const MapVec&       m_gamma;
    const double        m_l2penx;

public:
    ObjectiveX(
        MapVec XCi, const MapMat& YZ, int k,
        DataLoader* loader, int row_ind, const Distribution* distr,
        double si, const MapVec& gamma, double l2penx
    ) :
        m_XCi(XCi), m_YZ(YZ), m_k(k), m_loader(loader), m_row_ind(row_ind),
        m_distr(distr), m_si(si), m_gamma(gamma), m_l2penx(l2penx)
    {}

    double objfn(NumericVector x) override
    {
        // Write x to the XCi vector
        std::copy(x.begin(), x.end(), m_XCi.data());
        return objfn_Xi(m_XCi, m_YZ, m_k, m_loader, m_row_ind,
            m_distr, m_si, m_gamma, m_l2penx);
    }

    NumericVector grad(NumericVector x)
    {
        // Write x to the XCi vector
        std::copy(x.begin(), x.end(), m_XCi.data());
        return grad_Xi(m_XCi, m_YZ, m_k, m_loader, m_row_ind,
            m_distr, m_si, m_gamma, m_l2penx);
    }

    NumericMatrix hessian(NumericVector x)
    {
        // Write x to the XCi vector
        std::copy(x.begin(), x.end(), m_XCi.data());
        return hessian_Xi(m_XCi, m_YZ, m_k, m_loader, m_row_ind,
            m_distr, m_si, m_gamma, m_l2penx);
    }

    List direction(NumericVector x)
    {
        // Write x to the XCi vector
        std::copy(x.begin(), x.end(), m_XCi.data());
        return direction_Xi(m_XCi, m_YZ, m_k, m_loader, m_row_ind,
            m_distr, m_si, m_gamma, m_l2penx);
    }

    bool feas(NumericVector x)
    {
        // Write x to the XCi vector
        std::copy(x.begin(), x.end(), m_XCi.data());
        return feas_Xi(m_XCi, m_YZ, m_distr);
    }
};

// Objective function, gradient, and Hessian for Yj
class ObjectiveYZ: public Objective
{
private:
    const MapMat&       m_XC;
    MapVec              m_YZj;
    const int           m_k;
    IntegerVector       m_YZind;
    DataLoader*         m_loader;
    const int           m_col_ind;
    const Distribution* m_distr;
    const MapVec&       m_s;
    double              m_gammaj;
    const double        m_l2peny;
    const double        m_l2penz;

public:
    ObjectiveYZ(
        const MapMat& XC, MapVec YZj, int k, IntegerVector YZind,
        DataLoader* loader, int col_ind, const Distribution* distr,
        const MapVec& s, double gammaj, double l2peny, double l2penz
    ) :
        m_XC(XC), m_YZj(YZj), m_k(k), m_YZind(YZind),
        m_loader(loader), m_col_ind(col_ind), m_distr(distr),
        m_s(s), m_gammaj(gammaj), m_l2peny(l2peny), m_l2penz(l2penz)
    {}

    double objfn(NumericVector yz)
    {
        // Fill in the given vector
        subset_vector_assign(m_YZj, yz, m_YZind);
        return objfn_YZj(
            m_XC, m_YZj, m_k, m_YZind, m_loader, m_col_ind,
            m_distr, m_s, m_gammaj, m_l2peny, m_l2penz);
    }

    NumericVector grad(NumericVector yz)
    {
        // Fill in the given vector
        subset_vector_assign(m_YZj, yz, m_YZind);
        return grad_YZj(
            m_XC, m_YZj, m_k, m_YZind, m_loader, m_col_ind,
            m_distr, m_s, m_gammaj, m_l2peny, m_l2penz);
    }

    NumericMatrix hessian(NumericVector yz)
    {
        // Fill in the given vector
        subset_vector_assign(m_YZj, yz, m_YZind);
        return hessian_YZj(
            m_XC, m_YZj, m_k, m_YZind, m_loader, m_col_ind,
            m_distr, m_s, m_gammaj, m_l2peny, m_l2penz);
    }

    List direction(NumericVector yz)
    {
        // Fill in the given vector
        subset_vector_assign(m_YZj, yz, m_YZind);
        return direction_YZj(
            m_XC, m_YZj, m_k, m_YZind, m_loader, m_col_ind,
            m_distr, m_s, m_gammaj, m_l2peny, m_l2penz);
    }

    bool feas(NumericVector yz)
    {
        // Fill in the given vector
        subset_vector_assign(m_YZj, yz, m_YZind);
        return feas_YZj(m_XC, m_YZj, m_distr);
    }
};

// [[Rcpp::export]]
NumericMatrix opt_x_cpp(
    NumericMatrix XC0, MapMat YZ, int k, SEXP loader, List family,
    NumericVector s, NumericVector gamma,
    NumericVector l2penx, int verbose = 0, bool inplace = true
)
{
    MapMat init = Rcpp::as<MapMat>(XC0);
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr = family["internal"];
    MapVec gammav = Rcpp::as<MapVec>(gamma);

    // Get dimensions
    const int n = XC0.nrow();
    const int kr = XC0.ncol();

    // Make a copy of the XC0 matrix and transpose it
    MatrixXd XC = init.transpose();

    // Optimize each column of XC [(k+r) x n]
    for(int i = 0; i < n; i++)
    {
        if(verbose >= 2)
            Rcpp::Rcout << "===== Optimizing Row " << i + 1 << " of X =====" << std::endl;

        // XCi vector
        MapVec XCi(&XC.coeffRef(0, i), kr);

        // Objective function, gradient, Hessian, and feasibility
        ObjectiveX obj(XCi, YZ, k, data_loader, i, distr, s[i], gammav, l2penx[i]);

        // Run optimizer
        VectorXd Xi_init = XCi.head(k);
        List opt = constr_newton(Xi_init, obj, 100, 30, 0.001, (verbose >= 3));

        // Extract result
        NumericVector optx = opt["x"];
        std::copy(optx.begin(), optx.end(), XCi.data());

        if(verbose >= 3)
            Rcpp::Rcout << "==========" << std::endl << std::endl;
    }

    // Directly modify initial values
    if(inplace)
    {
        init.noalias() = XC.transpose();
        return XC0;
    }
    // Otherwise, create a new matrix
    NumericMatrix res = NumericMatrix(Rcpp::no_init_matrix(n, kr));
    MapMat resm = Rcpp::as<MapMat>(res);
    resm.noalias() = XC.transpose();
    return res;
}

// [[Rcpp::export]]
NumericMatrix opt_yb_cpp(
    NumericMatrix YZ0, MapMat XC, int k, IntegerVector YZind,
    SEXP loader, List family, NumericVector s, NumericVector gamma,
    NumericVector l2peny, NumericVector l2penz,
    int verbose = 0, bool inplace = true
)
{
    MapMat init = Rcpp::as<MapMat>(YZ0);
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr = family["internal"];
    MapVec sv = Rcpp::as<MapVec>(s);

    // Get dimensions
    const int p = YZ0.rows();
    const int kr = YZ0.cols();

    // Make a copy of the YZ0 matrix and transpose it
    MatrixXd YZ = init.transpose();

    // Optimize each row of YZ [(k+r) x p]
    for(int j = 0; j < p; j++)
    {
        if(verbose >= 2)
            Rcpp::Rcout << "===== Optimizing Row " << j + 1 << " of Y =====" << std::endl;

        // YZj vector
        MapVec YZj(&YZ.coeffRef(0, j), kr);

        // Objective function, gradient, Hessian, and feasibility
        ObjectiveYZ obj(XC, YZj, k, YZind, data_loader, j, distr,
                        sv, gamma[j], l2peny[j], l2penz[j]);

        // Run optimizer
        VectorXd YZj_init = subset_vector(YZj, YZind);
        List opt = constr_newton(YZj_init, obj, 100, 30, 0.001, (verbose >= 3));

        // Extract result
        NumericVector optx = opt["x"];
        subset_vector_assign(YZj, optx, YZind);

        if(verbose >= 3)
            Rcpp::Rcout << "==========" << std::endl << std::endl;
    }

    // Directly modify initial values
    if(inplace)
    {
        init.noalias() = YZ.transpose();
        return YZ0;
    }
    // Otherwise, create a new matrix
    NumericMatrix res = NumericMatrix(Rcpp::no_init_matrix(p, kr));
    MapMat resm = Rcpp::as<MapMat>(res);
    resm.noalias() = YZ.transpose();
    return res;
}
