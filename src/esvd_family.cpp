#include "distribution.h"
#include "family2.h"

// [[Rcpp::export(.esvd_family)]]
SEXP esvd_family(std::string family)
{
    Distribution* distr;

    // Look up family string
    if(family == "poisson")
        distr = get_poisson();
    else if(family == "neg_binom")
        distr = get_neg_binom();
    else if(family == "neg_binom2")
        distr = get_neg_binom2();
    else
        Rcpp::stop("unimplemented family");

    return Rcpp::XPtr<Distribution>(distr, true);
}



// X [n x k], C [n x r], Y [p x k], Z [p x r]
// [[Rcpp::export]]
double objfn_all_r(
    MapMat XC, MapMat YZ, int k, SEXP loader, SEXP family,
    NumericVector s, NumericVector gamma,
    NumericVector l2penx, NumericVector l2peny, NumericVector l2penz
)
{
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr(family);
    MapVec l2penxv = Rcpp::as<MapVec>(l2penx);
    MapVec l2penyv = Rcpp::as<MapVec>(l2peny);
    MapVec l2penzv = Rcpp::as<MapVec>(l2penz);

    const int n = XC.rows();
    const int p = YZ.rows();
    const int r = YZ.cols() - k;
    MatrixXd theta = XC * YZ.transpose();  // [n x p]

    // Log-likelihood part
    VectorXd log_prob(n);
    std::size_t total_non_na = 0;
    double loglik = 0.0;
    for(int j = 0; j < p; j++)
    {
        int non_na = distr->log_prob_col(data_loader, j, &theta.coeffRef(0, j),
                                         s.begin(), gamma[j], log_prob.data());
        if(non_na < 1)
            Rcpp::stop("all elements in Aj are NA");
        total_non_na += static_cast<std::size_t>(non_na);
        loglik += log_prob.sum();
    }

    // Penalty part
    double penx = XC.leftCols(k).rowwise().squaredNorm().cwiseProduct(l2penxv).sum();
    double peny = YZ.leftCols(k).rowwise().squaredNorm().cwiseProduct(l2penyv).sum();
    double penz = 0.0;
    if(r > 0)
        penz += YZ.rightCols(r).rowwise().squaredNorm().cwiseProduct(l2penzv).sum();

    return (-loglik + penx + peny + penz) / double(total_non_na);
}

// [[Rcpp::export]]
double objfn_Xi_r(
    MapVec XCi, MapMat YZ, int k,
    SEXP loader, int row_ind, SEXP family,
    double si, NumericVector gamma, double l2penx
)
{
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr(family);
    MapVec gammav = Rcpp::as<MapVec>(gamma);
    return objfn_Xi(XCi, YZ, k, data_loader, row_ind, distr, si, gammav, l2penx);
}

// [[Rcpp::export]]
NumericVector grad_Xi_r(
    MapVec XCi, MapMat YZ, int k,
    SEXP loader, int row_ind, SEXP family,
    double si, NumericVector gamma, double l2penx
)
{
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr(family);
    MapVec gammav = Rcpp::as<MapVec>(gamma);
    return grad_Xi(XCi, YZ, k, data_loader, row_ind, distr, si, gammav, l2penx);
}

// [[Rcpp::export]]
NumericMatrix hessian_Xi_r(
    MapVec XCi, MapMat YZ, int k,
    SEXP loader, int row_ind, SEXP family,
    double si, NumericVector gamma, double l2penx
)
{
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr(family);
    MapVec gammav = Rcpp::as<MapVec>(gamma);
    return hessian_Xi(XCi, YZ, k, data_loader, row_ind, distr, si, gammav, l2penx);
}

// [[Rcpp::export]]
bool feas_Xi_r(MapVec XCi, MapMat YZ, SEXP family)
{
    Rcpp::XPtr<Distribution> distr(family);
    return feas_Xi(XCi, YZ, distr);
}

// [[Rcpp::export]]
double objfn_YZj_r(
    MapMat XC, MapVec YZj, int k, IntegerVector YZind,
    SEXP loader, int col_ind, SEXP family,
    NumericVector s, double gammaj, double l2peny, double l2penz
)
{
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr(family);
    MapVec sv = Rcpp::as<MapVec>(s);
    return objfn_YZj(XC, YZj, k, YZind, data_loader, col_ind, distr, sv, gammaj, l2peny, l2penz);
}

// [[Rcpp::export]]
NumericVector grad_YZj_r(
    MapMat XC, MapVec YZj, int k, IntegerVector YZind,
    SEXP loader, int col_ind, SEXP family,
    NumericVector s, double gammaj, double l2peny, double l2penz
)
{
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr(family);
    MapVec sv = Rcpp::as<MapVec>(s);
    return grad_YZj(XC, YZj, k, YZind, data_loader, col_ind, distr, sv, gammaj, l2peny, l2penz);
}

// [[Rcpp::export]]
NumericMatrix hessian_YZj_r(
    MapMat XC, MapVec YZj, int k, IntegerVector YZind,
    SEXP loader, int col_ind, SEXP family,
    NumericVector s, double gammaj, double l2peny, double l2penz
)
{
    Rcpp::XPtr<DataLoader> data_loader(loader);
    Rcpp::XPtr<Distribution> distr(family);
    MapVec sv = Rcpp::as<MapVec>(s);
    return hessian_YZj(XC, YZj, k, YZind, data_loader, col_ind, distr, sv, gammaj, l2peny, l2penz);
}

// [[Rcpp::export]]
bool feas_YZj_r(MapMat XC, MapVec YZj, SEXP family)
{
    Rcpp::XPtr<Distribution> distr(family);
    return feas_YZj(XC, YZj, distr);
}
