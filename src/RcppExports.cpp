// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/eSVD2.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// data_loader
SEXP data_loader(SEXP mat);
RcppExport SEXP _eSVD2_data_loader(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(data_loader(mat));
    return rcpp_result_gen;
END_RCPP
}
// data_loader_description
void data_loader_description(SEXP loader_);
RcppExport SEXP _eSVD2_data_loader_description(SEXP loader_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type loader_(loader_SEXP);
    data_loader_description(loader_);
    return R_NilValue;
END_RCPP
}
// test_data_loader
void test_data_loader(SEXP mat);
RcppExport SEXP _eSVD2_test_data_loader(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mat(matSEXP);
    test_data_loader(mat);
    return R_NilValue;
END_RCPP
}
// esvd_family
List esvd_family(std::string family);
RcppExport SEXP _eSVD2_esvd_family(SEXP familySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    rcpp_result_gen = Rcpp::wrap(esvd_family(family));
    return rcpp_result_gen;
END_RCPP
}
// objfn_all_r
double objfn_all_r(MapMat XC, MapMat YZ, int k, SEXP loader, List family, NumericVector s, NumericVector gamma, double l2penx, double l2peny, double l2penz);
RcppExport SEXP _eSVD2_objfn_all_r(SEXP XCSEXP, SEXP YZSEXP, SEXP kSEXP, SEXP loaderSEXP, SEXP familySEXP, SEXP sSEXP, SEXP gammaSEXP, SEXP l2penxSEXP, SEXP l2penySEXP, SEXP l2penzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapMat >::type XC(XCSEXP);
    Rcpp::traits::input_parameter< MapMat >::type YZ(YZSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type l2penx(l2penxSEXP);
    Rcpp::traits::input_parameter< double >::type l2peny(l2penySEXP);
    Rcpp::traits::input_parameter< double >::type l2penz(l2penzSEXP);
    rcpp_result_gen = Rcpp::wrap(objfn_all_r(XC, YZ, k, loader, family, s, gamma, l2penx, l2peny, l2penz));
    return rcpp_result_gen;
END_RCPP
}
// objfn_Xi_r
double objfn_Xi_r(MapVec XCi, MapMat YZ, int k, SEXP loader, int row_ind, List family, double si, NumericVector gamma, double l2penx);
RcppExport SEXP _eSVD2_objfn_Xi_r(SEXP XCiSEXP, SEXP YZSEXP, SEXP kSEXP, SEXP loaderSEXP, SEXP row_indSEXP, SEXP familySEXP, SEXP siSEXP, SEXP gammaSEXP, SEXP l2penxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapVec >::type XCi(XCiSEXP);
    Rcpp::traits::input_parameter< MapMat >::type YZ(YZSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< int >::type row_ind(row_indSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< double >::type si(siSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type l2penx(l2penxSEXP);
    rcpp_result_gen = Rcpp::wrap(objfn_Xi_r(XCi, YZ, k, loader, row_ind, family, si, gamma, l2penx));
    return rcpp_result_gen;
END_RCPP
}
// grad_Xi_r
NumericVector grad_Xi_r(MapVec XCi, MapMat YZ, int k, SEXP loader, int row_ind, List family, double si, NumericVector gamma, double l2penx);
RcppExport SEXP _eSVD2_grad_Xi_r(SEXP XCiSEXP, SEXP YZSEXP, SEXP kSEXP, SEXP loaderSEXP, SEXP row_indSEXP, SEXP familySEXP, SEXP siSEXP, SEXP gammaSEXP, SEXP l2penxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapVec >::type XCi(XCiSEXP);
    Rcpp::traits::input_parameter< MapMat >::type YZ(YZSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< int >::type row_ind(row_indSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< double >::type si(siSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type l2penx(l2penxSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_Xi_r(XCi, YZ, k, loader, row_ind, family, si, gamma, l2penx));
    return rcpp_result_gen;
END_RCPP
}
// hessian_Xi_r
NumericMatrix hessian_Xi_r(MapVec XCi, MapMat YZ, int k, SEXP loader, int row_ind, List family, double si, NumericVector gamma, double l2penx);
RcppExport SEXP _eSVD2_hessian_Xi_r(SEXP XCiSEXP, SEXP YZSEXP, SEXP kSEXP, SEXP loaderSEXP, SEXP row_indSEXP, SEXP familySEXP, SEXP siSEXP, SEXP gammaSEXP, SEXP l2penxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapVec >::type XCi(XCiSEXP);
    Rcpp::traits::input_parameter< MapMat >::type YZ(YZSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< int >::type row_ind(row_indSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< double >::type si(siSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type l2penx(l2penxSEXP);
    rcpp_result_gen = Rcpp::wrap(hessian_Xi_r(XCi, YZ, k, loader, row_ind, family, si, gamma, l2penx));
    return rcpp_result_gen;
END_RCPP
}
// feas_Xi_r
bool feas_Xi_r(MapVec XCi, MapMat YZ, List family);
RcppExport SEXP _eSVD2_feas_Xi_r(SEXP XCiSEXP, SEXP YZSEXP, SEXP familySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapVec >::type XCi(XCiSEXP);
    Rcpp::traits::input_parameter< MapMat >::type YZ(YZSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    rcpp_result_gen = Rcpp::wrap(feas_Xi_r(XCi, YZ, family));
    return rcpp_result_gen;
END_RCPP
}
// objfn_YZj_r
double objfn_YZj_r(MapMat XC, MapVec YZj, int k, IntegerVector YZind, SEXP loader, int col_ind, List family, NumericVector s, double gammaj, double l2peny, double l2penz);
RcppExport SEXP _eSVD2_objfn_YZj_r(SEXP XCSEXP, SEXP YZjSEXP, SEXP kSEXP, SEXP YZindSEXP, SEXP loaderSEXP, SEXP col_indSEXP, SEXP familySEXP, SEXP sSEXP, SEXP gammajSEXP, SEXP l2penySEXP, SEXP l2penzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapMat >::type XC(XCSEXP);
    Rcpp::traits::input_parameter< MapVec >::type YZj(YZjSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YZind(YZindSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< int >::type col_ind(col_indSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type gammaj(gammajSEXP);
    Rcpp::traits::input_parameter< double >::type l2peny(l2penySEXP);
    Rcpp::traits::input_parameter< double >::type l2penz(l2penzSEXP);
    rcpp_result_gen = Rcpp::wrap(objfn_YZj_r(XC, YZj, k, YZind, loader, col_ind, family, s, gammaj, l2peny, l2penz));
    return rcpp_result_gen;
END_RCPP
}
// grad_YZj_r
NumericVector grad_YZj_r(MapMat XC, MapVec YZj, int k, IntegerVector YZind, SEXP loader, int col_ind, List family, NumericVector s, double gammaj, double l2peny, double l2penz);
RcppExport SEXP _eSVD2_grad_YZj_r(SEXP XCSEXP, SEXP YZjSEXP, SEXP kSEXP, SEXP YZindSEXP, SEXP loaderSEXP, SEXP col_indSEXP, SEXP familySEXP, SEXP sSEXP, SEXP gammajSEXP, SEXP l2penySEXP, SEXP l2penzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapMat >::type XC(XCSEXP);
    Rcpp::traits::input_parameter< MapVec >::type YZj(YZjSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YZind(YZindSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< int >::type col_ind(col_indSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type gammaj(gammajSEXP);
    Rcpp::traits::input_parameter< double >::type l2peny(l2penySEXP);
    Rcpp::traits::input_parameter< double >::type l2penz(l2penzSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_YZj_r(XC, YZj, k, YZind, loader, col_ind, family, s, gammaj, l2peny, l2penz));
    return rcpp_result_gen;
END_RCPP
}
// hessian_YZj_r
NumericMatrix hessian_YZj_r(MapMat XC, MapVec YZj, int k, IntegerVector YZind, SEXP loader, int col_ind, List family, NumericVector s, double gammaj, double l2peny, double l2penz);
RcppExport SEXP _eSVD2_hessian_YZj_r(SEXP XCSEXP, SEXP YZjSEXP, SEXP kSEXP, SEXP YZindSEXP, SEXP loaderSEXP, SEXP col_indSEXP, SEXP familySEXP, SEXP sSEXP, SEXP gammajSEXP, SEXP l2penySEXP, SEXP l2penzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapMat >::type XC(XCSEXP);
    Rcpp::traits::input_parameter< MapVec >::type YZj(YZjSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YZind(YZindSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< int >::type col_ind(col_indSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type gammaj(gammajSEXP);
    Rcpp::traits::input_parameter< double >::type l2peny(l2penySEXP);
    Rcpp::traits::input_parameter< double >::type l2penz(l2penzSEXP);
    rcpp_result_gen = Rcpp::wrap(hessian_YZj_r(XC, YZj, k, YZind, loader, col_ind, family, s, gammaj, l2peny, l2penz));
    return rcpp_result_gen;
END_RCPP
}
// feas_YZj_r
bool feas_YZj_r(MapMat XC, MapVec YZj, List family);
RcppExport SEXP _eSVD2_feas_YZj_r(SEXP XCSEXP, SEXP YZjSEXP, SEXP familySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MapMat >::type XC(XCSEXP);
    Rcpp::traits::input_parameter< MapVec >::type YZj(YZjSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    rcpp_result_gen = Rcpp::wrap(feas_YZj_r(XC, YZj, family));
    return rcpp_result_gen;
END_RCPP
}
// gamma_rate
double gamma_rate(NumericVector x, NumericVector mu, NumericVector s);
RcppExport SEXP _eSVD2_gamma_rate(SEXP xSEXP, SEXP muSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(gamma_rate(x, mu, s));
    return rcpp_result_gen;
END_RCPP
}
// log_gamma_rate
double log_gamma_rate(NumericVector x, NumericVector mu, NumericVector s, double lower, double upper);
RcppExport SEXP _eSVD2_log_gamma_rate(SEXP xSEXP, SEXP muSEXP, SEXP sSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(log_gamma_rate(x, mu, s, lower, upper));
    return rcpp_result_gen;
END_RCPP
}
// opt_x
NumericMatrix opt_x(NumericMatrix XC0, MapMat YZ, int k, SEXP loader, List family, NumericVector s, NumericVector gamma, double l2penx, int verbose, bool inplace);
RcppExport SEXP _eSVD2_opt_x(SEXP XC0SEXP, SEXP YZSEXP, SEXP kSEXP, SEXP loaderSEXP, SEXP familySEXP, SEXP sSEXP, SEXP gammaSEXP, SEXP l2penxSEXP, SEXP verboseSEXP, SEXP inplaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type XC0(XC0SEXP);
    Rcpp::traits::input_parameter< MapMat >::type YZ(YZSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type l2penx(l2penxSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type inplace(inplaceSEXP);
    rcpp_result_gen = Rcpp::wrap(opt_x(XC0, YZ, k, loader, family, s, gamma, l2penx, verbose, inplace));
    return rcpp_result_gen;
END_RCPP
}
// opt_yz
NumericMatrix opt_yz(NumericMatrix YZ0, MapMat XC, int k, IntegerVector YZind, SEXP loader, List family, NumericVector s, NumericVector gamma, double l2peny, double l2penz, int verbose, bool inplace);
RcppExport SEXP _eSVD2_opt_yz(SEXP YZ0SEXP, SEXP XCSEXP, SEXP kSEXP, SEXP YZindSEXP, SEXP loaderSEXP, SEXP familySEXP, SEXP sSEXP, SEXP gammaSEXP, SEXP l2penySEXP, SEXP l2penzSEXP, SEXP verboseSEXP, SEXP inplaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type YZ0(YZ0SEXP);
    Rcpp::traits::input_parameter< MapMat >::type XC(XCSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YZind(YZindSEXP);
    Rcpp::traits::input_parameter< SEXP >::type loader(loaderSEXP);
    Rcpp::traits::input_parameter< List >::type family(familySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type l2peny(l2penySEXP);
    Rcpp::traits::input_parameter< double >::type l2penz(l2penzSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type inplace(inplaceSEXP);
    rcpp_result_gen = Rcpp::wrap(opt_yz(YZ0, XC, k, YZind, loader, family, s, gamma, l2peny, l2penz, verbose, inplace));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_eSVD2_data_loader", (DL_FUNC) &_eSVD2_data_loader, 1},
    {"_eSVD2_data_loader_description", (DL_FUNC) &_eSVD2_data_loader_description, 1},
    {"_eSVD2_test_data_loader", (DL_FUNC) &_eSVD2_test_data_loader, 1},
    {"_eSVD2_esvd_family", (DL_FUNC) &_eSVD2_esvd_family, 1},
    {"_eSVD2_objfn_all_r", (DL_FUNC) &_eSVD2_objfn_all_r, 10},
    {"_eSVD2_objfn_Xi_r", (DL_FUNC) &_eSVD2_objfn_Xi_r, 9},
    {"_eSVD2_grad_Xi_r", (DL_FUNC) &_eSVD2_grad_Xi_r, 9},
    {"_eSVD2_hessian_Xi_r", (DL_FUNC) &_eSVD2_hessian_Xi_r, 9},
    {"_eSVD2_feas_Xi_r", (DL_FUNC) &_eSVD2_feas_Xi_r, 3},
    {"_eSVD2_objfn_YZj_r", (DL_FUNC) &_eSVD2_objfn_YZj_r, 11},
    {"_eSVD2_grad_YZj_r", (DL_FUNC) &_eSVD2_grad_YZj_r, 11},
    {"_eSVD2_hessian_YZj_r", (DL_FUNC) &_eSVD2_hessian_YZj_r, 11},
    {"_eSVD2_feas_YZj_r", (DL_FUNC) &_eSVD2_feas_YZj_r, 3},
    {"_eSVD2_gamma_rate", (DL_FUNC) &_eSVD2_gamma_rate, 3},
    {"_eSVD2_log_gamma_rate", (DL_FUNC) &_eSVD2_log_gamma_rate, 5},
    {"_eSVD2_opt_x", (DL_FUNC) &_eSVD2_opt_x, 10},
    {"_eSVD2_opt_yz", (DL_FUNC) &_eSVD2_opt_yz, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_eSVD2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
