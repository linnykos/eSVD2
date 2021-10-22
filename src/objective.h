#ifndef ESVD2_OBJECTIVE_H
#define ESVD2_OBJECTIVE_H

#include <RcppEigen.h>
#include "distribution.h"
#include "family.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Environment;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;

class Objective
{
public:
    virtual ~Objective() {}
    virtual double        objfn(NumericVector x) = 0;
    virtual NumericVector grad(NumericVector x) = 0;
    virtual NumericMatrix hessian(NumericVector x) = 0;
    virtual List          direction(NumericVector x) = 0;
    virtual bool          feas(NumericVector x) = 0;
};

// Objective function, gradient, and Hessian for Xi
class ObjectiveX: public Objective
{
private:
    MapMat       m_Y;
    SEXP         m_B;
    SEXP         m_Zi;
    MapVec       m_Ai;
    Environment  m_family;
    const double m_si;
    MapVec       m_gamma;
    const double m_offseti;
    const double m_l2pen;

public:
    ObjectiveX(MapMat Y, SEXP B, SEXP Zi, MapVec Ai, Environment family, double si, MapVec gamma, double offseti, double l2pen) :
        m_Y(Y), m_B(B), m_Zi(Zi), m_Ai(Ai), m_family(family), m_si(si), m_gamma(gamma), m_offseti(offseti), m_l2pen(l2pen)
    {}

    double objfn(NumericVector x)
    {
        MapVec Xi = Rcpp::as<MapVec>(x);
        return objfn_Xi_impl(Xi, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma, m_offseti, m_l2pen);
    }

    NumericVector grad(NumericVector x)
    {
        MapVec Xi = Rcpp::as<MapVec>(x);
        return grad_Xi_impl(Xi, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma, m_offseti, m_l2pen);
    }

    NumericMatrix hessian(NumericVector x)
    {
        MapVec Xi = Rcpp::as<MapVec>(x);
        return hessian_Xi_impl(Xi, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma, m_offseti, m_l2pen);
    }

    List direction(NumericVector x)
    {
        MapVec Xi = Rcpp::as<MapVec>(x);
        return direction_Xi_impl(Xi, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma, m_offseti, m_l2pen);
    }

    bool feas(NumericVector x)
    {
        MapVec Xi = Rcpp::as<MapVec>(x);
        return feas_Xi_impl(Xi, m_Y, m_B, m_Zi, m_family, m_offseti);
    }
};

// Objective function, gradient, and Hessian for Yj
class ObjectiveY: public Objective
{
private:
    MapMat       m_X;
    MapVec       m_Aj;
    Environment  m_family;
    MapVec       m_s;
    const double m_gammaj;
    MapVec       m_offset;
    const double m_l2pen;

public:
    ObjectiveY(MapMat X, MapVec Aj, Environment family, MapVec s, double gammaj, MapVec offset, double l2pen) :
        m_X(X), m_Aj(Aj), m_family(family), m_s(s), m_gammaj(gammaj), m_offset(offset), m_l2pen(l2pen)
    {}

    double objfn(NumericVector y)
    {
        MapVec Yj = Rcpp::as<MapVec>(y);
        return objfn_Yj_impl(Yj, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj, m_offset, m_l2pen);
    }

    NumericVector grad(NumericVector y)
    {
        MapVec Yj = Rcpp::as<MapVec>(y);
        return grad_Yj_impl(Yj, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj, m_offset, m_l2pen);
    }

    NumericMatrix hessian(NumericVector y)
    {
        MapVec Yj = Rcpp::as<MapVec>(y);
        return hessian_Yj_impl(Yj, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj, m_offset, m_l2pen);
    }

    List direction(NumericVector y)
    {
        MapVec Yj = Rcpp::as<MapVec>(y);
        return direction_Yj_impl(Yj, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj, m_offset, m_l2pen);
    }

    bool feas(NumericVector y)
    {
        MapVec Yj = Rcpp::as<MapVec>(y);
        return feas_Yj_impl(Yj, m_X, R_NilValue, R_NilValue, m_family, m_offset);
    }
};

#endif  // ESVD2_OBJECTIVE_H
