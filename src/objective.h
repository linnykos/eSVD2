#ifndef ESVD2_OBJECTIVE_H
#define ESVD2_OBJECTIVE_H

#include <Rcpp.h>
#include "distribution.h"
#include "family.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Environment;

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
    NumericMatrix m_Y;
    SEXP          m_B;
    SEXP          m_Zi;
    NumericVector m_Ai;
    Environment   m_family;
    const double  m_si;
    NumericVector m_gamma;

public:
    ObjectiveX(SEXP Y, SEXP B, SEXP Zi, SEXP Ai, SEXP family, double si, SEXP gamma) :
        m_Y(Y), m_B(B), m_Zi(Zi), m_Ai(Ai), m_family(family), m_si(si), m_gamma(gamma)
    {}

    double objfn(NumericVector x)
    {
        return objfn_Xi_impl(x, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma);
    }

    NumericVector grad(NumericVector x)
    {
        return grad_Xi_impl(x, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma);
    }

    NumericMatrix hessian(NumericVector x)
    {
        return hessian_Xi_impl(x, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma);
    }

    List direction(NumericVector x)
    {
        return direction_Xi_impl(x, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma);
    }

    bool feas(NumericVector x)
    {
        return feas_Xi_impl(x, m_Y, m_B, m_Zi, m_family);
    }
};

// Objective function, gradient, and Hessian for Yj
class ObjectiveY: public Objective
{
private:
    NumericMatrix m_X;
    NumericVector m_Aj;
    Environment   m_family;
    NumericVector m_s;
    const double  m_gammaj;

public:
    ObjectiveY(SEXP X, SEXP Aj, SEXP family, SEXP s, double gammaj) :
        m_X(X), m_Aj(Aj), m_family(family), m_s(s), m_gammaj(gammaj)
    {}

    double objfn(NumericVector y)
    {
        return objfn_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj);
    }

    NumericVector grad(NumericVector y)
    {
        return grad_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj);
    }

    NumericMatrix hessian(NumericVector y)
    {
        return hessian_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj);
    }

    List direction(NumericVector y)
    {
        return direction_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj);
    }

    bool feas(NumericVector y)
    {
        return feas_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_family);
    }
};

#endif  // ESVD2_OBJECTIVE_H
