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
    const double  m_offseti;

public:
    ObjectiveX(SEXP Y, SEXP B, SEXP Zi, SEXP Ai, SEXP family, double si, SEXP gamma, double offseti) :
        m_Y(Y), m_B(B), m_Zi(Zi), m_Ai(Ai), m_family(family), m_si(si), m_gamma(gamma), m_offseti(offseti)
    {}

    double objfn(NumericVector x)
    {
        return objfn_Xi_impl(x, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma, m_offseti);
    }

    NumericVector grad(NumericVector x)
    {
        return grad_Xi_impl(x, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma, m_offseti);
    }

    NumericMatrix hessian(NumericVector x)
    {
        return hessian_Xi_impl(x, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma, m_offseti);
    }

    List direction(NumericVector x)
    {
        return direction_Xi_impl(x, m_Y, m_B, m_Zi, m_Ai, m_family, m_si, m_gamma, m_offseti);
    }

    bool feas(NumericVector x)
    {
        return feas_Xi_impl(x, m_Y, m_B, m_Zi, m_family, m_offseti);
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
    NumericVector m_offset;

public:
    ObjectiveY(SEXP X, SEXP Aj, SEXP family, SEXP s, double gammaj, SEXP offset) :
        m_X(X), m_Aj(Aj), m_family(family), m_s(s), m_gammaj(gammaj), m_offset(offset)
    {}

    double objfn(NumericVector y)
    {
        return objfn_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj, m_offset);
    }

    NumericVector grad(NumericVector y)
    {
        return grad_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj, m_offset);
    }

    NumericMatrix hessian(NumericVector y)
    {
        return hessian_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj, m_offset);
    }

    List direction(NumericVector y)
    {
        return direction_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_Aj, m_family, m_s, m_gammaj, m_offset);
    }

    bool feas(NumericVector y)
    {
        return feas_Yj_impl(y, m_X, R_NilValue, R_NilValue, m_family, m_offset);
    }
};

#endif  // ESVD2_OBJECTIVE_H
