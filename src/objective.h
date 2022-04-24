#ifndef ESVD2_OBJECTIVE_H
#define ESVD2_OBJECTIVE_H

#include <Rcpp.h>

class Objective
{
public:
    virtual ~Objective() {}
    virtual double              objfn(Rcpp::NumericVector x) = 0;
    virtual Rcpp::NumericVector grad(Rcpp::NumericVector x) = 0;
    virtual Rcpp::NumericMatrix hessian(Rcpp::NumericVector x) = 0;
    virtual Rcpp::List          direction(Rcpp::NumericVector x) = 0;
    virtual bool                feas(Rcpp::NumericVector x) = 0;
};


#endif  // ESVD2_OBJECTIVE_H
