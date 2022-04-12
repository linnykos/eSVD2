#ifndef ESVD2_OBJECTIVE_H
#define ESVD2_OBJECTIVE_H

#include <Rcpp.h>

class Objective
{
private:
    using NumericVector = Rcpp::NumericVector;
    using NumericMatrix = Rcpp::NumericMatrix;
    using List = Rcpp::List;

public:
    virtual ~Objective() {}
    virtual double        objfn(NumericVector x) = 0;
    virtual NumericVector grad(NumericVector x) = 0;
    virtual NumericMatrix hessian(NumericVector x) = 0;
    virtual List          direction(NumericVector x) = 0;
    virtual bool          feas(NumericVector x) = 0;
};


#endif  // ESVD2_OBJECTIVE_H
