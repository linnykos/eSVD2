#include <Rcpp.h>
#include "distribution.h"

// Distribution: exponential
// Canonical parameter: lambda (mean = lambda)
// Natural parameter: theta = -1/lambda
// Log-likelihood: A * theta + s * log(-theta)

class Exponential: public Distribution
{
public:
    // Log-density for a single data point Aij
    inline double log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        return Aij * thetaij + si * std::log(-thetaij);
    }

    // Special case of Aij==0
    inline double log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        return si * std::log(-thetaij);
    }

    // 1st and 2nd derivatives of the log-density function w.r.t. thetaij
    inline void d12log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        if(compute_d1) *res1 = Aij + si / thetaij;
        if(compute_d2) *res2 = -si / (thetaij * thetaij);
    }

    // Special case of Aij==0
    inline void d12log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        if(compute_d1) *res1 = si / thetaij;
        if(compute_d2) *res2 = -si / (thetaij * thetaij);
    }

    // Feasibility of the natural parameter
    inline bool feas_always() const override { return false; }
    inline bool feasibility(int n, const double* theta) const override
    {
        for(int i = 0; i < n; i++)
            if(theta[i] >= 0.0)
                return false;

        return true;
    }
    Rcpp::NumericVector domain() const override
    {
        return Rcpp::NumericVector::create(R_NegInf, 0.0);
    }
};

Distribution* get_exponential()
{
    return new Exponential();
}
