#include <Rcpp.h>
#include "distribution.h"

// Distribution: Poisson
// Canonical parameter: lambda
// Natural parameter: theta = log(lambda)
// Log-likelihood: A * theta - exp(log(s) + theta)

class Poisson: public Distribution
{
public:
    // fn_si = log(si)
    inline void compute_fn_si(double si, double& fn_si) const override
    {
        fn_si = std::log(si);
    }

    // Log-density for a single data point Aij
    // fn_si = log(si)
    inline double log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        return Aij * thetaij - std::exp(fn_si + thetaij);
    }

    // Special case of Aij==0
    // fn_si = log(si)
    inline double log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        return -std::exp(fn_si + thetaij);
    }

    // 1st and 2nd derivatives of the log-density function w.r.t. thetaij
    // fn_si = log(si)
    inline void d12log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        const double d2 = -std::exp(fn_si + thetaij);
        if(compute_d1) *res1 = Aij + d2;
        if(compute_d2) *res2 = d2;
    }

    // Special case of Aij==0
    // fn_si = log(si)
    inline void d12log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        const double d2 = -std::exp(fn_si + thetaij);
        if(compute_d1) *res1 = d2;
        if(compute_d2) *res2 = d2;
    }

    // Feasibility of the natural parameter
    inline bool feas_always() const override { return true; }
    inline bool feasibility(int n, const double* theta) const override { return true; }
};

Distribution* get_poisson()
{
    return new Poisson();
}

// [[Rcpp::export]]
SEXP distribution_poisson()
{
    Poisson* distr = new Poisson();
    return Rcpp::XPtr<Distribution>(distr, true);
}
