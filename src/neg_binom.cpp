#include <Rcpp.h>
#include "distribution.h"

class NegBinom: public Distribution
{
public:
    // Log-density for a single data point Aij
    inline double log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj
    ) const override
    {
        return Aij * thetaij + si * gammaj * std::log(1.0 - std::exp(thetaij));
    }

    // Special case of Aij==0
    inline double log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj
    ) const override
    {
        return si * gammaj * std::log(1.0 - std::exp(thetaij));
    }

    // 1st and 2nd derivatives of the log-density function w.r.t. thetaij
    inline void d12log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2
    ) const override
    {
        const double exp_theta = std::exp(thetaij);
        const double one_minus_exp_theta = 1.0 - exp_theta;
        const double common = -si * gammaj * exp_theta / one_minus_exp_theta;
        if(compute_d1) *res1 = Aij + common;
        if(compute_d2) *res2 = common / one_minus_exp_theta;
    }

    // Special case of Aij==0
    inline void d12log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2
    ) const override
    {
        const double exp_theta = std::exp(thetaij);
        const double one_minus_exp_theta = 1.0 - exp_theta;
        const double common = -si * gammaj * exp_theta / one_minus_exp_theta;
        if(compute_d1) *res1 = common;
        if(compute_d2) *res2 = common / one_minus_exp_theta;
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
};

// [[Rcpp::export]]
SEXP distribution_neg_binom()
{
    NegBinom* distr = new NegBinom();
    return Rcpp::XPtr<Distribution>(distr, true);
}
