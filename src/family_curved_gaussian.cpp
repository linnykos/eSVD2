#include <Rcpp.h>
#include "distribution.h"

// Distribution: one-parameter Gaussian where sd = mean/scalar
// Canonical parameter: mu
// Natural parameter: theta = 1/mu
// Log-likelihood: ga * theta - 0.5 * ga * Ai * theta^2 / s + log(theta), ga = gamma^2 * A

class CurvedGaussian: public Distribution
{
public:
    // fn_gammaj = gamma^2
    inline void compute_fn_gammaj(double gammaj, double& fn_gammaj) const override
    {
        fn_gammaj = gammaj * gammaj;
    }

    // Log-density for a single data point Aij
    // fn_gammaj = gamma^2
    inline double log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        const double ga = fn_gammaj * Aij;
        return ga * thetaij - 0.5 * ga * Aij * thetaij * thetaij / si + std::log(thetaij);
    }

    // Special case of Aij==0
    inline double log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        return std::log(thetaij);
    }

    // 1st and 2nd derivatives of the log-density function w.r.t. thetaij
    // fn_gammaj = gamma^2
    inline void d12log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        const double ga = fn_gammaj * Aij;
        if(compute_d1) *res1 = ga - ga * Aij * thetaij / si + 1.0 / thetaij;
        if(compute_d2)
        {
            const double u = gammaj * Aij;
            *res2 = -u * u / si - 1.0 / (thetaij * thetaij);
        }
    }

    // Special case of Aij==0
    inline void d12log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        const double d1 = 1.0 / thetaij;
        if(compute_d1) *res1 = d1;
        if(compute_d2) *res2 = -d1 / thetaij;
    }

    // Feasibility of the natural parameter
    inline bool feas_always() const override { return false; }
    inline bool feasibility(int n, const double* theta) const override
    {
        for(int i = 0; i < n; i++)
            if(theta[i] <= 0.0)
                return false;

        return true;
    }
    Rcpp::NumericVector domain() const override
    {
        return Rcpp::NumericVector::create(0.0, R_PosInf);
    }
};

Distribution* get_curved_gaussian()
{
    return new CurvedGaussian();
}
