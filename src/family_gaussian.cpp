#include <Rcpp.h>
#include "distribution.h"

// Distribution: Gaussian
// Canonical parameter: mu
// Natural parameter: theta = mu
// Log-likelihood: -0.5 * (theta - A / s)^2 * s / gamma^2

class Gaussian: public Distribution
{
public:
    // fn_gammaj = 1 / gammaj^2
    inline void compute_fn_gammaj(double gammaj, double& fn_gammaj) const override
    {
        fn_gammaj = 1.0 / gammaj / gammaj;
    }

    // Log-density for a single data point Aij
    // fn_gammaj = 1 / gammaj^2
    inline double log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        const double u = thetaij - Aij / si;
        return -0.5 * u * u * si * fn_gammaj;
    }

    // Special case of Aij==0
    // fn_gammaj = 1 / gammaj^2
    inline double log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        return -0.5 * thetaij * thetaij * si * fn_gammaj;
    }

    // 1st and 2nd derivatives of the log-density function w.r.t. thetaij
    // fn_gammaj = 1 / gammaj^2
    inline void d12log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        if(compute_d1) *res1 = -(si * thetaij - Aij) * fn_gammaj;
        if(compute_d2) *res2 = -si * fn_gammaj;
    }

    // Special case of Aij==0
    // fn_si = log(si)
    inline void d12log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        const double d2 = -si * fn_gammaj;
        if(compute_d1) *res1 = thetaij * d2;
        if(compute_d2) *res2 = d2;
    }
};

Distribution* get_gaussian()
{
    return new Gaussian();
}
