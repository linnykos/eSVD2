#include <Rcpp.h>
#include "distribution.h"

// Distribution: Bernoulli
// Canonical parameter: p
// Natural parameter: theta = log(p/(1-p))
// Log-likelihood: A * theta - log(1 + exp(theta))

// Softplus function
// log(1 + exp(x)) = log(1 + exp(-|x|)) + max(x, 0)
// Can directly use R::log1pexp
inline double softplus(double x)
{
    return std::log(1.0 + std::exp(-std::abs(x))) + std::max(x, 0.0);
}

class Bernoulli: public Distribution
{
public:
    // Log-density for a single data point Aij
    inline double log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        return Aij * thetaij - R::log1pexp(thetaij);
    }

    // Special case of Aij==0
    inline double log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj) const override
    {
        return -R::log1pexp(thetaij);
    }

    // 1st and 2nd derivatives of the log-density function w.r.t. thetaij
    inline void d12log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        const double sigmoid = 1.0 / (1.0 + std::exp(-thetaij));
        if(compute_d1) *res1 = Aij - sigmoid;
        if(compute_d2) *res2 = -sigmoid * (1.0 - sigmoid);
    }

    // Special case of Aij==0
    inline void d12log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2) const override
    {
        const double sigmoid = 1.0 / (1.0 + std::exp(-thetaij));
        if(compute_d1) *res1 = -sigmoid;
        if(compute_d2) *res2 = -sigmoid * (1.0 - sigmoid);
    }
};

Distribution* get_bernoulli()
{
    return new Bernoulli();
}
