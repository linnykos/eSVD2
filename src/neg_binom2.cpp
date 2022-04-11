#include <Rcpp.h>
#include "distribution.h"

// Compute log(r + exp(x)), r > 0 with numerical stability
// log(r + exp(x)) = log(r) + log(1 + exp(x - log(r))) = log(r) + softplus(x - log(r))
// softplus(x) = log(1 + exp(x)) = log(1 + exp(-|x|)) + max(x, 0)
inline double softplusr(double x, double logr)
{
    const double xr = x - logr;
    const double softplus = std::log(1.0 + std::exp(-std::abs(xr))) + std::max(xr, 0.0);
    return logr + softplus;
}

// Compute exp(x) / (r + exp(x)) with numerical stability
// exp(x) / (r + exp(x)) = exp(x - log(r)) / (1 + exp(x - log(r))) = sigmoid(x - log(r))
inline double sigmoidr(double x, double logr)
{
    const double xr = x - logr;
    return 1.0 / (1.0 + std::exp(-xr));
}

class NegBinom2: public Distribution
{
public:
    // fn_gammaj = log(gammaj)
    inline void compute_fn_gammaj(double gammaj, double& fn_gammaj) const override
    {
        fn_gammaj = std::log(gammaj);
    }

    // Log-density for a single data point Aij
    // fn_gammaj = log(gammaj)
    inline double log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj
    ) const override
    {
        return Aij * thetaij - (Aij + gammaj) * softplusr(thetaij, fn_gammaj);
    }

    // Special case of Aij==0
    // fn_gammaj = log(gammaj)
    inline double log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj
    ) const override
    {
        return -gammaj * softplusr(thetaij, fn_gammaj);
    }

    // 1st and 2nd derivatives of the log-density function w.r.t. thetaij
    // fn_gammaj = log(gammaj)
    inline void d12log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2
    ) const override
    {
        const double common = -(Aij + gammaj) * sigmoidr(thetaij, fn_gammaj);
        if(compute_d1) *res1 = Aij + common;
        if(compute_d2) *res2 = common * gammaj / (gammaj + std::exp(thetaij));
    }

    // Special case of Aij==0
    // fn_gammaj = log(gammaj)
    inline void d12log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2
    ) const override
    {
        const double common = -gammaj * sigmoidr(thetaij, fn_gammaj);
        if(compute_d1) *res1 = common;
        if(compute_d2) *res2 = common * gammaj / (gammaj + std::exp(thetaij));
    }

    // Feasibility of the natural parameter
    inline bool feas_always() const override { return true; }
    inline bool feasibility(int n, const double* theta) const override { return true; }
};

Distribution* get_neg_binom2()
{
    return new NegBinom2();
}

// [[Rcpp::export]]
SEXP distribution_neg_binom2()
{
    NegBinom2* distr = new NegBinom2();
    return Rcpp::XPtr<Distribution>(distr, true);
}
