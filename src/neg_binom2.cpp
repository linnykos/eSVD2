#include <Rcpp.h>
#include "distribution.h"

using Rcpp::NumericVector;

// Compute log(r + exp(x)), r > 0 with numerical stability
// log(r + exp(x)) = log(r) + log(1 + exp(x - log(r))) = log(r) + softplus(x - log(r))
// softplus(x) = log(1 + exp(x)) = log(1 + exp(-|x|)) + max(x, 0)
inline double softplusr(double x, double r)
{
    const double logr = std::log(r);
    const double xr = x - logr;
    const double softplus = std::log(1.0 + std::exp(-std::abs(xr))) + std::max(xr, 0.0);
    return logr + softplus;
}

// Compute exp(x) / (r + exp(x)) with numerical stability
// exp(x) / (r + exp(x)) = exp(x - log(r)) / (1 + exp(x - log(r))) = sigmoid(x - log(r))
inline double sigmoidr(double x, double r)
{
    const double xr = x - std::log(r);
    return 1.0 / (1.0 + std::exp(-xr));
}

class NegBinom2: public Distribution
{
public:
    // Log-density for the i-th row of the data matrix, res[p x 1]
    // res[j] <- 0 if Ai[j] is NA
    // Returns the number of non-NA elements in Ai
    int log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    )
    {
        int non_na = 0;
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                res[j] = 0.0;
            } else {
                non_na++;
                res[j] = Ai[j] * thetai[j] - (Ai[j] + gamma[j]) * softplusr(thetai[j], gamma[j]);
            }
        }
        return non_na;
    }

    // Log-density for the j-th column of the data matrix, res [n x 1]
    // res[i] <- 0 if Aj[i] is NA
    // Returns the number of non-NA elements in Aj
    int log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res
    )
    {
        int non_na = 0;
        for(int i = 0; i < n; i++)
        {
            if(NumericVector::is_na(Aj[i]))
            {
                res[i] = 0.0;
            } else {
                non_na++;
                res[i] = Aj[i] * thetaj[i] - (Aj[i] + gammaj) * softplusr(thetaj[i], gammaj);
            }
        }
        return non_na;
    }

    // 1st and 2nd derivatives of log-density w.r.t. the i-th row of theta, res [p x 1]
    // res[j] <- 0 if Ai[j] is NA
    // res1 or res2 can be NULL, in which case the corresponding derivatives are not computed
    // Returns the number of non-NA elements in Ai
    int d12log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res1, double* res2
    )
    {
        const bool compute_d1 = (res1 != NULL);
        const bool compute_d2 = (res2 != NULL);
        int non_na = 0;
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                if(compute_d1)
                    res1[j] = 0.0;
                if(compute_d2)
                    res2[j] = 0.0;
            } else {
                non_na++;
                const double sig = sigmoidr(thetai[j], gamma[j]);
                const double common = -(Ai[j] + gamma[j]) * sig;
                if(compute_d1)
                    res1[j] = Ai[j] + common;
                if(compute_d2)
                    res2[j] = common * gamma[j] / (gamma[j] + std::exp(thetai[j]));
            }
        }
        return non_na;
    }

    // 1st and 2nd derivatives of log-density w.r.t. the j-th column of theta, res [n x 1]
    // res[i] <- 0 if Aj[i] is NA
    // res1 or res2 can be NULL, in which case the corresponding derivatives are not computed
    // Returns the number of non-NA elements in Aj
    int d12log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res1, double* res2
    )
    {
        const bool compute_d1 = (res1 != NULL);
        const bool compute_d2 = (res2 != NULL);
        int non_na = 0;
        for(int i = 0; i < n; i++)
        {
            if(NumericVector::is_na(Aj[i]))
            {
                if(compute_d1)
                    res1[i] = 0.0;
                if(compute_d2)
                    res2[i] = 0.0;
            } else {
                non_na++;
                const double sig = sigmoidr(thetaj[i], gammaj);
                const double common = -(Aj[i] + gammaj) * sig;
                if(compute_d1)
                    res1[i] = Aj[i] + common;
                if(compute_d2)
                    res2[i] = common * gammaj / (gammaj + std::exp(thetaj[i]));
            }
        }
        return non_na;
    }

    // Feasibility of the natural parameter
    bool feas_always()
    {
        return true;
    };

    bool feasibility(int n, const double* theta)
    {
        return true;
    }
};

// [[Rcpp::export]]
SEXP distribution_neg_binom2()
{
    NegBinom2* distr = new NegBinom2();
    return Rcpp::XPtr<Distribution>(distr, true);
}
