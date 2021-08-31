#include <Rcpp.h>
#include "distribution.h"

using Rcpp::NumericVector;

class NegBinom: public Distribution
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
                res[j] = Ai[j] * thetai[j] + si * gamma[j] * std::log(1.0 - std::exp(thetai[j]));
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
                res[i] = Aj[i] * thetaj[i] + s[i] * gammaj * std::log(1.0 - std::exp(thetaj[i]));
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
                const double exp_theta = std::exp(thetai[j]);
                const double one_minus_exp_theta = 1.0 - exp_theta;
                const double common = -si * gamma[j] * exp_theta / one_minus_exp_theta;
                if(compute_d1)
                    res1[j] = Ai[j] + common;
                if(compute_d2)
                    res2[j] = common / one_minus_exp_theta;
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
                const double exp_theta = std::exp(thetaj[i]);
                const double one_minus_exp_theta = 1.0 - exp_theta;
                const double common = -s[i] * gammaj * exp_theta / one_minus_exp_theta;
                if(compute_d1)
                    res1[i] = Aj[i] + common;
                if(compute_d2)
                    res2[i] = common / one_minus_exp_theta;
            }
        }
        return non_na;
    }
};

// [[Rcpp::export]]
SEXP distribution_neg_binom()
{
    NegBinom* distr = new NegBinom();
    return Rcpp::XPtr<Distribution>(distr, true);
}
