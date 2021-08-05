#include <Rcpp.h>
#include "distribution.h"

using Rcpp::NumericVector;

class NegBinom: public Distribution
{
public:
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

    int dlog_prob_row(
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
                const double exp_theta = std::exp(thetai[j]);
                res[j] = Ai[j] - si * gamma[j] * exp_theta / (1.0 - exp_theta);
            }
        }
        return non_na;
    }

    int dlog_prob_col(
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
                const double exp_theta = std::exp(thetaj[i]);
                res[i] = Aj[i] - s[i] * gammaj * exp_theta / (1.0 - exp_theta);
            }
        }
        return non_na;
    }

    int d2log_prob_row(
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
                const double exp_theta = std::exp(thetai[j]);
                const double one_minus_exp_theta = 1.0 - exp_theta;
                res[j] = -si * gamma[j] * exp_theta / one_minus_exp_theta / one_minus_exp_theta;
            }
        }
        return non_na;
    }

    int d2log_prob_col(
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
                const double exp_theta = std::exp(thetaj[i]);
                const double one_minus_exp_theta = 1.0 - exp_theta;
                res[i] = -s[i] * gammaj * exp_theta / one_minus_exp_theta / one_minus_exp_theta;
            }
        }
        return non_na;
    }

    int d12log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* resd1, double* resd2
    )
    {
        int non_na = 0;
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                resd1[j] = resd2[j] = 0.0;
            } else {
                non_na++;
                const double exp_theta = std::exp(thetai[j]);
                const double one_minus_exp_theta = 1.0 - exp_theta;
                const double common = -si * gamma[j] * exp_theta / one_minus_exp_theta;
                resd1[j] = Ai[j] + common;
                resd2[j] = common / one_minus_exp_theta;
            }
        }
        return non_na;
    }

    int d12log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* resd1, double* resd2
    )
    {
        int non_na = 0;
        for(int i = 0; i < n; i++)
        {
            if(NumericVector::is_na(Aj[i]))
            {
                resd1[i] = resd2[i] = 0.0;
            } else {
                non_na++;
                const double exp_theta = std::exp(thetaj[i]);
                const double one_minus_exp_theta = 1.0 - exp_theta;
                const double common = -s[i] * gammaj * exp_theta / one_minus_exp_theta;
                resd1[i] = Aj[i] + common;
                resd2[i] = common / one_minus_exp_theta;
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
