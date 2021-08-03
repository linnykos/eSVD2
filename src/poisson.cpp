#include <Rcpp.h>
#include "distribution.h"

using Rcpp::NumericVector;

class Poisson: public Distribution
{
public:
    int log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    )
    {
        const double log_si = std::log(si);
        int non_na = 0;
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                res[j] = 0.0;
            } else {
                non_na++;
                res[j] = Ai[j] * thetai[j] - std::exp(log_si + thetai[j]);
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
                res[i] = Aj[i] * thetaj[i] - std::exp(std::log(s[i]) + thetaj[i]);
            }
        }
        return non_na;
    }

    int dlog_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    )
    {
        const double log_si = std::log(si);
        int non_na = 0;
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                res[j] = 0.0;
            } else {
                non_na++;
                res[j] = Ai[j] - std::exp(log_si + thetai[j]);
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
                res[i] = Aj[i] - std::exp(std::log(s[i]) + thetaj[i]);
            }
        }
        return non_na;
    }

    int d2log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    )
    {
        const double log_si = std::log(si);
        int non_na = 0;
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                res[j] = 0.0;
            } else {
                non_na++;
                res[j] = -std::exp(log_si + thetai[j]);
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
                res[i] = -std::exp(std::log(s[i]) + thetaj[i]);
            }
        }
        return non_na;
    }

    int d12log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* resd1, double* resd2
    )
    {
        const double log_si = std::log(si);
        int non_na = 0;
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                resd1[j] = resd2[j] = 0.0;
            } else {
                non_na++;
                resd2[j] = -std::exp(log_si + thetai[j]);
                resd1[j] = Ai[j] + resd2[j];
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
                resd2[i] = -std::exp(std::log(s[i]) + thetaj[i]);
                resd1[i] = Aj[i] + resd2[i];
            }
        }
        return non_na;
    }
};

// [[Rcpp::export]]
SEXP distribution_poisson()
{
    Poisson* distr = new Poisson();
    return Rcpp::XPtr<Distribution>(distr, true);
}
