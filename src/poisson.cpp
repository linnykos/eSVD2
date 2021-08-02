#include <Rcpp.h>
#include "distribution.h"

using Rcpp::NumericVector;

class Poisson: public Distribution
{
public:
    void log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    )
    {
        const double log_si = std::log(si);
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                res[j] = NumericVector::get_na();
            } else {
                res[j] = Ai[j] * thetai[j] - std::exp(log_si + thetai[j]);
            }
        }
    }

    void log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res
    )
    {
        for(int i = 0; i < n; i++)
        {
            if(NumericVector::is_na(Aj[i]))
            {
                res[i] = NumericVector::get_na();
            } else {
                res[i] = Aj[i] * thetaj[i] - std::exp(std::log(s[i]) + thetaj[i]);
            }
        }
    }

    void dlog_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    )
    {
        const double log_si = std::log(si);
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                res[j] = NumericVector::get_na();
            } else {
                res[j] = Ai[j] - std::exp(log_si + thetai[j]);
            }
        }
    }

    void dlog_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res
    )
    {
        for(int i = 0; i < n; i++)
        {
            if(NumericVector::is_na(Aj[i]))
            {
                res[i] = NumericVector::get_na();
            } else {
                res[i] = Aj[i] - std::exp(std::log(s[i]) + thetaj[i]);
            }
        }
    }

    void d2log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    )
    {
        const double log_si = std::log(si);
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                res[j] = NumericVector::get_na();
            } else {
                res[j] = -std::exp(log_si + thetai[j]);
            }
        }
    }

    void d2log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res
    )
    {
        for(int i = 0; i < n; i++)
        {
            if(NumericVector::is_na(Aj[i]))
            {
                res[i] = NumericVector::get_na();
            } else {
                res[i] = -std::exp(std::log(s[i]) + thetaj[i]);
            }
        }
    }

    void d12log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* resd1, double* resd2
    )
    {
        const double log_si = std::log(si);
        for(int j = 0; j < p; j++)
        {
            if(NumericVector::is_na(Ai[j]))
            {
                resd1[j] = resd2[j] = NumericVector::get_na();
            } else {
                resd2[j] = -std::exp(log_si + thetai[j]);
                resd1[j] = Ai[j] + resd2[j];
            }
        }
    }

    void d12log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* resd1, double* resd2
    )
    {
        for(int i = 0; i < n; i++)
        {
            if(NumericVector::is_na(Aj[i]))
            {
                resd1[i] = resd2[i] = NumericVector::get_na();
            } else {
                resd2[i] = -std::exp(std::log(s[i]) + thetaj[i]);
                resd1[i] = Aj[i] + resd2[i];
            }
        }
    }
};

// [[Rcpp::export]]
SEXP distribution_poisson()
{
    Poisson* distr = new Poisson();
    return Rcpp::XPtr<Distribution>(distr, true);
}
