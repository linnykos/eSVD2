#include <Rcpp.h>

class Distribution
{
public:
    virtual ~Distribution() { Rcpp::Rcout << "** destructor **" << std::endl; }

    virtual int log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    ) = 0;

    virtual int log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res
    ) = 0;

    virtual int dlog_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    ) = 0;

    virtual int dlog_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res
    ) = 0;

    virtual int d2log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    ) = 0;

    virtual int d2log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res
    ) = 0;

    virtual int d12log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* resd1, double* resd2
    ) = 0;

    virtual int d12log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* resd1, double* resd2
    ) = 0;
};
