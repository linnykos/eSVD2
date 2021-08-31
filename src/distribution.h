#include <Rcpp.h>

// See eSVD2_writing/writeup/2021-05-20-covariates.pdf
class Distribution
{
public:
    virtual ~Distribution() { Rcpp::Rcout << "** destructor **" << std::endl; }

    // Log-density for the i-th row of the data matrix, res[p x 1]
    // res[j] <- 0 if Ai[j] is NA
    // Returns the number of non-NA elements in Ai
    virtual int log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res
    ) = 0;

    // Log-density for the j-th column of the data matrix, res [n x 1]
    // res[i] <- 0 if Aj[i] is NA
    // Returns the number of non-NA elements in Aj
    virtual int log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res
    ) = 0;

    // 1st and 2nd derivatives of log-density w.r.t. the i-th row of theta, res [p x 1]
    // res[j] <- 0 if Ai[j] is NA
    // res1 or res2 can be NULL, in which case the corresponding derivatives are not computed
    // Returns the number of non-NA elements in Ai
    virtual int d12log_prob_row(
        int p, const double* Ai, const double* thetai,
        double si, const double* gamma, double* res1, double* res2
    ) = 0;

    // 1st and 2nd derivatives of log-density w.r.t. the j-th column of theta, res [n x 1]
    // res[i] <- 0 if Aj[i] is NA
    // res1 or res2 can be NULL, in which case the corresponding derivatives are not computed
    // Returns the number of non-NA elements in Aj
    virtual int d12log_prob_col(
        int n, const double* Aj, const double* thetaj,
        const double* s, double gammaj, double* res1, double* res2
    ) = 0;
};
