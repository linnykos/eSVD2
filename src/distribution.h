#ifndef ESVD2_DISTRIBUTION_H
#define ESVD2_DISTRIBUTION_H

#include <Rcpp.h>
#include "data_loader.h"

// See eSVD2_writing/writeup/2021-05-20-covariates.pdf
class Distribution
{
public:
    virtual ~Distribution() { Rcpp::Rcout << "** destructor **" << std::endl; }

    // Log-density for a single data point Aij
    // fn_si is a function of si that can be shared by A[i, ]
    // fn_gammaj is a function of gammaj that can be shared by A[, j]
    virtual double log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj
    ) const = 0;
    // Special case of Aij==0
    virtual double log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj
    ) const
    {
        return log_prob_single(0.0, thetaij, si, gammaj, fn_si, fn_gammaj);
    }

    // 1st and 2nd derivatives of the log-density function w.r.t. thetaij
    virtual void d12log_prob_single(
        double Aij, double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2
    ) const = 0;
    // Special case of Aij==0
    virtual void d12log_prob_single(
        double thetaij, double si, double gammaj,
        double fn_si, double fn_gammaj, bool compute_d1, bool compute_d2,
        double* res1, double* res2
    ) const
    {
        d12log_prob_single(
            0.0, thetaij, si, gammaj, fn_si, fn_gammaj,
            compute_d1, compute_d2, res1, res2
        );
    }

    // Functions to compute fn_si and fn_gammaj
    virtual void compute_fn_si(double si, double& fn_si) const {}
    virtual void compute_fn_gammaj(double gammaj, double& fn_gammaj) const {}



    // Log-density for the i-th row of the data matrix, res[p x 1]
    // res[j] <- 0 if Aij is NA
    // Returns the number of non-NA elements in A[i, ]
    int log_prob_row(
        DataLoader* loader, std::size_t row_ind, const double* thetai,
        double si, const double* gamma, double* res
    ) const
    {
        // Function of si that can be shared by A[i, ]
        double fn_si = 0.0, fn_gammaj = 0.0;
        compute_fn_si(si, fn_si);

        // Iterate on A[i, ]
        VecIterator& iter = loader->row_iter(row_ind);
        int num_non_na = 0;
        for(; iter; ++iter, ++thetai, ++gamma, ++res)
        {
            // Compute fn_gammaj
            compute_fn_gammaj(*gamma, fn_gammaj);
            // Get the flag and the value
            double val;
            Flag flag = iter.value(val);
            switch(flag)
            {
            case Flag::regular:
                num_non_na++;
                *res = log_prob_single(val, *thetai, si, *gamma, fn_si, fn_gammaj);
                break;
            case Flag::zero:
                num_non_na++;
                *res = log_prob_single(*thetai, si, *gamma, fn_si, fn_gammaj);
                break;
            default:
                // Set result to be zero if Aij is NA
                *res = 0.0;
            }
        }
        return num_non_na;
    }

    // 1st and 2nd derivatives of log-density w.r.t. the i-th row of theta, res [p x 1]
    // res[j] <- 0 if Aij is NA
    // res1 or res2 can be nullptr, in which case the corresponding derivatives are not computed
    // Returns the number of non-NA elements in A[i, ]
    int d12log_prob_row(
        DataLoader* loader, std::size_t row_ind, const double* thetai,
        double si, const double* gamma, double* res1, double* res2
    ) const
    {
        // Function of si that can be shared by A[i, ]
        double fn_si = 0.0, fn_gammaj = 0.0;
        compute_fn_si(si, fn_si);
        // Whether to compute d1 and d2
        const bool compute_d1 = (res1 != nullptr);
        const bool compute_d2 = (res2 != nullptr);

        // Iterate on A[i, ]
        VecIterator& iter = loader->row_iter(row_ind);
        int num_non_na = 0;
        for(; iter; ++iter, ++thetai, ++gamma, ++res1, ++res2)
        {
            // Compute fn_gammaj
            compute_fn_gammaj(*gamma, fn_gammaj);
            // Get the flag and the value
            double val;
            Flag flag = iter.value(val);
            switch(flag)
            {
            case Flag::regular:
                num_non_na++;
                d12log_prob_single(
                    val, *thetai, si, *gamma, fn_si, fn_gammaj,
                    compute_d1, compute_d2, res1, res2
                );
                break;
            case Flag::zero:
                num_non_na++;
                d12log_prob_single(
                    *thetai, si, *gamma, fn_si, fn_gammaj,
                    compute_d1, compute_d2, res1, res2
                );
                break;
            default:
                // Set result to be zero if Aij is NA
                if(compute_d1) *res1 = 0.0;
                if(compute_d2) *res2 = 0.0;
            }
        }
        return num_non_na;
    }

    // Log-density for the j-th column of the data matrix, res [n x 1]
    // res[i] <- 0 if Aij is NA
    // Returns the number of non-NA elements in A[, j]
    int log_prob_col(
        DataLoader* loader, std::size_t col_ind, const double* thetaj,
        const double* s, double gammaj, double* res
    ) const
    {
        // Function of gammaj that can be shared by A[, j]
        double fn_si = 0.0, fn_gammaj = 0.0;
        compute_fn_gammaj(gammaj, fn_gammaj);

        // Iterate on A[, j]
        VecIterator& iter = loader->col_iter(col_ind);
        int num_non_na = 0;
        for(; iter; ++iter, ++thetaj, ++s, ++res)
        {
            // Compute fn_si
            compute_fn_si(*s, fn_si);
            // Get the flag and the value
            double val;
            Flag flag = iter.value(val);
            switch(flag)
            {
            case Flag::regular:
                num_non_na++;
                *res = log_prob_single(val, *thetaj, *s, gammaj, fn_si, fn_gammaj);
                break;
            case Flag::zero:
                num_non_na++;
                *res = log_prob_single(*thetaj, *s, gammaj, fn_si, fn_gammaj);
                break;
            default:
                // Set result to be zero if Aij is NA
                *res = 0.0;
            }
        }
        return num_non_na;
    }

    // 1st and 2nd derivatives of log-density w.r.t. the j-th column of theta, res [n x 1]
    // res[i] <- 0 if Aij is NA
    // res1 or res2 can be nullptr, in which case the corresponding derivatives are not computed
    // Returns the number of non-NA elements in A[, j]
    int d12log_prob_col(
        DataLoader* loader, std::size_t col_ind, const double* thetaj,
        const double* s, double gammaj, double* res1, double* res2
    ) const
    {
        // Function of gammaj that can be shared by A[, j]
        double fn_si = 0.0, fn_gammaj = 0.0;
        compute_fn_gammaj(gammaj, fn_gammaj);
        // Whether to compute d1 and d2
        const bool compute_d1 = (res1 != nullptr);
        const bool compute_d2 = (res2 != nullptr);

        // Iterate on A[, j]
        VecIterator& iter = loader->col_iter(col_ind);
        int num_non_na = 0;
        for(; iter; ++iter, ++thetaj, ++s, ++res1, ++res2)
        {
            // Compute fn_si
            compute_fn_si(*s, fn_si);
            // Get the flag and the value
            double val;
            Flag flag = iter.value(val);
            switch(flag)
            {
            case Flag::regular:
                num_non_na++;
                d12log_prob_single(
                    val, *thetaj, *s, gammaj, fn_si, fn_gammaj,
                    compute_d1, compute_d2, res1, res2
                );
                break;
            case Flag::zero:
                num_non_na++;
                d12log_prob_single(
                    *thetaj, *s, gammaj, fn_si, fn_gammaj,
                    compute_d1, compute_d2, res1, res2
                );
                break;
            default:
                // Set result to be zero if Aij is NA
                if(compute_d1) *res1 = 0.0;
                if(compute_d2) *res2 = 0.0;
            }
        }
        return num_non_na;
    }

    // Feasibility of the natural parameter
    virtual bool feas_always() const = 0;
    virtual bool feasibility(int n, const double* theta) const = 0;
};


#endif  // ESVD2_DISTRIBUTION_H
