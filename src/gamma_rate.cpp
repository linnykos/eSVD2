#include <Rcpp.h>
#include <boost/math/tools/roots.hpp>

using Rcpp::NumericVector;
using boost::math::tools::newton_raphson_iterate;

// Also see gamma_rate.R
//
// l(b) = xi * log(si) - log(xi!) + mui * b * log(b) - lgamma(mui * b)
//        + lgamma(mui * b + xi) - (mui * b + xi) * log(si + b)
//
// [l(b)]' = mui * (log(b) + 1) - mui * digamma(mui * b)
//           + mui * digamma(mui * b + xi)
//           - mui * log(si + b) - (mui * b + xi) / (si + b)
//         = mui * log(b / (si + b))
//           + mui * [digamma(mui * b + xi) - digamma(mui * b)]
//           + (mui * si - xi) / (si + b)
//
// [l(b)]'' = mui / b - mui^2 * trigamma(mui * b)
//            + mui^2 * trigamma(mui * b + xi)
//            - mui / (si + b) - (mui * si - xi) / (si + b)^2
//          = mui / b - mui / (si + b)
//            + mui^2 * [trigamma(mui * b + xi) - trigamma(mui * b)]
//            - (mui * si - xi) / (si + b)^2
//
// Want to find the root of [l(b)]'

class LogLikDeriv
{
private:
    const int     m_n;
    const double* m_x;
    const double* m_mu;
    const double* m_s;

public:
    LogLikDeriv(NumericVector x, NumericVector mu, NumericVector s) :
        m_n(x.length()), m_x(x.begin()), m_mu(mu.begin()), m_s(s.begin())
    {}

    // Returns [l(b)]' and [l(b)]''
    std::pair<double, double> operator()(const double& b)
    {
        double dl = 0.0, d2l = 0.0;
        for(int i = 0; i < m_n; i++)
        {
            const double xi = m_x[i];
            const double mui = m_mu[i];
            const double mui2 = mui * mui;
            const double si = m_s[i];
            const double mb = mui * b;
            const double mbx = mb + xi;
            const double bs = b + si;
            const double dli = mui * (std::log(b) + 1.0) +
                mui * (R::digamma(mbx) - R::digamma(mb)) -
                mui * std::log(bs) - mbx / bs;
            const double d2li = mui / b + mui2 * (R::trigamma(mbx) - R::trigamma(mb)) -
                mui / bs - (mui * si - xi) / bs / bs;
            dl += dli;
            d2l += d2li;
        }
        return std::make_pair(dl / m_n, d2l / m_n);
    }
};

// Note that [l(b)]' -> 0 when b -> +Inf, so we need to initialize b
// from a small value
//
// We first compute the result x* within [lb, ub] given initial value x0
// If [l(x*)]'' > 0, it means it is not a local maximum
// In this case, we set x0 <- x0/2, lb <- lb/2, and recompute the root

// [[Rcpp::export]]
double gamma_rate(NumericVector x, NumericVector mu, NumericVector s)
{
    // Function object
    LogLikDeriv deriv(x, mu, s);
    // Initial guess
    double guess = 1.0;
    // Lower bound
    double lb = 1e-6;
    // Upper bound
    double ub = Rcpp::max(s);
    ub = std::max(ub, 1.1 * guess);
    // Precision parameter
    const int digits = std::numeric_limits<double>::digits;
    int get_digits = static_cast<int>(digits * 0.5);
    // Maximum number of iterations
    const std::uintmax_t max_iter = 10;
    // Maximum number of tries
    const int max_try = 10;

    double res = 0.0;
    for(int i = 0; i < max_try; i++)
    {
        std::uintmax_t iter = max_iter;
        res = newton_raphson_iterate(deriv, guess, lb, ub, get_digits, iter);
        std::pair<double, double> dvals = deriv(res);
        // If d2l > 0, make guess smaller and try again
        if(dvals.second > 0.0)
        {
            guess *= 0.5;
            lb *= 0.5;
        } else {
            return res;
        }
    }

    return res;
}
