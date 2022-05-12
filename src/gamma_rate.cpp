#include <Rcpp.h>
#include <boost/math/tools/roots.hpp>

using Rcpp::NumericVector;
using boost::math::tools::newton_raphson_iterate;

/*

Gamma distribution with rate parameter beta

  Y ~ Gamma(alpha, beta)

p(y) = b^a / Gamma(a) * y^(a-1) * exp(-b * y)
log[p(y)] = a * log(b) - logGamma(a) + (a-1) * log(y) - b * y
mean = a / b

We assume Xi|lambdai ~ Pois(si * lambdai) and
  lambdai ~ Gamma with mean=mui and rate=beta, i.e.,
  lambdai ~ Gamma(mui * beta, beta)
Then
  log[p(xi|lambdai)] = xi * log(si * lambdai) - log(xi!) - si * lambdai * xi
  log[p(lambdai)] = mui * b * log(b) - logGamma(mui * b)
                    + (mui * b - 1) * log(lambdai) - b * lambdai

Marginally, Xi ~ NB(ri, p), where
  ri = mui * b, p = 1 / (1 + b / si) = si / (si + b)
  log[p(xi)] = logGamma(ri + xi) - log(xi!) - logGamma(ri)
               + xi * log(p) + ri * log(1 - p)
             = logGamma(mui * b + xi) - log(xi!)
               - logGamma(mui * b)
               + xi * log(si / (si + b))
               + mui * b * log(b / (si + b))
             = xi * log(si) - log(xi!)
               + mui * b * log(b) - logGamma(mui * b)
               + logGamma(mui * b + xi)
               - (mui * b + xi) * log(si + b)
So our target is to find the MLE of b given X and s

Let l(b) = log[p(xi; b)], and note that
  [b * log(b)]' = log(b) + 1
  [logGamma(x)]' = digamma(x)
  [digamma(x)]' = trigamma(x)
So
  [l(b)]' = mui * (log(b) + 1) - mui * digamma(mui * b)
            + mui * digamma(mui * b + xi)
            - mui * log(si + b) - (mui * b + xi) / (si + b)
  [l(b)]'' = mui / b - mui^2 * trigamma(mui * b)
             + mui^2 * trigamma(mui * b + xi)
             - mui / (si + b) - (mui * si - xi) / (si + b)^2

 */

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
            const double xi = m_x[i], mui = m_mu[i], si = m_s[i];
            const double mb = mui * b;
            const double mbx = mb + xi;
            const double msx = mui * si - xi;
            const double bs = b + si;
            // Two log() calls for numerical stability
            const double dli = mui * (std::log(b) - std::log(bs)) +
                mui * (R::digamma(mbx) - R::digamma(mb)) + msx / bs;
            const double d2li = mui / b - mui / bs +
                mui * mui * (R::trigamma(mbx) - R::trigamma(mb)) -
                msx / bs / bs;
            dl += dli;
            d2l += d2li;
        }
        // Rcpp::Rcout << "b = " << b << ", dl = " << dl << ", d2l = " << d2l << std::endl;
        return std::make_pair(dl / m_n, d2l / m_n);
    }
};

// A local maximum of l(b) should satisfy [l(b)]'=0 and [l(b)]''<0
// Note that [l(b)]' -> -0 and [l(b)]'' -> +0 when b -> +Inf,
// so we need to initialize b from a small value
//
// First we choose a proper upper bound for b*, the (local) maximum
// Let ub=max(s), if [l(ub)]''<0, then use ub;
// otherwise repeatedly multiply ub by a constant 0<gamma<1, until
// [l(ub)]''>0 and [l(gamma*ub)]''<0
//
// We then compute the result b* within [lb, ub] given initial value b0
// If [l(b*)]'' > 0, it means it is not a local maximum
// In this case, we set x0 <- gamma*x0, lb <- gamma*lb, and recompute the root

// [[Rcpp::export]]
double gamma_rate(NumericVector x, NumericVector mu, NumericVector s)
{
    // Function object
    LogLikDeriv deriv(x, mu, s);

    // Find a proper upper bound
    // gamma is the decrease factor, max_try is the maximum number of tries
    double ub = Rcpp::max(s);
    const double gamma = 0.5;
    const int max_try = 10;
    std::pair<double, double> dvals = deriv(ub);
    if(dvals.second > 0.0)
    {
        for(int i = 0; i < max_try; i++)
        {
            const double new_ub = gamma * ub;
            std::pair<double, double> new_dvals = deriv(new_ub);
            if(new_dvals.second <= 0.0)
                break;
            ub = new_ub;
        }
    }

    // Lower bound
    double lb = 1e-6;
    // Initial guess
    double guess = std::min(1.0, 0.5 * (lb + ub));
    // Precision parameter
    const int digits = std::numeric_limits<double>::digits;
    int get_digits = static_cast<int>(digits * 0.5);
    // Maximum number of iterations for Newton's method
    const std::uintmax_t max_iter = 10;

    double res = 0.0;
    for(int i = 0; i < max_try; i++)
    {
        std::uintmax_t iter = max_iter;
        res = newton_raphson_iterate(deriv, guess, lb, ub, get_digits, iter);
        std::pair<double, double> dvals = deriv(res);
        // If [l(b*)]'' > 0, make guess smaller and try again
        if(dvals.second > 0.0)
        {
            guess *= gamma;
            lb *= gamma;
        } else {
            return res;
        }
    }

    return res;
}
