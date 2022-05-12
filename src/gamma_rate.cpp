#include <Rcpp.h>
#include <boost/math/tools/roots.hpp>

using Rcpp::NumericVector;
using boost::math::tools::newton_raphson_iterate;

/*************************************************************************

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

*************************************************************************/

/*************************************************************************

Testing R code

# Examples
objfn <- function(b, xi, mui, si) {
  loglik <- xi * log(si) - lgamma(xi + 1) +
    mui * b * log(b) - lgamma(mui * b) + lgamma(mui * b + xi) -
    (mui * b + xi) * log(si + b)
  mean(loglik)
}
gradfn <- function(b, xi, mui, si) {
  grad <- mui * (log(b) + 1) - mui * digamma(mui * b) +
    mui * digamma(mui * b + xi) -
    mui * log(si + b) - (mui * b + xi) / (si + b)
  mean(grad)
}
hessfn <- function(b, xi, mui, si) {
  hess <- mui / b - mui^2 * trigamma(mui * b) +
    mui^2 * trigamma(mui * b + xi) -
    mui / (si + b) - (mui * si - xi) / (si + b)^2
  mean(hess)
}

library(MASS)
mod <- glm.nb(Days ~ .^2, data = quine)
xi <- quine$Days
mui <- fitted(mod)
si <- rep(1, length(xi))

curve(sapply(x, objfn, xi = xi, mui = mui, si = si), 0.01, 10)
curve(sapply(x, objfn, xi = xi, mui = mui, si = si), 0.01, 1)
curve(sapply(x, gradfn, xi = xi, mui = mui, si = si), 0.01, 1)
abline(h = 0, col = "red")
curve(sapply(x, hessfn, xi = xi, mui = mui, si = si), 0.2, 5)
# Log-scale
curve(sapply(exp(x), objfn, xi = xi, mui = mui, si = si), -5, 5)
curve(exp(x) * sapply(exp(x), gradfn, xi = xi, mui = mui, si = si), -5, 5)
abline(h = 0, col = "red")

yeast <- data.frame(cbind(numbers = 0:5, fr = c(213, 128, 37, 18, 3, 1)))
fit <- glm.nb(numbers ~ 1, weights = fr, data = yeast)
summary(fit)
mui <- fitted(fit)
xi <- yeast$numbers
si <- yeast$fr

*************************************************************************/

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



/*************************************************************************

Now we do the reparametrization, θ = 1 / β

  Xi|lambdai ~ Pois(si * lambdai)
  lambdai ~ Gamma(mui / θ, 1 / θ)

In the limit case θ = 0, we have lambdai = mui, a.s.

Marginally, Xi ~ NB(ri, p), where
  ri = mui / θ, p = si / (si + 1 / θ) = (si * θ) / (si * θ + 1)
  log[p(xi)] = logGamma(ri + xi) - log(xi!) - logGamma(ri)
               + xi * log(p) + ri * log(1 - p)

Further let ρ = log(θ) = -log(β), and note that
  log(p) = log(si) + log(θ) - log(si * θ + 1)
         = log(si) + ρ - log(1 + exp(ρ + log(si)))
         = log(si) + ρ - softplus(ρ + log(si))
  log(1 - p) = -log(si * θ + 1)
  ri * log(1 - p) = -mui * log(si * θ + 1) / θ = -mui * si * R(ρ + log(si))
where R(x) = softplus(x) / exp(x)

softplus(x) and R(x) are both numerically stable functions
  softplus(x) = log(1 + exp(x)) = max(0, x) + log(1 + exp(-|x|))
  If x > 0, direct compute R(x)
  If x < 0, compute y = exp(x), and then
  R(x) = log(1 + y) / y ~ 1 - y/2 + y^2/3 - y^3/4 + y^4/5
  Use the cutoff x < -5

Define the objective function as
  l(ρ) = logGamma(ri + xi) - log(xi!) - logGamma(ri)
         + xi * [log(si) + ρ - softplus(ρ + log(si))]
         - mui * si * R(ρ + log(si))
and our target is to find ρ that maximizes l(ρ)

Note that
  ri = mui / θ = mui * exp(-ρ), ri' = -ri
  [logGamma(x)]' = digamma(x)
  [digamma(x)]' = trigamma(x)
  [softplus(x)]' = sigmoid(x)
  [softplus(x)]'' = sigmoid(x) * [1 - sigmoid(x)]
  [R(x)]' = 1/(1+exp(x)) - R(x) = sigmoid(-x) - R(x) = 1 - sigmoid(x) - R(x)
  [R(x)]'' = -sigmoid(x) * [1 - sigmoid(x)] - 1 + sigmoid(x) + R(x)
           = [sigmoid(x)]^2 - 1 + R(x)
We have
  [l(ρ)]' = ri * [digamma(ri) - digamma(ri + xi)]
            + xi * [1 - sigmoid(ρ + log(si))]
            - mui * si * [1 - sigmoid(ρ + log(si)) - R(ρ + log(si))]
  [l(ρ)]'' = -ri * [digamma(ri) - digamma(ri + xi)]
             -ri^2 * [trigamma(ri) - trigamma(ri + xi)]
             - xi * sig * (1 - sig)
             - mui * si * [sig^2 - 1 + R(ρ + log(si))]

*************************************************************************/

// log(1 + exp(x))
inline double softplus(const double& x) { return R::log1pexp(x); }
// 1 / (1 + exp(-x))
inline double sigmoid(const double& x) { return R::plogis(x, 0.0, 1.0, 1, 0); }
// log(1 + exp(x)) / exp(x)
inline double log1pox(const double& x)
{
    const double y = std::exp(x);
    if(x >= -5)  return softplus(x) / y;

    const double y2 = y * y, y3 = y2 * y, y4 = y3 * y;
    return 1.0 - 0.5 * y + y2 / 3.0 - 0.25 * y3 + 0.2 * y4;
}

class LogLikDerivRho
{
private:
    const int     m_n;
    const double* m_x;
    const double* m_mus;
    const double* m_logmu;
    const double* m_logs;

public:
    LogLikDerivRho(NumericVector x, NumericVector mus, NumericVector logmu, NumericVector logs) :
        m_n(x.length()), m_x(x.begin()), m_mus(mus.begin()), m_logmu(logmu.begin()), m_logs(logs.begin())
    {}

    // Returns [l(ρ)]' and [l(ρ)]''
    std::pair<double, double> operator()(const double& rho)
    {
        double dl = 0.0, d2l = 0.0;
        for(int i = 0; i < m_n; i++)
        {
            const double xi = m_x[i], musi = m_mus[i], logmui = m_logmu[i], logsi = m_logs[i];
            const double ri = std::exp(logmui - rho), rs = rho + logsi;
            const double sig = sigmoid(rs), x1s = xi * (1.0 - sig), R = log1pox(rs);
            const double rd = ri * (R::digamma(ri) - R::digamma(ri + xi));
            const double dli = rd + x1s - musi * (1.0 - sig - R);
            const double d2li = -rd - ri * ri * (R::trigamma(ri) - R::trigamma(ri + xi)) -
                                x1s * sig - musi * (sig * sig - 1.0 + R);
            dl += dli;
            d2l += d2li;
        }
        // Rcpp::Rcout << "b = " << b << ", dl = " << dl << ", d2l = " << d2l << std::endl;
        return std::make_pair(dl / m_n, d2l / m_n);
    }
};

// [[Rcpp::export]]
double log_gamma_rate(NumericVector x, NumericVector mu, NumericVector s,
                      double lower = -10.0, double upper = 10.0)
{
    NumericVector logmu = Rcpp::log(mu);
    NumericVector logs = Rcpp::log(s);
    NumericVector mus = mu * s;

    // Function object
    LogLikDerivRho deriv(x, mus, logmu, logs);

    // Bounds on rho = -log(beta)
    const double lb = -upper, ub = -lower;

    // Test [l(ρ)]' on lower bound
    // If [l(lb)]' <= 0, l(ρ) is very likely to be monotonically decreasing
    // So we set ρ=lb, and log(β)=-ρ
    std::pair<double, double> dvals = deriv(lb);
    if(dvals.first <= 0.0)
        return -lb;

    // Now we have [l(lb)]' > 0
    // Test [l(ρ)]' on upper bound
    // If [l(ub)]' >= 0, then the optimal rho will be larger than ub
    // So we set ρ=ub, and log(β)=-ρ
    dvals = deriv(ub);
    if(dvals.first >= 0.0)
        return -ub;

    // Now we have bracketed the solution of [l(rho)]' = 0
    // [l(lb)]' > 0, [l(ub)]' < 0, lb < rho < ub
    // rho should satisfy [l(rho)]' = 0 and [l(rho)]'' < 0
    //
    // It is guaranteed that [l(ub)]'' < 0, and we want to find
    // some u such that [l(u)]' > 0 and [l(u)]'' < 0
    // We do a bisection search
    const int max_try = 10;
    double ulb = lb, uub = ub, u = 0.5 * (ulb + uub);
    for(int i = 0; i < max_try; i++)
    {
        std::pair<double, double> dvals = deriv(u);
        // If [l(u)]' < 0, u must be in the left interval
        if(dvals.first < 0.0)
        {
            uub = u;
            u = 0.5 * (ulb + uub);
        } else if(dvals.second < 0.0) {
            // If [l(u)]' >= 0 and [l(u)]'' < 0, interval found
            break;
        } else {
            // If [l(u)]' >= 0 and [l(u)]'' >= 0,
            // u must be in the right interval
            ulb = u;
            u = 0.5 * (ulb + uub);
        }
    }

    // Initial guess
    double guess = 0.5 * (u + ub);
    // Precision parameter
    const int digits = std::numeric_limits<double>::digits;
    int get_digits = static_cast<int>(digits * 0.5);
    // Maximum number of iterations for Newton's method
    const std::uintmax_t max_iter = 10;

    std::uintmax_t iter = max_iter;
    double res = newton_raphson_iterate(deriv, guess, u, ub, get_digits, iter);

    return -res;
}
