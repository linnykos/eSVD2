#include <RcppEigen.h>
#include "objective.h"

using Rcpp::NumericVector;
using Rcpp::List;

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;

List line_search(
    double alpha0, NumericVector x, double fx, NumericVector direction,
    Objective& objective, int max_linesearch, double scaling = 0.5
)
{
    double alpha = alpha0;
    NumericVector newx;
    double newfx = fx;

    for(int i = 0; i < max_linesearch; i++)
    {
        newx = x + alpha * direction;
        // Test feasibility
        const bool feasible = objective.feas(newx);
        if(!feasible)
        {
            alpha *= scaling;
            continue;
        }
        // Test function value
        newfx = objective.objfn(newx);
        if(newfx < fx)
            return List::create(
                Rcpp::Named("step") = alpha,
                Rcpp::Named("newx") = newx,
                Rcpp::Named("newfx") = newfx
            );

        alpha *= scaling;
    }

    // This function will early return if a proper step size is found
    // If no suitable alpha is obtained, return the initial x
    Rcpp::warning("line search failed, returning the initial x");
    return List::create(
        Rcpp::Named("step") = 0.0,
        Rcpp::Named("newx") = x,
        Rcpp::Named("newfx") = fx
    );
}

inline double vec_norm(NumericVector x)
{
    return std::sqrt(Rcpp::sum(x * x));
}

// -inv(H) * g
inline NumericVector compute_direction(NumericMatrix H, NumericVector g)
{
    const int n = H.nrow();
    NumericVector res(n);
    MapMat HH = Rcpp::as<MapMat>(H);
    MapVec gg = Rcpp::as<MapVec>(g);
    MapVec dd = Rcpp::as<MapVec>(res);

    Eigen::LLT<MatrixXd> solver(HH);
    if(solver.info() != Eigen::Success)
        Rcpp::stop("the Hessian matrix is singular");
    dd.noalias() = -solver.solve(gg);
    return res;
}

// Constrained Newton method
List constr_newton(
    NumericVector x0, Objective& objective, int max_iter = 100,
    int max_linesearch = 30, double eps_rel = 1e-5, bool verbose = false
)
{
    NumericVector x = Rcpp::clone(x0);
    double fx = objective.objfn(x);
    NumericVector grad = objective.grad(x);
    double xnorm = vec_norm(x);
    double xgrad = vec_norm(grad);
    if(verbose)
        Rcpp::Rcout << "Newton iter = 0, fx = " << fx << ", ||grad|| = " << xgrad << std::endl;

    // If gradient is close to zero, then early exit
    if(xgrad <= eps_rel * std::max(1.0, xnorm))
        return List::create(
            Rcpp::Named("x") = x,
            Rcpp::Named("fn") = fx,
            Rcpp::Named("grad") = grad
        );

    NumericMatrix hess = objective.hessian(x);
    NumericVector direction = compute_direction(hess, grad);

    for(int i = 0; i < max_iter; i++)
    {
        List lns = line_search(1.0, x, fx, direction, objective, max_linesearch, 0.5);
        NumericVector newx = lns["newx"];
        double newfx = Rcpp::as<double>(lns["newfx"]);

        const double oldxnorm = xnorm;
        const double xdiff = vec_norm(newx - x);
        x = newx;
        fx = newfx;
        List gd = objective.direction(x);
        grad = gd["grad"];
        direction = gd["direction"];
        xnorm = vec_norm(x);
        xgrad = vec_norm(grad);

        if(verbose)
            Rcpp::Rcout << "Newton iter = " << i + 1 << ", fx = " << fx << ", ||grad|| = " << xgrad << std::endl;

        if(xdiff <= eps_rel * oldxnorm || xgrad <= eps_rel * std::max(1.0, xnorm))
            break;
    }

    return List::create(
            Rcpp::Named("x") = x,
            Rcpp::Named("fn") = fx,
            Rcpp::Named("grad") = grad
        );
}

// [[Rcpp::export]]
NumericMatrix opt_x_cpp(
    NumericMatrix X0, NumericMatrix Y, SEXP B, SEXP Z,
    NumericMatrix A, Environment family,
    NumericVector s, NumericVector gamma,
    int verbose = 0
)
{
    // Get dimensions
    const int n = A.nrow();
    const int p = A.ncol();
    const int k = X0.ncol();
    NumericMatrix X(n, k);
    NumericVector Xi(k);
    NumericVector Ai(p);

    // Optimize each row of X
    int r = 0;
    if(Z != R_NilValue)
    {
        NumericMatrix Z_data(Z);
        r = Z_data.ncol();
    }
    NumericVector Zi_data(r);

    for(int i = 0; i < n; i++)
    {
        if(verbose >= 2)
            Rcpp::Rcout << "===== Optimizing Row " << i + 1 << " of X =====" << std::endl;

        // Determine Zi
        SEXP Zi = R_NilValue;
        if(Z != R_NilValue)
        {
            NumericMatrix Z_data(Z);
            for(int j = 0; j < r; j++)
                Zi_data[j] = Z_data(i, j);
            Zi = Zi_data;
        }

        // Prepare Xi
        for(int j = 0; j < k; j++)
        {
            Xi[j] = X0(i, j);
        }
        // Prepare Ai
        for(int j = 0; j < p; j++)
        {
            Ai[j] = A(i, j);
        }
        // Objective function, gradient, Hessian, and feasibility
        ObjectiveX obj(Y, B, Zi, Ai, family, s[i], gamma);

        // Run optimizer
        List opt = constr_newton(Xi, obj, 100, 30, 0.001, (verbose >= 3));

        // Extract result
        NumericVector optx = opt["x"];
        for(int j = 0; j < k; j++)
            X(i, j) = optx[j];

        if(verbose >= 3)
            Rcpp::Rcout << "==========" << std::endl << std::endl;
    }

    return X;
}
