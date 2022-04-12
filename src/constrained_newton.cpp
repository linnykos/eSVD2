#include <RcppEigen.h>
#include "objective.h"

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::List;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using MapMat = Eigen::Map<MatrixXd>;
using MapVec = Eigen::Map<VectorXd>;

// C++ translation of line_search() in constrained_newton_lbfgs.R
List line_search(
    double alpha0, NumericVector x, double fx, NumericVector direction,
    Objective& objective, int max_linesearch, double scaling = 0.5
)
{
    double alpha = alpha0;
    const int n = x.length();
    NumericVector newx = NumericVector(Rcpp::no_init_vector(n));
    double newfx = fx;

    MapVec newx_ = Rcpp::as<MapVec>(newx);
    MapVec x_ = Rcpp::as<MapVec>(x);
    MapVec direction_ = Rcpp::as<MapVec>(direction);

    for(int i = 0; i < max_linesearch; i++)
    {
        newx_.noalias() = x_ + alpha * direction_;
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

inline double vec_norm(NumericVector x_)
{
    MapVec x = Rcpp::as<MapVec>(x_);
    return x.norm();
}

inline double vec_dist(NumericVector x_, NumericVector y_)
{
    MapVec x = Rcpp::as<MapVec>(x_);
    MapVec y = Rcpp::as<MapVec>(y_);
    return (x - y).norm();
}

// -inv(H) * g
inline NumericVector compute_direction(NumericMatrix H_, NumericVector g_)
{
    const int n = H_.nrow();
    MapMat H = Rcpp::as<MapMat>(H_);
    MapVec g = Rcpp::as<MapVec>(g_);

    NumericVector res = NumericVector(Rcpp::no_init_vector(n));
    MapVec direc = Rcpp::as<MapVec>(res);

    Eigen::LLT<MatrixXd> solver(H);
    if(solver.info() != Eigen::Success)
    {
        // Fall back to gradient direction if Hessian is singular
        direc.noalias() = -g;
    } else {
        direc.noalias() = -solver.solve(g);
    }
    return res;
}

// C++ translation of constr_newton() in constrained_newton_lbfgs.R
List constr_newton(
    const VectorXd& x0, Objective& objective, int max_iter = 100,
    int max_linesearch = 30, double eps_rel = 1e-5, bool verbose = false
)
{
    NumericVector x(x0.data(), x0.data() + x0.size());
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
        const double xdiff = vec_dist(newx, x);
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
