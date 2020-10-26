# Minimize f(x) subject to Ax>b
# x: [p x 1]
# A: [m x p]
# b: [m x 1]
# Ax>b <=> a[i]'x>b[i], i=1,2,...,m
# f(x) is unbounded if any a[i]'x<=b[i]
#
# gr(x) is the gradient of f(x)
# hn(x) is the Hessian matrix of f(x)
#
# We use Newton's method to solve the optimization problem
# x[k+1]=x[k]-alpha*inv(H)*g, where g=gr(x), H=hn(x)
# alpha is chosen to satisfy the feasibility condition Ax>b

# Convenience function, ||x||
vnorm <- function(x) sqrt(sum(x^2))

# Line search function
#
# Find a step size such that the feasibility condition Ax>b holds and
# the function value decreases
#
# alpha0:         initial step size
# x:              current x value
# fx:             current objective function value, f(x)
# direction:      newx = x + alpha * direction
# f:              objective function
# feas:           function to test feasibility, returning TRUE if Ax>b
# max_linesearch: maximum number of line search tries
# scaling:        decrease factor of alpha
# ...:            additional arguments passed to f
line_search <- function(alpha0, x, fx, direction, f, feas,
                        max_linesearch, scaling = 0.5, ...)
{
    alpha <- alpha0
    for(i in seq_len(max_linesearch))
    {
        newx <- x + alpha * direction
        # Test feasibility
        feasible <- feas(newx, ...)
        if(!feasible)
        {
            alpha <- alpha * scaling
            next
        }
        # Test function value
        newfx = f(newx, ...)
        if(newfx < fx)
            return(list(step = alpha, newx = newx, newfx = newfx))

        alpha <- alpha * scaling
    }
    # This function will early return if a proper step size is found
    # Give an error if the function reaches here
    stop("line search failed")
}

constr_newton <- function(x0, f, gr, hn, feas,
                          max_iter = 100, max_linesearch = 30,
                          eps_rel = 1e-5, verbose = FALSE, ...)
{
    x <- x0
    fx <- f(x, ...)
    grad <- gr(x, ...)
    xnorm <- vnorm(x)
    xgrad <- vnorm(grad)
    if(verbose)
        cat(sprintf("Newton iter = 0, fx = %f, ||grad|| = %f\n", fx, xgrad))

    # If gradient is close to zero, then early exit
    if(xgrad <= eps_rel * max(1, xnorm))
        return(list(x = x, fn = fx, grad = grad))

    for(i in seq_len(max_iter))
    {
        Hess <- hn(x, ...)
        direction <- -solve(Hess, grad)
        lns <- line_search(1, x, fx, direction, f, feas,
                           max_linesearch, scaling = 0.5, ...)
        step <- lns$step
        newx <- lns$newx
        newfx <- lns$newfx

        oldxnorm <- xnorm
        xdiff <- vnorm(newx - x)
        x <- newx
        fx <- newfx
        grad <- gr(x, ...)
        xnorm <- vnorm(x)
        xgrad <- vnorm(grad)

        if(verbose)
            cat(sprintf("Newton iter = %d, fx = %f, ||grad|| = %f\n", i, fx, xgrad))

        if(xdiff < eps_rel * oldxnorm || xgrad <= eps_rel * max(1, xnorm))
            break
    }

    list(x = x, fn = fx, grad = grad)
}



# Optimization for eSVD

# Optimize u given v
opt_u_given_v <- function(u0, v_mat, dat, family_funcs, nuisance_param_vec,
                          library_size_vec, verbose = 0, ...)
{
    n <- nrow(dat)
    p <- ncol(dat)
    u_mat <- u0
    for(i in 1:n)
    {
        if(verbose >= 2)
            cat("===== Optimizing Row ", i, " of U =====\n", sep = "")
        opt <- constr_newton(
            u0[i, ], family_funcs$objfn, family_funcs$grad, family_funcs$hessian, family_funcs$feas,
            eps_rel = 1e-3, verbose = (verbose >= 3),
            other_mat = v_mat, dat_vec = dat[i, ],
            nuisance_param_vec = nuisance_param_vec,
            library_size = library_size_vec[i], ...
        )
        u_mat[i, ] <- opt$x
        if(verbose >= 3)
            cat("==========\n\n")
    }
    u_mat
}

# Optimize v given u
opt_v_given_u <- function(v0, u_mat, dat, family_funcs, nuisance_param_vec,
                          library_size_vec, verbose = 0, ...)
{
    n <- nrow(dat)
    p <- ncol(dat)
    v_mat <- v0
    for(j in 1:p)
    {
        if(verbose >= 2)
            cat("===== Optimizing Row ", j, " of V =====\n", sep = "")
        opt <- constr_newton(
            v0[j, ], family_funcs$objfn, family_funcs$grad, family_funcs$hessian, family_funcs$feas,
            eps_rel = 1e-3, verbose = (verbose >= 3),
            other_mat = u_mat, dat_vec = dat[, j],
            nuisance_param_vec = nuisance_param_vec,
            library_size = library_size_vec, ...
        )
        v_mat[j, ] <- opt$x
        if(verbose >= 3)
            cat("==========\n\n")
    }
    v_mat
}

#######################

#' Main optimization function for eSVD
#'
#' @param x_init initial estimate of a \eqn{n \times k} matrix for the \eqn{k}-dimensional
#' embedding for the \eqn{n} cells
#' @param y_init initial estimate of a \eqn{p \times k} matrix for the \eqn{k}-dimensional
#' embedding for the \eqn{p} genes
#' @param dat dataset where the \eqn{n} rows represent cells and \eqn{p} columns represent genes
#' @param family A character string, one of \code{"gaussian"}, \code{"exponential"},
#'               \code{"poisson"}, \code{"neg_binom"}, \code{"curved_gaussian"},
#'               and \code{"bernoulli"}.
#' @param nuisance_param_vec either \code{NA} or a single numeric or a length-\eqn{p}
#' vector of numerics representing nuisance parameters (for \code{family = "neg_binom"} or
#' \code{family = "curved_gausian"}).
#' @param library_size_vec either \code{NA} or a single numeric (default is \code{1}) or
#' a length-\eqn{n} vector of numerics.
#' If \code{NA}, the library size will be estimated.
#' @param max_iter positive integer
#' @param verbose no-negtaive integer
#' @param ... additional parameters, not used
#'
#' @return a list with elements \code{x_mat} and \code{y_mat} (and others), representing the two
#' latent matrices
#' @export
opt_esvd <- function(x_init, y_init, dat, family, nuisance_param_vec = NA,
                     library_size_vec = 1, max_iter = 100, verbose = 0, ...)
{
    stopifnot(nrow(x_init) == nrow(dat), nrow(y_init) == ncol(dat), ncol(x_init) == ncol(y_init))
    stopifnot(is.character(family))

    n <- nrow(dat); p <- nrow(dat); k <- ncol(x_init)
    family_funcs <- .string_to_distr_funcs(family)
    library_size_vec <- .parse_library_size(dat, library_size_vec)

    if(all(!is.na(nuisance_param_vec)) & length(nuisance_param_vec) == 1) {
        nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
    }

    u_mat <- x_init
    v_mat <- y_init
    losses <- c()
    for(i in 1:max_iter)
    {
        if(verbose >= 1)
            cat("========== eSVD Iter ", i, " ==========\n\n", sep = "")
        # Optimize u given v
        u_mat <- opt_u_given_v(u_mat, v_mat, dat, family_funcs,
                               nuisance_param_vec,
                               library_size_vec, verbose, ...)
        # Orthogonalize u
        # u_mat = sqrt(n) * svd(u_mat)$u
        # Optimize v given u
        v_mat <- opt_v_given_u(v_mat, u_mat, dat, family_funcs,
                               nuisance_param_vec,
                               library_size_vec, verbose, ...)
        # Loss function
        loss <- family_funcs$objfn_all(dat, u_mat, v_mat, nuisance_param_vec,
                                 library_size_vec, ...)
        losses <- c(losses, loss)
        if(verbose >= 1)
            cat("========== eSVD Iter ", i, ", loss = ", loss, " ==========\n\n", sep = "")
        # Orthogonalize v
        # v_mat = sqrt(p) * svd(v_mat)$u
    }

    list(x_mat = u_mat, y_mat = v_mat, loss = losses,
         nuisance_param_vec = nuisance_param_vec, library_size_vec = library_size_vec)
}
