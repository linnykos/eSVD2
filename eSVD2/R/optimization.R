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
