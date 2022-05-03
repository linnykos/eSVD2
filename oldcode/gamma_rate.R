# Gamma distribution with rate parameter beta
# Y ~ Gamma(alpha, beta)
# p(y) = b^a / Gamma(a) * y^(a-1) * exp(-b * y)
# log[p(y)] = a * log(b) - logGamma(a) + (a-1) * log(y) - b * y
# mean = a / b
#
# We assume Xi|lambdai ~ Pois(si * lambdai) and
# lambdai ~ Gamma with mean=mui and rate=beta, i.e.,
# lambdai ~ Gamma(mui * beta, beta)
# log[p(xi|lambdai)] = xi * log(si * lambdai) - log(xi!) - si * lambdai * xi
# log[p(lambdai)] = mui * b * log(b) - logGamma(mui * b)
#                   + (mui * b - 1) * log(lambdai) - b * lambdai
#
# Marginally, Xi ~ NB(ri, p), where
# ri = mui * b, p = 1 / (1 + b / si) = si / (si + b)
# log[p(xi)] = logGamma(ri + xi) - log(xi!) - logGamma(ri)
#              + xi * log(p) + ri * log(1 - p)
#            = logGamma(mui * b + xi) - log(xi!)
#              - logGamma(mui * b)
#              + xi * log(si / (si + b))
#              + mui * b * log(b / (si + b))
#            = xi * log(si) - log(xi!)
#              + mui * b * log(b) - logGamma(mui * b)
#              + logGamma(mui * b + xi)
#              - (mui * b + xi) * log(si + b)
#
# l(b) = log[p(xi; b)]
# [b * log(b)]' = log(b) + 1
# [logGamma(x)]' = digamma(x)
# [digamma(x)]' = trigamma(x)
#
# [l(b)]' = mui * (log(b) + 1) - mui * digamma(mui * b)
#           + mui * digamma(mui * b + xi)
#           - mui * log(si + b) - (mui * b + xi) / (si + b)
# [l(b)]'' = mui / b - mui^2 * trigamma(mui * b)
#            + mui^2 * trigamma(mui * b + xi)
#            - mui / (si + b) - (mui * si - xi) / (si + b)^2



# # Examples
# objfn <- function(b, xi, mui, si) {
#   loglik <- xi * log(si) - lgamma(xi + 1) +
#     mui * b * log(b) - lgamma(mui * b) + lgamma(mui * b + xi) -
#     (mui * b + xi) * log(si + b)
#   mean(loglik)
# }
# gradfn <- function(b, xi, mui, si) {
#   grad <- mui * (log(b) + 1) - mui * digamma(mui * b) +
#     mui * digamma(mui * b + xi) -
#     mui * log(si + b) - (mui * b + xi) / (si + b)
#   mean(grad)
# }
# hessfn <- function(b, xi, mui, si) {
#   hess <- mui / b - mui^2 * trigamma(mui * b) +
#     mui^2 * trigamma(mui * b + xi) -
#     mui / (si + b) - (mui * si - xi) / (si + b)^2
#   mean(hess)
# }
#
# library(MASS)
# mod <- glm.nb(Days ~ .^2, data = quine)
# xi <- quine$Days
# mui <- fitted(mod)
# si <- rep(1, length(xi))
#
# curve(sapply(x, objfn, xi = xi, mui = mui, si = si), 0.01, 10)
# curve(sapply(x, objfn, xi = xi, mui = mui, si = si), 0.01, 1)
# curve(sapply(x, gradfn, xi = xi, mui = mui, si = si), 0.01, 1)
# abline(h = 0, col = "red")
# curve(sapply(x, hessfn, xi = xi, mui = mui, si = si), 0.2, 5)
# # Log-scale
# curve(sapply(exp(x), objfn, xi = xi, mui = mui, si = si), -5, 5)
# curve(exp(x) * sapply(exp(x), gradfn, xi = xi, mui = mui, si = si), -5, 5)
# abline(h = 0, col = "red")
#
# yeast <- data.frame(cbind(numbers = 0:5, fr = c(213, 128, 37, 18, 3, 1)))
# fit <- glm.nb(numbers ~ 1, weights = fr, data = yeast)
# summary(fit)
# mui <- fitted(fit)
# xi <- yeast$numbers
# si <- yeast$fr
