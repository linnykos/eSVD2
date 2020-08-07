rm(list=ls())

devtools::install_github("linnykos/eSVD2", subdir = "eSVD2")

library(eSVD2)
session_info()

# poisson example
## generate matrix of natural parameters
set.seed(10)
n <- 100; p <- 150; k <- 5
x_mat <- matrix(abs(stats::rnorm(n*k)), nrow = n, ncol = k)
y_mat <- matrix(abs(stats::rnorm(p*k)), nrow = p, ncol = k)
nat_mat <- (x_mat %*% t(y_mat))/10

## generate data
dat <- eSVD2::generate_data(nat_mat, family = "poisson", nuisance_param_vec = NA, library_size_vec = NA)
quantile(dat)
image(eSVD2:::.rotate(dat), asp = T)

## determine initialization
init_res <- eSVD2::initialize_esvd(dat, k = k, family = "poisson", nuisance_param_vec = NA, library_size_vec = NA,
                                   config = eSVD2::initalization_default())

names(init_res)
dim(init_res$x_mat)
dim(init_res$y_mat)
range(init_res$x_mat %*% t(init_res$y_mat))
eSVD2:::.check_domain(nat_mat, init_res$domain)

######################

# curved gaussian example
## generate matrix of natural parameters
set.seed(15)
n <- 100; p <- 150; k <- 5
x_mat <- matrix(abs(stats::rnorm(n*k)), nrow = n, ncol = k)
y_mat <- matrix(abs(stats::rnorm(p*k)), nrow = p, ncol = k)
nat_mat <- x_mat %*% t(y_mat)

## generate data
dat <- eSVD2::generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2, library_size_vec = NA,
                            tol = 1e-3)
quantile(dat)
image(eSVD2:::.rotate(dat), asp = T)

## determine initialization
init_res <- eSVD2::initialize_esvd(dat, k = k, family = "curved_gaussian", nuisance_param_vec = 2, library_size_vec = NA,
                                   config = eSVD2::initalization_default())

names(init_res)
dim(init_res$x_mat)
dim(init_res$y_mat)
range(init_res$x_mat %*% t(init_res$y_mat))
eSVD2:::.check_domain(nat_mat, init_res$domain)


