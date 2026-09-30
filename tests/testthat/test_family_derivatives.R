# UNIT_TEST_PLAN.md section 3.2 -- analytic derivatives for all seven families.
#
# This is the single largest block in the plan and the cheapest per unit of risk
# removed: an analytic-derivative error currently produces a silently
# mis-converged fit with NO symptom -- the loss still decreases and the
# optimizer still terminates.
#
# `numDeriv` has been in `Suggests` since the beginning and has never been used.

test_that("T-CPP-FAM-00: all seven families construct and report their name", {
  for(family_name in .all_families()){
    family_obj <- esvd_family(family_name)
    expect_equal(family_obj$name, family_name, info = family_name)
    expect_true(is.function(family_obj$feasibility), info = family_name)
  }
})

test_that("T-CPP-FAM-01: esvd_family rejects an unknown family, listing the valid ones", {
  expect_error(esvd_family("not_a_family"))
})

test_that("T-CPP-FAM-02: the feasible points really are feasible", {
  for(family_name in .all_families()){
    point <- .feasible_point(family_name)
    family_obj <- esvd_family(family_name)

    # If this fails, every gradient test below is testing an infeasible point
    # and its result means nothing.
    expect_true(family_obj$feasibility(point$nat_mat), info = family_name)
  }
})

## The core of section 3.2: the analytic gradients and Hessians against
## `numDeriv`, for every family. A mismatch here is invisible in normal use --
## the loss still decreases and the optimizer still terminates.
##
## `XC` is the combined [X, C] matrix and `YZ` the combined [Y, Z]; with no
## covariates they are just X and Y, and `k` is the full column count.
test_that("T-CPP-FAM-03: grad_Xi_r matches a numerical gradient, all seven families", {
  skip_if_not_installed("numDeriv")

  for(family_name in .all_families()){
    point <- .feasible_point(family_name)
    family_obj <- esvd_family(family_name)
    loader <- data_loader(point$dat)
    cell_idx <- 1
    num_factors <- ncol(point$x_mat)

    objective_fn <- function(x_vec){
      objfn_Xi_r(XCi = x_vec, YZ = point$y_mat, k = num_factors,
                 loader = loader, row_ind = cell_idx - 1L,
                 family = family_obj, si = point$s_vec[cell_idx],
                 gamma = point$gamma_vec, l2penx = 0)
    }

    analytic_grad <- grad_Xi_r(XCi = point$x_mat[cell_idx, ],
                               YZ = point$y_mat, k = num_factors,
                               loader = loader, row_ind = cell_idx - 1L,
                               family = family_obj,
                               si = point$s_vec[cell_idx],
                               gamma = point$gamma_vec, l2penx = 0)
    numeric_grad <- numDeriv::grad(objective_fn, point$x_mat[cell_idx, ])

    expect_equal(as.numeric(analytic_grad), as.numeric(numeric_grad),
                 tolerance = 1e-5, info = family_name)
  }
})

test_that("T-CPP-FAM-04: hessian_Xi_r matches a numerical Hessian, all seven families", {
  skip_if_not_installed("numDeriv")

  for(family_name in .all_families()){
    point <- .feasible_point(family_name)
    family_obj <- esvd_family(family_name)
    loader <- data_loader(point$dat)
    cell_idx <- 1
    num_factors <- ncol(point$x_mat)

    objective_fn <- function(x_vec){
      objfn_Xi_r(XCi = x_vec, YZ = point$y_mat, k = num_factors,
                 loader = loader, row_ind = cell_idx - 1L,
                 family = family_obj, si = point$s_vec[cell_idx],
                 gamma = point$gamma_vec, l2penx = 0)
    }

    analytic_hessian <- hessian_Xi_r(XCi = point$x_mat[cell_idx, ],
                                     YZ = point$y_mat, k = num_factors,
                                     loader = loader, row_ind = cell_idx - 1L,
                                     family = family_obj,
                                     si = point$s_vec[cell_idx],
                                     gamma = point$gamma_vec, l2penx = 0)
    numeric_hessian <- numDeriv::hessian(objective_fn, point$x_mat[cell_idx, ])

    expect_equal(as.numeric(as.matrix(analytic_hessian)),
                 as.numeric(numeric_hessian),
                 tolerance = 1e-4, info = family_name)
  }
})

test_that("T-CPP-FAM-05: grad_YZj_r matches a numerical gradient, all seven families", {
  skip_if_not_installed("numDeriv")

  for(family_name in .all_families()){
    point <- .feasible_point(family_name)
    family_obj <- esvd_family(family_name)
    loader <- data_loader(point$dat)
    gene_idx <- 1
    num_factors <- ncol(point$y_mat)
    yz_ind <- seq_len(num_factors) - 1L

    objective_fn <- function(yz_vec){
      objfn_YZj_r(XC = point$x_mat, YZj = yz_vec, k = num_factors,
                  YZind = yz_ind, loader = loader, col_ind = gene_idx - 1L,
                  family = family_obj, s = point$s_vec,
                  gammaj = point$gamma_vec[gene_idx],
                  l2peny = 0, l2penz = 0)
    }

    analytic_grad <- grad_YZj_r(XC = point$x_mat,
                                YZj = point$y_mat[gene_idx, ],
                                k = num_factors, YZind = yz_ind,
                                loader = loader, col_ind = gene_idx - 1L,
                                family = family_obj, s = point$s_vec,
                                gammaj = point$gamma_vec[gene_idx],
                                l2peny = 0, l2penz = 0)
    numeric_grad <- numDeriv::grad(objective_fn, point$y_mat[gene_idx, ])

    expect_equal(as.numeric(analytic_grad), as.numeric(numeric_grad),
                 tolerance = 1e-5, info = family_name)
  }
})

test_that("T-CPP-FAM-06: hessian_YZj_r matches a numerical Hessian, all seven families", {
  skip_if_not_installed("numDeriv")

  for(family_name in .all_families()){
    point <- .feasible_point(family_name)
    family_obj <- esvd_family(family_name)
    loader <- data_loader(point$dat)
    gene_idx <- 1
    num_factors <- ncol(point$y_mat)
    yz_ind <- seq_len(num_factors) - 1L

    objective_fn <- function(yz_vec){
      objfn_YZj_r(XC = point$x_mat, YZj = yz_vec, k = num_factors,
                  YZind = yz_ind, loader = loader, col_ind = gene_idx - 1L,
                  family = family_obj, s = point$s_vec,
                  gammaj = point$gamma_vec[gene_idx],
                  l2peny = 0, l2penz = 0)
    }

    analytic_hessian <- hessian_YZj_r(XC = point$x_mat,
                                      YZj = point$y_mat[gene_idx, ],
                                      k = num_factors, YZind = yz_ind,
                                      loader = loader, col_ind = gene_idx - 1L,
                                      family = family_obj, s = point$s_vec,
                                      gammaj = point$gamma_vec[gene_idx],
                                      l2peny = 0, l2penz = 0)
    numeric_hessian <- numDeriv::hessian(objective_fn, point$y_mat[gene_idx, ])

    expect_equal(as.numeric(as.matrix(analytic_hessian)),
                 as.numeric(numeric_hessian),
                 tolerance = 1e-4, info = family_name)
  }
})

test_that("T-CPP-FAM-07: the l2 penalty enters the gradient as 2*l2pen*x", {
  point <- .feasible_point("poisson")
  family_obj <- esvd_family("poisson")
  loader <- data_loader(point$dat)
  cell_idx <- 1
  num_factors <- ncol(point$x_mat)
  l2pen_val <- 0.35

  grad_unpenalized <- grad_Xi_r(XCi = point$x_mat[cell_idx, ],
                                YZ = point$y_mat, k = num_factors,
                                loader = loader, row_ind = cell_idx - 1L,
                                family = family_obj,
                                si = point$s_vec[cell_idx],
                                gamma = point$gamma_vec, l2penx = 0)
  grad_penalized <- grad_Xi_r(XCi = point$x_mat[cell_idx, ],
                              YZ = point$y_mat, k = num_factors,
                              loader = loader, row_ind = cell_idx - 1L,
                              family = family_obj,
                              si = point$s_vec[cell_idx],
                              gamma = point$gamma_vec, l2penx = l2pen_val)

  # The penalty is what keeps the fit identifiable, so its contribution being
  # exactly right matters as much as the likelihood's. Note the normalization:
  # the objective is a MEAN over genes for a fixed cell, so the penalty
  # gradient carries the same 1/p factor. T-CPP-FAM-08 pins the corresponding
  # mean over cells.
  num_genes <- nrow(point$y_mat)
  expect_equal(as.numeric(grad_penalized - grad_unpenalized),
               as.numeric(2 * l2pen_val * point$x_mat[cell_idx, ] / num_genes),
               tolerance = 1e-8)
})

test_that("T-CPP-FAM-08: objfn_all_r agrees with the sum over cells of objfn_Xi_r", {
  point <- .feasible_point("poisson")
  family_obj <- esvd_family("poisson")
  loader <- data_loader(point$dat)
  num_factors <- ncol(point$x_mat)

  total_objective <- objfn_all_r(XC = point$x_mat, YZ = point$y_mat,
                                 k = num_factors, loader = loader,
                                 family = family_obj, s = point$s_vec,
                                 gamma = point$gamma_vec,
                                 l2penx = 0, l2peny = 0, l2penz = 0)

  per_cell_vec <- sapply(seq_len(nrow(point$x_mat)), function(cell_idx){
    objfn_Xi_r(XCi = point$x_mat[cell_idx, ], YZ = point$y_mat,
               k = num_factors, loader = loader, row_ind = cell_idx - 1L,
               family = family_obj, si = point$s_vec[cell_idx],
               gamma = point$gamma_vec, l2penx = 0)
  })

  # Two routes to one number: the all-at-once objective and the per-cell one
  # the optimizer actually minimizes. If the normalization differed between
  # them, the loss would still decrease monotonically and nothing would notice.
  expect_equal(total_objective, mean(per_cell_vec), tolerance = 1e-8)
})

test_that("T-CPP-FAM-12: bernoulli's dat_to_nat is deliberately lossy", {
  # Section 8 declines to write a round-trip test for `bernoulli`, because the
  # conversion maps to +/-1 regardless of the input magnitude. Assert THAT,
  # rather than asserting something unrelated and calling it covered.
  nat_one <- .dat_to_nat.bernoulli(1)
  nat_many <- .dat_to_nat.bernoulli(5)
  nat_zero <- .dat_to_nat.bernoulli(0)

  expect_equal(nat_one, nat_many)
  expect_false(isTRUE(all.equal(nat_one, nat_zero)))
})
