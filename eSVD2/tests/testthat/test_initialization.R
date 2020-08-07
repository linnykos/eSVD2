context("Test initialization")

## .matrix_completion is correct

test_that(".matrix_completion works", {
 set.seed(10)
 dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
 dat[sample(1:prod(dim(dat)), 10)] <- NA
 res <- .matrix_completion(dat, k = 2)
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(dat)))
})

########################

## .determine_initial_matrix is correct

test_that(".determine_initial_matrix works", {
 set.seed(10)
 dat <- matrix(1:40, nrow = 10, ncol = 4)
 
 res <- .determine_initial_matrix(dat, family = "poisson")
 
 expect_true(all(sort(names(res)) == sort(c("nat_mat", "domain"))))
 expect_true(length(res$domain) == 2)
 expect_true(res$domain[1] <= res$domain[2])
 expect_true(all(dim(res$nat_mat) == dim(dat)))
})


test_that(".determine_initial_matrix respects max_val", {
 set.seed(10)
 dat <- matrix(1:40, nrow = 10, ncol = 4)
         
 max_val <- 2
 res <- .determine_initial_matrix(dat, family = "poisson", max_val = max_val)
  
 expect_true(all(res$nat_mat >= 0))
 expect_true(all(res$nat_mat <= max_val))
})

################################

## .fix_rank_defficiency is correct

test_that(".fix_rank_defficiency works", {
 set.seed(10)
 x_mat <- abs(matrix(stats::rnorm(12), nrow = 4, ncol = 3))
 y_mat <- abs(matrix(stats::rnorm(15), nrow = 5, ncol = 3))
 
 res <- .fix_rank_defficiency(x_mat, y_mat, domain = c(0.1, 100))
 
 expect_true(length(res) == 2)
 expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat"))))
 
 nat_mat <- res$x_mat %*% t(res$y_mat)
 expect_true(all(nat_mat > 0))
 expect_true(Matrix::rankMatrix(nat_mat) == 3)
})

test_that(".fix_rank_defficiency can actually fix the rank", {
 set.seed(10)
 x_mat <- abs(matrix(stats::rnorm(12), nrow = 4, ncol = 3))
 y_mat <- abs(matrix(stats::rnorm(15), nrow = 5, ncol = 3))
 y_mat[,3] <- 0
 
 nat_mat <- x_mat %*% t(y_mat)
 
 res <- .fix_rank_defficiency(x_mat, y_mat, domain = c(min(nat_mat)/2, 100))
 
 nat_mat <- res$x_mat %*% t(res$y_mat)
 expect_true(all(nat_mat > 0))
 expect_true(Matrix::rankMatrix(nat_mat) == 3)
})

#############################

## initialize_esvd is correct

test_that("initialize_esvd works", {
 set.seed(10)
 dat <- matrix(1:40, nrow = 10, ncol = 4)
 
 res <- initialize_esvd(dat, k = 2, family = "poisson")
 
 expect_true(is.list(res))
 expect_true(class(res) == "eSVD")
 expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat"))))
 expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
 expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
})