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
 res <- .determine_initial_ matrix(dat, family = "poisson", max_val = max_val)
 
 expect_true(all(res$nat_mat >= 0))
 expect_true(all(res$nat_mat <= max_val))
})