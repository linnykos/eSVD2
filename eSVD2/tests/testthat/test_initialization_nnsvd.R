context("Test initialization method - NNDSVD")

## .initialization_nnsvd is correct

test_that(".initialization_nndsvd works", {
 set.seed(10)
 dat <- matrix(sample(1:40), 10, 4)
 res <- .initialization_nnsvd(dat, k = 2)
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(dat)))
})
