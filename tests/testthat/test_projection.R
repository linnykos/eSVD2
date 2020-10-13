context("Test projection")

## .projection_kmeans is correct

test_that(".projection_kmeans works", {
 set.seed(10)
 mat <- matrix(1:30, 6, 5)
 res <- .projection_kmeans(mat, k = 2, domain = NA, row = T)
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(mat)))
})