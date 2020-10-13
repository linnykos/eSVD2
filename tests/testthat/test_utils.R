context("Test utility functions")

## .intersect_intervals is correct

test_that(".intersect_intervals works", {
 vec1 <- c(1,10)
 vec2 <- c(5,15)
 res <- .intersect_intervals(vec1, vec2)
 
 expect_true(all(res == c(5,10)))
})

test_that(".intersect_intervals works correctly for randomly generated instances", {
 trials <- 100
 
 bool_vec <- sapply(1:trials, function(trial){
  set.seed(trial)
  vals <- stats::rnorm(4)
  idx <- sample(1:4, 2, replace = F)
  
  vec1 <- sort(vals[idx])
  vec2 <- sort(vals[-idx])
  
  res <- .intersect_intervals(vec1, vec2)
  
  res[1] >= vec1[1] & res[1] >= vec2[1] & res[2] <= vec1[2] & res[2] <= vec2[2]
 })
 
 expect_true(all(bool_vec))
})
