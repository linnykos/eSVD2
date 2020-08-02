.diag_matrix <- function(vec){
 k <- length(vec)
 if(k == 1) {
  matrix(vec, 1, 1)
 } else {
  diag(vec)
 }
}

.l2norm <- function(vec){
 sqrt(sum(vec^2))
}