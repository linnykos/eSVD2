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

.intersect_intervals <- function(vec1, vec2){
  stopifnot(length(vec1) == 2, length(vec2) == 2, vec1[1] <= vec1[2], vec2[1] <= vec2[2])
  
  if(vec1[1] >= vec2[1] & vec1[2] <= vec2[2]){
    return(vec1)
  } else if(vec2[1] >= vec1[1] & vec2[2] <= vec1[2]){
    return(vec2)
  } else {
    new_vec <- rep(NA, 2)
    new_vec[1] <- max(vec1[1], vec2[1])
    new_vec[2] <- min(vec1[2], vec2[2])
  }
  
  new_vec
}