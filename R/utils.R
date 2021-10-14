.l2norm <- function(x){sqrt(sum(x^2))}

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

.rotate = function(a) { t(a[nrow(a):1,]) }

.colorRamp_custom <- function(vec1, vec2, length, luminosity){
  mat <- matrix(0, nrow = length, ncol = 3)
  for(i in 1:3){
    mat[,i] <- seq(vec1[i], vec2[i], length.out = length)
  }

  if(luminosity){
    luminosity_vec <- apply(mat, 1, function(x){
      0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
    })

    target_luminosity <- mean(c(luminosity_vec[1], luminosity_vec[length]))

    mat <- t(sapply(1:nrow(mat), function(x){
      factor <- min(c(target_luminosity/luminosity_vec[x], 1/mat[x,]))
      mat[x,] * factor
    }))
  }

  apply(mat, 1, function(x){
    grDevices::rgb(x[1], x[2], x[3])
  })
}

# for diag(vec) %*% mat
.mult_vec_mat <- function(vec, mat){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == nrow(mat))
  vec * mat
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == ncol(mat))
  mat * rep(vec, rep(nrow(mat), length(vec)))
}
