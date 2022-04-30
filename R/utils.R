.l2norm <- function(x){sqrt(sum(x^2))}

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

## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))

  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]

  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}


## [[note to self: allow target_obj to be a Seurat obj]]
.append_rowcolnames <- function(bool_colnames,
                                bool_rownames,
                                source_obj,
                                target_obj){
  stopifnot(bool_colnames | bool_rownames)

  if(bool_rownames){
    if(inherits(source_obj, c("matrix", "dgCMatrix"))){
      rowname_vec <- rownames(source_obj)
    } else {
      stop("Not valid class for source_obj")
    }

    if(inherits(target_obj, c("matrix", "dgCMatrix")) & length(rowname_vec) > 0){
      rownames(target_obj) <- rowname_vec
    }
  }

  if(bool_colnames){
    if(inherits(source_obj, c("matrix", "dgCMatrix"))){
      colname_vec <- colnames(source_obj)
    } else {
      stop("Cannot identify class for target_obj")
    }

    if(inherits(target_obj, c("matrix", "dgCMatrix")) & length(colname_vec) > 0){
      colnames(target_obj) <- colname_vec
    }
  }

  target_obj
}

