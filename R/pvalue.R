pvalue_pairwise <- function(mat, membership_vec, verbose = T){
  stopifnot(is.factor(membership_vec))
  uniq_val <- unique(membership_vec)

  list_idx <- lapply(uniq_val, function(x){
    which(membership_vec == x)
  })

  res <- matrix(NA, nrow = length(unique(membership_vec)), ncol = ncol(mat))
  for(j in 1:ncol(mat)){
    for(i in 1:length(list_idx)){
      res[i,j] <- stats::wilcox.test(x = mat[list_idx[[i]],j], y = mat[-list_idx[[i]],j],
                         alternative = "greater")$p.value
    }
  }

  colnames(res) <- colnames(mat)
  rownames(res) <- uniq_val

  res
}
