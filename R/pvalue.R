pvalue_pairwise <- function(mat, membership_vec, celltype_vec,
                            verbose = T){
  stopifnot(is.factor(membership_vec), all(celltype_vec %in% levels(membership_vec)))

  k <- length(celltype_vec)
  list_idx <- lapply(celltype_vec, function(x){
    which(membership_vec == x)
  })

  permn_mat <- gtools::permutations(k, 2)
  res <- matrix(NA, nrow = nrow(permn_mat), ncol = ncol(mat))
  colnames(res) <- colnames(mat)
  rownames(res) <- apply(permn_mat, 1, function(x){paste0(x, collapse = "-")})

  lookup <- data.frame(idx1 = permn_mat[,1], idx2 = permn_mat[,2],
                       type1 = celltype_vec[permn_mat[,1]],
                       type2 = celltype_vec[permn_mat[,2]])

  for(j in 1:ncol(mat)){
    print(j)
    if(verbose && ncol(mat) > 10 && j %% floor(ncol(mat)/10) == 0) cat('*')
    for(i in 1:nrow(permn_mat)){
      res[i,j] <- stats::wilcox.test(x = mat[list_idx[[permn_mat[i,1]]],j],
                                     y = mat[list_idx[[permn_mat[i,2]]],j],
                                     alternative = "greater")$p.value
    }
  }

  list(pval_mat = res, lookup_mat = lookup)
}

pvalue_combine_de <- function(pval_obj){
  uniq_val <- sort(unique(as.numeric(pval_obj$lookup_mat[,1])))
  p <- ncol(pval_obj$pval_mat)

  # intersection-union
  max_pval_mat <- t(sapply(uniq_val, function(i){
    idx <- which(pval_obj$lookup[,1] == i)
    apply(pval_obj$pval_mat[idx,], 2, max)
  }))

  # select the smallest pvalue
  max_pval_vec <- pmin(apply(max_pval_mat*length(uniq_val), 2, min), 1)

  # [[note to self: should I compare this to a beta?]]
}
