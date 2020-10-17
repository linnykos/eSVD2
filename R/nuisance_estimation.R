initialize_nuisance_param <- function(dat, pred_mat, family, method = "mom"){
  stopifnot(all(dim(dat) == dim(pred_mat)), family %in% c("curved_gaussian", "neg_binom"))


}
