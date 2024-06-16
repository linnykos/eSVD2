#' Generate null data
#'
#' This function doesn't allow of much flexbility, and is meant primarily for simple null simulations.
#' Here, `n` (the number of cells) is equal to `cell_per_person*num_individuals`.
#'
#' @param cell_per_person The number of cells per individual
#' @param num_genes       The number of genes in the simulation
#' @param num_individuals The number of individuals
#'
#' @return a list with \code{covariates} (a `matrix` with `n` cells and
#' 4 columns, named `"Intercept"`, `"Log_UMI"`, `"Sex"`, and `"Age"`,
#' which will be used in `eSVD2::initialize_esvd`),
#' \code{df} (a `data.frame` with `num_individuals` rows
#' and 5 columns, named named `"Intercept"`, `"Log_UMI"`, `"Sex"`, `"Age"`, and `"Individual"`),
#' \code{metadata_individual} (a `matrix` with `n` rows and 1 column called `"Individual"` for
#' which cell originates from which individual)
#' \code{nat_mat} (a `n` by `p` `matrix` denoting the natural parameter for each cell's gene expression), and
#' \code{obs_mat} (a `n` by `p` `dgCMatrix` denoting the observed count for each cell's gene expression),
#' @export
generate_null <- function(cell_per_person = 100,
                          num_genes = 1000,
                          num_individuals = 20){
  n_each <- cell_per_person # number of cells per person
  s <- num_individuals # number of people
  p <- num_genes # number of genes

  n <- n_each*s # total number of cells

  # form covariates
  # first form the table
  df <- cbind(1,
              0,
              rep(c(0,1), each = s/2),
              rep(c(0,1), times = s/2),
              scale(round(rnorm(s, mean = 30, sd = 5)), center = F, scale = T),
              1:s)
  colnames(df) <- c("Intercept", "Log_UMI", "CC", "Sex", "Age", "Individual")
  # expand to covariate matrix
  covariates <- do.call(rbind, lapply(1:nrow(df), function(i){
    matrix(rep(df[i,], each = n_each), nrow = n_each, ncol = ncol(df))
  }))
  colnames(covariates) <- colnames(df)[1:ncol(df)]
  z_mat <- cbind(rep(0, p), # intercept
                 rep(0.1, p), # library
                 c(rep(0.8,10),rep(0, p-10)), # cc
                 rnorm(p, mean = 0, sd = 0.2), # sex
                 rnorm(p, mean = 0, sd = 0.5)) # age
  colnames(z_mat) <- colnames(df)[1:(ncol(df)-1)]

  # form nuisance
  dispersion_vec <- sample(rep(c(10, 1, 0.1), each = ceiling(p/3))[1:p])
  dispersion_vec[1:10] <- 1

  # construct natural matrix
  # first label what type of gene it'll be
  gene_type <- rep(NA, p)
  gene_type[1:10] <- "true_DE"
  gene_type[seq(11, p, by=2)] <- "null_large_var"
  gene_type[seq(12, p, by=2)] <- "null_interleaved"
  gene_type <- as.factor(gene_type)

  nat_mat <- matrix(NA, nrow = n, ncol = p)
  rownames(nat_mat) <- paste0("cell_", 1:n)
  colnames(nat_mat) <- paste0("gene_", 1:p)
  ind_idx <- lapply(sort(unique(covariates[,"Individual"])), function(indiv){
    which(covariates[,"Individual"] == indiv)
  })
  for(i in 1:length(ind_idx)){
    rownames(nat_mat)[ind_idx[[i]]] <- paste0("c", i, "_", ind_idx[[i]])
  }
  case_indiv <- df[which(df[,"CC"] == 1), "Individual"]
  control_indiv <- df[which(df[,"CC"] == 0), "Individual"]

  # set true DE
  gene_idx <- which(gene_type == "true_DE")
  cc_diff <- 1
  indiv_sd <- 0.1
  group_sd <- 0.1
  min_mean <- 0.5
  max_mean <- 3
  for(j in gene_idx){
    case_mean <- stats::runif(1, min = min_mean, max = max_mean)
    control_mean <- case_mean + sample(c(-1,1), size = 1)*cc_diff

    for(indiv in case_indiv){
      mean_val <- stats::rnorm(1, mean = case_mean, sd = group_sd)
      nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                  mean = mean_val,
                                                  sd = indiv_sd)
    }

    for(indiv in control_indiv){
      mean_val <- stats::rnorm(1, mean = control_mean, sd = group_sd)
      nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                  mean = mean_val,
                                                  sd = indiv_sd)
    }
  }

  # set null large variance
  gene_idx <- which(gene_type == "null_large_var")
  indiv_sd_low <- 0.1
  indiv_sd_high <- 0.75
  group_sd <- 0.1
  min_mean <- -0.5
  max_mean <- 1
  for(j in gene_idx){
    mean_val <- stats::runif(1, min = min_mean, max = max_mean)
    sign_val <- sample(c(-1,1), size = 1)

    for(indiv in case_indiv){
      mean_indiv <- stats::rnorm(1, mean = mean_val, sd = group_sd)
      if(sign_val == 1){
        nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                    mean = mean_indiv,
                                                    sd = indiv_sd_low)
      } else {
        nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                    mean = mean_indiv,
                                                    sd = indiv_sd_high)
      }
    }

    for(indiv in control_indiv){
      mean_indiv <- stats::rnorm(1, mean = mean_val, sd = group_sd)
      if(sign_val == 1){
        nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                    mean = mean_indiv,
                                                    sd = indiv_sd_high)
      } else {
        nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                    mean = mean_indiv,
                                                    sd = indiv_sd_low)
      }
    }
  }

  # set null interleaved
  gene_idx <- which(gene_type == "null_interleaved")
  indiv_sd_min <- 0.1
  indiv_sd_max <- 0.2
  group_sd_min <- 0.1
  group_sd_max <- 0.2
  min_mean <- 0.5
  max_mean <- 2
  for(j in gene_idx){
    mean_val <- stats::runif(1, min = min_mean, max = max_mean)
    sd_val <- stats::runif(1, min = group_sd_min, max = group_sd_max)

    for(indiv in 1:length(ind_idx)){
      mean_indiv <- stats::rnorm(1, mean = mean_val, sd = sd_val)
      sd_indiv <- stats::runif(1, min = indiv_sd_min, max = indiv_sd_max)

      nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                  mean = mean_indiv,
                                                  sd = sd_indiv)
    }
  }

  nat_mat <- pmin(nat_mat, log(100))
  nat_mat <- nat_mat - 1

  gamma_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    gamma_mat[,j] <- stats::rgamma(
      n = n,
      shape = exp(nat_mat[,j])*dispersion_vec[j],
      rate = dispersion_vec[j])
  }
  gamma_mat <- pmin(gamma_mat, 50)

  lib_mat <- tcrossprod(covariates[,c("Intercept", "Log_UMI", "Sex", "Age")],
                        z_mat[,c("Intercept", "Log_UMI", "Sex", "Age")])
  lib_mat <- exp(lib_mat)
  obs_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    obs_mat[,j] <- stats::rpois(n = n,
                                lambda = lib_mat[,j]*gamma_mat[,j])
  }
  covariates[,"Log_UMI"] <- log1p(Matrix::rowSums(obs_mat))
  length(which(obs_mat == 0))/prod(dim(obs_mat))

  rownames(obs_mat) <- rownames(nat_mat)
  colnames(obs_mat) <- paste0("gene_", 1:ncol(obs_mat))
  rownames(covariates) <- rownames(obs_mat)

  df <- data.frame(df)
  df$Individual <- paste0("indiv_", df$Individual)

  tmp <- paste0("indiv_", covariates[,"Individual"])
  metadata_individual <- factor(tmp)

  list(covariates = covariates[,c("Intercept", "Log_UMI", "CC", "Sex", "Age")],
       df = df,
       metadata_individual = metadata_individual,
       nat_mat = nat_mat,
       obs_mat = Matrix::Matrix(obs_mat, sparse = TRUE))
}
