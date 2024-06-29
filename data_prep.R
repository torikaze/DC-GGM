#' Data preparation
#'
#' Description:
#' Load two types of synthetic data
#'

par(family= "HiraKakuProN-W3")


# Library and data -----
library("dplyr")
library("tidyverse")
library("mvtnorm")
library("ShrinkCovMat")


# Synthetic data ----
check_Symmetric_PD = function(mat) {
  if (!isSymmetric(mat)) {
    stop('not symmetric')
  }
  eigen_min = min(eigen(mat)$values)
  if (eigen_min < 0) {
    cat('NG \nthe minimum eigen value is ', eigen_min)
  } else {
    cat('ok \nminimum eigen value is ', eigen_min, '\n')
  }
}
setting_n_p = data.frame(
  p = c(rep(50, 3), rep(100, 3), rep(200, 3), rep(400, 3)),
  n = c(c(25, 50, 100), c(50, 100, 200), c(100, 200, 400), c(200, 400, 800))
)

# Type 1: Yuan (2007)'s method
generate_true_omega_yuan = function(dim, num_nonzero = NULL, seed = 0) {
  trueOmega = diag(dim)
  for (i in 2:dim) {
    trueOmega[i, i-1] = 0.5
    trueOmega[i-1, i] = 0.5
    if (i > 2) {
      trueOmega[i, i-2] = 0.25
      trueOmega[i-2, i] = 0.25
    }
  }
  set.seed(seed)
  if (!is.null(num_nonzero)) {
    index_nonzero = which(trueOmega[lower.tri(trueOmega)] != 0)
    index_nonzero_new = sample(index_nonzero, size=num_nonzero, replace=FALSE)
    trueOmega[lower.tri(trueOmega)][-index_nonzero_new] = 0
    trueOmega = t(trueOmega)
    trueOmega[lower.tri(trueOmega)][-index_nonzero_new] = 0
  }
  return (trueOmega)
}

# Type 2: Mazumder (2012)'s method
n_eigen_binary_search = function(left, right, mat, n = 1) {
  if (left >= right) {
    stop('`left` must be smaller than `right`')
  }
  mid = (left + right) / 2
  for (i in 1:1000) {
    min_eigen = min(eigen(mat + mid*diag(x=1, nrow=nrow(mat)))$values)
    if (round(min_eigen, 7) == n) {
      break
    }
    if (min_eigen < n) {
      left = mid
      mid = (left + right) / 2
    } else {
      right = mid
      mid = (left + right) / 2
    }
  }
  res = list(
    opt = mid,
    eigen_value = min_eigen
  )
  return(res)
}
generate_true_omega_mazumder = function(dim, num_nonzero = NULL, seed = 0) {
  set.seed(0)
  B = matrix(rnorm(dim*dim), nrow=dim)
  B_symmetry = 0.5 * (B + t(B))

  if (is.null(num_nonzero)) {
    num_zero = ceiling(dim/2)
  }
  num_lower_tri_element = dim*(dim - 1) / 2
  index_lower_tri = seq(1, num_lower_tri_element)
  set.seed(seed)
  index_nonzero = sample(index_lower_tri, size=num_nonzero, replace=FALSE)
  index_zero = index_lower_tri[-which(index_lower_tri %in% index_nonzero)]

  B_symmetry[lower.tri(B_symmetry)][index_zero] = 0
  B_symmetry = t(B_symmetry)
  B_symmetry[lower.tri(B_symmetry)][index_zero] = 0
  isSymmetric(B_symmetry)
  eta = n_eigen_binary_search(left=0.1, right=50, mat=B_symmetry)
  trueOmega = B_symmetry + eta$opt*diag(x=1, nrow=nrow(B_symmetry))
  return(trueOmega)
}

# data generating function
generate_synthetic_data = function(p, n, method = 'chain', shrinkage = TRUE, num_nonzero = 30, seed = 0) {
  #' Generate synthetic data based on specified method and parameters
  #'
  #' This function generates synthetic data using either the 'chain' or 'random' method.
  #' It can also apply shrinkage to the covariance matrix of the generated data.
  #'
  #' @param p The dimension of the data.
  #' @param n The number of observations.
  #' @param method The method to generate the synthetic data. Can be either 'chain' or 'random'.
  #' @param shrinkage A logical value indicating whether to apply shrinkage to the covariance matrix.
  #' @param num_nonzero The number of non-zero elements in the true covariance matrix.
  #' @param seed The seed for random number generation.
  #'
  #' @return A list containing the generated data, the covariance matrix, the number of true edges, and the true graph.
  #'

  if (method != 'chain' && method != 'random') {
    stop('method must be "chain" or "random"')
  }
  if (method == 'chain') {
    trueOmega = generate_true_omega_yuan(dim=p, num_nonzero=num_nonzero, seed=seed)
  } else {
    trueOmega = generate_true_omega_mazumder(dim=p, num_nonzero=num_nonzero, seed=seed)
  }
  # check_Symmetric_PD(mat=trueOmega)
  trueS = solve(trueOmega)
  set.seed(0)
  data = rmvnorm(n, mean=rep(0,p), sigma=trueS)
  if (isTRUE(shrinkage)) {
    res_shrink = shrinkcovmat.unequal(t(data), centered=FALSE)
    if (res_shrink$lambdahat == 1) {
      data_cov = 0.9*diag(diag(cov(data))) + (1-0.9)*cov(data)
    } else {
      data_cov = res_shrink$Sigmahat
    }
  } else {
    data_cov = cov(data)
  }
  # check_Symmetric_PD(mat=data_cov)
  num_edge_true = sum(trueOmega[lower.tri(trueOmega)] != 0)
  true_graph = ifelse(trueOmega!=0 & row(trueOmega)!=col(trueOmega),1,0)
  result = list(
    data = data,
    cov = data_cov,
    num_edge_true = num_edge_true,
    true_graph = true_graph
  )
  return(result)
}
