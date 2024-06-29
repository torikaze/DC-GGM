# Library ----
library(glasso)
library(dplyr)
library(tidyverse)

subdiff_largest_K = function(vec, K = NULL) {
  if (length(vec) < K) {
    stop("K must be smaller than the length of 'vec'")
  }
  df = data.frame(
    index = 1:length(vec),
    value = abs(vec),
    s = vec)
  df_sort = df[order(-df$value), ]
  if (K != nrow(df_sort)) {
    df_sort[(K+1):nrow(df_sort), 's'] = 0
  }
  df = df_sort[order(df_sort$index),]
  return(sign(df['s'])[,1])
}

make_PDmat_from_large_lambda = function(S, S_sign, decrease_ratio = 0.5) {
  if (!identical(dim(S), dim(S_sign))) {
    stop('mat and mat_sign must be the same shape')
  } else if (any(eigen(S, only.values = TRUE)$values < 0)) {
    stop('S is not positive (semi-) definite')
  }
  # initial lambda is determined by the diagonal elements of sample coariance matrix
  lambda = min(diag(S))
  num_c = 0 # to check the number of iteration
  while (TRUE) {
    # S_modified = S - lambda * S_sign + lambda * diag(nrow(S))
    S_modified = S - lambda * S_sign
    eig_min = min(eigen(S_modified, only.values = TRUE)$values)
    num_c = num_c + 1
    if (eig_min > 0) {
      res = list(lambda = lambda, num_c = num_c)
      return(res)
    }
    lambda = lambda*decrease_ratio
  }
}

dc_ggm = function(S, K = NULL, iter_dc = 100) {
  #' `K` includes diagonal component, so K should be larger than the number of rows of S.
  #'

  p = nrow(S)
  if (K > (p*(p+1)/2) || K < 0) {
    stop('"K" must be in the range of {0, ..., p*(p+1)/2}')
  } else if (!isSymmetric(S)) {
    stop('S must be symmetric')
  }

  make_PDmat_iter = c()
  Omega = solve(S + diag(x=1, nrow=nrow(S)))
  for (iter in seq_len(iter_dc)) {
    Omega_old = Omega
    # order: from top left to bottom left -> top of next column
    Omega_flat = Omega[lower.tri(Omega, diag=TRUE)]
    # dual step: gradient of largest K components
    s_sign = subdiff_largest_K(vec=Omega_flat, K=K)
    mat_s_sign = matrix(0, nrow=p, ncol=p)
    mat_s_sign[lower.tri(mat_s_sign, diag=TRUE)] = s_sign
    mat_s_sign = mat_s_sign + t(mat_s_sign)
    diag(mat_s_sign) = diag(mat_s_sign) / 2
    # primal step: minimize penalized likelihood function
    res_make_PDmat = make_PDmat_from_large_lambda(S=S, S_sign=mat_s_sign)
    lambda = res_make_PDmat$lambda
    make_PDmat_iter = c(make_PDmat_iter, res_make_PDmat$num_c)
    result_glasso_in_dca = glasso(S - mat_s_sign*lambda, rho=lambda)
    Omega = result_glasso_in_dca$wi
    Sigma = result_glasso_in_dca$w

    diff = sum(abs(Omega - Omega_old)^2)
    if (diff < 1e-4) {
      break
    }
  }
  result = list(
    Omega = Omega,
    Sigma = Sigma,
    dc_iter = iter,
    make_PDmat_iter = make_PDmat_iter
  )
  return(result)
}
