# Library ----
library(ShrinkCovMat)
library(rsample)
library(GGMncv)
library(glasso)
library(parallel)
library(tictoc)

source('dc_ggm.R')

NUM_CORES = detectCores()
# DC
cv_dc = function(data, S, n_folds = 5, K_max = NULL) {
  p = ncol(data)
  if (is.null(K_max)) {
    K_max = p*(p+1)/2
  }

  K_range = seq(p+1, K_max, length.out=100) |> floor() |> unique()
  # Split into training and test data and compute covariance matrix for each fold
  set.seed(0)
  data_cv = data |> vfold_cv(v=n_folds)
  result_cv = data.frame(K=integer(), loglik_mean=double())
  train_list = lapply(data_cv$splits, analysis)
  val_list = lapply(data_cv$splits, assessment)
  cov_train_list = lapply(train_list, function(train) {
    res_shrink = shrinkcovmat.unequal(t(train), centered=FALSE)
    if (res_shrink$lambdahat == 1) {
      0.9*diag(diag(cov(train))) + (1-0.9)*cov(train)
    } else {
      res_shrink$Sigmahat
    }
  })
  cov_val_list = lapply(val_list, cov)

  res_list = mclapply(K_range, function(K) {
    loglik = numeric(n_folds)
    num_edge = numeric(n_folds)
    for (i in seq_len(n_folds)) {
      temp = dc_ggm(S=cov_train_list[[i]], K=K)
      loglik[i] = log(det(temp$Omega)) - sum(diag(temp$Omega %*% cov_val_list[[i]]))
      diag(temp$Omega) = 0
      num_edge[i] = sum(temp$Omega[lower.tri(temp$Omega)] != 0)
    }
    return(list(K = K, loglik_mean = mean(loglik), num_edge = mean(num_edge)))
  }, mc.cores = NUM_CORES)
  result_cv = do.call(rbind, lapply(res_list, data.frame))

  bestK = result_cv |>
    filter(loglik_mean == max(loglik_mean)) |> dplyr::select(K)
  model_opt = dc_ggm(S=S, K=bestK[1,1])
  result = list(
    model_opt = model_opt,
    result_cv = result_cv,
    bestK = bestK
  )
  return(result)
}


get_estimated_edge_num = function(lambda, S, n, method) {
    if (method == "glasso") {
        precision_mat = glasso(s=S, rho=lambda)$wi
    } else if (method == "scad") {
        precision_mat = ggmncv(R=S, n=n, penalty="scad", lambda=lambda, progress=FALSE)$Theta
    } else if (method == "adapt") {
        precision_mat = ggmncv(R=S, n=n, penalty="adapt", lambda=lambda, progress=FALSE)$Theta
    } else {
        stop("method ", method, "is invalid.")
    }
    diag(precision_mat) = 0
    num_edge = sum(precision_mat[lower.tri(precision_mat)] != 0)
    return(num_edge)
}

decide_lambda_range = function(S, n, method = 'glasso') {
    #' Decide the penalty parameter so that the model can select
    #' from a minimum to maximum number of edges.

    lambda = median(abs(S[lower.tri(S)]))
    # estimate the minimum number of edges.
    while (TRUE) {
        num_edge = get_estimated_edge_num(lambda=lambda, S=S, n=n, method=method)
        if (num_edge > 1) {
            lambda = lambda*1.5
        } else {
            lambda_max = lambda
            break
        }
    }

    lambda = min(abs(S[lower.tri(S)]))
    target_maximum = nrow(S) * (nrow(S) - 1) / 2
    while (TRUE) {
        num_edge = get_estimated_edge_num(lambda=lambda, S=S, n=n, method=method)
        if (num_edge < target_maximum) {
            lambda = lambda*0.1
        } else {
            lambda_min = lambda
            break
        }
    }
    return(seq(lambda_min, lambda_max, length.out = 100))
}

# glasso and Non-convex penalties
cv_ncv = function(data, S, n_folds = 5, lambda_range = NULL) {
  if (is.null(lambda_range)) {
    lam_min = min(abs(S[lower.tri(S)]))
    lam_max = max(abs(S[lower.tri(S)]))
    lambda_range = seq(lam_min, lam_max, length.out = 100)
  }

  lambda_range = decide_lambda_range(S=S, n=nrow(data), method="glasso")
  res_list = mclapply(lambda_range, function(lam) {
    res = obj_cvglasso(x = lam, n_folds = n_folds, data = data)
    return(list(lambda = lam, loglik_mean = res$loglik, num_edge = res$num_edge))
  }, mc.cores = NUM_CORES)
  result_cv_glasso = do.call(rbind, lapply(res_list, data.frame))
  bestlambda_glasso = result_cv_glasso |>
    filter(loglik_mean == max(loglik_mean)) |> dplyr::select(lambda)

  lambda_range = decide_lambda_range(S=S, n=nrow(data), method="scad")
  res_list = mclapply(lambda_range, function(lam) {
    res = obj_cvscad(x = lam, n_folds = n_folds, data = data)
    return(list(lambda = lam, loglik_mean = res$loglik, num_edge = res$num_edge))
  }, mc.cores = NUM_CORES)
  result_cv_scad = do.call(rbind, lapply(res_list, data.frame))
  bestlambda_scad = result_cv_scad |>
    filter(loglik_mean == max(loglik_mean)) |> dplyr::select(lambda)

  lambda_range = decide_lambda_range(S=S, n=nrow(data), method="adapt")
  res_list = mclapply(lambda_range, function(lam) {
    res = obj_cvadapt(x = lam, n_folds = n_folds, data = data)
    return(list(lambda = lam, loglik_mean = res$loglik, num_edge = res$num_edge))
  }, mc.cores = NUM_CORES)
  result_cv_adapt = do.call(rbind, lapply(res_list, data.frame))
  bestlambda_adapt = result_cv_adapt |>
    filter(loglik_mean == max(loglik_mean)) |> dplyr::select(lambda)

  model_glasso = glasso(S, rho=bestlambda_glasso[1,1])
  model_scad = ggmncv(R=S, n=nrow(data), penalty='scad', lambda = bestlambda_scad[1,1], progress=FALSE)
  model_adapt = ggmncv(R=S, n=nrow(data), penalty='adapt', lambda = bestlambda_adapt[1,1], progress=FALSE)

  result = list(
    model_glasso = model_glasso,
    model_scad = model_scad,
    model_adapt = model_adapt,
    result_cv_glasso = result_cv_glasso,
    result_cv_scad = result_cv_scad,
    result_cv_adapt = result_cv_adapt,
    lambda = data.frame(
      glasso = bestlambda_glasso[1,1],
      scad = bestlambda_scad[1,1],
      adapt = bestlambda_adapt[1,1]
    )
  )
  return(result)
}

obj_cvglasso = function(x, n_folds, data) {
  # Split into training and test data and compute covariance matrix for each fold
  set.seed(0)
  data_cv = data |> vfold_cv(v=n_folds)
  train_list = lapply(data_cv$splits, analysis)
  val_list = lapply(data_cv$splits, assessment)
  cov_train_list = lapply(train_list, function(train) {
    res_shrink = shrinkcovmat.unequal(t(train), centered=FALSE)
    if (res_shrink$lambdahat == 1) {
      0.9*diag(diag(cov(train))) + (1-0.9)*cov(train)
    } else {
      res_shrink$Sigmahat
    }
  })
  cov_val_list = lapply(val_list, cov)

  loglik = numeric(n_folds)
  num_edge = numeric(n_folds)
  for (i in seq_len(n_folds)) {
    temp = glasso(s=cov_train_list[[i]], rho=x)
    loglik[i] = log(det(temp$wi)) - sum(diag(temp$wi %*% cov_val_list[[i]]))
    diag(temp$wi) = 0
    num_edge[i] = sum(temp$wi[lower.tri(temp$wi)] != 0)
  }
  return(list(loglik = mean(loglik), num_edge = mean(num_edge)))
}

obj_cvscad = function(x, n_folds, data) {
  # Split into training and test data and compute covariance matrix for each fold
  set.seed(0)
  data_cv = data |> vfold_cv(v=n_folds)
  train_list = lapply(data_cv$splits, analysis)
  val_list = lapply(data_cv$splits, assessment)
  cov_train_list = lapply(train_list, function(train) {
    res_shrink = shrinkcovmat.unequal(t(train), centered=FALSE)
    if (res_shrink$lambdahat == 1) {
      0.9*diag(diag(cov(train))) + (1-0.9)*cov(train)
    } else {
      res_shrink$Sigmahat
    }
  })
  cov_val_list = lapply(val_list, cov)

  loglik = numeric(n_folds)
  num_edge = numeric(n_folds)
  for (i in seq_len(n_folds)) {
    temp = ggmncv(R=cov_train_list[[i]], n=nrow(train_list[[i]]), penalty='scad', lambda=x, progress=FALSE)
    loglik[i] = log(det(temp$Theta)) - sum(diag(temp$Theta %*% cov_val_list[[i]]))
    diag(temp$Theta) = 0
    num_edge[i] = sum(temp$Theta[lower.tri(temp$Theta)] != 0)
  }
  return(list(loglik = mean(loglik), num_edge = mean(num_edge)))
}

obj_cvadapt = function(x, n_folds, data) {
  # Split into training and test data and compute covariance matrix for each fold
  set.seed(0)
  data_cv = data |> vfold_cv(v=n_folds)
  train_list = lapply(data_cv$splits, analysis)
  val_list = lapply(data_cv$splits, assessment)
  cov_train_list = lapply(train_list, function(train) {
    res_shrink = shrinkcovmat.unequal(t(train), centered=FALSE)
    if (res_shrink$lambdahat == 1) {
      0.9*diag(diag(cov(train))) + (1-0.9)*cov(train)
    } else {
      res_shrink$Sigmahat
    }
  })
  cov_val_list = lapply(val_list, cov)

  loglik = numeric(n_folds)
  num_edge = numeric(n_folds)
  for (i in seq_len(n_folds)) {
    temp = ggmncv(R=cov_train_list[[i]], n=nrow(train_list[[i]]), penalty='adapt', lambda=x, progress=FALSE)
    loglik[i] = log(det(temp$Theta)) - sum(diag(temp$Theta %*% cov_val_list[[i]]))
    diag(temp$Theta) = 0
    num_edge[i] = sum(temp$Theta[lower.tri(temp$Theta)] != 0)
  }
  return(list(loglik = mean(loglik), num_edge = mean(num_edge)))
}


# integrate the results of cross-validation
cv_comparison = function(data, S, true_graph) {
  # DC
  tic(quiet = TRUE)
  res_cv_dc = cv_dc(data=data, S=S)
  a = toc(quiet = TRUE)
  time_dc = as.numeric(a$toc - a$tic)

  omega_dc = res_cv_dc$model_opt$Omega
  diag(omega_dc) = 0
  omega_dc[abs(omega_dc) > 0] = 1
  metrics_dc = graph_metrics(true_graph=true_graph, estimated_graph=omega_dc)
  metrics_dc$num_edge = sum(omega_dc[lower.tri(omega_dc)] != 0)

  tic(quiet = TRUE)
  res_cv_ncv = cv_ncv(data=data, S=S)
  a = toc(quiet = TRUE)
  time_ncv = as.numeric(a$toc - a$tic)

  # glasso
  omega_glasso = res_cv_ncv$model_glasso$wi
  diag(omega_glasso) = 0
  omega_glasso[abs(omega_glasso) > 0] = 1
  metrics_glasso = graph_metrics(true_graph=true_graph, estimated_graph=omega_glasso)
  metrics_glasso$num_edge = sum(omega_glasso[lower.tri(omega_glasso)] != 0)
  # SCAD
  omega_scad = res_cv_ncv$model_scad$Theta
  diag(omega_scad) = 0
  omega_scad[abs(omega_scad) > 0] = 1
  metrics_scad = graph_metrics(true_graph=true_graph, estimated_graph=omega_scad)
  metrics_scad$num_edge = sum(omega_scad[lower.tri(omega_scad)] != 0)
  # Adaptive
  omega_adapt = res_cv_ncv$model_adapt$Theta
  diag(omega_adapt) = 0
  omega_adapt[abs(omega_adapt) > 0] = 1
  metrics_adapt = graph_metrics(true_graph=true_graph, estimated_graph=omega_adapt)
  metrics_adapt$num_edge = sum(omega_adapt[lower.tri(omega_adapt)] != 0)

  time_all = data.frame(time_dc = time_dc, time_ncv = time_ncv)
  result = list(
    dc = metrics_dc,
    glasso = metrics_glasso,
    scad = metrics_scad,
    adapt = metrics_adapt,
    time = time_all,
    result_cv_dc = res_cv_dc$result_cv,
    result_cv_glasso = res_cv_ncv$result_cv_glasso,
    result_cv_scad = res_cv_ncv$result_cv_scad,
    result_cv_adapt = res_cv_ncv$result_cv_adapt
  )
  return(result)
}
