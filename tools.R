# Utility functions ----
graph_metrics = function(true_graph, estimated_graph) {
  value_true = true_graph[lower.tri(true_graph)]
  value_estimated = estimated_graph[lower.tri(estimated_graph)]
  TP = sum(value_true == 1 & value_estimated == 1)
  FP = sum(value_true == 0 & value_estimated == 1)
  TN = sum(value_true == 0 & value_estimated == 0)
  FN = sum(value_true == 1 & value_estimated == 0)
  accuracy = (TP + TN) / (TP + FP + TN + FN)
  precision = TP / (TP + FP)
  recall = TP / (TP + FN)
  F_measure = 2*precision*recall / (precision + recall)

  result = data.frame(
    accuracy = accuracy,
    recall = recall,
    precision = precision,
    F_measure = F_measure
    )
  return(result)
}

precisionmat_to_graphmetrics = function(precision_mat, true_graph) {
    diag(precision_mat) = 0
    num_edge = sum(precision_mat[lower.tri(precision_mat)] != 0)
    precision_mat[abs(precision_mat) > 0] = 1
    met = graph_metrics(true_graph=true_graph, estimated_graph=precision_mat)
    result = list(
        num_edge = num_edge,
        metrics = met
    )
    return(result)
}

num_edge_from_precisionmat = function(Omega) {
  diag(Omega) = 0
  num_edge = sum(Omega[lower.tri(Omega)] != 0)
  return(num_edge)
}

calc_loglikelihood = function(Omega, S, original=FALSE) {
  loglik = NULL
  if (original) {
    p = nrow(S)
    loglik = -p*log(2*pi) + log(det(Omega)) - sum(diag(Omega %*% S))
  } else {
    loglik = log(det(Omega)) - sum(diag(Omega %*% S))
  }
  return(loglik)
}

save_table_safely = function(data, path) {
  if (!file.exists(path)) {
    write_csv(data, path)
  } else {
    cat("the file already exists")
  }
}

save_pdf_safely = function(fig, path, width = 8, height = 10) {
  if (!file.exists(path)) {
    ggsave(file=path, plot=fig, device = cairo_pdf, width = width, height = height)
  } else {
    cat("the file already exists")
  }
}
