par(family= "HiraKakuProN-W3")
path_data = './data/'
path_export = './output/'
NUM_CORES = detectCores()

# Library ----
library(dplyr)
library(tidyverse)
library(cowplot)

source('dc_ggm.R')
source('data_prep.R')
source('tools.R')

# Exp2: estimate fixed number of nonzeo-elements ----
exp2_estimate_same_number_of_edges = function(S, num_edge_target, true_graph, n) {
  # check the estimated edge number `num_edge_dc` of DC at first.
  # if `num_edge_dc` is not equal to `num_edge_target`, terminate the function

  p = nrow(S)

  # DC
  K_min = max(p + num_edge_target - 30, p + 1)
  K_range = seq(K_min, p + num_edge_target + 100)
  res_list = mclapply(K_range, function(K) {
    omega_estimated = dc_ggm(S=S, K=K)$Omega
    diag(omega_estimated) = 0
    num_edge_dc = sum(omega_estimated[lower.tri(omega_estimated)] != 0)
    return(list(K = K, num_edge = num_edge_dc))
  }, mc.cores = NUM_CORES)
  result_dc = do.call(rbind, lapply(res_list, data.frame)) |> filter(num_edge == num_edge_target)
  if (dim(result_dc)[1] == 0) {
    return(FALSE)
  }
  omega_estimated = dc_ggm(S=S, K=result_dc$K[1])$Omega
  diag(omega_estimated) = 0
  omega_estimated[omega_estimated != 0] = 1
  num_edge_dc = sum(omega_estimated[lower.tri(omega_estimated)] != 0)
  if (num_edge_dc != num_edge_target) {
    return(FALSE)
  }
  metrics_dc = graph_metrics(true_graph=true_graph, estimated_graph = omega_estimated)

  result_dc = data.frame(
    num_edge = num_edge_dc,
    metrics = metrics_dc)

  # glasso
  lam_min = min(abs(S[lower.tri(S)]))*0.1
  lam_max = max(abs(S[lower.tri(S)]))*10
  num_edge_glasso = c()
  metrics_glasso = c()
  left = lam_min
  right = lam_max
  mid = (left + right) / 2
  for (i in 1:100) {
    omega_estimated = glasso(s=S, rho=mid)$wi
    diag(omega_estimated) = 0
    num_edge_glasso = sum(omega_estimated[lower.tri(omega_estimated)] != 0)
    if (num_edge_glasso == num_edge_target) {
      break
    } else if (num_edge_glasso > num_edge_target) {
      left = mid
      mid = (left + right) / 2
    } else if (num_edge_glasso < num_edge_target) {
      right = mid
      mid = (left + right) / 2
    }
  }
  omega_estimated[omega_estimated != 0] = 1
  metrics_glasso = graph_metrics(true_graph=true_graph, estimated_graph = omega_estimated)

  result_glasso = data.frame(
    num_edge = num_edge_glasso,
    metrics = metrics_glasso)

  # SCAD
  num_edge_scad = c()
  metrics_scad = c()
  left = lam_min
  right = lam_max
  mid = (left + right) / 2
  for (i in 1:100) {
    omega_estimated = ggmncv(R=S, n=n, penalty='scad',lambda = mid, progress=FALSE)$adj
    diag(omega_estimated) = 0
    num_edge_scad = sum(omega_estimated[lower.tri(omega_estimated)] != 0)
    if (num_edge_scad == num_edge_target) {
      break
    } else if (num_edge_scad > num_edge_target) {
      left = mid
      mid = (left + right) / 2
    } else if (num_edge_scad < num_edge_target) {
      right = mid
      mid = (left + right) / 2
    }
  }
  metrics_scad = graph_metrics(true_graph=true_graph, estimated_graph=omega_estimated)

  result_scad = data.frame(
    num_edge = num_edge_scad,
    metrics = metrics_scad)

  # adapt
  num_edge_adapt = c()
  metrics_adapt = c()
  left = lam_min
  right = lam_max
  mid = (left + right) / 2
  for (i in 1:100) {
    omega_estimated = ggmncv(R=S, n=n, penalty='adapt',lambda = mid, progress=FALSE)$adj
    diag(omega_estimated) = 0
    num_edge_adapt = sum(omega_estimated[lower.tri(omega_estimated)] != 0)
    if (num_edge_adapt == num_edge_target) {
      break
    } else if (num_edge_adapt > num_edge_target) {
      left = mid
      mid = (left + right) / 2
    } else if (num_edge_adapt < num_edge_target) {
      right = mid
      mid = (left + right) / 2
    }
  }
  metrics_adapt = graph_metrics(true_graph=true_graph, estimated_graph=omega_estimated)

  result_adapt = data.frame(
    num_edge = num_edge_adapt,
    metrics = metrics_adapt)

  if (num_edge_dc != num_edge_target | num_edge_glasso != num_edge_target | num_edge_scad != num_edge_target) {
    browser()
  }

  result_all = list(
    result_dc = result_dc,
    result_glasso = result_glasso,
    result_scad = result_scad,
    result_adapt = result_adapt
  )
  return(result_all)
}

METHOD = 'random'
n_simulation = 30
num_nonzero = 30
num_edge_target_list = c(num_nonzero - 10, num_nonzero, num_nonzero + 10)

compare_dc_gg_all = c()
compare_glasso_gg_all = c()
compare_scad_gg_all = c()
compare_adapt_gg_all = c()
for (i in seq_len(nrow(setting_n_p))) {
  n = setting_n_p$n[i]
  p = setting_n_p$p[i]
  res_dc_all = c()
  res_glasso_all = c()
  res_scad_all = c()
  res_adapt_all = c()

  cat("\n(n, p) = (",n,",", p,")")
  for (num_edge_target in num_edge_target_list) {
    num_valid_iteration = 0
    for (seed in seq_len(100)) {
      data_info = generate_synthetic_data(p, n, method=METHOD, shrinkage=TRUE, num_nonzero=num_nonzero, seed=seed)
      # sometimes DC failed to estimate `num_edge_target` edges.
      # So experiment is repeated untill there are 30 successful results of DC with `num_edge_target`.
      res_compare = exp2_estimate_same_number_of_edges(S=data_info$cov, num_edge_target=num_edge_target, true_graph=data_info$true_graph, n=n)
      if (isFALSE(res_compare)) {
        next
      } else {
        res_dc = res_compare$result_dc
        res_dc_all = rbind(res_dc_all, res_dc)
        res_glasso = res_compare$result_glasso
        res_glasso_all = rbind(res_glasso_all, res_glasso)
        res_scad = res_compare$result_scad
        res_scad_all = rbind(res_scad_all, res_scad)
        res_adapt = res_compare$result_adapt
        res_adapt_all = rbind(res_adapt_all, res_adapt)
        num_valid_iteration = num_valid_iteration + 1

      }
      if (num_valid_iteration == n_simulation) {
        break
      }
    }
    if (num_valid_iteration != n_simulation) {
      stop('dc fialed to estimate `num_edge_target` edges')
    }
    cat(which(num_edge_target == num_edge_target_list), '/', length(num_edge_target_list), '')
  }
  colnames(res_dc_all) = c('num_edge', 'accuracy', 'recall', 'precision', 'F_measure')
  colnames(res_glasso_all) = c('num_edge', 'accuracy', 'recall', 'precision', 'F_measure')
  colnames(res_scad_all) = c('num_edge', 'accuracy', 'recall', 'precision', 'F_measure')
  colnames(res_adapt_all) = c('num_edge', 'accuracy', 'recall', 'precision', 'F_measure')

  compare_dc_gg = res_dc_all |>
    mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
    group_by(num_edge) |>
    summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), num=n()) |>
    mutate(model = "DC", n = n, p = p, n_p = str_c(n, '_', p))
  compare_glasso_gg = res_glasso_all |>
    mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
    group_by(num_edge) |>
    summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), num=n()) |>
    mutate(model = 'glasso', n = n, p = p, n_p = str_c(n, '_', p))
  compare_scad_gg = res_scad_all |>
    mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
    group_by(num_edge) |>
    summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), num=n()) |>
    mutate(model = 'SCAD', n = n, p = p, n_p = str_c(n, '_', p))
  compare_adapt_gg = res_adapt_all |>
    mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
    group_by(num_edge) |>
    summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), num=n()) |>
    mutate(model = 'adapt', n = n, p = p, n_p = str_c(n, '_', p))

  compare_dc_gg_all = rbind(compare_dc_gg_all, compare_dc_gg)
  compare_glasso_gg_all = rbind(compare_glasso_gg_all, compare_glasso_gg)
  compare_scad_gg_all = rbind(compare_scad_gg_all, compare_scad_gg)
  compare_adapt_gg_all = rbind(compare_adapt_gg_all, compare_adapt_gg)
  cat('\n', round(i / nrow(setting_n_p) * 100, 0),"%\n")
}

n_p_pairs2 = c('25_50', '50_50', '100_50', '50_100', '100_100', '200_100', '100_200', '200_200', '400_200', '200_400', '400_400', '800_400')
result_exp2 = rbind(compare_glasso_gg_all, compare_dc_gg_all, compare_scad_gg_all, compare_adapt_gg_all)

result_exp2 %>%
  group_by(n_p, num) %>%
  filter(mean == max(mean)) |>
  ungroup() |>
  summarize(DC_top = sum(model == 'DC'), glasso_top = sum(model == 'glasso'), SCAD_top = sum(model == 'SCAD'), adapt_top = sum(model == 'adapt'))

today_date = format(Sys.Date(), "%y%m%d")
save_table_safely(data = arrange(result_exp2, p, n, num_edge), path = str_c(path_export, 'table/result_exp2_', METHOD, '_', today_date, '.csv'))


# Visualization ----
levels_n = c("25", "50", "100", "200", "400", "800")
levels_p = str_c('p=', c("50", "100", "200", "400"))
levels_model = c('DC', 'glasso', 'SCAD', 'adapt')
result_exp2 = result_exp2 |> mutate(n_p = str_c('p=', p, ', N=', n))
levels_pn = unique(result_exp2$n_p)
result_exp2
exp2_fig = result_exp2  |>
  ggplot(aes(x=as.factor(num_edge), y=mean)) +
  facet_wrap(~factor(n_p, levels=levels_pn), scale='free_x', ncol=3) +
  geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
  geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
  geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity') +
  theme_minimal(base_family="HiraKakuPro-W3") +
  theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(x='# of edges', y='F-measure', col='', fill='')
exp2_fig
save_pdf_safely(fig=exp2_fig, path = str_c(path_export, 'figure/exp2_result_', METHOD, '_', today_date, '.pdf'), width = 8, height = 10)
