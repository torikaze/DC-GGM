par(family= "HiraKakuProN-W3")
path_data = './data/'
path_export = './output/'

# Library ----
library(dplyr)
library(tidyverse)
library(cowplot)

source('cv_tools_wide_lamrange.R')
source('data_prep.R')
source('tools.R')

# Exp 6: log-likelihood during cross validation ----
METHOD = 'random'
num_nonzero = 30
# CV in various setting
edge_loglik_dc_all = data.frame()
edge_loglik_glasso_all = data.frame()
edge_loglik_scad_all = data.frame()
edge_loglik_adapt_all = data.frame()

for (i in seq_len(nrow(setting_n_p))) {
  n = setting_n_p$n[i]
  p = setting_n_p$p[i]
  cat("\n(n, p) = (",n,",", p,")\n")
  # use only one seed
  data_info = generate_synthetic_data(p=p, n=n, method=METHOD, num_nonzero=num_nonzero, seed=0)
  cv_by_seed = cv_comparison(data=data_info$data, S=data_info$cov, true_graph=data_info$true_graph)

  edge_loglik_dc = cv_by_seed$result_cv_dc |> mutate(n_p=str_c(n, '_', p), n=n, p=p, model = 'DC')
  edge_loglik_glasso = cv_by_seed$result_cv_glasso |> mutate(n_p=str_c(n, '_', p), n=n, p=p, model = 'glasso')
  edge_loglik_scad = cv_by_seed$result_cv_scad |> mutate(n_p=str_c(n, '_', p), n=n, p=p, model = 'SCAD')
  edge_loglik_adapt = cv_by_seed$result_cv_adapt |> mutate(n_p=str_c(n, '_', p), n=n, p=p, model = 'adapt')

  edge_loglik_dc_all = rbind(edge_loglik_dc_all, edge_loglik_dc)
  edge_loglik_glasso_all = rbind(edge_loglik_glasso_all, edge_loglik_glasso)
  edge_loglik_scad_all = rbind(edge_loglik_scad_all, edge_loglik_scad)
  edge_loglik_adapt_all = rbind(edge_loglik_adapt_all, edge_loglik_adapt)
  cat("\n", i / nrow(setting_n_p) * 100,"%\n")
}

n_p_pairs = c('25_50', '50_50', '100_50', '50_100', '100_100', '200_100', '100_200', '200_200', '400_200', '200_400', '400_400', '800_400')
result_exp6 = rbind(edge_loglik_dc_all[,-1], edge_loglik_glasso_all[,-1], edge_loglik_scad_all[,-1], edge_loglik_adapt_all[,-1])

today_date = format(Sys.Date(), "%y%m%d")
save_table_safely(data = result_exp6, path = str_c(path_export, 'table/result_exp6_', METHOD, '_', today_date, '.csv'))


# Visualization ----
levels_n = c("25", "50", "100", "200", "400", "800")
levels_p = str_c('p=', c("50", "100", "200", "400"))
levels_model = c('DC', 'glasso', 'SCAD', 'adapt')
levels_pn = unique(result_exp6$n_p)

# Start of Selection
edge_loglik_fig = result_exp6 |>
  ggplot(aes(x=num_edge, y=loglik_mean)) +
  facet_wrap(~factor(n_p, levels=levels_pn), scale='free', nrow=4) +
  geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
  geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
  geom_vline(xintercept=30, linetype="dashed", color="black") +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
  theme_minimal(base_family="HiraKakuPro-W3") +
  theme(
    text = element_text(size = 15),
    legend.key.size = unit(1.0, "cm"),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(x="# of edges", y='log-likelihood', color='', fill='')

save_pdf_safely(fig=edge_loglik_fig, path = str_c(path_export, 'figure/exp6_log_result_', METHOD, '_', today_date, '.pdf'), width = 8, height = 9)

