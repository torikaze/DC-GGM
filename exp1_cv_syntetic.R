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

# Exp1: cross validation ----
# decide penalty parameters by K-folds cross-validation.
METHOD = 'random'
n_simulation = 30
num_nonzero = 30
seed_list = seq(1, n_simulation)

result_cv_glasso_all = data.frame()
result_cv_dc_all = data.frame()
result_cv_scad_all = data.frame()
result_cv_adapt_all = data.frame()
time_cv_all = data.frame()
for (i in seq_len(nrow(setting_n_p))) {
  result_cv_glasso = data.frame()
  result_cv_dc = data.frame()
  result_cv_scad = data.frame()
  result_cv_adapt = data.frame()
  time_cv = data.frame()

  n = setting_n_p$n[i]
  p = setting_n_p$p[i]
  cat("\n(n, p) = (",n,",", p,")\n")
  cat('cv process: ')
  for (seed in seed_list) {
    data_info = generate_synthetic_data(p=p, n=n, method=METHOD, num_nonzero=num_nonzero, seed=seed)
    cv_by_seed = cv_comparison(data=data_info$data, S=data_info$cov, true_graph=data_info$true_graph)
    result_cv_glasso = rbind(result_cv_glasso, cv_by_seed$glasso)
    result_cv_dc = rbind(result_cv_dc, cv_by_seed$dc)
    result_cv_scad = rbind(result_cv_scad, cv_by_seed$scad)
    result_cv_adapt = rbind(result_cv_adapt, cv_by_seed$adapt)
    time_cv = rbind(time_cv, cv_by_seed$time)
    cat(which(seed == seed_list), '')
  }
  result_cv_glasso = result_cv_glasso |> mutate(n_p=str_c(n, '_', p), n=n, p=p)
  result_cv_dc = result_cv_dc |> mutate(n_p=str_c(n, '_', p), n=n, p=p)
  result_cv_scad = result_cv_scad |> mutate(n_p=str_c(n, '_', p), n=n, p=p)
  result_cv_adapt = result_cv_adapt |> mutate(n_p=str_c(n, '_', p), n=n, p=p)
  time_cv = time_cv |> mutate(n_p=str_c(n, '_', p), n=n, p=p)

  result_cv_glasso_all = rbind(result_cv_glasso_all, result_cv_glasso)
  result_cv_dc_all = rbind(result_cv_dc_all, result_cv_dc)
  result_cv_scad_all = rbind(result_cv_scad_all, result_cv_scad)
  result_cv_adapt_all = rbind(result_cv_adapt_all, result_cv_adapt)
  time_cv_all = rbind(time_cv_all, time_cv)
  cat("\n", i / nrow(setting_n_p) * 100,"%\n")
  cat('time (DC, ncv): ', round(apply(time_cv[,c('time_dc', 'time_ncv')], 2, mean), 1), '\n')
}

# calculate mean and standard error by the groups with the same settings (n and p)
summary_cv_glasso = result_cv_glasso_all |>
  mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
  group_by(n_p) |>
  summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), n_mean=mean(num_edge), n_std=sd(num_edge)/sqrt(n_simulation)) |> mutate(model = 'glasso')
summary_cv_dc = result_cv_dc_all |>
  mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
  group_by(n_p) |>
  summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), n_mean=mean(num_edge), n_std=sd(num_edge)/sqrt(n_simulation)) |> mutate(model = "DC")
summary_cv_scad = result_cv_scad_all |>
  mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
  group_by(n_p) |>
  summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), n_mean=mean(num_edge), n_std=sd(num_edge)/sqrt(n_simulation)) |> mutate(model = 'SCAD')
summary_cv_adapt = result_cv_adapt_all |>
  mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
  group_by(n_p) |>
  summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), n_mean=mean(num_edge), n_std=sd(num_edge)/sqrt(n_simulation)) |> mutate(model = 'adapt')

n_p_pairs = c('25_50', '50_50', '100_50', '50_100', '100_100', '200_100', '100_200', '200_200', '400_200', '200_400', '400_400', '800_400')
# save results
result_cv_glasso_all = result_cv_glasso_all |> mutate(model="glasso")
result_cv_dc_all = result_cv_dc_all |> mutate(model="DC")
result_cv_scad_all = result_cv_scad_all |> mutate(model="SCAD")
result_cv_adapt_all = result_cv_adapt_all |> mutate(model="adapt")

result_exp1 = rbind(result_cv_glasso_all, result_cv_dc_all, result_cv_scad_all, result_cv_adapt_all)

today_date = format(Sys.Date(), "%y%m%d")
# save_table_safely(data = result_exp1, path = str_c(path_export, 'table/result_exp1_', METHOD, '_', today_date, '.csv'))
save_table_safely(data = result_exp1, path = str_c('AAA', METHOD, '_', today_date, '.csv'))

# METHOD = 'chain'
# result_exp1 = read_csv(str_c(path_export, 'table/result_exp1_', METHOD, '_240614.csv'))
# result_exp1 = result_exp1 |>
#   mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
#   group_by(n_p, model) |>
#   summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), n_mean=mean(num_edge), n_std=sd(num_edge)/sqrt(n_simulation))

result_exp1 = rbind(summary_cv_glasso, summary_cv_dc, summary_cv_scad, summary_cv_adapt)
exp1_np = result_exp1$n_p |> str_split("_", simplify=TRUE)
result_exp1['n'] = exp1_np[,1]
result_exp1['p'] = str_c('p=', exp1_np[,2])

levels_n = c("25", "50", "100", "200", "400", "800")
levels_p = str_c('p=', c("50", "100", "200", "400"))
levels_model = c('DC', 'glasso', 'SCAD', 'adapt')

a = result_exp1 |>
  ggplot(aes(x=factor(n, levels=levels_n), y=mean)) +
  facet_wrap(~factor(p, levels=levels_p), scale='free_x') +
  geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
  geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
  geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity') +
  theme_minimal(base_family="HiraKakuPro-W3") +
  theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(x='N', y='F-measure', col='', fill='')

num_edge_gg = result_exp1 |>
  ggplot(aes(x=factor(n, levels=levels_n), y=n_mean)) +
  facet_wrap(~factor(p, levels=levels_p), scale='free') +
  geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
  geom_errorbar(aes(ymin = n_mean - 2*n_std, ymax = n_mean+2*n_std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
  geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity') +
  theme_minimal(base_family="HiraKakuPro-W3") +
  theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(x='N', y='# of edges', col='', fill='')

a = a + theme(plot.margin = margin(5.5, 20, 5.5, 20))
num_edge_gg = num_edge_gg + theme(plot.margin = margin(5.5, 20, 5.5, 20))
combined_plot = cowplot::plot_grid(a, num_edge_gg, ncol = 1, align = 'v')

save_pdf_safely(fig=combined_plot, path = str_c(path_export, 'figure/combined_exp1_result_', METHOD, '_', today_date, '.pdf'), width = 8, height = 10)

