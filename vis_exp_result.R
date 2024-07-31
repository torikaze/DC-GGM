par(family= "HiraKakuProN-W3")
path_export = './output/'

# Library ----
library(dplyr)
library(gridExtra)
library(tidyverse)
library(grid)

source('tools.R')

# Visualization for each experiment ----
# Exp1: cross validation ----
METHOD = 'random'
result_exp1 = read_csv(str_c(path_export, 'table/result_exp1_', METHOD, '_240615.csv'))
result_exp1 = result_exp1 |>
  mutate_at("F_measure", ~replace(., is.na(.), 0)) |>
  group_by(n_p, model) |>
  summarise(mean=mean(F_measure), std=sd(F_measure)/sqrt(n_simulation), n_mean=mean(num_edge), n_std=sd(num_edge)/sqrt(n_simulation))
result_exp1
exp1_np = result_exp1$n_p |> str_split("_", simplify=TRUE)
result_exp1['n'] = exp1_np[,1]
result_exp1['p'] = str_c('p=', exp1_np[,2])
result_exp1
levels_n = c("25", "50", "100", "200", "400", "800")
levels_p = str_c('p=', c("50", "100", "200", "400"))
levels_model = c('DC', 'glasso', 'SCAD', 'adapt')
levels_p
result_exp1

# F1 score
y_axis_max = result_exp1$mean |> max()
p1 = result_exp1 |>
    filter(p == 'p=50') |>
    ggplot(aes(x=factor(n, levels=levels_n), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(a) p=50') +
    ylim(0, y_axis_max)
p2 = result_exp1 |>
    filter(p == 'p=100') |>
    ggplot(aes(x=factor(n, levels=levels_n), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(b) p=100') +
    ylim(0, y_axis_max)
p3 = result_exp1 |>
    filter(p == 'p=200') |>
    ggplot(aes(x=factor(n, levels=levels_n), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(c) p=200') +
    ylim(0, y_axis_max)
p4 = result_exp1 |>
    filter(p == 'p=400') |>
    ggplot(aes(x=factor(n, levels=levels_n), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(d) p=400') +
    ylim(0, y_axis_max)
x_label = textGrob("Sample size", gp=gpar(fontsize=14), vjust=0.3)
y_label = textGrob("F1 score", gp=gpar(fontsize=14), rot=90, vjust=1.5)

legend_exp1 = get_legend(p1 + theme(legend.position = "right"))

fig_exp1_n_vs_f1 = grid.arrange(
    arrangeGrob(
        p1 + theme(legend.position = "none"),
        p2 + theme(legend.position = "none"),
        p3 + theme(legend.position = "none"),
        p4 + theme(legend.position = "none"),
        ncol=2),
    legend_exp1,
    ncol=2, widths=c(4, 1),
    bottom = x_label,
    left = y_label)
save_pdf_safely(fig=fig_exp1_n_vs_f1, path = str_c(path_export, 'figure/exp1_n_vs_f1_', METHOD, '.pdf'), width = 6, height = 4)

# Number of estimated edges
p1 = result_exp1 |>
    filter(p == 'p=50') |>
    ggplot(aes(x=factor(n, levels=levels_n), y=n_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = n_mean - 2*n_std, ymax = n_mean+2*n_std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(a) p=50')
p2 = result_exp1 |>
    filter(p == 'p=100') |>
    ggplot(aes(x=factor(n, levels=levels_n), y=n_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = n_mean - 2*n_std, ymax = n_mean+2*n_std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(b) p=100')
p3 = result_exp1 |>
    filter(p == 'p=200') |>
    ggplot(aes(x=factor(n, levels=levels_n), y=n_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = n_mean - 2*n_std, ymax = n_mean+2*n_std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(c) p=200')
p4 = result_exp1 |>
    filter(p == 'p=400') |>
    ggplot(aes(x=factor(n, levels=levels_n), y=n_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = n_mean - 2*n_std, ymax = n_mean+2*n_std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(d) p=400')

x_label = textGrob("Sample size", gp=gpar(fontsize=14), vjust=0.3)
y_label = textGrob("Number of estimated edges", gp=gpar(fontsize=14), rot=90, vjust=1.5)

legend_exp1 = get_legend(p1 + theme(legend.position = "right"))

fig_exp1_n_vs_numedge = grid.arrange(
    arrangeGrob(
        p1 + theme(legend.position = "none"),
        p2 + theme(legend.position = "none"),
        p3 + theme(legend.position = "none"),
        p4 + theme(legend.position = "none"),
        ncol=2),
    legend_exp1,
    ncol=2, widths=c(4, 1),
    bottom = x_label,
    left = y_label)

save_pdf_safely(fig=fig_exp1_n_vs_numedge, path = str_c(path_export, 'figure/exp1_n_vs_numedge_', METHOD, '.pdf'), width = 6, height = 4)


# Exp 6: log-likelihood during cross validation ----
METHOD = 'chain'
result_exp6 = read_csv(str_c(path_export, 'table/result_exp6_', METHOD, '_240615.csv'))

levels_n = c("25", "50", "100", "200", "400", "800")
levels_p = str_c('p=', c("50", "`100", "200", "400"))
levels_model = c('DC', 'glasso', 'SCAD', 'adapt')
levels_pn = unique(result_exp6$n_p)
p1 = result_exp6 |>
    filter(n_p == '25_50') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(a) p=50, n=25')
p2 = result_exp6 |>
    filter(n_p == '50_50') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(b) p=50, n=50')
p3 = result_exp6 |>
    filter(n_p == '100_50') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(c) p=50, n=100')
p4 = result_exp6 |>
    filter(n_p == '50_100') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(d) p=100, n=50')
p5 = result_exp6 |>
    filter(n_p == '100_100') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(e) p=100, n=100')
p6 = result_exp6 |>
    filter(n_p == '200_100') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(f) p=100, n=200')
p7 = result_exp6 |>
    filter(n_p == '100_200') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(g) p=200, n=100')
p8 = result_exp6 |>
    filter(n_p == '200_200') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(h) p=200, n=200')
p9 = result_exp6 |>
    filter(n_p == '400_200') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(i) p=200, n=400')
p10 = result_exp6 |>
    filter(n_p == '200_400') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(j) p=400, n=200')
p11 = result_exp6 |>
    filter(n_p == '400_400') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(k) p=400, n=400')
p12 = result_exp6 |>
    filter(n_p == '800_400') |>
    ggplot(aes(x=num_edge, y=loglik_mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), size=0.3) +
    geom_line(aes(color = factor(model, levels=levels_model)), size=1.2, alpha=0.7) +
    geom_vline(xintercept=30, linetype="dashed", color="black") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x, format = function(x) sprintf("%.0f", x)))) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(1.3, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(l) p=400, n=800')
x_label = textGrob("Number of estimated edges", gp=gpar(fontsize=14), vjust=0.3)
y_label = textGrob("Log-likelihood", gp=gpar(fontsize=14), rot=90, vjust=1.5)

legend_exp6 = get_legend(p1 + theme(legend.position = "right"))

fig_exp6_edgenum_vs_loglik = grid.arrange(
    arrangeGrob(
        p1 + theme(legend.position = "none"),
        p2 + theme(legend.position = "none"),
        p3 + theme(legend.position = "none"),
        p4 + theme(legend.position = "none"),
        p5 + theme(legend.position = "none"),
        p6 + theme(legend.position = "none"),
        p7 + theme(legend.position = "none"),
        p8 + theme(legend.position = "none"),
        p9 + theme(legend.position = "none"),
        p10 + theme(legend.position = "none"),
        p11 + theme(legend.position = "none"),
        p12 + theme(legend.position = "none"),
        ncol=3),
    legend_exp6,
    ncol=3, widths=c(4, 1, 0),
    bottom = x_label,
    left = y_label)

save_pdf_safely(fig=fig_exp6_edgenum_vs_loglik, path = str_c(path_export, 'figure/exp6_edgenum_vs_loglik_', METHOD, '.pdf'), width = 9, height = 9)


# Exp2: estimate fixed number of nonzeo-elements ----
METHOD = 'random'
result_exp2 = read_csv(str_c(path_export, 'table/result_exp2_', METHOD, '_240512.csv'))

levels_n = c("25", "50", "100", "200", "400", "800")
levels_p = str_c('p=', c("50", "`100", "200", "400"))
levels_model = c('DC', 'glasso', 'SCAD', 'adapt')
result_exp2 = result_exp2 |> mutate(n_p = str_c(n, '_', p))
levels_pn = unique(result_exp2$n_p)
y_axis_max = result_exp2$mean |> max()
p1 = result_exp2 |>
    filter(n_p == '25_50') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.3, margin = margin(t = -10))) +
    labs(x='', y='', col='', fill='', caption='(a) p=50, n=25') +
    ylim(0, y_axis_max)
p2 = result_exp2 |>
    filter(n_p == '50_50') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(b) p=50, n=50') +
    ylim(0, y_axis_max)
p3 = result_exp2 |>
    filter(n_p == '100_50') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(c) p=50, n=100') +
    ylim(0, y_axis_max)
p4 = result_exp2 |>
    filter(n_p == '50_100') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(d) p=100, n=50') +
    ylim(0, y_axis_max)
p5 = result_exp2 |>
    filter(n_p == '100_100') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(e) p=100, n=100') +
    ylim(0, y_axis_max)
p6 = result_exp2 |>
    filter(n_p == '200_100') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(f) p=100, n=200') +
    ylim(0, y_axis_max)
p7 = result_exp2 |>
    filter(n_p == '100_200') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(g) p=200, n=100') +
    ylim(0, y_axis_max)
p8 = result_exp2 |>
    filter(n_p == '200_200') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(h) p=200, n=200') +
    ylim(0, y_axis_max)
p9 = result_exp2 |>
    filter(n_p == '400_200') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(i) p=200, n=400') +
    ylim(0, y_axis_max)
p10 = result_exp2 |>
    filter(n_p == '200_400') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(j) p=400, n=200') +
    ylim(0, y_axis_max)
p11 = result_exp2 |>
    filter(n_p == '400_400') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(k) p=400, n=400') +
    ylim(0, y_axis_max)
p12 = result_exp2 |>
    filter(n_p == '800_400') |>
    ggplot(aes(x=as.factor(num_edge), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model)), position=position_dodge(.9), show.legend=FALSE) +
    geom_errorbar(aes(ymin = mean - 2*std, ymax = mean+2*std, color=factor(model, levels=levels_model)), position=position_dodge(.9), width = 0.4, alpha = 1, show.legend = FALSE) +
    geom_bar(aes(col=factor(model, levels=levels_model), fill=factor(model, levels=levels_model)), color='gray', alpha=0.5, position='dodge', stat='identity', show.legend = TRUE) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 15), legend.key.size = unit(0.7, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = -10))) +
    labs(x='', y='', color='', fill='', caption='(l) p=400, n=800') +
    ylim(0, y_axis_max)
x_label = textGrob("Number of estimated edges", gp=gpar(fontsize=14), vjust=0.3)
y_label = textGrob("F1 score", gp=gpar(fontsize=14), rot=90, vjust=1.5)

legend_exp2 = get_legend(p1 + theme(legend.position = "right"))

fig_exp2_same_edge_num = grid.arrange(
    arrangeGrob(
        p1 + theme(legend.position = "none"),
        p2 + theme(legend.position = "none"),
        p3 + theme(legend.position = "none"),
        p4 + theme(legend.position = "none"),
        p5 + theme(legend.position = "none"),
        p6 + theme(legend.position = "none"),
        p7 + theme(legend.position = "none"),
        p8 + theme(legend.position = "none"),
        p9 + theme(legend.position = "none"),
        p10 + theme(legend.position = "none"),
        p11 + theme(legend.position = "none"),
        p12 + theme(legend.position = "none"),
        ncol=3),
    legend_exp2,
    ncol=3, widths=c(4, 1, 0),
    bottom = x_label,
    left = y_label)

save_pdf_safely(fig=fig_exp2_same_edge_num, path = str_c(path_export, 'figure/exp2_same_edge_num_', METHOD, '.pdf'), width = 9, height = 9)


# Exp4: measure exeution time for each mode ----
METHOD = 'chain'
result_exp4 = read_csv(str_c(path_export, 'table/exp4_time_comparison_', METHOD, '_240728.csv'))

levels_n = str_c("50", "100", "200", "400")
levels_model = c('DC', 'glasso', 'SCAD', 'adapt')
time_summary = result_exp4 |> dplyr::select(-n, -p) |> pivot_longer(-n_p) |> group_by(n_p, name) |>
    summarise(mean = mean(value), std = sd(value) / sqrt(n_simulation)) |> rename(model = name)
time_summary_np = time_summary$n_p |> str_split("_", simplify=TRUE)
time_summary['p'] = as.numeric(time_summary_np[,2])
time_summary['n'] = str_c(time_summary_np[,1])
time_summary = time_summary |> filter(n_p %in% c('100_50', '100_100', '100_200', '400_200', '400_400', '400_800', '50_100', '100_100', '200_100', '200_400', '400_400', '800_400'))
# p vs time
p1 = time_summary |> filter(model %in% c('DC', 'glasso')) |>
    filter(n == 100, p %in% c(50, 100, 200)) |>
    ggplot(aes(x=p, y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model), shape=factor(model, levels=levels_model)), size=3) +
    geom_line(aes(color=factor(model, levels=levels_model))) +
    scale_shape_manual(values=c('DC'=19, 'glasso'=17)) +
    scale_x_continuous(breaks = unique(time_summary$p)) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = 1)), axis.title.x = element_text(size = 13)) +
    labs(x='Number of variables', y='', color=NULL, shape=NULL, caption='(a) n=100')
p2 = time_summary |> filter(model %in% c('DC', 'glasso')) |>
    filter(n == 400, p %in% c(200, 400, 800)) |>
    ggplot(aes(x=p, y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model), shape=factor(model, levels=levels_model)), size=3) +
    geom_line(aes(color=factor(model, levels=levels_model))) +
    scale_shape_manual(values=c('DC'=19, 'glasso'=17)) +
    scale_x_continuous(breaks = unique(time_summary$p)) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = 1)), axis.title.x = element_text(size = 13)) +
    labs(x='Number of variables', y='', color=NULL, shape=NULL, caption='(b) n=400')
x_label = textGrob(NULL, gp=gpar(fontsize=14), vjust=0.3)
y_label = textGrob("Computation time (s)", gp=gpar(fontsize=14), rot=90, vjust=1.5)

legend_exp4 = get_legend(p1 + theme(legend.position = "right"))

fig_exp4 = grid.arrange(
    arrangeGrob(
        p1 + theme(legend.position = "none"),
        p2 + theme(legend.position = "none"),
        ncol=2),
    legend_exp4,
    ncol=2, widths=c(4, 1),
    bottom = x_label,
    left = y_label)


save_pdf_safely(fig=fig_exp4, path = str_c(path_export, 'figure/exp4_p_vs_time_', METHOD, '.pdf'), width = 8, height = 4)

# n vs time
p1 = time_summary |> filter(model %in% c('DC', 'glasso')) |>
    filter(p == 100, n %in% c(50, 100, 200)) |>
    ggplot(aes(x=factor(n, levels=c(50, 100, 200)), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model), shape=factor(model, levels=levels_model)), size=3) +
    geom_line(aes(color=factor(model, levels=levels_model), group=model)) +
    scale_shape_manual(values=c('DC'=19, 'glasso'=17)) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = 1)), axis.title.x = element_text(size = 13)) +
    labs(x="Sample size", y='', color=NULL, shape=NULL, caption='(a) p=100')
p2 = time_summary |> filter(model %in% c('DC', 'glasso')) |>
    filter(p == 400, n %in% c(200, 400, 800)) |>
    ggplot(aes(x=factor(n, levels=c(200, 400, 800)), y=mean)) +
    geom_point(aes(color=factor(model, levels=levels_model), shape=factor(model, levels=levels_model)), size=3) +
    geom_line(aes(color=factor(model, levels=levels_model), group=model)) +
    scale_shape_manual(values=c('DC'=19, 'glasso'=17)) +
    theme_minimal(base_family="HiraKakuPro-W3") +
    theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), plot.caption = element_text(hjust = -0.1, margin = margin(t = 1)), axis.title.x = element_text(size = 13)) +
    labs(x="Sample size", y='', color=NULL, shape=NULL, caption='(b) p=400')

time_summary |> filter(model %in% c('DC', 'glasso')) |>
    filter(p == 100, n %in% c(50, 100, 200))

x_label = textGrob(NULL, gp=gpar(fontsize=14), vjust=0.3)
y_label = textGrob("Computation time (s)", gp=gpar(fontsize=14), rot=90, vjust=1.5)

legend_exp4 = get_legend(p1 + theme(legend.position = "right"))
fig_exp4 = grid.arrange(
    arrangeGrob(
        p1 + theme(legend.position = "none"),
        p2 + theme(legend.position = "none"),
        ncol=2),
    legend_exp4,
    ncol=2, widths=c(4, 1),
    bottom = x_label,
    left = y_label)

save_pdf_safely(fig=fig_exp4, path = str_c(path_export, 'figure/exp4_n_vs_time_', METHOD, '.pdf'), width = 8, height = 4)
