par(family= "HiraKakuProN-W3")
path_data = './data/'
path_export = './output/'

# Library ----
library(dplyr)
library(tidyverse)
library(tictoc)

source('dc_ggm.R')
source('data_prep.R')

METHOD = 'random'
n_simulation = 10
num_nonzero = 30
seed_list = seq(1, n_simulation)

# Exp4: measure exeution time for each mode ----
n_list = c(100, 400)
time_result1 = data.frame()
for (n in n_list) {
    p_list = as.integer(c(n / 2, n, 2*n))
    for (p in p_list) {
        for (seed in seed_list) {
            data_info = generate_synthetic_data(p=p, n=n, method=METHOD, num_nonzero=num_nonzero, seed=seed)
            K = ceiling(p*(p-1) / 2 / 2)
            tic()
            dc_result = dc_ggm(S=data_info$cov, K=K)
            a = toc(quiet = TRUE)
            time_dc = as.numeric(a$toc - a$tic)

            lambda = abs(data_info$cov[lower.tri(data_info$cov)]) |> median()
            tic()
            glasso_result = glasso(s=data_info$cov, rho=lambda)
            a = toc(quiet = TRUE)
            time_glasso = as.numeric(a$toc - a$tic)

            tic()
            temp = ggmncv(R=data_info$cov, n=n, penalty='scad', lambda=lambda, progress=FALSE)
            a = toc(quiet = TRUE)
            time_scad = as.numeric(a$toc - a$tic)

            tic()
            temp = ggmncv(R=data_info$cov, n=n, penalty='adapt', lambda=lambda, progress=FALSE)
            a = toc(quiet = TRUE)
            time_adapt = as.numeric(a$toc - a$tic)

            time_result1 = rbind(time_result1, data.frame(
                p = p,
                n = n,
                time_dc = time_dc,
                time_glasso = time_glasso,
                time_scad = time_scad,
                time_adapt = time_adapt
            ))
        }
        cat("(n, p)=(",n,",",p,") ")
    }
}

p_list = c(100, 400)
time_result2 = data.frame()
for (p in p_list) {
    n_list = as.integer(c(p / 2, p, 2*p))
    for (n in n_list) {
        for (seed in seed_list) {
            data_info = generate_synthetic_data(p=p, n=n, method=METHOD, num_nonzero=num_nonzero, seed=seed)
            K = ceiling(p*(p-1) / 2 / 2)
            tic()
            dc_result = dc_ggm(S=data_info$cov, K=K)
            a = toc(quiet = TRUE)
            time_dc = as.numeric(a$toc - a$tic)

            lambda = abs(data_info$cov[lower.tri(data_info$cov)]) |> median()
            tic()
            glasso_result = glasso(s=data_info$cov, rho=lambda)
            a = toc(quiet = TRUE)
            time_glasso = as.numeric(a$toc - a$tic)

            tic()
            temp = ggmncv(R=data_info$cov, n=n, penalty='scad', lambda=lambda, progress=FALSE)
            a = toc(quiet = TRUE)
            time_scad = as.numeric(a$toc - a$tic)

            tic()
            temp = ggmncv(R=data_info$cov, n=n, penalty='adapt', lambda=lambda, progress=FALSE)
            a = toc(quiet = TRUE)
            time_adapt = as.numeric(a$toc - a$tic)

            time_result2= rbind(time_result2, data.frame(
                p = p,
                n = n,
                time_dc = time_dc,
                time_glasso = time_glasso,
                time_scad = time_scad,
                time_adapt = time_adapt
            ))
        }
        cat("(n, p)=(",n,",",p,") ")
    }
}

time_result = rbind(time_result1, time_result2)
time_result
colnames(time_result)[3:6] = c('DC', 'glasso', 'SCAD', 'adapt')
time_result = time_result |> mutate(n_p = str_c(n, '_', p))
today_date = format(Sys.Date(), "%y%m%d")
save_table_safely(data=time_result, path = str_c(path_export, 'table/exp4_time_comparison_', METHOD, '_', today_date, '.csv'))


# Visualization ----
levels_n = str_c('N=', c("50", "100", "200", "400"))
levels_model = c('DC', 'glasso', 'SCAD', 'adapt')

time_summary_original = time_result |> dplyr::select(-n, -p) |> pivot_longer(-n_p) |> group_by(n_p, name) |>
    summarise(mean = mean(value), std = sd(value) / sqrt(n_simulation)) |> rename(model = name)
time_summary_np = time_summary_original$n_p |> str_split("_", simplify=TRUE)
time_summary_original['p'] = as.numeric(time_summary_np[,2])
time_summary_original['n'] = str_c('N=', time_summary_np[,1])

time_original_gg = time_summary_original |> filter(model %in% c('DC', 'glasso')) |>
  ggplot(aes(x=p, y=mean)) +
  facet_wrap(~factor(n, levels=levels_n), scale='free', ncol=2) +
  geom_point(aes(color=factor(model, levels=levels_model), shape=factor(model, levels=levels_model)), size=3) +
  geom_line(aes(color=factor(model, levels=levels_model))) +
  scale_shape_manual(values=c('DC'=19, 'glasso'=17)) +
  theme_minimal(base_family="HiraKakuPro-W3") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"), panel.spacing = unit(2, "lines")) +
  labs(x='p', y='sec', color=NULL, shape=NULL)
save_pdf_safely(fig=time_original_gg, path=str_c(path_export, 'figure/exp4_time_comparison_', METHOD, '_', today_date, '.pdf'), width=9, height=8)

