# this script is to reproduce the main results figures of the paper -- figs 3 and 5 -- with different start points

#script DOB: 9/5/2025
#author: Kaleigh Davis, Postdoc UoG

#### packages and referencing #####
#load necessary pkgs
library(tidyverse)
library(janitor)
library(MCMCpack)
library(bayesplot)
library(MCMCvis)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)
library(colorspace)
library(purrr)
library(viridis)
library(beepr)

# get referencing set up for macarthur temp dependence function
source("R-scripts/kd-analysis-scripts/02-temp-dep-macarthur-KD.R") #this contains the macarthur translation function, with all parameters flexibly defined in the function for assigning at time of use, and the arrhenius function.

#load in distributions for parameter values.
# these are continuous distributions generated from empirical data using MCMC regression, in kd_analysis_repro.csv
param_vals <- read_csv(file = "data/processed-data/param_post_dists.csv")

#split these into dfs for each parameter
param_vals %>%
  mutate(parameter = str_replace(parameter, "resource_growth_rate", "rgr"),
         parameter = str_replace(parameter, "carrying_capacity", "k"),
         parameter = str_replace(parameter, "conversion_efficiency", "v"),
         parameter = str_replace(parameter, "mortality_rate", "m"),
         parameter = str_replace(parameter, "consumption rate", "c")) %>% 
  group_by(parameter) %>%
  group_split() %>%
  set_names(unique(param_vals$parameter)) %>%  # Set the names based on unique category values
  walk(~ assign(paste0(.x$parameter[1], "_post_dist"), .x, envir = .GlobalEnv))

#get summary stats for all parameters ########
param_sum <- param_vals %>%
  group_by(parameter) %>% 
  summarize(
    across(
      intercept,
      list(
        Mean = mean,
        Q1 = ~quantile(., 0.25),
        Median = median,
        Q3 = ~quantile(., 0.75),
        Min = min,
        Max = max,
        sd = sd
      ),
      .names = "{.fn}" )
  ) 

#make longform
param_sum1 <- param_sum %>% 
  pivot_longer(cols = c(Mean:sd), 
               names_to = "summary_stat",
               values_to = "value")

#split df up into dfs for each summary statistic
param_sum1 %>% 
  group_by(summary_stat) %>% 
  group_split() %>% 
  purrr::walk(~ assign(paste0(.x$summary_stat[1]), .x, envir = .GlobalEnv))

#### ONE-BY-ONE; rrc, equal growth rates ####
# c ####
c_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  c_var <- bind_rows(c_var, hold) 
}

#log plot
c_var_plot <- 
  ggplot() +
  geom_path(data = c_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.55, label = expression("Consumption \nrate," ~ italic(c)[italic(ik)]), size = 6)


# r ####
r_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept,
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  r_var <- bind_rows(r_var, hold) 
}

#log plot
r_var_plot <- 
  ggplot() +
  geom_path(data = r_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Resource \ngrowth rate," ~ italic(r)[italic(k)]), size = 6)


# K ####
k_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept,
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  k_var <- bind_rows(k_var, hold) 
}

#log plot
k_var_plot <- 
  ggplot() +
  geom_path(data = k_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.55, label = expression("Resource \ncarrying \ncapacity," ~ italic(K)[italic(k)]), size = 6)

# v ####
v_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  v_var <- bind_rows(v_var, hold) 
}

#log plot
v_var_plot <- 
  ggplot() +
  geom_path(data = v_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(v_var, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Conversion \nefficiency," ~ italic(v)[italic(k)]), size = 6)

# m ####
m_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  m_var <- bind_rows(m_var, hold) 
}

#log plot
m_var_plot <- 
  ggplot() +
  geom_path(data = m_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(m_var, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.7, label = expression("Consumer \nmortality rate," ~ italic(m)[italic(i)]), size = 6)

# rrc equal base growth plot ###
rrce_plots <- c_var_plot + r_var_plot + k_var_plot + v_var_plot + m_var_plot
# ggsave(plot = rrce_plots, filename = "figures/kd-figs/param_var_startpoint1.pdf", width = 14, height = 10)

#### POMPOM; rrc, equal growth rates ####
rrce_all <- data.frame()
for(f in 1:500){ #was 200
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), #was by 0.1
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  rrce_all <- bind_rows(rrce_all, hold) 
}

#get average change in position after 5, 10, 20C warming
rrce_all_avg_new <- rrce_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_rrce_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrce_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrce_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrce_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrce_all_avg_new, aes(x = new_median_stab_pot, y = new_median_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrce_all_avg_new, aes(x = new_median_stab_pot, y = new_median_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrce_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

### rrc equal start two plots #####
startpoint1 <- rrce_plots / log_pom_rrce_all
# ggsave(plot = startpoint1, filename = "figures/kd-figs/supp_startpoint1.pdf", height = 20, width = 26)

#### ONE-BY-ONE; uneven reciprocal preference, equal growth rates ####
# c ####
c_var1 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.2, c1P_b = 1.8, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b = 2000, K_P_b = 2000, 
                      v1N_b = 0.2, v1P_b = 1.8, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  c_var1 <- bind_rows(c_var1, hold) 
}

#log plot
c_var1_plot <- 
  ggplot() +
  geom_path(data = c_var1, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var1, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.55, label = expression("Consumption \nrate," ~ italic(c)[italic(ik)]), size = 6)


# r ####
r_var1 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept,
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.2, c1P_b = 1.8, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.2, v1P_b = 1.8, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  r_var1 <- bind_rows(r_var1, hold) 
}

#log plot
r_var1_plot <- 
  ggplot() +
  geom_path(data = r_var1, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var1, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Resource \ngrowth rate," ~ italic(r)[italic(k)]), size = 6)


# K ####
k_var1 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept,
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.2, c1P_b = 1.8, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.2, v1P_b = 1.8, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  k_var1 <- bind_rows(k_var1, hold) 
}

k_var1col <- k_var1 %>% 
  mutate(inequality = ifelse(K_EaP < K_EaN, "P more negative",
                             ifelse(K_EaN < K_EaP, "N more negative", "other")))

#log plot
k_var1_plot <-
  ggplot() +
  geom_path(data = k_var1, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var1, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  # scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  # theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.55, label = expression("Resource \ncarrying \ncapacity," ~ italic(K)[italic(k)]), size = 6)

# v ####
v_var1 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.2, c1P_b = 1.8, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.2, v1P_b = 1.8, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  v_var1 <- bind_rows(v_var1, hold) 
}

#log plot
v_var1_plot <- 
  ggplot() +
  geom_path(data = v_var1, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(v_var1, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Conversion \nefficiency," ~ italic(v)[italic(k)]), size = 6)

# m ####
m_var1 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.2, c1P_b = 1.8, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.2, v1P_b = 1.8, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  m_var1 <- bind_rows(m_var1, hold) 
}

#log plot
m_var1_plot <- 
  ggplot() +
  geom_path(data = m_var1, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(m_var1, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.7, label = expression("Consumer \nmortality rate," ~ italic(m)[italic(i)]), size = 6)

# uneven reciprocal, equal base growth plot ###
urrce_plots <- c_var1_plot + r_var1_plot + k_var1_plot + v_var1_plot + m_var1_plot
# ggsave(plot = urrce_plots, filename = "figures/kd-figs/param_var1_startpoint2.pdf", width = 14, height = 10)

#### POMPOM; uneven reciprocal preference, equal growth rates ####
urrce_all <- data.frame()
for(f in 1:500){ #was 200
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), #was by 0.1
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.2, v1P_b = 1.8, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce_all <- bind_rows(urrce_all, hold) 
}

#get average change in position after 5, 10, 20C warming
urrce_all_avg_new <- urrce_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

### uneven reciprocal preference, equal start two plots #####
startpoint2 <- urrce_plots / log_pom_urrce_all
ggsave(plot = startpoint2, filename = "figures/kd-figs/supp_startpoint2.pdf", height = 20, width = 26)

#### ONE-BY-ONE; generalist-specialist ####
# c ####
c_var2 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.4, c1P_b = 0.8, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.8, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  c_var2 <- bind_rows(c_var2, hold) 
}

#log plot
c_var2_plot <- 
  ggplot() +
  geom_path(data = c_var2, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var2, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.7, label = expression("Consumption \nrate," ~ italic(c)[italic(ik)]), size = 6)


# r ####
r_var2 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept,
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.4, c1P_b = 0.8, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.8, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  r_var2 <- bind_rows(r_var2, hold) 
}

#log plot
r_var2_plot <- 
  ggplot() +
  geom_path(data = r_var2, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var2, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Resource \ngrowth rate," ~ italic(r)[italic(k)]), size = 6)


# K ####
k_var2 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept,
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.4, c1P_b = 0.8, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.8, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  k_var2 <- bind_rows(k_var2, hold) 
}

#log plot
k_var2_plot <- 
  ggplot() +
  geom_path(data = k_var2, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var2, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.45, y = 0.55, label = expression("Resource \ncarrying \ncapacity," ~ italic(K)[italic(k)]), size = 6)

# v ####
v_var2 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.4, c1P_b = 0.8, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.8, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  v_var2 <- bind_rows(v_var2, hold) 
}

#log plot
v_var2_plot <- 
  ggplot() +
  geom_path(data = v_var2, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(v_var2, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = 0.7, label = expression("Conversion \nefficiency," ~ italic(v)[italic(k)]), size = 6)

# m ####
m_var2 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.4, c1P_b = 0.8, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.8, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  m_var2 <- bind_rows(m_var2, hold) 
}

#log plot
m_var2_plot <- 
  ggplot() +
  geom_path(data = m_var2, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(m_var2, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.7, label = expression("Consumer \nmortality rate," ~ italic(m)[italic(i)]), size = 6)

# generalist-specialist growth plot ###
gs_plots <- c_var2_plot + r_var2_plot + k_var2_plot + v_var2_plot + m_var2_plot
ggsave(plot = gs_plots, filename = "figures/kd-figs/param_var2_startpoint3.pdf", width = 14, height = 10)

#### POMPOM; gs - strong ####
gs_all <- data.frame()
for(f in 1:500){ #was 200
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), #was by 0.1
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.4, c1P_b = 0.8, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.8, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  gs_all <- bind_rows(gs_all, hold) 
}

#get average change in position after 5, 10, 20C warming
gs_all_avg_new <- gs_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>%
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))
  

gs_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  ggplot() + 
  geom_histogram(aes(x = new_fit_ratio)) + 
  geom_vline(aes(xintercept = mean(new_fit_ratio))) +
  geom_vline(data = gs_all_avg_new, aes(xintercept = new_mean_fit_rat),colour = "red")

#pompom
log_pom_gs_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = gs_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(gs_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(gs_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = gs_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 4) +
  geom_point(data = gs_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 2) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = gs_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

### gs two plots #####
startpoint3 <- gs_plots / log_pom_gs_all
ggsave(plot = startpoint3, filename = "figures/kd-figs/supp_startpoint3.pdf", height = 20, width = 26)

#### ONE-BY-ONE; generalist-specialist - mild ####
# c ####
c_var3 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.4, c1P_b = 0.6, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.4, v1P_b = 0.6, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  c_var3 <- bind_rows(c_var3, hold) 
}

#log plot
c_var3_plot <- 
  ggplot() +
  geom_path(data = c_var3, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var3, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.7, label = expression("Consumption \nrate," ~ italic(c)[italic(ik)]), size = 6)


# r ####
r_var3 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept,
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.4, c1P_b = 0.6, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.4, v1P_b = 0.6, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  r_var3 <- bind_rows(r_var3, hold) 
}

#log plot
r_var3_plot <- 
  ggplot() +
  geom_path(data = r_var3, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var3, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Resource \ngrowth rate," ~ italic(r)[italic(k)]), size = 6)


# K ####
k_var3 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept,
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.4, c1P_b = 0.6, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.4, v1P_b = 0.6, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  k_var3 <- bind_rows(k_var3, hold) 
}

#log plot
k_var3_plot <- 
  ggplot() +
  geom_path(data = k_var3, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var3, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.45, y = 0.55, label = expression("Resource \ncarrying \ncapacity," ~ italic(K)[italic(k)]), size = 6)

# v ####
v_var3 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept,
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.4, c1P_b = 0.6, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.4, v1P_b = 0.6, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  v_var3 <- bind_rows(v_var3, hold) 
}

#log plot
v_var3_plot <- 
  ggplot() +
  geom_path(data = v_var3, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(v_var3, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = 0.7, label = expression("Conversion \nefficiency," ~ italic(v)[italic(k)]), size = 6)

# m ####
m_var3 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = 0,
                      c_Ea2N = 0,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.4, c1P_b = 0.6, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.4, v1P_b = 0.6, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  m_var3 <- bind_rows(m_var3, hold) 
}

#log plot
m_var3_plot <- 
  ggplot() +
  geom_path(data = m_var3, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(m_var3, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  # coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  #dims from full pompom plot figure
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #back to this plot code
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.7, label = expression("Consumer \nmortality rate," ~ italic(m)[italic(i)]), size = 6)

# generalist-specialist growth plot ###
gsm_plots <- c_var3_plot + r_var3_plot + k_var3_plot + v_var3_plot + m_var3_plot
# ggsave(plot = gs_plots, filename = "figures/kd-figs/param_var3_startpoint3.pdf", width = 14, height = 10)

#### POMPOM; gs-mild ####
gsm_all <- data.frame()
for(f in 1:500){ #was 200
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, 
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.3, c1P_b = 0.7, 
                      c2N_b = 0.5, c2P_b = 0.5, 
                      r_N_b = 0.1, r_P_b = 0.05, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.3, v1P_b = 0.7, 
                      v2N_b = 0.5, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  gsm_all <- bind_rows(gsm_all, hold) 
}

#get average change in position after 5, 10, 20C warming
gsm_all_avg_new <- gsm_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_gsm_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = gsm_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(gsm_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(gsm_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = gsm_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = gsm_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = gsm_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

### gs two plots #####
startpoint4 <- gsm_plots / log_pom_gsm_all
ggsave(plot = startpoint4, filename = "figures/kd-figs/supp_startpoint4.pdf", height = 20, width = 26)
