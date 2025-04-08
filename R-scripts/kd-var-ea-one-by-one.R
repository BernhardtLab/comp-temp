#this script is to explore the effects of temperature dependence of different parameters on changes in niche and fitness differences

#script DOB: April 1, 2025
#author: Kaleigh Davis, University of Guelph postdoc

# Objective 1: plot out effects of temp dep on each parameter while each other is held at mean -- look at total effects and relative effects on niche and fitness differences -- DONE 4/1/2025
# Objective 2: Using the logged definitions of niche and fitness differences, calculate euclidean distances to 20C warming point and see how TAs affect this. 

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

# get referencing set up for macarthur temp dependence function
source("R-scripts/temp-dep-macarthur-KD.R") #this contains the macarthur translation function, with all parameters flexibly defined in the function for assigning at time of use, and the arrhenius function.

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
        Max = max
      ),
      .names = "{.fn}" )
  ) 

#make longform
param_sum1 <- param_sum %>% 
  pivot_longer(cols = c(Mean:Max), 
               names_to = "summary_stat",
               values_to = "value")

#split df up into dfs for each summary statistic
param_sum1 %>% 
  group_by(summary_stat) %>% 
  group_split() %>% 
  purrr::walk(~ assign(paste0(.x$summary_stat[1]), .x, envir = .GlobalEnv))

# basic simulation setup -- here all param EAs are drawn from distribution, consumers have reciprocal resource use, and N grows faster than P at ref temp
r_var <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  r_var <- bind_rows(r_var, hold) 
}

# plot -- I've got the axes and shaded region converted correctly, I think, but I cannot get the start point to appear on the line
# holding off on the log-log plotting for now
# ggplot() +
#   geom_path(data = r_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 2) +
#   geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
#               aes(x = x,
#                   y = NULL,
#                   ymin = exp(-x),
#                   ymax = 1/(exp(-x))),
#               fill = "grey", color = "black", alpha = 0.2) +
#   geom_point(data = filter(r_var, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 4) +
#   geom_hline(yintercept = 1, linetype = 5) +
#   # scale_colour_continuous_diverging() +
#   # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) +
#   xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
#   ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))")))

#regular plot
ggplot() +
  geom_path(data = r_var, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #strong effects along fitness differences axis

# calculate euclidean distances at 20C for each iteration ####
r_var_e <- r_var %>% 
  filter(T %in% c(25, 45)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(r_ta = abs(r_EaN - r_EaP),
         dist20 = sqrt((T45_new_stabil_potential - T25_new_stabil_potential)^2 + (T45_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T45_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T45_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist20, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

r_var_e %>% 
  # dplyr::select(iteration, r_ta, dist20, shift_fitrat, shift_nichediffs) %>% 
  pivot_longer(cols = c(dist20, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") %>% 
  ggplot(aes(x = r_ta, y = value, colour = )) + 
  geom_point() + 
  facet_wrap(~response_var) + 
  labs(x = "Abs diff in r_Ea", y = "Response variable value") #need to think about this more
  


### vary r_EaN, not r_EaP ####
r_var1 <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  r_var1 <- bind_rows(r_var1, hold) 
}

#plot
ggplot() +
  geom_path(data = r_var1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #same trend in effect as both drawn from dist, but smaller range of effects, supporting that thermal asymmetries generate those

### vary all c_EAs ####
c_var <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  c_var <- bind_rows(c_var, hold) 
}

#plot
ggplot() +
  geom_path(data = c_var, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

### vary all c_EaPs ####
c_var1 <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  c_var1 <- bind_rows(c_var1, hold) 
}

#plot
ggplot() +
  geom_path(data = c_var1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #same pattern, shrunk in a bit.

### vary c_Eas for just one consumer####
c_var2 <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  c_var2 <- bind_rows(c_var2, hold) 
}

#plot
ggplot() +
  geom_path(data = c_var2, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var2, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #down to one line; still varying mostly along x axis

### vary K_Eas ####
k_var <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  k_var <- bind_rows(k_var, hold) 
}

#plot
ggplot() +
  geom_path(data = k_var, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #interesting! Only moves points along x axis & only down

### vary only one K_Ea ####
k_var1 <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  k_var1 <- bind_rows(k_var1, hold) 
}

#plot
ggplot() +
  geom_path(data = k_var1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #same, but constrained

### vary v_Eas ####
v_var <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept, 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  v_var <- bind_rows(v_var, hold) 
}

#plot
ggplot() +
  geom_path(data = v_var, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(v_var, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #interesting! Only moves points along x axis & only down

### vary v_EaN only ####
v_var1 <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  v_var1 <- bind_rows(v_var1, hold) 
}

#plot
ggplot() +
  geom_path(data = v_var1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(v_var1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #interesting! Only moves points along x axis & only down

### vary m_Eas ####
m_var <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.10, m2_b = 0.10) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  m_var <- bind_rows(m_var, hold) 
}

#plot
ggplot() +
  geom_path(data = m_var, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(m_var, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #pretty much does not move the point

### vary m_Ea1 only ####
m_var1 <- data.frame()
for(f in 1:200){ 
  hold = temp_dep_mac(T = seq(25, 50, by = 1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.10, m2_b = 0.10) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  m_var1 <- bind_rows(m_var1, hold) 
}

#plot
ggplot() +
  geom_path(data = m_var1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(m_var1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #pretty much does not move the point

### trying to force conditions where m_EA has an impact ####
#This seems to only be possible by really increasing the temperature range over which we run the simulation. The lowest where you can see any effect at all is around 75C warming. Better seen around 75C warming.
m_var2 <- data.frame()
for(f in 1:40){ 
  hold = temp_dep_mac(T = seq(25, 125, by = 1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.10, m2_b = 0.10) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  m_var2 <- bind_rows(m_var2, hold) 
}

#plot
ggplot() +
  geom_path(data = m_var2, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(m_var2, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
  facet_wrap(~iteration) #only moves up and down -- still just a couple scenarios that have this effect

m_var2 %>% 
  filter(iteration == 39) %>% 
  dplyr::select(c(m_Ea1, m_Ea2)) 

m_var_tas <- m_var2 %>% 
  group_by(iteration) %>% 
  mutate(m_ta = abs(m_Ea1 - m_Ea2)) %>% 
  distinct(iteration, m_ta, .keep_all = T)

m_var_tas %>% 
  ungroup() %>% 
  mutate(iteration = as.factor(iteration)) %>% 
  mutate(iteration = fct_reorder(iteration, m_ta)) %>%
  ggplot(aes(x = iteration, y = m_ta)) + 
  geom_point() #cross-checking this against the facet plot above confirms that it is only when temps can go very high and thermal asymmetries are very large that we see m_EA affect niche/fitness differences



