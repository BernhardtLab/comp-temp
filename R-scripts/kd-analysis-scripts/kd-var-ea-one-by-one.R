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
library(viridis)
library(beepr)
library(furrr) #for running large simulations in parallel

# get referencing set up for macarthur temp dependence function
source("R-scripts/kd-analysis-scripts/temp-dep-macarthur-KD.R") #this contains the macarthur translation function, with all parameters flexibly defined in the function for assigning at time of use, and the arrhenius function.

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


# r_Ea varies --------------------------------------------------------
# basic simulation setup -- here all param EAs are drawn from distribution, consumers have reciprocal resource use, and N grows faster than P at ref temp
r_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
                      # r_EaN = sample_n(k_post_dist, size = 1)$intercept, #testing whether a dist centered on a negative value will drive the same pattern -- it does
                      # r_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      # r_EaN = sample_n(v_post_dist, size = 1)$intercept, #testing whether a dist centered on 0 will drive the same pattern -- this is extreme! 
                      # r_EaP = sample_n(v_post_dist, size = 1)$intercept,
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  r_var <- bind_rows(r_var, hold) 
}
beep(2)

r_var_plot <-
  ggplot() +
  geom_path(data = r_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.70)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C \nWarming") +
  annotate("text", x = 0.3, y = 0.85, label = expression("Resource growth rate," ~ italic(r)), size = 6) +
  theme_cowplot(font_size = 20) +
  theme(legend.position = "none")

#get legend for composite plot -- need to generate r_var_plot WITH the legend first, then add in the legend.position = "none" for the actual composite plot
# rvar_legend  <- get_legend(r_var_plot)

#regular plot
ggplot() +
  geom_path(data = r_var, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  # ymin = log(1 - x),
                  ymax = 1/(1-x)),
                  # ymax = log(1/(1-x))),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) +
  xlab(expression(paste("Niche differences (1-", rho, ")"))) +
  ylab(expression(paste("Fitness differences (", f[2], "/", f[1], ")"))) #strong effects along fitness differences axis


## calculate euclidean distances at 20C for each iteration #####
r_var_e <- r_var %>% 
  filter(T %in% c(25, 40)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(r_ta = abs(r_EaN - r_EaP),
         dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

r_var_e %>% 
  # ggplot(aes(x = r_ta, y = value)) + 
  ggplot(aes(x = r_EaN, y = value, colour = r_EaP)) +
  geom_point() + 
  facet_wrap(~response_var) 
  # labs(x = "Abs diff in r_Ea", y = "Response variable value") #need to think about this more

# thermal asymmetry - euclidean distance plot
# r_var_plot_e <-
#   r_var_e %>%
#   filter(response_var == "dist15") %>% 
#   ggplot(aes(x = scale(r_ta), y = value)) +
#   # ggplot(aes(x = scale(r_ta), y = value, colour = r_EaN < r_EaP)) +
#   geom_point(size = 3) + 
#   labs(x = "Scaled absolute value \nof thermal asymmetry", y = "Displacement of species pair \n(Euclidean distance) after 15C warming") + 
#   coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(0, 0.5)) + 
#   annotate("text", x = 0, y = 0.45, label = "Resource \ngrowth rate, r", size = 5.5) + 
#   theme_cowplot(font_size = 20)
  
#unscales TA - euclidean distance plot
r_var_plot_e2 <-
  r_var_e %>% 
    filter(response_var == "dist15") %>% 
    # ggplot(aes(x = r_ta, y = value)) +
    ggplot(aes(x = r_ta, y = value, colour = r_EaN < r_EaP)) +
    geom_point(size = 3) + 
    labs(x = "Absolute value \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
    coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
    annotate("text", x = 0.88, y = 0.35, label = "Resource \ngrowth rate, r", size = 5.5) + 
    theme_cowplot(font_size = 20) 

r_var_plot_e3 <-
  r_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = r_EaN, y = r_EaP, colour = value)) + 
  labs(x = "Temperature dependence \nof resource 1 growth rate", y = "Temperature dependence \nof resource 2 growth rate", colour = "ED") + 
  geom_point(size = 3) + 
  annotate("text", x = 1.07, y = 0.5, label = "r", size = 7, fontface = "italic")
  
# calculate change in alpha in addition to thermal asymmetry ####
r_var_alpha <- r_var %>% 
  filter(T %in% c(25, 40)) %>% 
  dplyr::select(-c(g1, g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(a11, a22, a12, a21, new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(r_ta = abs(r_EaN - r_EaP),
         a12_shift_abs = abs(T40_a12 - T25_a12),
         a21_shift_abs = abs(T40_a21 - T25_a21),
         a11_shift_abs = abs(T40_a11 - T25_a11),
         a22_shift_abs = abs(T40_a22 - T25_a22),
         dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 
  
r_var_alpha %>% 
  ggplot(aes(x = r_ta, y = a21_shift_abs, colour = a12_shift_abs)) + 
  geom_point() + 
  facet_wrap(~response_var) + 
  scale_colour_viridis_c(option = "inferno") +
  labs(x = "Abs change in a21", y = "Response variable value") #need to think about this more

r_var_alpha %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = a21_shift_abs, y = a12_shift_abs, colour = r_ta)) + 
  geom_point() + 
  facet_wrap(~response_var) + 
  scale_colour_viridis_c(option = "inferno") 
  # coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 0.5)) 

r_var_alpha %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = a11_shift_abs, y = a22_shift_abs, colour = r_ta)) + 
  geom_point() + 
  scale_colour_viridis_c(option = "inferno") + 
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 0.5)) 

### restrict draws of r_EaN relative to that of r_EaP ####
# basic simulation setup -- here all param EAs are drawn from distribution, consumers have reciprocal resource use, and N grows faster than P at ref temp
r_var_res <- data.frame()
for(f in 1:500){ 
  
  r_EaP_val <- sample_n(rgr_post_dist, size = 1)$intercept
  
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = sample_n(filter(rgr_post_dist, intercept < r_EaP_val), size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = r_EaP_val,
                      # r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 0.5, r_P_b = 1, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  r_var_res <- bind_rows(r_var_res, hold) 
}

# # r_var_res_plot <-
  ggplot() +
  geom_path(data = r_var_res, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(r_var_res, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "inferno") +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.70)) +
  xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) +
  # labs(colour = "Degrees C \nWarming") +
  annotate("text", x = 0.2, y = 0.7, label = "Resource \ngrowth rate, r", size = 6) +
  theme_cowplot(font_size = 20) +
  theme(legend.position = "none")

#get legend for composite plot -- need to generate r_var_res_plot WITH the legend first, then add in the legend.position = "none" for the actual composite plot
# rvar_legend  <- get_legend(r_var_res_plot)

#regular plot
# ggplot() +
#   geom_path(data = r_var_res, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
#   geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
#               aes(x = x,
#                   y = NULL,
#                   ymin = 1-x,
#                   # ymin = log(1 - x),
#                   ymax = 1/(1-x)),
#               # ymax = log(1/(1-x))),
#               fill = "grey", color = "black", alpha = 0.2) +
#   geom_point(data = filter(r_var_res, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
#   geom_hline(yintercept = 1, linetype=5) +
#   # scale_colour_continuous_diverging() +
#   # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) +
#   xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
#   ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #strong effects along fitness differences axis


## calculate euclidean distances at 20C for each iteration #####
r_var_res_e <- r_var_res %>% 
  filter(T %in% c(25, 40)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(r_ta = abs(r_EaN - r_EaP),
         dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

r_var_res_e %>% 
  # dplyr::select(iteration, r_ta, dist20, shift_fitrat, shift_nichediffs) %>% 
  # ggplot(aes(x = r_ta, y = value, colour = r_EaN > r_EaP)) + 
  ggplot(aes(x = r_ta, y = value)) + 
  geom_point() + 
  facet_wrap(~response_var) + 
  labs(x = "Abs diff in r_Ea", y = "Response variable value") #need to think about this more

r_var_res_plot_e <-
  r_var_res_e %>%
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = scale(r_ta), y = value)) +
  # ggplot(aes(x = scale(r_ta), y = value, colour = r_EaN < r_EaP)) +
  geom_point(size = 3) + 
  labs(x = "Scaled absolute value \nof thermal asymmetry", y = "Displacement of species pair \n(Euclidean distance) after 15C warming") + 
  coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(0, 0.5)) + 
  annotate("text", x = 0, y = 0.45, label = "Resource \ngrowth rate, r", size = 5.5) + 
  theme_cowplot(font_size = 20)

r_var_res_plot_e2 <-
  r_var_res_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = r_ta, y = value)) +
  # ggplot(aes(x = r_ta, y = value, colour = r_EaN < r_EaP)) +
  geom_point(size = 3) + 
  labs(x = "Absolute value \nof thermal asymmetry", y = "Displacement of species pair \n(Euclidean distance) after 15C warming") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 1.0, y = 0.05, label = "Resource \ngrowth rate, r", size = 5.5) + 
  theme_cowplot(font_size = 20) 

#### vary r_EaN, not r_EaP ####
r_var1 <- data.frame()
for(f in 1:500){ 
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


# c_Ea varies -------------------------------------------------------------

### vary all c_EAs ####
c_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  c_var <- bind_rows(c_var, hold) 
}

# #plot
# ggplot() +
#   geom_path(data = c_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 2) +
#   geom_ribbon(data = data.frame(x = seq(0, 0.7, 0.001)),
#               aes(x = x,
#                   y = NULL,
#                   ymin = exp(-x),
#                   ymax = 1/(exp(-x))),
#               fill = "grey", color = "black", alpha = 0.2) +
#   geom_point(data = filter(c_var, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 4) +
#   geom_hline(yintercept = 1, linetype = 5) +
#   scale_colour_viridis_c(option = "inferno") +
#   coord_cartesian(ylim=c(0, 1.2), xlim = c(0, 0.7)) +
#   xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
#   ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
#   labs(colour = "Degrees C Warming")

#log plot
c_var_plot <-
  ggplot() +
  geom_path(data = c_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") 
  # annotate("text", x = 0.25, y = 0.7, label = expression("Consumption rate," ~ italic(c)), size = 6)

### vary all c_EaPs ####
c_var1 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  c_var1 <- bind_rows(c_var1, hold) 
}

# #plot
# ggplot() +
#   geom_path(data = c_var1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
#   geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
#               aes(x = x,
#                   y = NULL,
#                   ymin = 1-x,
#                   ymax = 1/(1-x)),
#               fill = "grey", color = "black", alpha = 0.2) +
#   geom_point(data = filter(c_var1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
#   geom_hline(yintercept = 1, linetype=5) + 
#   # scale_colour_continuous_diverging() +
#   # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
#   xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
#   ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #same pattern, shrunk in a bit.

#log plot
c_var1_plot <-
  ggplot() +
  geom_path(data = c_var1, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var1, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") 
  # annotate("text", x = 0.25, y = 0.7, label = expression("Consumption rate," ~ italic(c)), size = 6)

### vary c_Eas for just one consumer - this is the one for the MS -- COME BACK TO THIS! ####
c_var2 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      K_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "carrying_capacity" & summary_stat == "Mean"), value)), 
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  c_var2 <- bind_rows(c_var2, hold) 
}
beep(2)

#reg plot
# ggplot() +
#   geom_path(data = c_var2, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
#   geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
#               aes(x = x,
#                   y = NULL,
#                   ymin = 1-x,
#                   ymax = 1/(1-x)),
#               fill = "grey", color = "black", alpha = 0.2) +
#   geom_point(data = filter(c_var2, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
#   geom_hline(yintercept = 1, linetype=5) + 
#   # scale_colour_continuous_diverging() +
#   # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
#   xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
#   ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) #same pattern, shrunk in a bit.

#log plot
c_var2_plot <-
  ggplot() +
  geom_path(data = c_var2, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var2, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.85, label = expression("Consumption rate," ~ italic(c)), size = 6)
  
#three types of consumption rate simulations
c_three <- c_var2_plot + c_var1_plot + c_var_plot + plot_annotation(tag_levels = "A")

#save
ggsave(plot = c_three, filename = "figures/kd-figs/c_drawtypes.pdf", width = 16, height = 10)


## calculate euclidean distances at 20C for each iteration #####
c_var_e <- c_var2 %>% 
  filter(T %in% c(25, 40)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(c_ta = abs(c_Ea2N - c_Ea2P),
         dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#effects on euclidean distance and niche and fitness diffs
c_var_e %>% 
  ggplot(aes(x = c_ta, y = value)) + 
  geom_point() + 
  facet_wrap(~response_var) + 
  labs(x = "Abs diff in r_Ea", y = "Response variable value") #need to think about this more

c_var_plot_e <-
  c_var_e %>% 
  filter(response_var == "dist15") %>% 
  # ggplot(aes(x = scale(c_ta), y = value)) + 
  ggplot(aes(x = scale(c_ta), y = value, colour = c_Ea2N < c_Ea2P)) +
  geom_point(size = 3) + 
  labs(x = "Scaled absolute value \nof thermal asymmetry", y = "Displacement of species pair \nafter 15°C warming \n(Euclidean distance)") + 
  coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(0, 0.5)) + 
  annotate("text", x = 0, y = 0.45, label = "Consumer 1 resource \nconsumption rate", size = 5.5) + 
  theme_cowplot(font_size = 20)

c_var_plot_e2 <-
  c_var_e %>% 
  filter(response_var == "dist15") %>% 
  # ggplot(aes(x = c_ta, y = value)) + 
  ggplot(aes(x = c_ta, y = value, colour = c_Ea2N < c_Ea2P)) +
  geom_point(size = 3) + 
  labs(x = "Absolute value \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 0.70, y = 0.35, label = "Consumer 2 \nconsumption rate, c", size = 5.5) + 
  theme_cowplot(font_size = 20)

c_var_plot_e3 <- c_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = c_Ea1N, y = c_Ea1P, colour = value)) + 
  geom_point(size = 3) +
  labs(x = "Temperature dependence \nof consumer 1 consumption rate \non resource 1", y = "Temperature dependence \nof consumer 1 consumption rate \non resource 2", colour = "ED") + 
  annotate("text", x = 1.07, y = 0.5, label = "c", size = 7, fontface = "italic")

# K_Eas vary --------------------------------------------------------------

### vary both K_Eas ####
k_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  k_var <- bind_rows(k_var, hold) 
}
beep(2)

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
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

#log plot
k_var_plot <- 
  ggplot() +
  geom_path(data = k_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Resource \ncarrying capacity," ~ italic(K)), size = 6)

## calculate euclidean distances at 20C for each iteration #####
k_var_e <- k_var %>% 
  filter(T %in% c(25, 40)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(k_ta = abs(K_EaN - K_EaP),
         dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

k_var_e %>% 
  ggplot(aes(x = k_ta, y = value, colour = K_EaN > K_EaP)) + 
  geom_point() + 
  facet_wrap(~response_var) + 
  labs(x = "Abs diff in K_Ea", y = "Response variable value") #need to think about this more

k_var_plot_e <-
  k_var_e %>% 
  filter(response_var == "dist15") %>% 
  # ggplot(aes(x = scale(k_ta), y = value)) + 
  ggplot(aes(x = scale(k_ta), y = value, colour = K_EaN < K_EaP)) +
  geom_point(size = 3) + 
  labs(x = "Scaled absolute value \nof thermal asymmetry", y = "Displacement of species pair \nafter 15°C warming \n(Euclidean distance) ") + 
  coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(0, 0.5)) + 
  annotate("text", x = 0.25, y = 0.45, label = "Resource \ncarrying capacity, K", size = 6) + 
  theme_cowplot(font_size = 20)

k_var_plot_e2 <-
  k_var_e %>% 
  filter(response_var == "dist15") %>% 
  # ggplot(aes(x = k_ta, y = value)) +
  ggplot(aes(x = k_ta, y = value, colour = K_EaN < K_EaP)) +
  geom_point(size = 3) + 
  labs(x = "Absolute value \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 0.7, y = 0.35, label = "Resource \ncarrying capacity, K", size = 6) + 
  theme_cowplot(font_size = 20)

k_var_plot_e3 <- k_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = K_EaN, y = K_EaP, colour = value)) + 
  labs(x = "Temperature dependence \nof resource 1 carrying capacity", y = "Temperature dependence \nof resource 2 carrying capacity", colour = "ED") + 
  annotate("text", x = 1.07, y = 0.5, label = "K", size = 7, fontface = "italic") +
  geom_point(size = 3)

### restrict so that K_EaN > K_EaP, or vice versa #####
k_var_res <- data.frame()
for(f in 1:500){ 
  
  #first draw K_EaN
  K_EaP_val <- sample_n(k_post_dist, size = 1)$intercept
  # K_EaN_val <- sample_n(k_post_dist, size = 1)$intercept
  
  #then draw all other params from dist, using 
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
                      ref_temp = 25,
                      r_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      r_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "resource_growth_rate" & summary_stat == "Mean"), value)), 
                      c_Ea1N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea1P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      c_Ea2N = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)),
                      c_Ea2P = unlist(dplyr::select(filter(param_sum1, parameter == "consumption rate" & summary_stat == "Mean"), value)), 
                      K_EaN = sample_n(filter(k_post_dist, intercept < K_EaP_val), size = 1)$intercept,
                      # K_EaN = sample_n(k_post_dist, size = 1)$intercept,
                      K_EaP = K_EaP_val,
                      # K_EaP = sample_n(k_post_dist, size = 1)$intercept,
                      v_EaN = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)),
                      v_EaP = unlist(dplyr::select(filter(param_sum1, parameter == "conversion_efficiency" & summary_stat == "Mean"), value)), 
                      m_Ea1 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)), 
                      m_Ea2 = unlist(dplyr::select(filter(param_sum1, parameter == "mortality_rate" & summary_stat == "Mean"), value)),
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 0.5, r_P_b = 1, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  k_var_res <- bind_rows(k_var_res, hold) 
}
beep(1)

#check var in K_EaN and K_eaP across iterations
hist(k_post_dist$intercept)
hist(k_var_res$K_EaN)
hist(k_var_res$K_EaP)


#plot
ggplot() +
  geom_path(data = k_var_res, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var_res, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  # scale_colour_continuous_diverging() +
  # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

#log plot
# k_var_res_plot <- 
  ggplot() +
  geom_path(data = k_var_res, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(k_var_res, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "inferno") +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = "Resource \ncarrying capacity, K", size = 6)

## calculate euclidean distances at 20C for each iteration #####
k_var_res_e <- k_var_res %>% 
  filter(T %in% c(25, 40)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(k_ta = abs(K_EaN - K_EaP),
         dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value")

k_var_res_e %>% 
  ggplot(aes(x = k_ta, y = value)) + #good, it worked
  geom_point() + 
  facet_wrap(~response_var) + 
  labs(x = "Abs diff in K_Ea", y = "Response variable value") #need to think about this more

k_var_res_plot_e <-
  k_var_res_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = scale(k_ta), y = value)) + 
  # ggplot(aes(x = scale(r_ta), y = value, colour = r_EaN < r_EaP)) +
  geom_point(size = 3) + 
  labs(x = "Scaled absolute value \nof thermal asymmetry", y = "Displacement of species pair \n(Euclidean distance) after 15C warming") + 
  coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(0, 0.5)) + 
  annotate("text", x = 0.25, y = 0.45, label = "Resource \ncarrying capacity, K", size = 6) + 
  theme_cowplot(font_size = 20)

k_var_res_plot_e2 <-
  k_var_res_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = k_ta, y = value)) + 
  geom_point(size = 3) + 
  labs(x = "Absolute value \nof thermal asymmetry", y = "Displacement of species pair \n(Euclidean distance) after 15C warming") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 0.95, y = 0.05, label = "Resource \ncarrying capacity, K", size = 6) + 
  theme_cowplot(font_size = 20)


### vary only one K_Ea ####
k_var1 <- data.frame()
for(f in 1:500){ 
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



# v_Eas vary --------------------------------------------------------------

### vary both v_Eas ####
v_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  v_var <- bind_rows(v_var, hold) 
}
beep(2)


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
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

v_var_plot <- 
  ggplot() +
  geom_path(data = v_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(v_var, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0.2, y = 0.7, label = expression("Conversion \nefficiency," ~ italic(v)), size = 6, hjust = 0.5)

## calculate euclidean distances at 20C for each iteration #####
v_var_e <- v_var %>% 
  filter(T %in% c(25, 40)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(v_ta = abs(v_EaN - v_EaP),
         dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

v_var_e %>% 
  ggplot(aes(x = v_ta, y = value, colour = v_EaN > v_EaP)) + 
  geom_point() + 
  facet_wrap(~response_var) + 
  labs(x = "Abs diff in v_Ea", y = "Response variable value") #need to think about this more

v_var_plot_e <-
  v_var_e %>% 
  filter(response_var == "dist15") %>% 
  # ggplot(aes(x = scale(v_ta), y = value)) + 
  ggplot(aes(x = scale(v_ta), y = value, colour = v_EaN < v_EaP)) +
  geom_point(size = 3) +
  labs(x = "Scaled absolute value \nof thermal asymmetry", y = "Displacement of species pair \n(Euclidean distance) after 15C warming") + 
  coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(0, 0.5)) + 
  annotate("text", x = 0, y = 0.45, label = "Conversion \nefficiency, v", size = 6) + 
  theme_cowplot(font_size = 20)

v_var_plot_e2 <-
  v_var_e %>% 
  filter(response_var == "dist15") %>% 
  # ggplot(aes(x = v_ta, y = value)) +
  ggplot(aes(x = v_ta, y = value, colour = v_EaN < v_EaP)) +
  geom_point(size = 3) + 
  labs(x = "Absolute value \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 0.5, y = 0.35, label = "Conversion \nefficiency, v", size = 6) + 
  theme_cowplot(font_size = 20)

v_var_plot_e3 <- v_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = v_EaN, y = v_EaP, colour = value)) + 
  labs(x = "Temperature dependence \nof resource 1 conversion efficiency", y = "Temperature dependence \nof resource 2 conversion efficiency", colour = "ED") + 
  annotate("text", x = 1.07, y = 0.5, label = "v", size = 7, fontface = "italic") +
  geom_point(size = 3)

### vary v_EaN only ####
v_var1 <- data.frame()
for(f in 1:500){ 
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


# m_Eas vary --------------------------------------------------------------

### vary both m_Eas ####
m_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), 
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  m_var <- bind_rows(m_var, hold) 
}
beep(2)

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

#log plot
m_var_plot <- 
  ggplot() +
  geom_path(data = m_var, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(m_var, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim=c(-1,1), xlim = c(0, 0.7)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  # labs(colour = "Degrees C \nWarming") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Consumer \nmortality rate," ~ italic(m)), size = 6)

## calculate euclidean distances at 20C for each iteration #####
m_var_e <- m_var %>% 
  filter(T %in% c(25, 40)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

m_var_e %>% 
  ggplot(aes(x = m_ta, y = value)) + 
  geom_point() + 
  facet_wrap(~response_var) + 
  labs(x = "Abs diff in m_Ea", y = "Response variable value") #need to think about this more

m_var_plot_e <-
  m_var_e %>% 
  filter(response_var == "dist15") %>% 
  # ggplot(aes(x = scale(m_ta), y = value)) + 
  ggplot(aes(x = scale(m_ta), y = value, colour = m_Ea1 < r_Ea2)) +
  geom_point(size = 3) + 
  labs(x = "Scaled absolute value \nof thermal asymmetry", y = "Displacement of species pair \n(Euclidean distance) after 15C warming") + 
  coord_cartesian(xlim = c(-1.5, 3.5), ylim = c(0, 0.5)) + 
  annotate("text", x = 0.25, y = 0.45, label = "Consumer \nmortality rate, m", size = 6) +
  theme_cowplot(font_size = 20)

m_var_plot_e2 <-
  m_var_e %>% 
  filter(response_var == "dist15") %>% 
  # ggplot(aes(x = m_ta, y = value)) +
  ggplot(aes(x = m_ta, y = value, colour = m_Ea1 < m_Ea2)) +
  geom_point(size = 3) + 
  labs(x = "Absolute value \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 0.75, y = 0.35, label = "Consumer \nmortality rate, m", size = 6) +
  theme_cowplot(font_size = 20)

m_var_plot_e3 <- m_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = m_Ea1, y = m_Ea2, colour = value)) + 
  labs(x = "Temperature dependence \nof consumer 1 mortality rate", y = "Temperature dependence \nof consumer 2 mortality rate", colour = "ED") +
  annotate("text", x = 1.07, y = 0.5, label = "m", size = 7, fontface = "italic") + 
  geom_point(size = 3)

### vary m_Ea1 only ####
m_var1 <- data.frame()
for(f in 1:500){ 
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


# Combining all param_var and euclidean displacement plots -------------------------------------------

param_var_plot <- c_var_plot + r_var_plot + k_var_plot + v_var_plot + m_var_plot + rvar_legend + 
  plot_annotation(tag_levels = "A")
# ggsave(plot = param_var_plot, filename = "figures/kd-figs/param_var_plots.pdf", width = 14, height = 10)

param_e_plot <- r_var_plot_e + c_var_plot_e + v_var_plot_e + k_var_plot_e + m_var_plot_e
# ggsave(plot = param_e_plot, filename = "figures/kd-figs/param_e_plots.pdf", width = 16, height = 12)

param_e_unscaled_plot <- c_var_plot_e2 + r_var_plot_e2 + k_var_plot_e2 + v_var_plot_e2 + m_var_plot_e2
# ggsave(plot = param_e_unscaled_plot, filename = "figures/kd-figs/param_e_plots_unscaled.pdf", width = 16, height = 12)
# ggsave(plot = param_e_unscaled_plot, filename = "figures/kd-figs/param_e_plots_unscaled_inequality.pdf", width = 18, height = 12)

param_e_plot3 <- r_var_plot_e3 + c_var_plot_e3 + v_var_plot_e3 + k_var_plot_e3 + m_var_plot_e3
ggsave(plot = param_e_plot3, filename = "figures/kd-figs/param_e_plots2.pdf", width = 18, height = 12)