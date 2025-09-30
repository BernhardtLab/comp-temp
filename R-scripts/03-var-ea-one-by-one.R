# This script is to explore the effects of temperature sensitivity of different MacArthur consumer-resource parameters on changes in niche and fitness differences. In this script, each MacArthur consumer-resource parameter is given by an Arrhenius function, with a temperature sensitivity (activation energy, slope) term and an intercept term, which determines the value of the function at ambient temperatures (Tref, ref temp). In each simulation, temperature sensitivities are defined as "{shorthand-parameter_EAik}", where ik captures the relevant consumer, resource, or both, and intercepts are defined as "{shorthand-parameter_b}". Consumers are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively, throughout the script. The script simulates the effects of warming when each parameter is given a temperature sensitivity, randomly drawn from the parameter's empirical distribution (generated in 01-param-dists), while all other parameters of the model are assigned a temperature sensitivity of 0. Script 04-full-temp-var-analysis simulates warming when all parameters have temperature sensitivities drawn from their empirical distributions simultaneously. 

# script DOB: April 1, 2025
# author: Kaleigh Davis, University of Guelph postdoc

#### packages and referencing #####
# load necessary pkgs
library(tidyverse)
library(janitor)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)
library(viridis)
library(beepr)

# get referencing set up for MacArthur temp sensitivity function
source("R-scripts/02-temp-dep-macarthur.R") #this script contains the MacArthur translation function, with all parameters flexibly defined in the function for assigning at time of use, and the Arrhenius function.

# load in distributions for parameter values
# these are distributions generated from empirical data using MCMC regression, in 01-param-dists.R
param_vals <- read_csv(file = "data/processed-data/param_post_dists.csv")

# split these into dfs for each parameter
param_vals %>%
  mutate(parameter = str_replace(parameter, "resource_growth_rate", "rgr"),
         parameter = str_replace(parameter, "carrying_capacity", "k"),
         parameter = str_replace(parameter, "conversion_efficiency", "v"),
         parameter = str_replace(parameter, "mortality_rate", "m"),
         parameter = str_replace(parameter, "consumption rate", "c")) %>% 
  group_by(parameter) %>%
  group_split() %>%
  set_names(unique(param_vals$parameter)) %>% 
  walk(~ assign(paste0(.x$parameter[1], "_post_dist"), .x, envir = .GlobalEnv))

# get summary stats for all parameters ########
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

################################################################################
##########################    MAIN ANALYSIS   ##################################
################################################################################

# Simulation setup for the analysis in the main text: here all param EAs (temperature sensitivities) are drawn from their estimated empirical distributions, consumers have reciprocal resource use (i.e. each consumer specializes on a different resource and has equal strength preference for its preferred resource), and resource N (== a) grows faster than resource P (== b) at the ambient temperature (== ref temp). The species pair starts on the boundary of coexistence.

####################### r_Ea varies ----------------------
r_var <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all r_EAs from empirical distributions
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b = 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  r_var <- bind_rows(r_var, hold) 
}

#plot - for fig 3
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C \nWarming") +
  annotate("text", x = 0.75, y = 0.05, label = expression("Resource \ngrowth rate,"~italic(r)[italic(k)]), size = 6) +
  theme_cowplot(font_size = 20) +
  theme(legend.position = "none")

# get legend for composite plot -- need to generate r_var_plot WITH the legend first, then add in the legend.position = "none" for the actual composite plot
rvar_legend  <- get_legend(r_var_plot)

# calculate euclidean distances at 25C for each iteration
r_var_e <- r_var %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(r_var_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#plot effect of thermal asymmetry on distance and NFD - for figure S4
r3p <- r_var_e %>% 
  ggplot(aes(x = r_ta, y = value)) +
  geom_point() + 
  facet_wrap(~response_var,
             labeller = labeller(response_var = c("dist15" = "Eucl Dist", "shift_fitrat" = "Change FD", "shift_nichediffs" = "Change ND"))) +
  labs(x = "Magnitude \n of thermal asymmetry", y = "Value") + 
  ggtitle("Resource growth rate, r_k") + 
  theme(legend.position = "none")

#unscaled TA - euclidean distance plot - for figure 4
r_var_plot_e2 <-
  r_var_e %>% 
    filter(response_var == "dist15") %>% 
    ggplot(aes(x = abs_r_ta, y = value, colour = r_EaP > r_EaN)) +
    geom_point(size = 3) + 
    labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
    coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
    annotate("text", x = 0.88, y = 0.35, label = "Resource \ngrowth rate, r", size = 5.5) + 
    theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[ra] * "- E"[rb]*"|"))

############### c_Ea varies -------------------------------------------------------------
#vary all four c_Eas
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
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  c_var <- bind_rows(c_var, hold) 
}

#log plot - for fig 3
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.6, y = 0.05, label = expression("Consumption rate," ~ italic(c)[italic(ik)]), size = 6)

# vary just c's for one consumer for thermal asymmetry plot
c_var_ta <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
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
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  c_var_ta <- bind_rows(c_var_ta, hold) 
}

#log plot
c_var_ta_plot <-
  ggplot() +
  geom_path(data = c_var_ta, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none")

# calculate euclidean distances at 25C for each iteration
c_var_e <- c_var_ta %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(abs_c_ta = abs(c_Ea1P - c_Ea1N),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(c_var_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#effects on euclidean distance (fig 4) and niche and fitness diffs (fig s4)
#for figure S4
c3p <-
  c_var_e %>% 
  ggplot() + 
  geom_point(aes(x = abs_c_ta, y = value)) + 
  facet_wrap(~response_var) + 
  labs(x = "Manitude of thermal asymmetries (aggregated)", y = "Response variable value") +
  facet_wrap(~response_var,
             labeller = labeller(response_var = c("dist15" = "Eucl Dist", "shift_fitrat" = "Change FD", "shift_nichediffs" = "Change ND")))+ 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Value") + 
  ggtitle("Consumption rate, c_ik")

#for figure 4
c_var_plot_e2 <-
  c_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_c_ta, y = value)) +
  geom_point(aes(colour = c_Ea1P > c_Ea1N)) + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) +
  annotate("text", x = 0.4, y = 0.35, label = "Consumption rate, c", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[c1a] * "- E"[c1b]*"|"))

############## K_Eas vary  --------------------------------------------------------------
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
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b = 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  k_var <- bind_rows(k_var, hold) 
}

#log plot - for fig 3
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.65, y = 0.05, label = expression("Resource \ncarrying capacity," ~ italic(K)[italic(k)]), size = 6)

# calculate euclidean distances at 25C for each iteration 
k_var_e <- k_var %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(k_ta = K_EaN - K_EaP,
         abs_k_ta = abs(K_EaN - K_EaP),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(k_var_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#figure s4
k3p <- 
  k_var_e %>% 
  ggplot(aes(x = k_ta, y = value)) + 
  geom_point() + 
  facet_wrap(~response_var,
             labeller = labeller(response_var = c("dist15" = "Eucl Dist", "shift_fitrat" = "Change FD", "shift_nichediffs" = "Change ND"))) +
  labs(x = "Magnitude \nof thermal asymmetry", y = "Value") + 
  ggtitle("Resource carrying capacity, K_k")

#figure 4
k_var_plot_e2 <-
  k_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_k_ta, y = value, colour = K_EaP > K_EaN)) +
  geom_point(size = 3) + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 0.7, y = 0.35, label = "Resource \ncarrying capacity, K", size = 6) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[Ka] * "- E"[Kb]*"|"))

# v_Eas vary --------------------------------------------------
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
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000,
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  v_var <- bind_rows(v_var, hold) 
}

#plot - for fig 3
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0.7, y = 0.05, label = expression("Conversion \nefficiency," ~ italic(v)[italic(ik)]), size = 6, hjust = 0.5)

# calculate euclidean distances at 25C for each iteration
v_var_e <- v_var %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(v_ta = v_EaN - v_EaP,
         abs_v_ta = abs(v_EaN - v_EaP),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 


#anywhere where ND and FD are exactly the same at end? No, great.
nrow(v_var_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#figure s4
v3p <- 
  v_var_e %>% 
  ggplot(aes(x = v_ta, y = value)) + 
  geom_point() + 
  facet_wrap(~response_var,
             labeller = labeller(response_var = c("dist15" = "Eucl Dist", "shift_fitrat" = "Change FD", "shift_nichediffs" = "Change ND"))) +
  labs(x = "Magnitude \nof thermal asymmetry", y = "Value") + 
  ggtitle("Conversion efficiency, v_k")

#figure 4
v_var_plot_e2 <-
  v_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_v_ta, y = value, colour = v_EaP > v_EaN)) +
  geom_point(size = 3) + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 0.5, y = 0.35, label = "Conversion \nefficiency, v", size = 6) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[va] * "- E"[vb]*"|"))

# m_Eas vary  --------------------------------------------------------------
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
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000,
                      v1N_b = 0.5, v1P_b = 1,
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  m_var <- bind_rows(m_var, hold) 
}

# plot - for fig 3
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.7, y = 0.05, label = expression("Consumer \nmortality rate," ~ italic(m)[italic(i)]), size = 6)

# calculate euclidean distances at 25C for each iteration
m_var_e <- m_var %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(m_ta = m_Ea1 - m_Ea2,
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(m_var_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#figure s4
m3p <-
  m_var_e %>% 
  ggplot(aes(x = m_ta, y = value)) + 
  geom_point() + 
  facet_wrap(~response_var,
             labeller = labeller(response_var = c("dist15" = "Eucl Dist", "shift_fitrat" = "Change FD", "shift_nichediffs" = "Change ND"))) +
  labs(x = "Magnitude of thermal asymmetry", y = "Value") + 
  ggtitle("Mortality, m_i")

#figure 4
m_var_plot_e2 <-
  m_var_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_m_ta, y = value, colour = m_Ea2 > m_Ea1)) +
  geom_point(size = 3) + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  annotate("text", x = 0.75, y = 0.35, label = "Consumer \nmortality rate, m", size = 6) +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[m1] * "- E"[m2]*"|"))

# Combine all param_var and euclidean displacement plots -------------------------------------------

# Figure 3 - effect trajectories --------------------------------------
param_var_plot <- c_var_plot + r_var_plot + k_var_plot + v_var_plot + m_var_plot + rvar_legend + 
  plot_annotation(tag_levels = "A")
# ggsave(plot = param_var_plot, filename = "figures/param_var_plots_EA0.pdf", width = 14, height = 10)

# Figure 4 - effect size of thermal asymmetries and effect size --------------------------------------
param_e_unscaled_plot <- c_var_plot_e2 + r_var_plot_e2 + k_var_plot_e2 + v_var_plot_e2 + m_var_plot_e2 + plot_annotation(tag_levels = "A")
# ggsave(plot = param_e_unscaled_plot, filename = "figures/param_e_plots_unscaled_inequality_EA0.pdf", width = 18, height = 12)

# Figure S4 - three panel NFD plots -------------------------------------------
threeps <- r3p + c3p + v3p + k3p + m3p
# ggsave(plot = threeps, filename = "figures/param_three_panels.pdf", width = 16, height = 12)

#############################################################################
############################# SUPPLEMENTAL ANALYSES #########################
#############################################################################

# vary different combinations of consumption rates - Figure S2----------------------------

#####  vary all c_EaPs ####
c_var1 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2N = 0,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0,
                      m_Ea2 = 0,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000,
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5,
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") 

### vary all c_EaNs  ####
c_var5 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P =  0,
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0,
                      m_Ea2 = 0,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000,
                      v1N_b = 0.5, v1P_b = 1,
                      v2N_b = 1, v2P_b = 0.5,
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  c_var5 <- bind_rows(c_var5, hold) 
}

#log plot
c_var5_plot <-
  ggplot() +
  geom_path(data = c_var5, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var5, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none")

### vary c_EAs for each species' preferred resource ####
c_var4 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = 0,
                      r_EaP = 0,
                      c_Ea1N = 0,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = 0,
                      K_EaN = 0,
                      K_EaP = 0,
                      v_EaN = 0,
                      v_EaP = 0,
                      m_Ea1 = 0,
                      m_Ea2 = 0,
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5,
                      r_N_b = 1, r_P_b = 0.5,
                      K_N_b= 2000, K_P_b = 2000,
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  c_var4 <- bind_rows(c_var4, hold) 
}

#log plot
c_var4_plot <-
  ggplot() +
  geom_path(data = c_var4, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(c_var4, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 5) +
  geom_hline(yintercept = 0, linetype = 5) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") 

# combined plot for different combinations of C draws
c_0 <- c_var4_plot + c_var5_plot + c_var1_plot + c_var_ta_plot  + c_var_plot + rvar_legend + plot_annotation(tag_levels = "A")

# ggsave(plot = c_0, filename = "figures/c_drawtypes_EAs0.pdf", width = 16, height = 10)
#in order, these have focal params: consumption rate of preferred resource only, consumption rate of N only, consumption rate of P only, consumption rates of consumer 1 N & P, all four consumption rates
