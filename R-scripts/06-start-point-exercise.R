# This script is to reproduce the main supplementary figures of the paper -- with different start points. We achieve these different start points by manipulating the resource preferences (analogous to niche widths) of each species at the models starting conditions, i.e. ambient temperature. Scenario 1 is given in the main text (script 04) and features two consumers that each specialize on one of two resources and they have equally strong preference for this resource. The two resources have uneven growth rates under ambient conditions, which places the species pair on the boundary of coexistence under ambient conditions. In scenario 2, consumer preferences are the same as in Scenario 1, but the resources grow at the same rate under ambient conditions, which moves the start point from the coexistence boundary to the middle of the coexistence region. In scenario 3, species have very different resource preferences, where one consumer has a strong preference for its preferred resource and the other consumer has a weak preference for its preferred resources. Both resources grow at the same rate under ambient conditions. Parameters defining each of these starting conditions are given in each simulation below, and in summary in Table S1. In each simulation, each MacArthur consumer-resource parameter is given by an Arrhenius function, with a temperature sensitivity (activation energy, slope) term and an intercept term, which determines the value of the function at ambient temperatures (Tref, ref temp). In each simulation, temperature sensitivities are defined as "{parameter_EAik}", where ik captures the relevant consumer, resource, or both, and intercepts are defined as "{parameter-ik_b}". Consumers are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively, throughout the script. The script simulates the effects of warming when each parameter is given a temperature sensitivity, randomly drawn from the parameter's empirical distribution (generated in 01-param-dists), simultaneously.

#script DOB: 9/5/2025
#author: Kaleigh Davis, Postdoc UoG with Joey Bernhardt

# sourced scripts:
#   R-scripts/02-temp-dep-macarthur.R — temperature-dependent MacArthur
#     consumer-resource function
#   R-scripts/03-arrhenius.R — Arrhenius temperature-dependence function

# inputs:
#   data/processed-data/param_post_dists.csv — posterior distributions of
#     temperature sensitivity for each model parameter (generated in script 01)

# outputs:
#   figures/S2-supp_startpoint_errce.pdf — supplementary figures S2
#   figures/S3-supp_startpoint_urrce.png — supplementary figure S3
#   figures/S4-extra-start-points.png — additional start point comparisons; Figure S4

#### packages and referencing #####
#load necessary pkgs
library(tidyverse)
library(janitor)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)
library(purrr)
library(viridis)
library(beepr)
library(car)

# get referencing set up for MacArthur temp dependence function
source("R-scripts/02-temp-dep-macarthur.R") #this calls the MacArthur translation function, with all parameters flexibly defined in the function for assigning at time of use
source("R-scripts/03-arrhenius.R") #this calls the arrhenius function.

#load in distributions for parameter values.
# these are continuous distributions generated from empirical data using MCMC regression, in 01-param-dists.R
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
        Median = median,
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

### rrc, equal growth rates -- scenario 2 - Fig S2 ####
#### ONE-BY-ONE ####
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences"))) + 
  theme_cowplot(font_size = 26) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.55, label = expression("Consumption \nrate," ~ italic(c)[italic(ik)]), size = 8)


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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences"))) + 
  theme_cowplot(font_size = 26) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.6, label = expression("Resource \ngrowth rate," ~ italic(r)[italic(k)]), size = 8)


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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences"))) + 
  theme_cowplot(font_size = 26) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.45, label = expression("Resource \ncarrying \ncapacity," ~ italic(K)[italic(k)]), size = 8)

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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences"))) + 
  theme_cowplot(font_size = 26) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.6, label = expression("Conversion \nefficiency," ~ italic(v)[italic(ik)]), size = 8)

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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences"))) + 
  theme_cowplot(font_size = 24) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.6, label = expression("Consumer \nmortality rate," ~ italic(m)[italic(i)]), size = 8)

# rrc equal base growth plot ###
rrce_plots <- c_var_plot + r_var_plot + k_var_plot + v_var_plot + m_var_plot

#### POMPOM; rrc, equal growth rates ####
rrce_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), 
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

#get average change in position after 15C warming
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
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
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
  geom_point(data = rrce_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrce_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrce_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 9, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 9, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 9, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 9, fontface = 2) +
  theme_cowplot(font_size = 28)

rrce_all %>% 
  group_by(T) %>% 
  summarise(median_ND = median(new_stabil_potential),
            median_FD = median(new_fit_ratio)) %>% 
  pivot_longer(cols = c(median_ND, median_FD), names_to = "NFD", values_to = "median") %>% 
  ggplot() +
  geom_point(aes(x = T, y = median, colour = NFD)) +
  facet_wrap(~NFD, scales = "free_y")

# rrc equal start two plots -- Figure S2
startpoint1 <- rrce_plots / log_pom_rrce_all + plot_annotation(tag_levels = "A")
# ggsave(plot = startpoint1, filename = "figures/S2-supp_startpoint_errce.pdf", height = 20, width = 26)

######### uneven reciprocal preference, equal growth rates -- scenario 3 - Fig S3 ########

#### ONE-BY-ONE; ####
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
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b = 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.9, 
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
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
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.9, 
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
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
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.9, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  k_var1 <- bind_rows(k_var1, hold) 
}


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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
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
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.9, 
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.25, y = 0.7, label = expression("Conversion \nefficiency," ~ italic(v)[italic(ik)]), size = 6)

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
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.9, 
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
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = 0.3, y = 0.7, label = expression("Consumer \nmortality rate," ~ italic(m)[italic(i)]), size = 6)

# uneven reciprocal, equal base growth plot ###
urrce_plots <- c_var1_plot + r_var1_plot + k_var1_plot + v_var1_plot + m_var1_plot

#### POMPOM; uneven reciprocal preference, equal growth rates ####
urrce_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 1),
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
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.9, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce_all <- bind_rows(urrce_all, hold) 
}

#get average change in position after 15C warming
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
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
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
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 9, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 9, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 9, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 9, fontface = 2) +
  theme_cowplot(font_size = 26)

urrce_parmT <- urrce_all %>% 
  group_by(T) %>% 
  summarise(ND = median(new_stabil_potential),
            FD = median(new_fit_ratio),
            a11 = median(a11),
            a12 = median(a12),
            a21 = median(a21),
            a22 = median(a22),
            KN = median(KN),
            KN_over_rN = median(KN/rN),
            KP = median(KP),
            KP_over_rP = median(KP/rP)) %>% 
  pivot_longer(cols = c(ND:KP_over_rP), names_to = "parameter", values_to = "median") %>% 
  ggplot() +
  geom_point(aes(x = T, y = median)) + 
  facet_wrap(~parameter, scales = "free_y")

# uneven reciprocal preference, equal start two plots
startpoint2 <- urrce_plots / log_pom_urrce_all + plot_annotation(tag_levels = "A")
# ggsave(plot = startpoint2, filename = "figures/S3-supp_startpoint_urrce.png", height = 20, width = 26, bg = "white")

#repeat pompom plot, but without annotations
log_pom_urrce_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
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
  theme_cowplot(font_size = 30) +
  ggtitle("Full model")

#for how many simulations does ND decrease and does FD decrease
urrce_start <- urrce_all %>% 
  filter(T == 10) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_sd_stab_pot = sd(new_stabil_potential),
            new_sd_fit_rat = sd(new_fit_ratio)) %>% 
  dplyr::select(new_mean_stab_pot, new_mean_fit_rat)

urrce_shifts <- urrce_all %>% 
  dplyr::select(iteration, T, new_stabil_potential, new_fit_ratio) %>% 
  mutate(stab_pot_start = urrce_start$new_mean_stab_pot,
         fit_rat_start = urrce_start$new_mean_fit_rat) %>% 
  mutate(fd_shift = ifelse(new_fit_ratio < fit_rat_start, "decrease", 
                           ifelse(new_fit_ratio > fit_rat_start, "increase",
                                  ifelse(new_fit_ratio == fit_rat_start, "no change", "potato"))),
         nd_shift = ifelse(new_stabil_potential < stab_pot_start, "decrease",
                           ifelse(new_stabil_potential > stab_pot_start, "increase",
                                  ifelse(new_stabil_potential == stab_pot_start, "no change","potato")))) %>% 
  filter(T == 10 | T == 25)

urrce_shifts %>% 
  group_by(T, fd_shift, nd_shift) %>% 
  tally()

#FD decrease in 253/500 cases; 2nd sim: 254/500; 3rd:279/500
#ND decrease in 250/500 cases; 2nd sim: 255/500; 3rd: 263/500

# calculate euclidean distances at 25C for each iteration
urrce_all_e <- urrce_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(urrce_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#which thermal asymmetries drive large shifts in competition?
urrce_all_ed <- urrce_all_e %>% 
  filter(response_var == "dist15")

#multiple regression to test for effects of TA on ED
model <- lm(value ~ abs_r_ta + abs_c1_ta + abs_c2_ta + abs_c3_ta + abs_c4_ta + abs_c5_ta + abs_c6_ta + abs_k_ta + abs_v_ta + abs_m_ta, data = urrce_all_ed)

# #visreg::visreg(model)
summary(model) #here for effect estimates
Anova(model, type = "II") #here for p-values

#### Simulate from different start points -- Figure S4 ##########
#### POMPOM 1; uneven reciprocal preference, unequal growth rates ####
urrce1_all <- data.frame()
for(f in 1:500){ 
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
                      c1N_b = 0.6, c1P_b = 0.9, 
                      c2N_b = 0.5, c2P_b = 0.2, 
                      r_N_b = 0.7, r_P_b = 0.2, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.3, v1P_b = 0.9, 
                      v2N_b = 0.3, v2P_b = 0.1, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce1_all <- bind_rows(urrce1_all, hold) 
}

#get average change in position after 15C warming
urrce1_all_avg_new <- urrce1_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce1_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce1_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce1_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce1_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce1_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce1_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce1_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = -0.8, 
           label = "c_1a = 0.6, c_1b = 0.9, 
c_2a = 0.5, c_2b = 0.2, 
r_a = 0.7, r_b = 0.2, 
K_a= 1000, K_b = 1000, 
v_1a = 0.3, v_1b = 0.9, 
v_2a = 0.3, v_2b = 0.1, 
m_1 = 0.01, m_2 = 0.01", 
           size = 9, fontface = 1, hjust = 0) 

#make the pompom again, but without annotations
#pompom
log_pom_urrce1_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce1_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce1_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce1_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce1_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce1_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce1_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 36) + 
  theme(legend.position = "none") 

# calculate euclidean distances at 25C for each iteration
urrce1_all_e <- urrce1_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end?
nrow(urrce1_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential)) #No, great.

urrce1_all_ed <- urrce1_all_e %>% 
  filter(response_var == "dist15")

model1 <- lm(value ~ abs_r_ta + abs_c1_ta + abs_c2_ta + abs_c3_ta + abs_c4_ta + abs_c5_ta + abs_c6_ta + abs_k_ta + abs_v_ta + abs_m_ta, data = urrce1_all_ed)

#visreg::visreg(model1)
summary(model1)
Anova(model1, type = "II")

#### POMPOM 2; uneven reciprocal preference, unequal growth rates ####
urrce2_all <- data.frame()
for(f in 1:500){ 
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
                      c1N_b = 0.5, c1P_b = 0.2, 
                      c2N_b = 0.9, c2P_b = 0.9, 
                      r_N_b = 0.2, r_P_b = 0.7, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.1, v1P_b = 0.6, 
                      v2N_b = 0.4, v2P_b = 0.7, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce2_all <- bind_rows(urrce2_all, hold) 
}

#get average change in position after 15C warming
urrce2_all_avg_new <- urrce2_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce2_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce2_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce2_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce2_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce2_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce2_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce2_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = 0.8, 
           label = "c_1a = 0.5, c_1b = 0.2, 
c_2a = 0.9, c_2b = 0.9, 
r_a = 0.2, r_b = 0.7, 
K_a= 1000, K_b = 1000, 
v_1a = 0.1, v_1b = 0.6, 
v_2a = 0.4, v_2b = 0.7, 
m_1 = 0.01, m_2 = 0.01", 
           size = 9, fontface = 1, hjust = 0) 

log_pom_urrce2_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce2_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce2_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce2_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce2_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce2_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce2_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 36) + 
  theme(legend.position = "none")

# calculate euclidean distances at 25C for each iteration
urrce2_all_e <- urrce2_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(urrce2_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#which thermal asymmetries drive large shifts in competition?
urrce2_all_ed <- urrce2_all_e %>% 
  filter(response_var == "dist15")

model2 <- lm(value ~ abs_r_ta + abs_c1_ta + abs_c2_ta + abs_c3_ta + abs_c4_ta + abs_c5_ta + abs_c6_ta + abs_k_ta + abs_v_ta + abs_m_ta, data = urrce2_all_ed)

#visreg::visreg(model2)
summary(model2)
anova(model2)

#### POMPOM 3; uneven reciprocal preference, unequal growth rates ####
urrce3_all <- data.frame()
for(f in 1:500){ 
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
                      c1N_b = 1.3, c1P_b = 0.9, 
                      c2N_b = 0.9, c2P_b = 0.2, 
                      r_N_b = 0.7, r_P_b = 2, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.9, 
                      v2N_b = 0.9, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce3_all <- bind_rows(urrce3_all, hold) 
}

#get average change in position after 15C warming
urrce3_all_avg_new <- urrce3_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce3_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce3_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce3_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce3_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce3_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce3_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce3_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = 0.8, 
           label = "c_1a = 1.3, c_1b = 0.9, 
c_2a = 0.9, c_2b = 0.2, 
r_a = 0.7, r_b = 2, 
K_a= 1000, K_b = 1000, 
v_1a = 0.9, v_1b = 0.9, 
v_2a = 0.9, v_2b = 0.9, 
m_1 = 0.1, m_2 = 0.1", 
           size = 9, fontface = 1, hjust = 0) 

#make the pompom again, but without annotations
#pompom
log_pom_urrce3_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce3_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce3_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce3_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce3_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce3_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce3_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 36) + 
  theme(legend.position = "none") 

# calculate euclidean distances at 25C for each iteration
urrce3_all_e <- urrce3_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end?
nrow(urrce3_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential)) #No, great.

#which thermal asymmetries drive large shifts in competition?
urrce3_all_ed <- urrce3_all_e %>% 
  filter(response_var == "dist15")

model3 <- lm(value ~ abs_r_ta + abs_c1_ta + abs_c2_ta + abs_c3_ta + abs_c4_ta + abs_c5_ta + abs_c6_ta + abs_k_ta + abs_v_ta + abs_m_ta, data = urrce3_all_ed)

#visreg::visreg(model3)
summary(model3)
Anova(model3, type = "II")

#### POMPOM 4; uneven reciprocal preference, unequal growth rates ####
urrce4_all <- data.frame()
for(f in 1:500){ 
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
                      c1N_b = 1.3, c1P_b = 0.4, 
                      c2N_b = 0.9, c2P_b = 0.4, 
                      r_N_b = 0.7, r_P_b = 1, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.4, 
                      v2N_b = 0.9, v2P_b = 0.4, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce4_all <- bind_rows(urrce4_all, hold) 
}

#get average change in position after 15C warming
urrce4_all_avg_new <- urrce4_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce4_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce4_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce4_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce4_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce4_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce4_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce4_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = -0.8, 
           label = "c_1a = 1.3, c_1b = 0.4, 
c_2a = 0.9, c_2b = 0.4, 
r_a = 0.7, r_b = 1, 
K_a= 1000, K_b = 1000, 
v_1a = 0.9, v_1b = 0.4, 
v_2a = 0.9, v_2b = 0.4, 
m_1 = 0.1, m_2 = 0.1", 
           size = 9, fontface = 1, hjust = 0) 

log_pom_urrce4_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce4_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce4_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce4_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce4_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce4_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce4_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 36) + 
  theme(legend.position = "none") 

# calculate euclidean distances at 25C for each iteration
urrce4_all_e <- urrce4_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(urrce4_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#which thermal asymmetries drive large shifts in competition?
urrce4_all_ed <- urrce4_all_e %>% 
  filter(response_var == "dist15")

model4 <- lm(value ~ abs_r_ta + abs_c1_ta + abs_c2_ta + abs_c3_ta + abs_c4_ta + abs_c5_ta + abs_c6_ta + abs_k_ta + abs_v_ta + abs_m_ta, data = urrce4_all_ed)

#visreg::visreg(model4)
summary(model4)
anova(model4)

#### POMPOM 5; uneven reciprocal preference, unequal growth rates ####
urrce5_all <- data.frame()
for(f in 1:500){ 
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
                      c1N_b = 1.3, c1P_b = 0.8, 
                      c2N_b = 0.3, c2P_b = 2, 
                      r_N_b = 0.8, r_P_b = 0.8, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.8, 
                      v2N_b = 0.1, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce5_all <- bind_rows(urrce5_all, hold) 
}

#get average change in position after 15C warming
urrce5_all_avg_new <- urrce5_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce5_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce5_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce5_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce5_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce5_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce5_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce5_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = 0.8, 
           label = "c_1a = 1.3, c_1b = 0.8, 
c_2a = 0.3, c_2b = 2, 
r_a = 0.8, r_b = 0.8, 
K_a= 1000, K_b = 1000, 
v_1a = 0.9, v_1b = 0.8, 
v_2a = 0.1, v_2b = 0.9, 
m_1 = 0.1, m_2 = 0.1", 
           size = 9, fontface = 1, hjust = 0)

log_pom_urrce5_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce5_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce5_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce5_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce5_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce5_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce5_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  theme_cowplot(font_size = 36) + 
  theme(legend.position = "none") 

# calculate euclidean distances at 25C for each iteration
urrce5_all_e <- urrce5_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(urrce5_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#which thermal asymmetries drive large shifts in competition?
urrce5_all_ed <- urrce5_all_e %>% 
  filter(response_var == "dist15")

model5 <- lm(value ~ abs_r_ta + abs_c1_ta + abs_c2_ta + abs_c3_ta + abs_c4_ta + abs_c5_ta + abs_c6_ta + abs_k_ta + abs_v_ta + abs_m_ta, data = urrce5_all_ed)

#visreg::visreg(model5)
summary(model5)


#### POMPOM 6; uneven reciprocal preference, unequal growth rates ####
urrce6_all <- data.frame()
for(f in 1:500){ 
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
                      c1N_b = 1.3, c1P_b = 0.6, 
                      c2N_b = 0.3, c2P_b = 2, 
                      r_N_b = 0.8, r_P_b = 1.5, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.8, 
                      v2N_b = 0.1, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce6_all <- bind_rows(urrce6_all, hold) 
}

#get average change in position after 15C warming
urrce6_all_avg_new <- urrce6_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce6_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce6_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce6_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce6_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce6_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce6_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce6_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 26) + 
  annotate("text", x = 0, y = -0.8, 
           label = "c_1a = 1.3, c_1b = 0.6,
           c_2a = 0.3, c_2b = 2,
           r_a = 0.8, r_b = 1.5,
           K_a = 1000, K_b = 1000,
           v_1a = 0.9, v_1b = 0.8,
           v_2a = 0.1, v_2b = 0.9,
           m_1 = 0.1, m_2 = 0.1", 
           size = 9, fontface = 1, hjust = 0) 

log_pom_urrce6_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce6_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce6_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce6_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce6_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce6_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce6_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 36)

# calculate euclidean distances at 25C for each iteration
urrce6_all_e <- urrce6_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(urrce6_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#which thermal asymmetries drive large shifts in competition?
urrce6_all_ed <- urrce6_all_e %>% 
  filter(response_var == "dist15")

model6 <- lm(value ~ abs_r_ta + abs_c1_ta + abs_c2_ta + abs_c3_ta + abs_c4_ta + abs_c5_ta + abs_c6_ta + abs_k_ta + abs_v_ta + abs_m_ta, data = urrce6_all_ed)

#visreg::visreg(model6)
summary(model6)
Anova(model6, type = "II")

### make supplemental figure -- pompoms from many start points #####
diff_start_points <- (log_pom_urrce2_all + log_pom_urrce3_all +  log_pom_urrce5_all) / (log_pom_urrce4_all + log_pom_urrce1_all + log_pom_urrce6_all) + plot_annotation(tag_levels = "A")

diff_start_points_noanno <- (log_pom_urrce2_all_noanno + log_pom_urrce3_all_noanno +  log_pom_urrce5_all_noanno) / (log_pom_urrce4_all_noanno + log_pom_urrce1_all_noanno + log_pom_urrce6_all_noanno) + plot_annotation(tag_levels = "A")

# ggsave(plot = diff_start_points_noanno, filename = "figures/S4-extra-start-points.png", height = 24, width = 30, device = "png")

### No thermal asymmetries --> no shift in competition from diverse start points - reviewer response #########
######## no TAs ##################
# pom base - urrce -------
urrce_nota_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
                      ref_temp = 10,
                      r_EaN = mean(rgr_post_dist$intercept),
                      r_EaP = mean(rgr_post_dist$intercept), 
                      c_Ea1N = mean(c_post_dist$intercept),
                      c_Ea1P = mean(c_post_dist$intercept), 
                      c_Ea2N = mean(c_post_dist$intercept),
                      c_Ea2P = mean(c_post_dist$intercept), 
                      K_EaN = mean(k_post_dist$intercept), 
                      K_EaP = mean(k_post_dist$intercept), 
                      v_EaN = mean(v_post_dist$intercept),
                      v_EaP = mean(v_post_dist$intercept), 
                      m_Ea1 = mean(m_post_dist$intercept), 
                      m_Ea2 = mean(m_post_dist$intercept),
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.9, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce_nota_all <- bind_rows(urrce_nota_all, hold) 
}

#get average change in position after 15C warming
urrce_nota_all_avg_new <- urrce_nota_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce_nota_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce_nota_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 12) +
  geom_point(data = filter(urrce_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 11) +
  # position after 15C warming
  geom_point(data = urrce_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce_nota_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  theme_cowplot(font_size = 26)

#### POMPOM 1; uneven reciprocal preference, unequal growth rates ####
urrce1_nota_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
                      ref_temp = 10,
                      r_EaN = mean(rgr_post_dist$intercept),
                      r_EaP = mean(rgr_post_dist$intercept), 
                      c_Ea1N = mean(c_post_dist$intercept),
                      c_Ea1P = mean(c_post_dist$intercept), 
                      c_Ea2N = mean(c_post_dist$intercept),
                      c_Ea2P = mean(c_post_dist$intercept), 
                      K_EaN = mean(k_post_dist$intercept), 
                      K_EaP = mean(k_post_dist$intercept), 
                      v_EaN = mean(v_post_dist$intercept),
                      v_EaP = mean(v_post_dist$intercept), 
                      m_Ea1 = mean(m_post_dist$intercept), 
                      m_Ea2 = mean(m_post_dist$intercept),
                      c1N_b = 0.6, c1P_b = 0.9, 
                      c2N_b = 0.5, c2P_b = 0.2, 
                      r_N_b = 0.7, r_P_b = 0.2, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.3, v1P_b = 0.9, 
                      v2N_b = 0.3, v2P_b = 0.1, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce1_nota_all <- bind_rows(urrce1_nota_all, hold) 
}

#get average change in position after 15C warming
urrce1_nota_all_avg_new <- urrce1_nota_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce1_nota_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce1_nota_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce1_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 12) +
  geom_point(data = filter(urrce1_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 11) +
  # position after 15C warming
  geom_point(data = urrce1_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce1_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce1_nota_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none")

#### POMPOM 2; uneven reciprocal preference, unequal growth rates ####
urrce2_nota_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
                      ref_temp = 10,
                      r_EaN = mean(rgr_post_dist$intercept),
                      r_EaP = mean(rgr_post_dist$intercept), 
                      c_Ea1N = mean(c_post_dist$intercept),
                      c_Ea1P = mean(c_post_dist$intercept), 
                      c_Ea2N = mean(c_post_dist$intercept),
                      c_Ea2P = mean(c_post_dist$intercept), 
                      K_EaN = mean(k_post_dist$intercept), 
                      K_EaP = mean(k_post_dist$intercept), 
                      v_EaN = mean(v_post_dist$intercept),
                      v_EaP = mean(v_post_dist$intercept), 
                      m_Ea1 = mean(m_post_dist$intercept), 
                      m_Ea2 = mean(m_post_dist$intercept),
                      c1N_b = 0.5, c1P_b = 0.2, 
                      c2N_b = 0.9, c2P_b = 0.9, 
                      r_N_b = 0.2, r_P_b = 0.7, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.1, v1P_b = 0.6, 
                      v2N_b = 0.4, v2P_b = 0.7, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce2_nota_all <- bind_rows(urrce2_nota_all, hold) 
}

#get average change in position after 15C warming
urrce2_nota_all_avg_new <- urrce2_nota_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce2_nota_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce2_nota_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce2_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 12) +
  geom_point(data = filter(urrce2_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 11) +
  # position after 15C warming
  geom_point(data = urrce2_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce2_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce2_nota_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") 

#### POMPOM 3; uneven reciprocal preference, unequal growth rates ####
urrce3_nota_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
                      ref_temp = 10,
                      r_EaN = mean(rgr_post_dist$intercept),
                      r_EaP = mean(rgr_post_dist$intercept), 
                      c_Ea1N = mean(c_post_dist$intercept),
                      c_Ea1P = mean(c_post_dist$intercept), 
                      c_Ea2N = mean(c_post_dist$intercept),
                      c_Ea2P = mean(c_post_dist$intercept), 
                      K_EaN = mean(k_post_dist$intercept), 
                      K_EaP = mean(k_post_dist$intercept), 
                      v_EaN = mean(v_post_dist$intercept),
                      v_EaP = mean(v_post_dist$intercept), 
                      m_Ea1 = mean(m_post_dist$intercept), 
                      m_Ea2 = mean(m_post_dist$intercept),
                      c1N_b = 1.3, c1P_b = 0.9, 
                      c2N_b = 0.9, c2P_b = 0.2, 
                      r_N_b = 0.7, r_P_b = 2, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.9, 
                      v2N_b = 0.9, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce3_nota_all <- bind_rows(urrce3_nota_all, hold) 
}

#get average change in position after 15C warming
urrce3_nota_all_avg_new <- urrce3_nota_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce3_nota_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce3_nota_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce3_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 12) +
  geom_point(data = filter(urrce3_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 11) +
  # position after 15C warming
  geom_point(data = urrce3_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce3_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce3_nota_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none")

#### POMPOM 4; uneven reciprocal preference, unequal growth rates ####
urrce4_nota_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
                      ref_temp = 10,
                      r_EaN = mean(rgr_post_dist$intercept),
                      r_EaP = mean(rgr_post_dist$intercept), 
                      c_Ea1N = mean(c_post_dist$intercept),
                      c_Ea1P = mean(c_post_dist$intercept), 
                      c_Ea2N = mean(c_post_dist$intercept),
                      c_Ea2P = mean(c_post_dist$intercept), 
                      K_EaN = mean(k_post_dist$intercept), 
                      K_EaP = mean(k_post_dist$intercept), 
                      v_EaN = mean(v_post_dist$intercept),
                      v_EaP = mean(v_post_dist$intercept), 
                      m_Ea1 = mean(m_post_dist$intercept), 
                      m_Ea2 = mean(m_post_dist$intercept),
                      c1N_b = 1.3, c1P_b = 0.4, 
                      c2N_b = 0.9, c2P_b = 0.4, 
                      r_N_b = 0.7, r_P_b = 1, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.4, 
                      v2N_b = 0.9, v2P_b = 0.4, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce4_nota_all <- bind_rows(urrce4_nota_all, hold) 
}

#get average change in position after 15C warming
urrce4_nota_all_avg_new <- urrce4_nota_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce4_nota_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce4_nota_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce4_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 12) +
  geom_point(data = filter(urrce4_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 11) +
  # position after 15C warming
  geom_point(data = urrce4_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce4_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce4_nota_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none")

#### POMPOM 5; uneven reciprocal preference, unequal growth rates ####
urrce5_nota_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
                      ref_temp = 10,
                      r_EaN = mean(rgr_post_dist$intercept),
                      r_EaP = mean(rgr_post_dist$intercept), 
                      c_Ea1N = mean(c_post_dist$intercept),
                      c_Ea1P = mean(c_post_dist$intercept), 
                      c_Ea2N = mean(c_post_dist$intercept),
                      c_Ea2P = mean(c_post_dist$intercept), 
                      K_EaN = mean(k_post_dist$intercept), 
                      K_EaP = mean(k_post_dist$intercept), 
                      v_EaN = mean(v_post_dist$intercept),
                      v_EaP = mean(v_post_dist$intercept), 
                      m_Ea1 = mean(m_post_dist$intercept), 
                      m_Ea2 = mean(m_post_dist$intercept),
                      c1N_b = 1.3, c1P_b = 0.8, 
                      c2N_b = 0.3, c2P_b = 2, 
                      r_N_b = 0.8, r_P_b = 0.8, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.8, 
                      v2N_b = 0.1, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce5_nota_all <- bind_rows(urrce5_nota_all, hold) 
}

#get average change in position after 15C warming
urrce5_nota_all_avg_new <- urrce5_nota_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce5_nota_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce5_nota_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce5_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 12) +
  geom_point(data = filter(urrce5_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 11) +
  # position after 15C warming
  geom_point(data = urrce5_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce5_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce5_nota_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") 


#### POMPOM 6; uneven reciprocal preference, unequal growth rates ####
urrce6_nota_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
                      ref_temp = 10,
                      r_EaN = mean(rgr_post_dist$intercept),
                      r_EaP = mean(rgr_post_dist$intercept), 
                      c_Ea1N = mean(c_post_dist$intercept),
                      c_Ea1P = mean(c_post_dist$intercept), 
                      c_Ea2N = mean(c_post_dist$intercept),
                      c_Ea2P = mean(c_post_dist$intercept), 
                      K_EaN = mean(k_post_dist$intercept), 
                      K_EaP = mean(k_post_dist$intercept), 
                      v_EaN = mean(v_post_dist$intercept),
                      v_EaP = mean(v_post_dist$intercept), 
                      m_Ea1 = mean(m_post_dist$intercept), 
                      m_Ea2 = mean(m_post_dist$intercept),
                      c1N_b = 1.3, c1P_b = 0.6, 
                      c2N_b = 0.3, c2P_b = 2, 
                      r_N_b = 0.8, r_P_b = 1.5, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.8, 
                      v2N_b = 0.1, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce6_nota_all <- bind_rows(urrce6_nota_all, hold) 
}

#get average change in position after 15C warming
urrce6_nota_all_avg_new <- urrce6_nota_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce6_nota_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce6_nota_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce6_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 12) +
  geom_point(data = filter(urrce6_nota_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 11) +
  # position after 15C warming
  geom_point(data = urrce6_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce6_nota_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce6_nota_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 26) 

# no TAs multiplot ------------
diff_start_points_notas <- (log_pom_urrce2_nota_all + log_pom_urrce3_nota_all +  log_pom_urrce5_nota_all) / (log_pom_urrce4_nota_all + log_pom_urrce1_nota_all + log_pom_urrce6_nota_all) + plot_annotation(tag_levels = "A")

