# This script is to reproduce the main results figures of the paper -- figs 3 and 5 -- with different start points. We achieve these different start points by manipulating the resource preferences (analogous to niche widths) of each species at the models starting conditions, i.e. ambient temperature. Scenario 1 is given in the main text (script 04) and features two consumers that each specialize on one of two resources and they have equally strong preference for this resource. The two resources have uneven growth rates under ambient conditions, which places the species pair on the boundary of coexistence under ambient conditions. In scenario 2, consumer preferences are the same as in Scenario 1, but the resources grow at the same rate under ambient conditions, which moves the start point from the coexistence boundary to the middle of the coexistence region. In scenario 3, species have very different resource preferences, where one consumer has a strong preference for its preferred resource and the other consumer has a weak preference for its preferred resources. Both resources grow at the same rate under ambient conditions. Parameters defining each of these starting conditions are given in each simulation below, and in summary in Table S1. In each simulation, each MacArthur consumer-resource parameter is given by an Arrhenius function, with a temperature sensitivity (activation energy, slope) term and an intercept term, which determines the value of the function at ambient temperatures (Tref, ref temp). In each simulation, temperature sensitivities are defined as "{parameter_EAik}", where ik captures the relevant consumer, resource, or both, and intercepts are defined as "{parameter-ik_b}". Consumers are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively, throughout the script. The script simulates the effects of warming when each parameter is given a temperature sensitivity, randomly drawn from the parameter's empirical distribution (generated in 01-param-dists), simultaneously.

#script DOB: 9/5/2025
#author: Kaleigh Davis, Postdoc UoG with Joey Bernhardt

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

# get referencing set up for MacArthur temp dependence function
source("R-scripts/02-temp-dep-macarthur.R") #this contains the MacArthur translation function, with all parameters flexibly defined in the function for assigning at time of use, and the arrhenius function.

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

#### ONE-BY-ONE; rrc, equal growth rates -- scenario 2####
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

### rrc equal start two plots #####
startpoint1 <- rrce_plots / log_pom_rrce_all + plot_annotation(tag_levels = "A")
# ggsave(plot = startpoint1, filename = "figures/supp_startpoint_errce.pdf", height = 20, width = 26)

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

#### POMPOM; uneven reciprocal preference, equal growth rates - scenario 3 ####
urrce_all <- data.frame()
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

### uneven reciprocal preference, equal start two plots #####
startpoint2 <- urrce_plots / log_pom_urrce_all + plot_annotation(tag_levels = "A")
# ggsave(plot = startpoint2, filename = "figures/supp_startpoint_urrce.pdf", height = 20, width = 26)
