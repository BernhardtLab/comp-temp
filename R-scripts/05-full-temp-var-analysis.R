# This script tests effects of all parameters varying simultaneously with temperature, and investigates drivers of large changes in competition with warming. The competing consumers in this analysis each specialize on one of two resources and they have equally strong preference for this resource. The two resources have uneven growth rates under ambient conditions, which places the species pair on the boundary of coexistence under ambient conditions. Parameters defining these starting conditions are given in the simulations below, and in summary in Table S1. Exploration of other starting conditions is conducted in script 06. In each simulation, each MacArthur consumer-resource parameter is given by an Arrhenius function, with a temperature sensitivity (activation energy, slope) term and an intercept term, which determines the value of the function at ambient temperatures (Tref, ref temp). In each simulation, temperature sensitivities are defined as "{parameter_EAik}", where ik captures the relevant consumer, resource, or both, and intercepts are defined as "{parameter-ik_b}". Consumers are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively, throughout the script. The script simulates the effects of warming when each parameter is given a temperature sensitivity, randomly drawn from the parameter's empirical distribution (generated in 01-param-dists), simultaneously.

#This script produces Figures 5 in the main text and all supplementary figures from the same set of starting conditions with four resources.

# Author: Kaleigh Davis, PDF University of Guelph
# Script DOB: 30 April 2025

# sourced scripts:
#   R-scripts/02-temp-dep-macarthur.R — temperature-dependent MacArthur
#     consumer-resource function
#   R-scripts/03-arrhenius.R — Arrhenius temperature-dependence function

# inputs:
#   data/processed-data/param_post_dists.csv — posterior distributions of
#     temperature sensitivity for each model parameter (generated in script 01)

# outputs:
#   figures/Fig5-pom_hist_nfd.pdf — main text figure 5
#   figures/S5-5C_warm_pom_hist_nfd.pdf — supplementary figure S5
#   figures/S12-50C_warm_pom_hist_nfd.pdf — supplementary figure S12
#   figures/S8-NOTA_pom_hist_nfd.pdf — supplementary figure S8
#   figures/S10-200C_param_trajectories.pdf — supplementary figure S10
#   figures/S11-extreme_ea_pom_hist_nfd.pdf - supplementary figure S11

#load packages
library(tidyverse)
library(janitor)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)
library(purrr)
library(viridis)
library(beepr)
library(see)
library(visreg)

# get referencing set up for macarthur temp dependence function
source("R-scripts/02-temp-dep-macarthur.R") #this calls the macarthur translation function, with all parameters flexibly defined in the function for assigning at time of use
source("R-scripts/03-arrhenius.R") # this calls the arrhenius function.

#load in distributions for parameter values.
# these are continuous distributions generated from empirical data using MCMC regression, in /R-scripts/01-param-dists.R
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

################################################################################
##############################  MAIN ANALYSIS   ################################
################################################################################

##### draw all param EAs at random -- for figure 5 ##############
rrc <- data.frame()
for(f in 1:500){ 
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  rrc <- bind_rows(rrc, hold) 
}

#get average change in position after 15C warming
rrc_avg_new <- rrc %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  # annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  # annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

#for how many species pairs does ND decrease and does FD decrease
rrc_start <- rrc %>% 
  filter(T == 10) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_sd_stab_pot = sd(new_stabil_potential),
            new_sd_fit_rat = sd(new_fit_ratio)) %>% 
  dplyr::select(new_mean_stab_pot, new_mean_fit_rat)

rrc_shifts <- rrc %>% 
  dplyr::select(iteration, T, new_stabil_potential, new_fit_ratio) %>% 
  mutate(stab_pot_start = rrc_start$new_mean_stab_pot,
         fit_rat_start = rrc_start$new_mean_fit_rat,
         shift_fd = new_fit_ratio - fit_rat_start,
         shift_nd = new_stabil_potential - stab_pot_start,
         shift_fd_prop = shift_fd/fit_rat_start,
         shift_nd_prop = shift_nd/stab_pot_start) %>% 
  mutate(fd_shift_cat = ifelse(new_fit_ratio < fit_rat_start, "decrease", 
                           ifelse(new_fit_ratio > fit_rat_start, "increase",
                                  ifelse(new_fit_ratio == fit_rat_start, "no change", "potato"))),
         nd_shift_cat = ifelse(new_stabil_potential < stab_pot_start, "decrease",
                           ifelse(new_stabil_potential > stab_pot_start, "increase",
                                  ifelse(new_stabil_potential == stab_pot_start, "no change","potato")))) %>% 
  filter(T == 10 | T == 25)

rrc_shifts %>% 
  group_by(T, fd_shift_cat, nd_shift_cat) %>% 
  tally()

## FD: 288/500 decrease (another sim - 275/500, 300/500)
## ND: 289/500 decrease (another sim - 291/500, 326/500)
## ND + FD: 185/500 decrease in both (another sim - 170/500, 204/500)

shift_sums <- rrc_shifts %>% 
  filter(T == 25) %>% 
  summarise(across(
    c(shift_fd, shift_fd_prop, shift_nd, shift_nd_prop),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      median = ~median(.x, na.rm = TRUE))))

# does euclidean distance relate to thermal asymmetry -----------
# get euclidean distances for each species pair
rrc_e <- rrc %>% 
  filter(T %in% c(10, 25)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  #thermal asymmetries
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2)) %>% 
         #shifts in competition
  mutate(dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(rrc_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#which thermal asymmetries drive large shifts in competition?
rrc_ed <- rrc_e %>% 
  filter(response_var == "dist15")

# run multiple regression to test for correlations between TAs and Euclidean distance from start point to warmed point - Table S4 ----------------------------
model <- lm(value ~ abs_r_ta + abs_c1_ta + abs_c2_ta + abs_c3_ta + abs_c4_ta + abs_c5_ta + abs_c6_ta + abs_k_ta + abs_v_ta + abs_m_ta, data = rrc_ed)

visreg::visreg(model)
summary(model) #1-3 consumption rates are significant, depending on the simulation run. Effect sizes are repeatable.
anova(model)
# r, c4 (E_c2b-E_c2a), k, v significant

# histogram plot of euclidean distances in the pom pom plot
pom_hist <- rrc_e %>% 
  filter(response_var == "dist15") %>% 
  ggplot() + 
  geom_histogram(aes(x = value), binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n15°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_p <- rrc %>% 
  filter(T %in% c(10, 25)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+15°C Warming"))

rrc_p_avg <- rrc_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            med_stabil_potential = median(new_stabil_potential),
            med_fitrat = median(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift <-
  ggplot() + 
  geom_jitter(data = filter(rrc_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_p_avg, aes(x = temp, y = med_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.75)) 

#shift in fitness ratio
fd_shift <-
  ggplot() + 
  geom_jitter(data = filter(rrc_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none")  +
  coord_cartesian(ylim = c(0, 0.75))

nd_shift + fd_shift

# Pompom with subplots - Figure 5 
bottom_patch <- pom_hist + nd_shift + fd_shift

comb_plot1 <- log_pom / bottom_patch + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

# ggsave(plot = comb_plot1, filename = "figures/Fig5-pom_hist_nfd.pdf", width = 12, height = 10)

################################################################################
#########################  SUPPLEMENTARY ANALYSES   ############################
################################################################################

# Repeat analysis with 50C warming - Figure S12 #####
rrc50 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 60, by = 0.5), 
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
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1, 
                      v2N_b = 1, v2P_b = 0.5, 
                      m1_b = 0.01, m2_b = 0.01)
  hold$iteration <- f
  rrc50 <- bind_rows(rrc50, hold) 
}

#get average change in position after 50C warming
rrc50_avg_new <- rrc50 %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 50) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_fit_rat = median(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential))

#pompom
log_pom50 <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.50, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc50, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc50, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc50, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 50C warming
  geom_point(data = rrc50_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc50_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc50_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.50)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# get euclidean distances
rrc50_e <- rrc50 %>% 
  filter(T %in% c(10, 60)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T60_new_stabil_potential - T10_new_stabil_potential)^2 + (T60_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T60_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T60_new_stabil_potential - T10_new_stabil_potential) 

hist(rrc50_e$dist15)

#histogram plot of euclidean dsistances in the pom pom plot
pom_hist50 <- rrc50_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n50°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc50_p <- rrc50 %>% 
  filter(T %in% c(10, 60)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+50°C Warming"))

rrc50_p_avg <- rrc50_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            med_stabil_potential = median(new_stabil_potential),
            med_fitrat = median(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift50 <-
  ggplot() + 
  geom_jitter(data = filter(rrc50_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc50_p_avg, aes(x = temp, y = med_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+50°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(-0.5, 0.65)) 

#shift in fitness ratio
fd_shift50 <-
  ggplot() + 
  geom_jitter(data = filter(rrc50_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc50_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+50°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  coord_cartesian(ylim = c(-0.5, 0.65))

nd_shift50 + fd_shift50

# Pompom with subplots
bottom_patch50 <- pom_hist50 + nd_shift50 + fd_shift50

comb_plot50 <- log_pom50 / bottom_patch50 + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

# ggsave(plot = comb_plot50, filename = "figures/S12-50C_warm_pom_hist_nfd.pdf", width = 12, height = 10)

# Repeat analysis with 5C warming - Figure S5 ####
rrc5 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 15, by = 0.1), 
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
                      c1N_b = 0.5, c1P_b = 1,
                      c2N_b = 1, c2P_b = 0.5,
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000,
                      v1N_b = 0.5, v1P_b = 1,
                      v2N_b = 1, v2P_b = 0.5,
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  rrc5 <- bind_rows(rrc5, hold) 
}

#get average change in position after 5C warming
rrc5_avg_new <- rrc5 %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 5) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom5 <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc5, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc5, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc5, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc5_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc5_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc5_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.022, 0.5), xlim = c(-0.022, 0.5)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5)) +
  theme_cowplot(font_size = 20)

# get euclidean distances
rrc5_e <- rrc5 %>% 
  filter(T %in% c(10, 15)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T15_new_stabil_potential - T10_new_stabil_potential)^2 + (T15_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T15_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T15_new_stabil_potential - T10_new_stabil_potential) 

hist(rrc5_e$dist15)

#histogram plot of euclidean dsistances in the pom pom plot
pom_hist5<- rrc5_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n5°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc5_p <- rrc5 %>% 
  filter(T %in% c(10, 15)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+5°C Warming"))

rrc5_p_avg <- rrc5_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            med_fitrat = median(new_fit_ratio),
            med_stab_pot = median(new_stabil_potential),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift5 <-
  ggplot() + 
  geom_jitter(data = filter(rrc5_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc5_p_avg, aes(x = temp, y = med_stab_pot, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+5°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.65)) 

#shift in fitness ratio
fd_shift5 <-
  ggplot() + 
  geom_jitter(data = filter(rrc5_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc5_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+5°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.65))

nd_shift5 + fd_shift5

# pompom with subplots
bottom_patch5 <- pom_hist5 + nd_shift5 + fd_shift5

comb_plot5 <- log_pom5 / bottom_patch5 + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

# ggsave(plot = comb_plot5, filename = "figures/S5-5C_warm_pom_hist_nfd.pdf", width = 12, height = 10)

# No thermal asymmetries -- Figure S8 ####
nota <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), # was by 0.1
                      ref_temp = 10,
                      r_EaN = mean(rgr_post_dist$intercept),
                      r_EaP = mean(rgr_post_dist$intercept),
                      # r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
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
                      c1N_b = 0.5, c1P_b = 1, 
                      c2N_b = 1, c2P_b = 0.5, 
                      r_N_b = 1, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.5, v1P_b = 1,
                      v2N_b = 1, v2P_b = 0.5,
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  nota <- bind_rows(nota, hold) 
}

#get average change in position after 5, 10, 20C warming
nota_avg_new <- nota %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#base pompom for comparison
log_pom_nota <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = nota, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(nota, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 9.5) +
  geom_point(data = filter(nota, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 8.5) +
  # position after 15C warming
  geom_point(data = nota_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7) +
  geom_point(data = nota_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 5.5) +
  geom_hline(yintercept = 0, linetype=5) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "Degrees \nC Warming") +
  annotate("text", x = 0.37, y = -0.05, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# get euclidean distances
nota_e <- nota %>% 
  filter(T %in% c(10, 25)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) 

hist(nota_e$dist15)

#histogram plot of euclidean distances in the pom pom plot
pom_hist_nota <- nota_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 1e-16, colour = "black") + 
  labs(x = "Euclidean distance with \n15C warming", y = "Count") + 
  scale_x_continuous(
    labels = scales::number_format(accuracy = 0.01)) +
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
nota_p <- nota %>% 
  filter(T %in% c(10, 25)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+15C Warming"))

nota_p_avg <- nota_p %>% 
  group_by(temp) %>% 
  summarise(med_stabil_potential = median(new_stabil_potential), 
            med_fitrat = median(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift_nota <-
  ggplot() + 
  geom_point(data = nota_p, aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = nota_p_avg, aes(x = temp, y = med_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.65)) 

#shift in fitness ratio
fd_shift_nota <-
  ggplot() + 
  geom_point(data = nota_p, aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = nota_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.65))

nd_shift_nota + fd_shift_nota

# NO TA Pompom with subplots 
# big pompom with small panels underneath
bottom_patch_nota <- pom_hist_nota + nd_shift_nota + fd_shift_nota

comb_plot_nota <- log_pom_nota / bottom_patch_nota + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

# ggsave(plot = comb_plot_nota, filename = "figures/S8-NOTA_pom_hist_nfd.pdf", width = 12, height = 10)

##### simulate extreme (200C) warming in order to observe model behaviour ##############
rrc_200 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 210, by = 1), #was by 0.1
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  rrc_200 <- bind_rows(rrc_200, hold) 
}

#get average change in position after 200C warming
rrc_200_avg_new <- rrc_200 %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 200) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom200 <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_200, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_200, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_200, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_200_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_200_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_200_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.75, 1.5), xlim = c(-0.022, 0.75)) +
  # scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  # annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  # annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# get euclidean distances for each species pair
rrc_200_e <- rrc_200 %>% 
  filter(T %in% c(10, 210)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist200 = sqrt((T210_new_stabil_potential - T10_new_stabil_potential)^2 + (T210_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T210_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T210_new_stabil_potential - T10_new_stabil_potential) 

hist(rrc_200_e$dist200)

#histogram plot of euclidean distances in the pom pom plot
pom_hist200 <- rrc_200_e %>% 
  ggplot(aes(x = dist200)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n200°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_200_p <- rrc_200 %>% 
  filter(T %in% c(10, 210)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+200°C Warming"))

rrc_200_p_avg <- rrc_200_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            med_stabil_potential = median(new_stabil_potential),
            med_fitrat = median(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift200 <-
  ggplot() + 
  geom_jitter(data = filter(rrc_200_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_200_p_avg, aes(x = temp, y = med_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+200°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(-1, 1)) 

#shift in fitness ratio
fd_shift200 <-
  ggplot() + 
  geom_jitter(data = filter(rrc_200_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_200_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+200°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none")  +
  coord_cartesian(ylim = c(-1, 1))

nd_shift200 + fd_shift200

# Pompom with subplots
bottom_patch200 <- pom_hist200 + nd_shift200 + fd_shift200

comb_plot200 <- log_pom200 / bottom_patch200 + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

### Plot vars driving neutrality -- Figure S10 ####
# develop relevant composite variables
rrc200a <-
  rrc_200 %>% 
  mutate(KN_over_rN = KN/rN,
         KP_over_rP = KP/rP,
         intra_inter1 = a11/a12,
         intra_inter2 = a22/a21,
         inter_intra1 = a12/a11,
         inter_intra2 = a21/a22) %>% 
  dplyr::select(-ref_temp, -coexist) %>% 
  dplyr::select(iteration, T, everything()) 

#plot KN/rN as a function of KN
rrc200a %>% 
  ggplot() + 
  geom_line(aes(x = KN, y = KN_over_rN, group = iteration, colour = K_EaN))

#because K has a negative temperature dependence, K approaches 0 exponentially as systems heat up. That means that K/r grows rapidly as the system heats up

#KN temp response relative to KN/rN
rrc200a %>% 
  ggplot() + 
  geom_line(aes(x = T, y = KN, group = iteration), colour = "blue") + #negative exponential
  geom_line(aes(x = T, y = KN_over_rN, group = iteration), colour = "red") #negative exponential divided by a positive exponential

#get summary trajectories for each parameter
rrc200a_sum <- rrc200a %>% 
  group_by(T) %>% 
  summarise(across(
    .cols = c(a11, a12, a21, a22, KN, KN_over_rN, new_fit_ratio, new_stabil_potential, rN, intra_inter1, intra_inter2, inter_intra1, inter_intra2, fit_ratio, rho),
    .fns = list(
      mean   = ~mean(.x, na.rm = TRUE),
      median = ~median(.x, na.rm = TRUE),
      sd     = ~sd(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ),
  .groups = "drop"
  )

rrc200_sums <- rrc200a_sum %>% 
  pivot_longer(cols = c(a11_mean:rho_sd), names_to = "parameter", values_to = "value") %>% 
  separate(parameter, sep = "_(?=[^_]+$)", into = c("parameter", "statistic")) %>%  
  filter(parameter %in% c("a11", "a12", "a21", "a22", "KN", "KN_over_rN", "new_fit_ratio", "new_stabil_potential", "rN", "intra_inter1", "intra_inter2", "inter_intra1", "inter_intra2", "fit_ratio", "rho") & statistic %in% c("median", "mean")) %>% 
  mutate(parameter = forcats::fct_relevel(parameter, "KN", "rN", "KN_over_rN", "a11", "a12", "a21", "a22", "intra_inter1", "intra_inter2", "inter_intra1", "inter_intra2", "fit_ratio", "rho","new_fit_ratio", "new_stabil_potential")) 

param_labels <- c(
  a11 = "alpha[11]",
  a12 = "alpha[12]",
  a21 = "alpha[21]",
  a22 = "alpha[22]",
  KN  = "K[a]",
  KN_over_rN = "K[a]/r[a]",
  new_fit_ratio = "log(Fitness~ratio)",
  new_stabil_potential = "-log(rho)",
  rN = "r[a]",
  intra_inter1 = "alpha[11]/alpha[12]",
  intra_inter2 = "alpha[22]/alpha[21]",
  inter_intra1 = "alpha[12]/alpha[11]",
  inter_intra2 = "alpha[21]/alpha[22]",
  fit_ratio = "Fitness~ratio",
  rho = "rho"
)

#divide temperatures
rrc200_sums_to25 <- rrc200_sums %>% 
  filter(T > 9, T < 26)

# trajectories of each param over 200C - Figure S10
warming_tendencies200 <-
  rrc200_sums %>% 
  filter(statistic == "median") %>% 
  filter(! parameter %in% c("inter_intra1", "inter_intra2", "fit_ratio", "rho")) %>% 
  # filter(T > 150) %>% 
  ggplot() + 
  geom_point(aes(x = T, y = value)) +
  # geom_hline(yintercept = 1, colour = "green4") + 
  facet_wrap(~parameter, 
             labeller = as_labeller(param_labels, label_parsed), 
             scales = "free") + 
  theme(strip.text = element_text(size = 16))

# ggsave(plot = warming_tendencies200, filename = "figures/S10-200C_param_trajectories.pdf", width = 12, height = 8)

#summaries of each param over the first 15C
warming_tendencies15 <- 
  rrc200_sums_to25 %>% 
  filter(statistic == "median") %>% 
  filter(! parameter %in% c("inter_intra1", "inter_intra2")) %>% 
  ggplot() + 
  # geom_hline(yintercept = 1, colour = "green4") + 
  geom_point(aes(x = T, y = value)) +
  facet_wrap(~parameter, 
             labeller = as_labeller(param_labels, label_parsed), 
             scales = "free")

#simulate effects of warming on species pairs with extreme thermal asymmetries - Figure S11 ----------------------------------
param_extremes <- param_vals %>% 
  group_by(parameter) %>% 
  mutate(q5 = quantile(intercept, 0.05),
         q95 = quantile(intercept, 0.95)) %>% 
  filter(intercept <= q5 | intercept >= q95) %>% 
  ungroup() %>% 
  dplyr::select(-q5, -q95)

hist(filter(param_extremes, parameter == "resource_growth_rate")$intercept)
hist(filter(param_extremes, parameter == "carrying_capacity")$intercept)

#split these into dfs for each parameter
param_extremes %>%
  mutate(parameter = str_replace(parameter, "resource_growth_rate", "rgr"),
         parameter = str_replace(parameter, "carrying_capacity", "k"),
         parameter = str_replace(parameter, "conversion_efficiency", "v"),
         parameter = str_replace(parameter, "mortality_rate", "m"),
         parameter = str_replace(parameter, "consumption rate", "c")) %>% 
  group_by(parameter) %>%
  group_split() %>%
  set_names(unique(param_vals$parameter)) %>%  # Set the names based on unique category values
  walk(~ assign(paste0(.x$parameter[1], "_extreme_pd"), .x, envir = .GlobalEnv))

# sim for all params drawn randomly from the extremes of the distribution --------------------
rrc_ex <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), #was by 0.1
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_extreme_pd, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_extreme_pd, size = 1)$intercept, 
                      c_Ea1N = sample_n(c_extreme_pd, size = 1)$intercept,
                      c_Ea1P = sample_n(c_extreme_pd, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_extreme_pd, size = 1)$intercept,
                      c_Ea2P = sample_n(c_extreme_pd, size = 1)$intercept, 
                      K_EaN = sample_n(k_extreme_pd, size = 1)$intercept, 
                      K_EaP = sample_n(k_extreme_pd, size = 1)$intercept, 
                      v_EaN = sample_n(v_extreme_pd, size = 1)$intercept,
                      v_EaP = sample_n(v_extreme_pd, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_extreme_pd, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_extreme_pd, size = 1)$intercept,
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  rrc_ex <- bind_rows(rrc_ex, hold) 
}

hist(rrc_ex$r_EaN)

#get average change in position after 15C warming
rrc_ex_avg_new <- rrc_ex %>% 
  mutate(rel_T = T - 10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_ex <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_ex, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_ex, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_ex, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_ex_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_ex_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_ex_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  # coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 24) 

#for how many species pairs does ND decrease and does FD decrease
rrc_ex_start <- rrc_ex %>% 
  filter(T == 10) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_sd_stab_pot = sd(new_stabil_potential),
            new_sd_fit_rat = sd(new_fit_ratio)) %>% 
  dplyr::select(new_mean_stab_pot, new_mean_fit_rat)

rrc_ex_shifts <- rrc_ex %>% 
  dplyr::select(iteration, T, new_stabil_potential, new_fit_ratio) %>% 
  mutate(stab_pot_start = rrc_ex_start$new_mean_stab_pot,
         fit_rat_start = rrc_ex_start$new_mean_fit_rat) %>% 
  mutate(fd_shift = ifelse(new_fit_ratio < fit_rat_start, "decrease", 
                           ifelse(new_fit_ratio > fit_rat_start, "increase",
                                  ifelse(new_fit_ratio == fit_rat_start, "no change", "potato"))),
         nd_shift = ifelse(new_stabil_potential < stab_pot_start, "decrease",
                           ifelse(new_stabil_potential > stab_pot_start, "increase",
                                  ifelse(new_stabil_potential == stab_pot_start, "no change","potato")))) %>% 
  filter(T == 10 | T == 25)

rrc_ex_shifts %>% 
  group_by(T, fd_shift, nd_shift) %>% 
  tally()

## FD: 316/500 decrease; 282/500
## ND: 299/500 decrease; 293/500
## ND + FD: 203/500 decrease in both; 177/500

# get euclidean distances for each species pair
rrc_ex_e <- rrc_ex %>% 
  filter(T %in% c(10, 25)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  #thermal asymmetries
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2)) %>% 
  #shifts in competition
  mutate(dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(rrc_ex_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#just euclidean distances
rrc_ex_ed <- rrc_ex_e %>% 
  filter(response_var == "dist15")

#histogram plot of euclidean distances in the pom pom plot
pom_hist_ex <- rrc_ex_ed %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n15°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_ex_p <- rrc_ex %>% 
  filter(T %in% c(10, 25)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+15°C Warming"))

rrc_ex_p_avg <- rrc_ex_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            med_stabil_potential = median(new_stabil_potential),
            med_fitrat = median(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift_ex <-
  ggplot() + 
  geom_jitter(data = filter(rrc_ex_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_ex_p_avg, aes(x = temp, y = med_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.75)) 

#shift in fitness ratio
fd_shift_ex <-
  ggplot() + 
  geom_jitter(data = filter(rrc_ex_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_ex_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none")  +
  coord_cartesian(ylim = c(0, 0.75))

nd_shift_ex + fd_shift_ex

# Pompom with subplots - Figure 5 
bottom_patch_ex <- pom_hist_ex + nd_shift_ex + fd_shift_ex

comb_plot_ex <- log_pom_ex / bottom_patch_ex + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

# ggsave(plot = comb_plot_ex, filename = "figures/S11-extreme_ea_pom_hist_nfd.pdf", width = 12, height = 10)

# prove results are robust to different ranges of temperatures - reviewer response #####
# ref temp 0
rrc0 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(0, 15, by = 0.1), #was by 0.1
                      ref_temp = 0,
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  rrc0 <- bind_rows(rrc0, hold) 
}

#get average change in position after 15C warming
rrc0_avg_new <- rrc0 %>% 
  mutate(rel_T = T) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom0 <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc0, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-0, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc0, T==0), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc0, T==0), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-0), size = 6) +
  # position after 15C warming
  geom_point(data = rrc0_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc0_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc0_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  # annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  # annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 24) + 
  ggtitle("0-15°C Warming")

# ref temp 0
rrc20 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(20, 35, by = 0.1), #was by 0.1
                      ref_temp = 20,
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  rrc20 <- bind_rows(rrc20, hold) 
}

#get average change in position after 15C warming
rrc20_avg_new <- rrc20 %>% 
  mutate(rel_T = T - 20) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom20 <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc20, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-20, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc20, T==20), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc20, T==20), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-20), size = 6) +
  # position after 15C warming
  geom_point(data = rrc20_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc20_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc20_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  # annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  # annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 24) + 
  ggtitle("20-35°C Warming")

ranges <- log_pom0 + log_pom20
#this was a reviewer requested proof of concept that did not make it into the manuscript or supplemental materials

# What if resources are all unicellular vs multicellular? -- reviewer response -----------------

#load in data with cellularity indicated
cparam_vals <- read_csv("data/processed-data/unimulti_param_post_dists.csv")

#split these into dfs for each parameter
cparam_vals %>%
  mutate(parameter = str_replace(parameter, "resource_growth_rate", "rgr"),
         parameter = str_replace(parameter, "carrying_capacity", "k"),
         parameter = str_replace(parameter, "conversion_efficiency", "v"),
         parameter = str_replace(parameter, "mortality_rate", "m"),
         parameter = str_replace(parameter, "consumption rate", "c")) %>% 
  group_by(parameter) %>%
  group_split() %>%
  set_names(unique(cparam_vals$parameter)) %>%  # Set the names based on unique category values
  walk(~ assign(paste0(.x$parameter[1], "_c_post_dist"), .x, envir = .GlobalEnv))

#get summary stats for all parameters ########
cparam_sum <- cparam_vals %>%
  group_by(parameter) %>% 
  summarize(
    across(
      intercept,
      list(
        bottom10 = ~quantile(., 0.10),
        top10 = ~quantile(., 0.90),
        Median   = ~median(., na.rm = TRUE),
        Min      = ~min(., na.rm = TRUE),
        Max      = ~max(., na.rm = TRUE),
        sd       = ~sd(., na.rm = TRUE)
      ),
      .names = "{.fn}" )
  ) 


#make longform
cparam_sum1 <- cparam_sum %>% 
  pivot_longer(cols = c(bottom10:sd), 
               names_to = "summary_stat",
               values_to = "value")

#split df up into dfs for each summary statistic
cparam_sum1 %>% 
  group_by(summary_stat) %>% 
  group_split() %>% 
  purrr::walk(~ assign(paste0(.x$summary_stat[1]), .x, envir = .GlobalEnv))

rrc_unic <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), #was by 0.1
                      ref_temp = 10,
                      r_EaN = sample_n(filter(rgr_c_post_dist, cellularity == "uni"), size = 1)$intercept,
                      r_EaP = sample_n(filter(rgr_c_post_dist, cellularity == "uni"), size = 1)$intercept, 
                      c_Ea1N = sample_n(c_c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_c_post_dist, size = 1)$intercept, 
                      K_EaN = sample_n(k_c_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_c_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_c_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_c_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_c_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_c_post_dist, size = 1)$intercept,
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  rrc_unic <- bind_rows(rrc_unic, hold) 
}

#get average change in position after 15C warming
rrc_unic_avg_new <- rrc_unic %>% 
  mutate(rel_T = T - 10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

# pompom
log_pom_unic <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_unic, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_unic, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_unic, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_unic_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_unic_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_unic_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  # annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  # annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 24) + 
  ggtitle("All unicellular rgr")

#repeat with multicellular only
rrc_multic <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), #was by 0.1
                      ref_temp = 10,
                      r_EaN = sample_n(filter(rgr_c_post_dist, cellularity == "multi"), size = 1)$intercept,
                      r_EaP = sample_n(filter(rgr_c_post_dist, cellularity == "multi"), size = 1)$intercept, 
                      c_Ea1N = sample_n(c_c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_c_post_dist, size = 1)$intercept, 
                      K_EaN = sample_n(k_c_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_c_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_c_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_c_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_c_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_c_post_dist, size = 1)$intercept,
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  rrc_multic <- bind_rows(rrc_multic, hold)
}

#get average change in position after 15C warming
rrc_multic_avg_new <- rrc_multic %>% 
  mutate(rel_T = T - 10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

# pompom
log_pom_multic <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_multic, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_multic, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_multic, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_multic_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_multic_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_multic_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  # annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  # annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 24) + 
  ggtitle("All multicellular rgr")

cellularity_plots <- (log_pom_unic + theme(legend.position = "none")) + (log_pom_multic + theme(legend.position = "none"))
#this was a reviewer requested analysis that did not make it into the manuscript or supplemental materials
