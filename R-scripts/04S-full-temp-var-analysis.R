# This is a daughter script of script 04 to explore the reviewers' analysis suggestions.

# Parent script info: This script is to test effects of all parameters varying simultaneously with temperature, and investigate drivers of large changes in competition with warming. The competing consumers in this analysis each specialize on one of two resources and they have equally strong preference for this resource. The two resources have uneven growth rates under ambient conditions, which places the species pair on the boundary of coexistence under ambient conditions. Parameters defining these starting conditions are given in the simulations below, and in summary in Table S1. Exploration of other starting conditions is conducted in script 05. In each simulation, each MacArthur consumer-resource parameter is given by an Arrhenius function, with a temperature sensitivity (activation energy, slope) term and an intercept term, which determines the value of the function at ambient temperatures (Tref, ref temp). In each simulation, temperature sensitivities are defined as "{parameter_EAik}", where ik captures the relevant consumer, resource, or both, and intercepts are defined as "{parameter-ik_b}". Consumers are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively, throughout the script. The script simulates the effects of warming when each parameter is given a temperature sensitivity, randomly drawn from the parameter's empirical distribution (generated in 01-param-dists), simultaneously.

# Author: Kaleigh Davis, PDF University of Guelph
# Script DOB: 15 December 2025

#load packages
library(tidyverse)
library(janitor)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)
library(purrr)
library(viridis)
library(beepr)

# get referencing set up for macarthur temp dependence function
source("R-scripts/02-temp-dep-macarthur.R") #base CR model - two resources and arrhenius-style temperature dependence
source("R-scripts/supplemental-analysis/02-temp-dep-macarthur-4r.R") #CR model with four shared resources, rather than two
source("R-scripts/supplemental-analysis/02S-temp-dep-macarthur-unimodal.R") #unimodal temperature dependence model and CR model adapted for this

#load in distributions for parameter values ####
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


##### draw all param EAs at random -- simulate extreme warming in order to observe model behaviour ##############
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

#get average change in position after 15C warming
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
  geom_ribbon(data = data.frame(x = seq(0, 2.5, 0.001)),
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
  coord_cartesian(ylim = c(-10, 10), xlim = c(-0.022, 1.5)) +
  # scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  # annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  # annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

#KD - pick up here adjusting to 200 degree analysis
# get euclidean distances for each species pair
rrc_200_e <- rrc_200 %>% 
  filter(T %in% c(10, 25)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) 

hist(rrc_200_e$dist15)

#histogram plot of euclidean distances in the pom pom plot
pom_hist <- rrc_200_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n15°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_200_p <- rrc_200 %>% 
  filter(T %in% c(10, 25)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+15°C Warming"))

rrc_200_p_avg <- rrc_200_p %>% 
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
  geom_jitter(data = filter(rrc_200_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_200_p_avg, aes(x = temp, y = med_stabil_potential, fill = temp), size = 5, pch = 21) +
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
  geom_jitter(data = filter(rrc_200_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_200_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
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

# ggsave(plot = comb_plot1, filename = "figures/pom_hist_nfd.pdf", width = 12, height = 10)

# Repeat analysis with higher degrees of warming #####
rrc50 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 210, by = 1), 
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
  filter(rel_T == 200) %>% # = degrees warmed
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_fit_rat = median(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential))

#pompom
# log_pom50 <-
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

# ggsave(plot = comb_plot50, filename = "figures/50C_warm_pom_hist_nfd.pdf", width = 12, height = 10)

# Repeat analysis more resources - 4 resources in favour of 2 ####
#two more resources - C and D, which species have the exact same preferences for as the original resources
#need to make a temp_dep_mac that can handle four resources -- NEXT UP, KD
rrc_4r <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac_4r(T = seq(10, 25, by = 0.1), 
                      ref_temp = 10,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept,
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
                      r_EaC = sample_n(rgr_post_dist, size = 1)$intercept,
                      r_EaD = sample_n(rgr_post_dist, size = 1)$intercept,
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea1C = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea1D = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2C = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2D = sample_n(c_post_dist, size = 1)$intercept, 
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaC = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaD = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaC = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaD = sample_n(v_post_dist, size = 1)$intercept,
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.5, c1P_b = 1,
                      c1C_b = 0.5, c1D_b = 1,
                      c2N_b = 1, c2P_b = 0.5,
                      c2C_b = 1, c2D_b = 0.5,
                      r_N_b = 1, r_P_b = 0.5, 
                      r_C_b = 1, r_D_b = 0.5,
                      K_N_b= 2000, K_P_b = 2000,
                      K_C_b = 2000, K_D_b = 2000,
                      v1N_b = 0.5, v1P_b = 1,
                      v1C_b = 0.5, v1D_b = 1,
                      v2N_b = 1, v2P_b = 0.5,
                      v2C_b = 1, v2D_b = 0.5,
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  rrc_4r <- bind_rows(rrc_4r, hold) 
}

#get average change in position after 15C warming
rrc_4r_avg_new <- rrc_4r %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
# log_pom4r <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_4r, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_4r, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_4r, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_4r_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_4r_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_4r_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.022, 0.5), xlim = c(-0.022, 0.5)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5)) +
  theme_cowplot(font_size = 20)

# get euclidean distances
rrc_4r_e <- rrc_4r %>% 
  filter(T %in% c(10, 15)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T15_new_stabil_potential - T10_new_stabil_potential)^2 + (T15_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T15_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T15_new_stabil_potential - T10_new_stabil_potential) 

hist(rrc_4r_e$dist15)

#histogram plot of euclidean dsistances in the pom pom plot
pom_hist5<- rrc_4r_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n5°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_4r_p <- rrc_4r %>% 
  filter(T %in% c(10, 15)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+5°C Warming"))

rrc_4r_p_avg <- rrc_4r_p %>% 
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
  geom_jitter(data = filter(rrc_4r_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_4r_p_avg, aes(x = temp, y = med_stab_pot, fill = temp), size = 5, pch = 21) +
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
  geom_jitter(data = filter(rrc_4r_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_4r_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
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

comb_plot5 <- log_pom4r / bottom_patch5 + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

# ggsave(plot = comb_plot5, filename = "figures/4r_pom_hist_nfd.pdf", width = 12, height = 10)

#can I track the euclidean distance from neutrality for each species pair with warming?
#get start position
start_new_stab_pot <- rrc_4r %>% 
  filter(T == 10) %>% 
  summarise(mean_start_nd = mean(new_stabil_potential),
            sd_stard_nd = sd(new_stabil_potential)) %>% 
  dplyr::select(mean_start_nd) %>% 
  unlist()

start_new_fit_rat <- rrc_4r %>% 
  filter(T == 10) %>% 
  summarise(mean_start_fd = mean(new_fit_ratio),
            sd_stard_fd = sd(new_fit_ratio)) %>% 
  dplyr::select(mean_start_fd) %>% 
  unlist()


track_4r <- rrc_4r%>% 
  dplyr::select(T, iteration, new_stabil_potential, new_fit_ratio, r_EaN:m_Ea2) %>%  #param combo identifies the simulation replicate
  mutate(start_nd = start_new_stab_pot,
         start_fd = start_new_fit_rat,
         ED = sqrt((new_stabil_potential - start_nd)^2 + (new_fit_ratio - start_fd)^2),
         dist_neut = sqrt(new_stabil_potential^2 + new_fit_ratio^2),
         iter_for_grouping = sprintf("%03d", iteration),
         iter_group = str_extract(as.character(iter_for_grouping), "^[0-9]", group = NULL))

track_4r %>% 
  ggplot() + 
  geom_line(aes(x = T, y = ED, colour = iteration, group = iteration)) 
# facet_wrap(~iter_group)

track_4r %>% 
  ggplot() + 
  geom_line(aes(x = T, y = dist_neut, colour = iteration, group = iteration)) + 
  ylab("Distance to neutrality") + 
  facet_wrap(~iter_group, scales = "free")

# Repeat analysis more resources - 4 resources more neutral ####
#two more resources - C and D, which species have the opposite preferences from the original resources
rrc_4r1 <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac_4r(T = seq(10, 25, by = 0.1), 
                         ref_temp = 10,
                         r_EaN = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaC = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaD = sample_n(rgr_post_dist, size = 1)$intercept,
                         c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea1C = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea1D = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2C = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2D = sample_n(c_post_dist, size = 1)$intercept, 
                         K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaC = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaD = sample_n(k_post_dist, size = 1)$intercept, 
                         v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaP = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaC = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaD = sample_n(v_post_dist, size = 1)$intercept,
                         m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                         m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                         c1N_b = 0.5, c1P_b = 1,
                         c1C_b = 1, c1D_b = 0.5,
                         c2N_b = 1, c2P_b = 0.5,
                         c2C_b = 0.5, c2D_b = 1,
                         r_N_b = 1, r_P_b = 0.5, 
                         r_C_b = 1, r_D_b = 0.5,
                         K_N_b= 2000, K_P_b = 2000,
                         K_C_b = 2000, K_D_b = 2000,
                         v1N_b = 0.5, v1P_b = 1,
                         v1C_b = 1, v1D_b = 0.5,
                         v2N_b = 1, v2P_b = 0.5,
                         v2C_b = 0.5, v2D_b = 1,
                         m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  rrc_4r1 <- bind_rows(rrc_4r1, hold) 
}

#get average change in position after 15C warming
rrc_4r1_avg_new <- rrc_4r1 %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
# log_pom4r <-
ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_4r1, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_4r1, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_4r1, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_4r1_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_4r1_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_4r1_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.5, 0.5), xlim = c(-0.022, 0.5)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5)) +
  theme_cowplot(font_size = 20)

# get euclidean distances
rrc_4r1_e <- rrc_4r1 %>% 
  filter(T %in% c(10, 15)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T15_new_stabil_potential - T10_new_stabil_potential)^2 + (T15_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T15_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T15_new_stabil_potential - T10_new_stabil_potential) 

hist(rrc_4r1_e$dist15)

#histogram plot of euclidean dsistances in the pom pom plot
pom_hist5<- rrc_4r1_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n5°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_4r1_p <- rrc_4r1 %>% 
  filter(T %in% c(10, 15)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+5°C Warming"))

rrc_4r1_p_avg <- rrc_4r1_p %>% 
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
  geom_jitter(data = filter(rrc_4r1_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_4r1_p_avg, aes(x = temp, y = med_stab_pot, fill = temp), size = 5, pch = 21) +
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
  geom_jitter(data = filter(rrc_4r1_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_4r1_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
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

comb_plot5 <- log_pom4r / bottom_patch5 + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

# ggsave(plot = comb_plot5, filename = "figures/4r_pom_hist_nfd.pdf", width = 12, height = 10)

# make temperature dependence unimodal ##############
rrc_uni <- data.frame()
for(f in 1:500){ 
  hold = uni_temp_dep_mac_spec_diffs(T = seq(10, 25, by = 0.1), 
                      ED = 3,
                      Topt_C1 = 25,
                      Topt_C2 = 25,
                      Topt_Cr = 25,
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
  rrc_uni <- bind_rows(rrc_uni, hold) 
}

#get average change in position after 15C warming
rrc_uni_avg_new <- rrc_uni %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
# log_pom_uni <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_uni, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_uni, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_uni, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_uni_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_uni_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_uni_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
    # coord_cartesian(ylim = c(-5, 5), xlim = c(-2.5, 5)) +
  # scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  # annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  # annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  # annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)
  
# ggsave(filename = "figures/unimodal-pompom.pdf", plot = last_plot())

  #can I track the euclidean distance from neutrality for each species pair with warming?
  #get start position
  start_new_stab_pot <- rrc_uni %>% 
    filter(T == 10) %>% 
    summarise(mean_start_nd = mean(new_stabil_potential),
              sd_stard_nd = sd(new_stabil_potential)) %>% 
    dplyr::select(mean_start_nd) %>% 
    unlist()
  
  start_new_fit_rat <- rrc_uni %>% 
    filter(T == 10) %>% 
    summarise(mean_start_fd = mean(new_fit_ratio),
              sd_stard_fd = sd(new_fit_ratio)) %>% 
    dplyr::select(mean_start_fd) %>% 
    unlist()
  
  
  track2 <- rrc_uni %>% 
    dplyr::select(T, iteration, new_stabil_potential, new_fit_ratio, r_EaN:m_Ea2) %>%  #param combo identifies the simulation replicate
    mutate(start_nd = start_new_stab_pot,
           start_fd = start_new_fit_rat,
           ED = sqrt((new_stabil_potential - start_nd)^2 + (new_fit_ratio - start_fd)^2),
           dist_neut = sqrt(new_stabil_potential^2 + new_fit_ratio^2),
           iter_for_grouping = sprintf("%03d", iteration),
           iter_group = str_extract(as.character(iter_for_grouping), "^[0-9]", group = NULL))
  
  track2 %>% 
    ggplot() + 
    geom_line(aes(x = T, y = ED, colour = iteration, group = iteration)) 
  # facet_wrap(~iter_group)
  
  track2 %>% 
    ggplot() + 
    geom_line(aes(x = T, y = dist_neut, colour = iteration, group = iteration)) + 
    scale_x_continuous(breaks = c(10, 5, 20, 25)) +
    labs(y = "Distance to neutrality", x = "Temperature") + 
    geom_vline(xintercept = 23, linetype = 5) 
  

# make temperature dependence unimodal with different Topts for different species##############
rrc_uni_diffs <- data.frame()
for(f in 1:500){ 
  hold = uni_temp_dep_mac_spec_diffs(T = seq(10, 25, by = 0.1), 
                          ED = 3,
                          Topt_Cr = 25,
                          Topt_C1 = 23,
                          Topt_C2 = 24,
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
  rrc_uni_diffs <- bind_rows(rrc_uni_diffs, hold) 
}

#get average change in position after 15C warming
rrc_uni_diffs_avg_new <- rrc_uni_diffs %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_uni_diffs <-
ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_uni_diffs, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_uni_diffs, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_uni_diffs, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_uni_diffs_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_uni_diffs_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_uni_diffs_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  # coord_cartesian(ylim = c(-5, 5), xlim = c(-2.5, 5)) +
  # scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# ggsave(filename = "figures/unimodal-vartopt-pompom.pdf", plot = last_plot(), width = 8, height = 6, units = "in", device = "pdf")

#can I track the euclidean distance from neutrality for each species pair with warming?
#get start position
start_new_stab_pot <- rrc_uni_diffs %>% 
  filter(T == 10) %>% 
  summarise(mean_start_nd = mean(new_stabil_potential),
            sd_stard_nd = sd(new_stabil_potential)) %>% 
  dplyr::select(mean_start_nd) %>% 
  unlist()

start_new_fit_rat <- rrc_uni_diffs %>% 
  filter(T == 10) %>% 
  summarise(mean_start_fd = mean(new_fit_ratio),
            sd_stard_fd = sd(new_fit_ratio)) %>% 
  dplyr::select(mean_start_fd) %>% 
  unlist()


track <- rrc_uni_diffs %>% 
  dplyr::select(T, iteration, new_stabil_potential, new_fit_ratio, r_EaN:m_Ea2) %>%  #param combo identifies the simulation replicate
  mutate(start_nd = start_new_stab_pot,
         start_fd = start_new_fit_rat,
         ED = sqrt((new_stabil_potential - start_nd)^2 + (new_fit_ratio - start_fd)^2),
         dist_neut = sqrt(new_stabil_potential^2 + new_fit_ratio^2),
         iter_for_grouping = sprintf("%03d", iteration),
         iter_group = str_extract(as.character(iter_for_grouping), "^[0-9]", group = NULL))

track %>% 
  ggplot() + 
  geom_line(aes(x = T, y = ED, colour = iteration, group = iteration)) 
  # facet_wrap(~iter_group)

track %>% 
  ggplot() + 
  geom_line(aes(x = T, y = dist_neut, colour = iteration, group = iteration)) + 
  scale_x_continuous(breaks = c(10, 5, 20, 23, 25)) +
  labs(y = "Distance to neutrality", x = "Temperature") + 
  geom_vline(xintercept = 23, linetype = 5) 
  
