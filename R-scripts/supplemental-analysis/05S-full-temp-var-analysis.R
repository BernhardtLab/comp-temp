# This is a daughter script of script 04 to explore the effects of 1) a higher degree of warming; 2) additional resources; 3) unimodal temperature dependence on competitive outcomes; and 4) temperature sensitivities only at the extremes of observed values

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
source("R-scripts/02-temp-dep-macarthur.R") #base CR model - two resources and arrhenius-form temperature dependence
source("R-scripts/03-arrhenius.R") #arrhenius function
source("R-scripts/supplemental-analysis/02S-temp-dep-macarthur-4r.R") #CR model with four shared resources, rather than two, arrhenius-form temperature dependence
source("R-scripts/supplemental-analysis/02S-temp-dep-macarthur-8r.R") #CR model with eight shared resources, rather than two, arrhenius-form temperature dependence
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

# ggsave(plot = comb_plot200, filename = "figures/pom_hist_nfd200.pdf", width = 12, height = 10)


# Repeat analysis more resources - 4 resources in favour of 2 ####
#two more resources - C and D, which species have the exact same preferences for as the original resources
#need to make a temp_dep_mac that can handle four resources
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

#track the euclidean distance from neutrality for each species pair with warming
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

# Repeat analysis more resources - 8 resources in favour of 2 ####
rrc_8r <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac_8r(T = seq(10, 25, by = 0.1), 
                         ref_temp = 10,
                         r_EaN = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaP = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaC = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaD = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaE = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaF = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaG = sample_n(rgr_post_dist, size = 1)$intercept,
                         r_EaH = sample_n(rgr_post_dist, size = 1)$intercept,
                         c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea1C = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea1D = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea1E = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea1F = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea1G = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea1H = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                         c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2C = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2D = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2E = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2F = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2G = sample_n(c_post_dist, size = 1)$intercept, 
                         c_Ea2H = sample_n(c_post_dist, size = 1)$intercept, 
                         K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaC = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaD = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaE = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaF = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaG = sample_n(k_post_dist, size = 1)$intercept, 
                         K_EaH = sample_n(k_post_dist, size = 1)$intercept, 
                         v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaP = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaC = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaD = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaE = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaF = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaG = sample_n(v_post_dist, size = 1)$intercept,
                         v_EaH = sample_n(v_post_dist, size = 1)$intercept,
                         m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                         m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                         c1N_b = 0.5, c1P_b = 1,
                         c1C_b = 0.5, c1D_b = 1,
                         c1E_b = 0.5, c1F_b = 1,
                         c1G_b = 0.5, c1H_b = 1,
                         c2N_b = 1, c2P_b = 0.5,
                         c2C_b = 1, c2D_b = 0.5,
                         c2E_b = 1, c2F_b = 0.5,
                         c2G_b = 1, c2H_b = 0.5,
                         r_N_b = 1, r_P_b = 0.5, 
                         r_C_b = 1, r_D_b = 0.5,
                         r_E_b = 1, r_F_b = 0.5,
                         r_G_b = 1, r_H_b = 0.5,
                         K_N_b= 2000, K_P_b = 2000,
                         K_C_b = 2000, K_D_b = 2000,
                         K_E_b = 2000, K_F_b = 2000,
                         K_G_b = 2000, K_H_b = 2000,
                         v1N_b = 0.5, v1P_b = 1,
                         v1C_b = 0.5, v1D_b = 1,
                         v1E_b = 0.5, v1F_b = 1,
                         v1G_b = 0.5, v1H_b = 1,
                         v2N_b = 1, v2P_b = 0.5,
                         v2C_b = 1, v2D_b = 0.5,
                         v2E_b = 1, v2F_b = 0.5,
                         v2G_b = 1, v2H_b = 0.5,
                         m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  rrc_8r <- bind_rows(rrc_8r, hold) 
}

#get average change in position after 15C warming
rrc_8r_avg_new <- rrc_8r %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom8r <-
ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_8r, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_8r, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_8r, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_8r_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_8r_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_8r_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  # coord_cartesian(ylim = c(-0.022, 0.5), xlim = c(-0.022, 0.5)) +
  # scale_y_continuous(breaks = c(0, 0.25, 0.5)) +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  theme_cowplot(font_size = 20)

# get euclidean distances
rrc_8r_e <- rrc_8r %>% 
  filter(T %in% c(10, 25)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) 

hist(rrc_8r_e$dist15)

#histogram plot of euclidean distances in the pom pom plot
pom_hist8r <- rrc_8r_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n15°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)

# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_8r_p <- rrc_8r %>% 
  filter(T %in% c(10, 25)) %>% 
  mutate(temp = ifelse(T == 10, "Ambient", "+15°C Warming"))

rrc_8r_p_avg <- rrc_8r_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            med_fitrat = median(new_fit_ratio),
            med_stab_pot = median(new_stabil_potential),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift8r <-
  ggplot() + 
  geom_jitter(data = filter(rrc_8r_p, T>10), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_8r_p_avg, aes(x = temp, y = med_stab_pot, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.65)) 

#shift in fitness ratio
fd_shift8r <-
  ggplot() + 
  geom_jitter(data = filter(rrc_8r_p, T>10), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_8r_p_avg, aes(x = temp, y = med_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.65))

nd_shift8r + fd_shift8r

# pompom with subplots
bottom_patch8r <- pom_hist8r + nd_shift8r + fd_shift8r

comb_plot8r <- log_pom8r / bottom_patch8r + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")

# ggsave(plot = comb_plot8r, filename = "figures/8r_pom_hist_nfd.pdf", width = 12, height = 10)

#track the euclidean distance from neutrality for each species pair with warming
#get start position
start_new_stab_pot <- rrc_8r %>% 
  filter(T == 10) %>% 
  summarise(mean_start_nd = mean(new_stabil_potential),
            sd_stard_nd = sd(new_stabil_potential)) %>% 
  dplyr::select(mean_start_nd) %>% 
  unlist()

start_new_fit_rat <- rrc_8r %>% 
  filter(T == 10) %>% 
  summarise(mean_start_fd = mean(new_fit_ratio),
            sd_stard_fd = sd(new_fit_ratio)) %>% 
  dplyr::select(mean_start_fd) %>% 
  unlist()


track_8r <- rrc_8r%>% 
  dplyr::select(T, iteration, new_stabil_potential, new_fit_ratio, r_EaN:m_Ea2) %>%  #param combo identifies the simulation replicate
  mutate(start_nd = start_new_stab_pot,
         start_fd = start_new_fit_rat,
         ED = sqrt((new_stabil_potential - start_nd)^2 + (new_fit_ratio - start_fd)^2),
         dist_neut = sqrt(new_stabil_potential^2 + new_fit_ratio^2),
         iter_for_grouping = sprintf("%03d", iteration),
         iter_group = str_extract(as.character(iter_for_grouping), "^[0-9]", group = NULL))

track_8r %>% 
  ggplot() + 
  geom_line(aes(x = T, y = ED, colour = iteration, group = iteration)) 
# facet_wrap(~iter_group)

track_8r %>% 
  ggplot() + 
  geom_line(aes(x = T, y = dist_neut, colour = iteration, group = iteration)) + 
  ylab("Distance to neutrality") + 
  facet_wrap(~iter_group, scales = "free")


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
  
