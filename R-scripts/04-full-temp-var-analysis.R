# This script is to test effects of all parameters varying simultaneously with temperature, and investigate drivers of large changes in competition with warming

# Author: Kaleigh Davis, PDF University of Guelph
# Script DOB: 30 April 2025

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
source("R-scripts/02-temp-dep-macarthur.R") 
#this contains the macarthur translation function, with all parameters flexibly defined in the function for assigning at time of use, and the arrhenius function.

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

##### draw all param EAs at random (rrc) ##############
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

#get average change in position after 5, 10, 20C warming
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
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# get euclidean distances for each species pair
rrc_e <- rrc %>% 
  filter(T %in% c(10, 25)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) 

hist(rrc_e$dist15)

#histogram plot of euclidean distances in the pom pom plot
pom_hist <- rrc_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
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

# ggsave(plot = comb_plot1, filename = "figures/pom_hist_nfd.pdf", width = 12, height = 10)

################################################################################
#########################  SUPPLEMENTARY ANALYSES   ############################
################################################################################

# Repeat analysis with 50C warming (rrc) - Figure S5 #####
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

# ggsave(plot = comb_plot50, filename = "figures/50C_warm_pom_hist_nfd.pdf", width = 12, height = 10)

# Repeat analysis with 5C warming (rrc) - Figure S1 ####
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

# ggsave(plot = comb_plot5, filename = "figures/5C_warm_pom_hist_nfd.pdf", width = 12, height = 10)

# No thermal asymmetries -- Figure S7 ####
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
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio))

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
  geom_point(data = nota_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat), colour = "black",  size = 7) +
  geom_point(data = nota_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat, colour = rel_T),  size = 5.5) +
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
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift_nota <-
  ggplot() + 
  geom_point(data = nota_p, aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = nota_p_avg, aes(x = temp, y = mean_stabil_potential, fill = temp), size = 5, pch = 21) +
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
  geom_point(data = nota_p_avg, aes(x = temp, y = mean_fitrat, fill = temp), size = 5, pch = 21) + 
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
# ggsave(plot = comb_plot_nota, filename = "figures/NOTA_pom_hist_nfd.pdf", width = 12, height = 10)

