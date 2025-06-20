# This script is to test effects of all parameters varying simultaneously with temperature, and investigate drivers of large changes in competition with warming

# Author: Kaleigh Davis, PDF University of Guelph
# Script DOB: 30 April 2025

#load packages
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

### NO THERMAL ASYMMETRIES SIM ####
### no tas  #####
nota <- data.frame()
for(f in 1:500){ #was 200
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), #was by 0.1
                      ref_temp = 25,
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.7, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  nota <- bind_rows(nota, hold) 
}

#get average change in position after 5, 10, 20C warming
nota_avg_new <- nota %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio))

#base pompom for comparison
# log_pom <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = nota, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(nota, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(nota, T==25), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-25), size = 6) +
  # position after 15C warming
  geom_point(data = nota_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = nota_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "Degrees \nC Warming") +
  # coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(0, 0.55)) +
  # coord_cartesian(ylim = c(-1, 1.5), xlim = c(-0.05, 1)) + 
  # scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.37, y = -0.05, label = "Co-existence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# get euclidean distances
nota_e <- nota %>% 
  filter(T %in% c(25, 125)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T125_new_stabil_potential - T25_new_stabil_potential)^2 + (T125_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T125_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T125_new_stabil_potential - T25_new_stabil_potential) 

hist(nota_e$dist15)

#how many paths cross 0 on the niche differences axis
nrow(nota_e %>% filter(T125_new_fit_ratio < 0) %>% distinct(iteration, .keep_all = T)) #150!

#histogram plot of euclidean dsistances in the pom pom plot
pom_hist <- nota_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n100C warming", y = "Count") + 
  theme_cowplot(font_size = 14)


# plot absolute shift in niche diffs and fitness diffs with warming #
nota_p <- nota %>% 
  filter(T %in% c(25, 125)) %>% 
  mutate(temp = ifelse(T == 25, "Ambient", "+100C Warming"))

nota_p_avg <- nota_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift <-
  ggplot() + 
  geom_point(data = nota_p, aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = nota_p_avg, aes(x = temp, y = mean_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+100C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.65)) 

#shift in fitness ratio
fd_shift <-
  ggplot() + 
  geom_point(data = nota_p, aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = nota_p_avg, aes(x = temp, y = mean_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+100C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.65))

nd_shift + fd_shift



##### POMPOM PLOT FOR MANUSCRIPT -- draw all param EAs at random ##############
### rrc  #####
rrc <- data.frame()
for(f in 1:500){ #was 200
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), #was by 0.1
                      ref_temp = 25,
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  rrc <- bind_rows(rrc, hold) 
}

#get average change in position after 5, 10, 20C warming
rrc_avg_new <- rrc %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio))

# # base pompom
# # pom <- 
#   ggplot() +
#   geom_path(data = rrc, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
#   geom_ribbon(data = data.frame(x = seq(0.01, 0.5, 0.001)),
#               aes(x = x,
#                   y = NULL,
#                   ymin = 1-x,
#                   ymax = 1/(1-x)),
#               fill = "grey", color = "black", alpha = 0.2) +
#   geom_point(data = filter(rrc, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "white", size = 5) +
#   geom_point(data = filter(rrc, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) + #N grows faster
#   geom_point(data = rrc_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 5) +
#   geom_point(data = rrc_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 4) +
#   geom_hline(yintercept = 1, linetype=5) +
#   scale_colour_viridis_c(option = "inferno") +
#   labs(colour = "Degrees warming") +
#   # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) + 
#   xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
#   ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

#base pompom for comparison
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
  geom_path(data = rrc, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc, T==25), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-25), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = rrc_avg_new, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 0.55)) +
  # coord_cartesian(ylim = c(-1, 1.5), xlim = c(-0.05, 1)) + 
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Co-existence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# ggsave(plot = log_pom, filename = "figures/kd-figs/log-pom1.pdf", width = 12, height = 10)

# get euclidean distances
rrc_e <- rrc %>% 
  filter(T %in% c(25, 40)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) 

hist(rrc_e$dist15)

#how many paths cross 0 on the niche differences axis
nrow(rrc_e %>% filter(T40_new_fit_ratio < 0) %>% distinct(iteration, .keep_all = T)) #2!

#histogram plot of euclidean dsistances in the pom pom plot
pom_hist <- rrc_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n15°C warming", y = "Count") + 
  theme_cowplot(font_size = 14)


# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_p <- rrc %>% 
  filter(T %in% c(25, 40)) %>% 
  mutate(temp = ifelse(T == 25, "Ambient", "+15°C Warming"))

rrc_p_avg <- rrc_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift <-
  ggplot() + 
  geom_jitter(data = filter(rrc_p, T>25), aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_p_avg, aes(x = temp, y = mean_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Niche differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.65)) 

#shift in fitness ratio
fd_shift <-
  ggplot() + 
  geom_jitter(data = filter(rrc_p, T>25), aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3, width = 0.03) +
  geom_point(data = rrc_p_avg, aes(x = temp, y = mean_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness differences"))) +
  scale_x_discrete(limits = c("Ambient", "+15°C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.65))

nd_shift + fd_shift

##### Pompom with subplots ######
# big pompom with small panels underneath
bottom_patch <- pom_hist + nd_shift + fd_shift

comb_plot1 <- log_pom / bottom_patch + 
  plot_layout(heights = c(2.25, 1)) + 
  plot_annotation(tag_levels = "A")
# ggsave(plot = comb_plot1, filename = "figures/kd-figs/pom_hist_nfd.pdf", width = 12, height = 10)
# ggsave(plot = comb_plot1, filename = "figures/kd-figs/100C_warm_pom_hist_nfd.pdf", width = 12, height = 10)


#try just a set of four equal panels in a square 
top_patch <- log_pom + pom_hist
bottom_patch1 <- nd_shift + fd_shift

top_patch/bottom_patch1

#use pompom plot to test for drivers of increase in ND #####
# get TAs for each parameter, at each iteration
rrc_parm_draws <- rrc %>% 
  group_by(iteration) %>% 
  distinct(r_EaN, r_EaP, c_Ea1N, c_Ea1P, c_Ea2N, c_Ea2P, K_EaN, K_EaP, v_EaN, v_EaP, m_Ea1, m_Ea2,
           c1N_b, c2P_b, c1P_b, c2N_b, r_N_b, r_P_b, K_N_b, K_P_b, v1N_b, v1P_b, v2N_b, v2P_b, m1_b,m2_b) %>% 
  ungroup()

rrc_tas <- rrc_parm_draws %>% 
  dplyr::select(iteration, r_EaN:m_Ea2) %>% 
  mutate(r_ta = r_EaN - r_EaP,
         c_ta1 = c_Ea1N - c_Ea1P,
         K_ta = K_EaN - K_EaP,
         v_ta = v_EaN - v_EaP,
         m_ta = m_Ea1 - m_Ea2) %>% 
  dplyr::select(iteration, r_ta:m_ta)

rrc_w_tas <- left_join(rrc, rrc_tas)

#plot
ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_w_tas, aes(x = new_stabil_potential, y = new_fit_ratio, color = as.factor(c_Ea2P > r_EaP), group = iteration), linewidth = 3, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype=5) +
  #aesthetic customization
  # scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(0, 0.55)) + 
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  facet_wrap(~new_stabil_potential >= 0.4054651) #this is just the start point

#have tried this with all of the temperature sensitive parameters (ta > 0) and some tests that make some intuitive sense to me as sources of categorical differences in outcomes (c_Ea > r_Ea)and there doesn't seem to be any clear patterns emerging.

# how many paths have increased ND after 15C warming? #
start_pos <- rrc_e %>% 
  dplyr::select(iteration, T25_new_fit_ratio, T25_new_stabil_potential) %>% 
  rename(start_fit_rat = T25_new_fit_ratio,
         start_stab_pot = T25_new_stabil_potential)

diffs_up_down <- rrc %>% 
  filter(T == 40) %>%
  dplyr::select(iteration, new_fit_ratio, new_stabil_potential) %>% 
  left_join(start_pos, .) %>% 
  mutate(nd_up_down = ifelse(new_stabil_potential > start_stab_pot, "up", "down"),
         fd_up_down = ifelse(new_fit_ratio > start_fit_rat, "up", "down"),
         both_up = ifelse(nd_up_down == "up" & fd_up_down == "up", "yes", "no"),
         both_down = ifelse(nd_up_down == "down" & fd_up_down == "down", "yes", "n"))
  
#possible pairs are: 
# FD UP, ND UP
# FD UP, ND DOWN
# FD DOWN, ND UP
# FD DOWN, ND DOWN

#both up
nrow(diffs_up_down %>% 
  filter(nd_up_down == "up" & fd_up_down == "up")) #88

#both down
nrow(diffs_up_down %>% 
       filter(nd_up_down == "down" & fd_up_down == "down")) #180

#nd up, fd down
nrow(diffs_up_down %>% 
       filter(nd_up_down == "up" & fd_up_down == "down")) #103

#nd down, fd up
nrow(diffs_up_down %>% 
       filter(nd_up_down == "down" & fd_up_down == "up")) #129

#these total 500 -- good

# can we get the same results if we just allow r and c to vary with temperature? #####
### rrc_rc  #####
rrc_rc <- data.frame()
for(f in 1:500){ #was 200
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), #was by 0.1
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 0.5, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  rrc_rc <- bind_rows(rrc_rc, hold) 
}

#get average change in position after 5, 10, 20C warming
rrc_rc_avg_new <- rrc_rc %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio))

# log_pom_rc <-
ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_rc, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_rc, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_rc, T==25), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-25), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_rc_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_rc_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "Degrees \nC Warming") +
  # coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(0, 0.55)) + 
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.37, y = -0.05, label = "Co-existence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# ggsave(plot = log_pom, filename = "figures/kd-figs/log-pom1.pdf", width = 12, height = 10)

# get euclidean distances
rrc_rc_e <- rrc_rc %>% 
  filter(T %in% c(25, 40)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) 

hist(rrc_rc_e$dist15)

#histogram plot of euclidean dsistances in the pom pom plot
pom_hist_rc <- rrc_rc_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n15C warming", y = "Count") + 
  theme_cowplot(font_size = 14)


# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_rc_p <- rrc_rc %>% 
  filter(T %in% c(25, 40)) %>% 
  mutate(temp = ifelse(T == 25, "Ambient", "+15C Warming"))

rrc_rc_p_avg <- rrc_rc_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift_rc <-
  ggplot() + 
  geom_point(data = rrc_rc_p, aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = rrc_rc_p_avg, aes(x = temp, y = mean_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Stabilization potential"))) +
  scale_x_discrete(limits = c("Ambient", "+15C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.65)) 

#shift in fitness ratio
fd_shift_rc <-
  ggplot() + 
  geom_point(data = rrc_rc_p, aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = rrc_rc_p_avg, aes(x = temp, y = mean_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness difference"))) +
  scale_x_discrete(limits = c("Ambient", "+15C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.65))


##### Pompom with subplots ######
# big pompom with small panels underneath
bottom_patch <- pom_hist_rc + nd_shift_rc + fd_shift_rc

comb_plot1 <- log_pom_rc / bottom_patch + 
  plot_layout(heights = c(2.25, 1))
# ggsave(plot = comb_plot1, filename = "figures/kd-figs/pom_hist_nfd.pdf", width = 12, height = 10)

#try just a set of four equal panels in a square 
top_patch <- log_pom + pom_hist
bottom_patch1 <- nd_shift + fd_shift

top_patch/bottom_patch1

# can we get the same results if we just allow r and c to vary with temperature? #####
### rrc_rc1  #####
rrc_rc1 <- data.frame()
for(f in 1:500){ #was 200
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), #was by 0.1
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept, 
                      c_Ea2N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2P = sample_n(c_post_dist, size = 1)$intercept, 
                      K_EaN = mean(k_post_dist$intercept), 
                      K_EaP = mean(k_post_dist$intercept), 
                      v_EaN = mean(v_post_dist$intercept),
                      v_EaP = mean(v_post_dist$intercept), 
                      m_Ea1 = mean(m_post_dist$intercept), 
                      m_Ea2 = mean(m_post_dist$intercept),
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 0.5, r_P_b = 0.5, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  rrc_rc1 <- bind_rows(rrc_rc1, hold) 
}

#get average change in position after 5, 10, 20C warming
rrc_rc1_avg_new <- rrc_rc1 %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio))

log_pom_rc1 <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrc_rc1, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrc_rc1, T==25), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrc_rc1, T==25), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-25), size = 6) +
  # position after 15C warming
  geom_point(data = rrc_rc1_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = rrc_rc1_avg_new, aes(x = new_mean_stab_pot, y = new_mean_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "Degrees \nC Warming") +
  # coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(0, 0.55)) + 
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.37, y = -0.05, label = "Co-existence", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 5, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 5, fontface = 2) +
  theme_cowplot(font_size = 20)

# ggsave(plot = log_pom, filename = "figures/kd-figs/log-pom1.pdf", width = 12, height = 10)

# get euclidean distances
rrc_rc1_e <- rrc_rc1 %>% 
  filter(T %in% c(25, 40)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist15 = sqrt((T40_new_stabil_potential - T25_new_stabil_potential)^2 + (T40_new_fit_ratio - T25_new_fit_ratio)^2),
         shift_fitrat = T40_new_fit_ratio - T25_new_fit_ratio,
         shift_nichediffs = T40_new_stabil_potential - T25_new_stabil_potential) 

hist(rrc_rc1_e$dist15)

#histogram plot of euclidean dsistances in the pom pom plot
pom_hist_rc1 <- rrc_rc1_e %>% 
  ggplot(aes(x = dist15)) + 
  geom_histogram(binwidth = 0.05, colour = "black") + 
  labs(x = "Euclidean distance with \n15C warming", y = "Count") + 
  theme_cowplot(font_size = 14)


# plot absolute shift in niche diffs and fitness diffs with warming #
rrc_rc1_p <- rrc_rc1 %>% 
  filter(T %in% c(25, 40)) %>% 
  mutate(temp = ifelse(T == 25, "Ambient", "+15C Warming"))

rrc_rc1_p_avg <- rrc_rc1_p %>% 
  group_by(temp) %>% 
  summarise(mean_stabil_potential = mean(new_stabil_potential), 
            mean_fitrat = mean(new_fit_ratio),
            sd_stabil_potential = sd(new_stabil_potential),
            sd_fitrat = sd(new_fit_ratio))

#shift in stabilization potential
nd_shift_rc1 <-
  ggplot() + 
  geom_point(data = rrc_rc1_p, aes(x = temp, y = new_stabil_potential), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = rrc_rc1_p_avg, aes(x = temp, y = mean_stabil_potential, fill = temp), size = 5, pch = 21) +
  labs(x = "Temperature", y = expression(paste("Stabilization potential"))) +
  scale_x_discrete(limits = c("Ambient", "+15C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.65)) 

#shift in fitness ratio
fd_shift_rc1 <-
  ggplot() + 
  geom_point(data = rrc_rc1_p, aes(x = temp, y = new_fit_ratio), colour = "lightgrey", alpha = 0.3) +
  geom_point(data = rrc_rc1_p_avg, aes(x = temp, y = mean_fitrat, fill = temp), size = 5, pch = 21) + 
  labs(x = "Temperature", y = expression(paste("Fitness difference"))) +
  scale_x_discrete(limits = c("Ambient", "+15C Warming")) + 
  scale_fill_manual(values = c("#C23A75", "#FBFCBE")) +
  theme_cowplot(font_size = 14) + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.65))


##### Pompom with subplots ######
# big pompom with small panels underneath
bottom_patch <- pom_hist_rc1 + nd_shift_rc1 + fd_shift_rc1

comb_plot1 <- log_pom_rc1 / bottom_patch + 
  plot_layout(heights = c(2.25, 1))
# ggsave(plot = comb_plot1, filename = "figures/kd-figs/pom_hist_nfd.pdf", width = 12, height = 10)

#try just a set of four equal panels in a square 
top_patch <- log_pom + pom_hist
bottom_patch1 <- nd_shift + fd_shift

top_patch/bottom_patch1


## Now I want to just save the rrc_avg_new dataset for each iterations and jam it onto the last iteration, so I can plot the start point and the distribution of average endpoints ####
### FIRST, WITH FOR LOOP ####
# Initialize an empty list to store results
results_list <- list()

# Loop through 15 iterations
for (i in 1:25) {
  print(paste("Iteration", i, Sys.time()))
  # Initialize an empty data frame to collect results for this iteration
  rrc <- data.frame()
  
  # Simulation loop (your existing code)
  for (f in 1:500) { 
    hold <- temp_dep_mac(
      T = seq(25, 40, by = 0.1),
      ref_temp = 25,
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
      r_N_b = 0.5, r_P_b = 1,
      K_N_b = 2000, K_P_b = 2000,
      v1N_b = 0.5, v1P_b = 1,
      v2N_b = 1, v2P_b = 0.5,
      m1_b = 0.01, m2_b = 0.01
    )
    
    hold$iteration <- f
    rrc <- bind_rows(rrc, hold)
  }
  
  # Calculate average change in position after 5, 10, 20C warming
  rrc_avg_new <- rrc %>%
    mutate(rel_T = T - 25) %>%
    filter(rel_T %in% c(0, 15)) %>%
    group_by(rel_T) %>%
    summarise(new_mean_stab_pot = mean(new_stabil_potential),
              new_mean_fit_rat = mean(new_fit_ratio))
  
  # Add iteration number to rrc_avg_new
  rrc_avg_new$iteration <- i
  
  # Store results in the list
  results_list[[i]] <- rrc_avg_new
  
  print(paste("Iteration", i, "complete", Sys.time()))
  print(length(results_list))
  
}

# Combine all results into a single data frame
final_results <- bind_rows(results_list)

fr <- final_results %>% 
  distinct(new_mean_stab_pot, new_mean_fit_rat, .keep_all = T)

fr %>% 
  ggplot() + 
  geom_ribbon(data = data.frame(x = seq(0, 0.6, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(aes(x = new_mean_stab_pot, y = new_mean_fit_rat, colour = as.factor(rel_T)), size = 3, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 5) +
  coord_cartesian(ylim=c(-0.5, 0), xlim = c(0.25, 0.5)) + 
  theme(legend.position = "none")

### NEXT, WITH PURRR, FURRR ####
# Define a function for running your simulation and calculating averages
run_simulation <- function(iteration) {
  rrc <- data.frame()
  
  for (f in 1:500) { 
    hold <- temp_dep_mac(
      T = seq(25, 40, by = 0.1),
      ref_temp = 25,
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
      r_N_b = 0.5, r_P_b = 1,
      K_N_b = 2000, K_P_b = 2000,
      v1N_b = 0.5, v1P_b = 1,
      v2N_b = 1, v2P_b = 0.5,
      m1_b = 0.01, m2_b = 0.01
    )
    
    hold$iteration <- f
    rrc <- bind_rows(rrc, hold)
  }
  
  # Calculate average change in position after 5, 10, 20C warming
  rrc_avg_new <- rrc %>%
    mutate(rel_T = T - 25) %>%
    filter(rel_T %in% c(0, 15)) %>%
    group_by(rel_T) %>% 
    summarise(new_mean_stab_pot = mean(new_stabil_potential),
              new_mean_fit_rat = mean(new_fit_ratio))
  
  # Add iteration number to rrc_avg_new
  rrc_avg_new$iteration <- iteration
  
  # Print timestamp and iteration info
  print(paste("Iteration", iteration, "complete at", Sys.time()))
  
  return(rrc_avg_new)
}

# # Use purrr::map to run simulations iteratively and combine results
# results_list <- map_df(1:4, run_simulation) #This took 7 minutes
# 
# # Print final results (optional)
# print(results_list)

# try running the sims in parallel across cores
plan(multisession, workers = 4)  # Set the number of workers for parallel processing

# Use furrr::future_map_df for parallel processing
results_list <- future_map(1:50, run_simulation) #takes just over 2 minutes

beep(2)

#switch back to sequential 
plan(sequential)

results_list1 <- bind_rows(results_list)

results_list1 %>% 
  group_by(new_mean_stab_pot, new_mean_fit_rat) %>% 
  tally() %>% 
  filter(n > 1) #just 50 copies of start point -- great

results_list2 <- results_list1 %>% 
  distinct(new_mean_stab_pot, new_mean_fit_rat, .keep_all = T)

avg_move <- results_list2 %>% 
  ggplot() + 
  geom_ribbon(data = data.frame(x = seq(0, 0.6, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(aes(x = new_mean_stab_pot, y = new_mean_fit_rat, colour = as.factor(rel_T)), size = 6, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 5) +
  coord_cartesian(ylim=c(-0.45, -0.25), xlim = c(0.30, 0.42)) + 
  theme(legend.position = "none") + 
  xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  labs(title = "Average movement of species pair with 15C warming, \nacross 50 simulations",
       caption = "Each blue point (n = 50) represents the average movement \nin MCT space of 500 species pairs, drawn from the empirical distribution.") + 
  theme_cowplot(font_size = 18)

# ggsave(avg_move, filename = "figures/kd-figs/avg-movement.pdf", width = 12, height = 10)

#just jack up the iterations and plot all the end points ###
rrc_ep <- data.frame()

for(f in 1:5000){ #was 200
  hold = temp_dep_mac(T = seq(25, 40, by = 0.1), #was by 0.1
                      ref_temp = 25,
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
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 0.2, 0.4
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 0.4, 0.2
                      r_N_b = 0.5, r_P_b = 1, #growth rate for each resource at ref temp 0.1, 0.1
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 2000, 2000
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 0.2, 0.4
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 0.4, 0.2
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m 0.1, 0.1
  hold$iteration <- f
  rrc_ep <- bind_rows(rrc_ep, hold) 
}

ribbon <- geom_ribbon(data = data.frame(x = seq(0, 0.6, 0.001)),
                      aes(x = x,
                          y = NULL,
                          ymin = -x,
                          ymax = x),
                      fill = "grey", color = "black", alpha = 0.2)

rrc_ep %>% 
  filter(T %in% c(25, 40)) %>%
  ggplot() + 
  # geom_path(data = rrc_ep, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-25, group = iteration), linewidth = 3) +
  geom_point(aes(x = new_stabil_potential, y = new_fit_ratio, colour = T), alpha = 0.2, size = 2) +
  ribbon +
  geom_hline(yintercept = 0, linetype = 5) +
  # coord_cartesian(ylim=c(-0.45, -0.25), xlim = c(0.30, 0.42)) + 
  theme(legend.position = "none") + 
  xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) + 
  theme_cowplot(font_size = 18)

