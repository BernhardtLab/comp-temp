#this script is to play around with some specific biological scenarios, and see how coex/compex might play out with warming in those scenarios

#script DOB: Feb 27, 2025
#author: Kaleigh Davis, University of Guelph postdoc
#parent script here is kd-comptemp-scenarios.R


### Perspective and setup ######
# "Starting point," or the coordinates at the reference temperature, before any warming or cooling occurs, makes a huge difference in competitive outcomes. As such, we will begin with a set of biological scenarios that impose different starting positions.
# Then, we will expose all of those to warming, both symmetrical warming, and asymmetrical warming, according to biologicall relevant warming scenarios, and see what effects this has on competitive outcomes. 

#From old script #
#Questions I'd like to answer for each scenario:
## 1. How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming
####### In order to answer this I need to calculate the coordinates of these lines? Or can I just set up the inequality and see who satisfies is at these different time points?
## 2. Do niche or fitness differences change more with warming?
####### Set up a scaled position axis and see what values are after 25 deg warming or cooling along each axis
## 3. Are there more lines that move net W (decrease in stab potential) or net E (increase stab potential) after 20 degrees of warming; net N vs net S
######## Count number of species pairs that end up more than and less than the starting stab potential and fitness ratio values

#Other questions to explore:
# max distance at max warming as a function of stab potential or as a function of start point euclidean dist from origin

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

#get quantiles for all parameters ########
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

param_sum1 <- param_sum %>% 
  pivot_longer(cols = c(Mean:Max), 
               names_to = "summary_stat",
               values_to = "value")

#split df up into dfs for each summary statistic
param_sum1 %>% 
  group_by(summary_stat) %>% 
  group_split() %>% 
  purrr::walk(~ assign(paste0(.x$summary_stat[1]), .x, envir = .GlobalEnv))

# establishing start point scenarios #######

##### change magnitude of difference between consumption vectors to see if I can move the starting point along the threshold ####
#there should be an analytical way of identifying pairs of c_b values that will fall along the barrier line -- I need to figure this out

# Species are exactly on borderline of coexistence/competitive exclusion #
#Only difference is N has faster growth rate at reference temperature
df1 <- temp_dep_mac(T = 25, #was by 0.1
                    ref_temp = 25,
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
                    m_Ea1 = 0, 
                    m_Ea2 = 0,
                    c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                    c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                    r_N_b = 0.1, r_P_b = 0.05, #growth rate for each resource at ref temp; changing the relationship between these two moves from upper line to center to lower line
                    K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                    v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                    v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                    m1_b = 0.01, m2_b = 0.01) %>% 
  mutate(scenario = "rrc",
         sub_scenario = "N > P")
 
#same, but switch which resource grows faster (Now P)
df2 <- temp_dep_mac(T = 25,
                    ref_temp = 25,
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
                    m_Ea1 = 0, 
                    m_Ea2 = 0,
                    c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                    c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                    r_N_b = 0.05, r_P_b = 0.1, 
                    K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                    v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                    v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                    m1_b = 0.01, m2_b = 0.01) %>% 
  mutate(scenario = "rrc",
         sub_scenario = "P > N")

#both resources grow at the same rate
df3 <- temp_dep_mac(T = 25,
                    ref_temp = 25,
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
                    m_Ea1 = 0, 
                    m_Ea2 = 0,
                    c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                    c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                    r_N_b = 0.05, r_P_b = 0.05, 
                    K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                    v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                    v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                    m1_b = 0.01, m2_b = 0.01) %>% 
  mutate(scenario = "rrc",
         sub_scenario = "N==P")

# generalist-specialist #
df4 <- temp_dep_mac(T = 25, #was by 0.1
                    ref_temp = 25,
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
                    m_Ea1 = 0, 
                    m_Ea2 = 0,
                    c1N_b = 0.6, c1P_b = 0.4, #spec 1 consumes more P
                    c2N_b = 0.5, c2P_b = 0.5, #spec 2 consumes more N
                    r_N_b = 0.05, r_P_b = 0.05, 
                    K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                    v1N_b = 0.6, v1P_b = 0.4, #sp 1 converts P more efficiently
                    v2N_b = 0.5, v2P_b = 0.5, #sp 2 converts N more efficiently
                    m1_b = 0.01, m2_b = 0.01) %>% 
  mutate(scenario = "gen-spec",
         sub_scenario = "modest")

#slightly more extreme generalist-specialist trade off
df5 <- temp_dep_mac(T = 25,
                    ref_temp = 25,
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
                    m_Ea1 = 0, 
                    m_Ea2 = 0,
                    c1N_b = 0.7, c1P_b = 0.3, #spec 1 consumes more N
                    c2N_b = 0.5, c2P_b = 0.5, #spec 2 consumes both, equally
                    r_N_b = 0.05, r_P_b = 0.05, 
                    K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                    v1N_b = 0.7, v1P_b = 0.3, #sp 1 converts N more efficiently
                    v2N_b = 0.5, v2P_b = 0.5, #sp 2 converts both, equally
                    m1_b = 0.01, m2_b = 0.01) %>% 
  mutate(scenario = "gen-spec",
         sub_scenario = "extreme")

#see if I can hack the starting point to the 0.3, 0.3 range (== the average across the dataset in Buche et al 2022)
#same, but switch which resource grows faster (Now P)
df6 <- temp_dep_mac(T = 25,
                    ref_temp = 25,
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
                    m_Ea1 = 0, 
                    m_Ea2 = 0,
                    c1N_b = 0.1, c1P_b = 0.22, # moves mostly left/right
                    c2N_b = 0.5, c2P_b = 0.4, #
                    r_N_b = 0.05, r_P_b = 0.2, #increasing this difference drives the point diagonally down and left or up and right
                    K_N_b= 5000, K_P_b = 2000, #
                    v1N_b = 0.1, v1P_b = 0.22, #
                    v2N_b = 0.5, v2P_b = 0.005, #
                    m1_b = 0.01, m2_b = 50) %>% # starting point is extremely insensitive to changes in mortality rate
  mutate(scenario = "average-sp-pair",
         sub_scenario = "none")

#glue them all together
starts <- bind_rows(df1, df2, df3, df4, df5, df6)

#distribution of points at ref temp
starts %>% 
  filter(!str_detect(scenario, "average")) %>% 
  ggplot(aes(x = stabil_potential, y = fit_ratio)) + 
  geom_ribbon(data = data.frame(x = seq(0.005, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)), #can log these in order to make the coexistence zone more symmetrical
              fill = "darkgrey", color = "black", alpha = 0.2) +
  geom_point(aes(shape = scenario, colour = sub_scenario), size = 4) +
  geom_hline(yintercept = 1, linetype=5) +
  geom_point(x = 0.39, y = 0.35, size = 4, shape = 21) + 
  # coord_cartesian(ylim = c(0.2, 2)) +
  # guides(colour = "none") +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
  scale_x_log10() +
  scale_y_log10()

#get euclidean distance of start points to "origin"
starts_e <- starts %>% 
  dplyr::select(scenario, sub_scenario, c1N_b:m2_b, stabil_potential, fit_ratio) %>%
  mutate(start_dist_origin = sqrt((stabil_potential - 0)^2 + (fit_ratio - 1)^2)) %>% 
  dplyr::select(scenario, stabil_potential, fit_ratio, sub_scenario, start_dist_origin) %>% 
  rename(start_stabil_potential = stabil_potential,
         start_fit_ratio = fit_ratio)
#need to bind these onto the right hand side of the final dataset? Or can cram all the simulations together into a massive df and make a column for scenario and then calculate start euclidean distance there.


############################################################################################
################## RECIPROCAL RESOURCE CONSUMPTION AT T_REF SCENARIOS ######################
############################################################################################

##### draw all params at random ##############
### rrc equal base rates #####
rrc_e <- data.frame()
for(f in 1:200){ #was 200
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), #was by 0.1
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
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.05, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  rrc_e <- bind_rows(rrc_e, hold) 
}

rrc_e %>% filter(rho == 1)
hist(rrc_e$rho)

rrc_e <- rrc_e %>% 
  mutate(scenario = "rrc",
         sub_scenario = "N==P")

#get average change in position after 5, 10, 20C warming
# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
rrc_e_avg <- rrc_e %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#base pompom for comparison
ggplot() +
  geom_path(data = rrc_e, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0.01, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = filter(rrc_e, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = filter(rrc_e, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4, shape = "triangle") + #N grows faster
  geom_point(data = filter(rrc_e, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "turquoise", size = 3, shape = "triangle") + #N grows faster
  geom_point(data = rrc_e_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = rrc_e_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 3) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
  scale_x_log10() +
  scale_y_log10()

# QUESTION: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
rrc_e %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(5, 10, 15, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally()

### rrc N faster ########
rrc_n <- data.frame()
for(f in 1:200){ #was 200
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), #was by 0.1
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
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  rrc_n <- bind_rows(rrc_n, hold) 
}

rrc_n %>% filter(rho == 1)
hist(rrc_n$rho)

rrc_n <- rrc_n %>% 
  mutate(scenario = "rrc",
         sub_scenario = "N > P")

#get average change in position after 5, 10, 20C warming
# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
rrc_n_avg <- rrc_n %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#base pompom for comparison
ggplot() +
  geom_path(data = rrc_n, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0.01, 0.6, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = filter(rrc_n, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = filter(rrc_n, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4, shape = "triangle") + #N grows faster
  geom_point(data = filter(rrc_n, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "turquoise", size = 3, shape = "triangle") + #N grows faster
  geom_point(data = rrc_n_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = rrc_n_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 3) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
  scale_y_log10() + 
  scale_x_log10()

# QUESTION: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
rrc_n %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(5, 10, 15, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally()

### rrc P faster ########
rrc_p <- data.frame()
for(f in 1:200){ #was 200
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), #was by 0.1
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
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.05, r_P_b = 0.1, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  rrc_p <- bind_rows(rrc_p, hold) 
}

rrc_p %>% filter(rho == 1)
hist(rrc_p$rho)

rrc_p <- rrc_p %>% 
  mutate(scenario = "rrc",
         sub_scenario = "P > N")

#get average change in position after 5, 10, 20C warming
# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
rrc_p_avg <- rrc_p %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#base pompom for comparison
ggplot() +
  geom_path(data = rrc_p, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0.01, 0.6, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = filter(rrc_p, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = filter(rrc_p, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4, shape = "triangle") + #N grows faster
  geom_point(data = filter(rrc_p, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "turquoise", size = 3, shape = "triangle") + #N grows faster
  geom_point(data = rrc_p_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = rrc_p_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 3) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
  scale_y_log10() + 
  scale_x_log10()

# QUESTION: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
rrc_p %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(5, 10, 15, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally()


##### all params highly temperature sensitive ##############
#I've implemented this by drawing all EAs from the top quantile of their empirical distributions
jb1 <- data.frame()
for(f in 1:200){ #was 200
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), #was by 0.1
                      ref_temp = 25,
                      r_EaN = sample_n(filter(rgr_post_dist, intercept > filter(Q3, parameter == "resource_growth_rate")$value), size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(filter(rgr_post_dist, intercept > filter(Q3, parameter == "resource_growth_rate")$value), size = 1)$intercept, 
                      c_Ea1N = sample_n(filter(c_post_dist, intercept > filter(Q3, parameter == "consumption rate")$value), size = 1)$intercept,
                      c_Ea1P = sample_n(filter(c_post_dist, intercept > filter(Q3, parameter == "consumption rate")$value), size = 1)$intercept, 
                      c_Ea2N = sample_n(filter(c_post_dist, intercept > filter(Q3, parameter == "consumption rate")$value), size = 1)$intercept,
                      c_Ea2P = sample_n(filter(c_post_dist, intercept > filter(Q3, parameter == "consumption rate")$value), size = 1)$intercept, 
                      K_EaN = sample_n(filter(k_post_dist, intercept > filter(Q3, parameter == "carrying_capacity")$value), size = 1)$intercept, 
                      K_EaP = sample_n(filter(k_post_dist, intercept > filter(Q3, parameter == "carrying_capacity")$value), size = 1)$intercept, 
                      v_EaN = sample_n(filter(v_post_dist, intercept > filter(Q3, parameter == "conversion_efficiency")$value), size = 1)$intercept,
                      v_EaP = sample_n(filter(v_post_dist, intercept > filter(Q3, parameter == "conversion_efficiency")$value), size = 1)$intercept, 
                      m_Ea1 = sample_n(filter(m_post_dist, intercept > filter(Q3, parameter == "mortality_rate")$value), size = 1)$intercept, 
                      m_Ea2 = sample_n(filter(m_post_dist, intercept > filter(Q3, parameter == "mortality_rate")$value), size = 1)$intercept,
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  jb1 <- bind_rows(jb1, hold) 
}

#tassel plot
ggplot() +
  geom_path(data = jb1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(jb1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

##### all params minimally temperature sensitive ##############
#I've implemented this by drawing all EAs from the bottom quantile of their empirical distributions
jb2 <- data.frame()
for(f in 1:200){ #was 200
  hold = temp_dep_mac(T = seq(0, 50, by = 0.1), #was by 0.1
                      ref_temp = 25,
                      r_EaN = sample_n(filter(rgr_post_dist, intercept < filter(Q1, parameter == "resource_growth_rate")$value), size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(filter(rgr_post_dist, intercept < filter(Q1, parameter == "resource_growth_rate")$value), size = 1)$intercept, 
                      c_Ea1N = sample_n(filter(c_post_dist, intercept < filter(Q1, parameter == "consumption rate")$value), size = 1)$intercept,
                      c_Ea1P = sample_n(filter(c_post_dist, intercept < filter(Q1, parameter == "consumption rate")$value), size = 1)$intercept, 
                      c_Ea2N = sample_n(filter(c_post_dist, intercept < filter(Q1, parameter == "consumption rate")$value), size = 1)$intercept,
                      c_Ea2P = sample_n(filter(c_post_dist, intercept < filter(Q1, parameter == "consumption rate")$value), size = 1)$intercept, 
                      K_EaN = sample_n(filter(k_post_dist, intercept < filter(Q1, parameter == "carrying_capacity")$value), size = 1)$intercept, 
                      K_EaP = sample_n(filter(k_post_dist, intercept < filter(Q1, parameter == "carrying_capacity")$value), size = 1)$intercept, 
                      v_EaN = sample_n(filter(v_post_dist, intercept < filter(Q1, parameter == "conversion_efficiency")$value), size = 1)$intercept,
                      v_EaP = sample_n(filter(v_post_dist, intercept < filter(Q1, parameter == "conversion_efficiency")$value), size = 1)$intercept, 
                      m_Ea1 = sample_n(filter(m_post_dist, intercept < filter(Q1, parameter == "mortality_rate")$value), size = 1)$intercept, 
                      m_Ea2 = sample_n(filter(m_post_dist, intercept < filter(Q1, parameter == "mortality_rate")$value), size = 1)$intercept,
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  jb2 <- bind_rows(jb2, hold) 
}

#tassel plot
ggplot() +
  geom_path(data = jb2, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(jb2, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

##### consumers have high EAs, resources have low EAs ##############
jb3 <- data.frame()
for(f in 1:200){ #was 200
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), #was by 0.1
                      ref_temp = 25,
                      r_EaN = sample_n(filter(rgr_post_dist, intercept < filter(Q1, parameter == "resource_growth_rate")$value), size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(filter(rgr_post_dist, intercept < filter(Q1, parameter == "resource_growth_rate")$value), size = 1)$intercept, 
                      c_Ea1N = sample_n(filter(c_post_dist, intercept > filter(Q3, parameter == "consumption rate")$value), size = 1)$intercept,
                      c_Ea1P = sample_n(filter(c_post_dist, intercept > filter(Q3, parameter == "consumption rate")$value), size = 1)$intercept, 
                      c_Ea2N = sample_n(filter(c_post_dist, intercept > filter(Q3, parameter == "consumption rate")$value), size = 1)$intercept,
                      c_Ea2P = sample_n(filter(c_post_dist, intercept > filter(Q3, parameter == "consumption rate")$value), size = 1)$intercept, 
                      K_EaN = sample_n(filter(k_post_dist, intercept < filter(Q1, parameter == "carrying_capacity")$value), size = 1)$intercept, 
                      K_EaP = sample_n(filter(k_post_dist, intercept < filter(Q1, parameter == "carrying_capacity")$value), size = 1)$intercept, 
                      v_EaN = sample_n(filter(v_post_dist, intercept < filter(Q1, parameter == "conversion_efficiency")$value), size = 1)$intercept,
                      v_EaP = sample_n(filter(v_post_dist, intercept < filter(Q1, parameter == "conversion_efficiency")$value), size = 1)$intercept, 
                      m_Ea1 = sample_n(filter(m_post_dist, intercept > filter(Q3, parameter == "mortality_rate")$value), size = 1)$intercept, 
                      m_Ea2 = sample_n(filter(m_post_dist, intercept > filter(Q3, parameter == "mortality_rate")$value), size = 1)$intercept,
                      c1N_b = 0.2, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.4, c2P_b = 0.2, #spec 2 consumes more N
                      r_N_b = 0.1, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.2, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.4, v2P_b = 0.2, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  jb3 <- bind_rows(jb3, hold) 
}

#tassel plot
ggplot() +
  geom_path(data = jb3, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(jb3, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

############################################################################################
################## GENERALIST-SPECIALIST TRADE-OFF AT T_REF SCENARIOS ######################
############################################################################################

##### Draw all param EAs from distribution ##############
# extreme generalist-specialist trade-off ####
gs_ext <- data.frame()
for(f in 1:200){ #was 200
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), #was by 0.1
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
                      c1N_b = 0.7, c1P_b = 0.3, #spec 1 consumes more P
                      c2N_b = 0.5, c2P_b = 0.5, #spec 2 consumes more N
                      r_N_b = 0.05, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.7, v1P_b = 0.3, #sp 1 converts P more efficiently
                      v2N_b = 0.5, v2P_b = 0.5, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  gs_ext <- bind_rows(gs_ext, hold) 
}

gs_ext <- gs_ext %>% 
  mutate(scenario = "gen-spec", 
         sub_scenario = "extreme")

#get average change in position after 5, 10, 20C warming
# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
gs_ext_avg <- gs_ext %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#pompom plot
ggplot() +
  geom_path(data = gs_ext, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0.0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = filter(gs_ext, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = filter(gs_ext, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4, shape = "square") +
  geom_point(data = filter(gs_ext, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "salmon", size = 3, shape = "square") + 
  geom_point(data = gs_ext_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = gs_ext_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 3) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
  scale_x_log10() + 
  scale_y_log10()

#moderate generalist-specialist trade-off #####
gs_mod <- data.frame()
for(f in 1:200){ #was 200
  hold = temp_dep_mac(T = seq(25, 50, by = 0.1), #was by 0.1
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
                      c1N_b = 0.6, c1P_b = 0.4, #spec 1 consumes more P
                      c2N_b = 0.5, c2P_b = 0.5, #spec 2 consumes more N
                      r_N_b = 0.05, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.6, v1P_b = 0.4, #sp 1 converts P more efficiently
                      v2N_b = 0.5, v2P_b = 0.5, #sp 2 converts N more efficiently
                      m1_b = 0.01, m2_b = 0.01) #same for both species; model v insensitive to changes in m
  hold$iteration <- f
  gs_mod <- bind_rows(gs_mod, hold) 
}

gs_mod <- gs_mod %>% 
  mutate(scenario = "gen-spec", 
         sub_scenario = "modest")

#get average change in position after 5, 10, 20C warming
# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
gs_mod_avg <- gs_mod %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#pompom plot
ggplot() +
  geom_path(data = gs_mod, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0.0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = filter(gs_mod, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = filter(gs_mod, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4, shape = "square") +
  geom_point(data = filter(gs_mod, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "salmon", size = 3, shape = "square") + 
  geom_point(data = gs_mod_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = gs_mod_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 3) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
  scale_x_log10() +
  scale_y_log10()

# QUESTION: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
gs_mod %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(5, 10, 15, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally()

#Calculate Euclidean distance at 5C, 10C, 20C #####
rrc_n_e <- rrc_n %>% 
  filter(T %in% c(25, 30, 35, 45)) %>% 
  dplyr::select(-c(new_stabil_potential, a11:g2, new_fit_ratio:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(stabil_potential, fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist5 = sqrt((T30_stabil_potential - T25_stabil_potential)^2 + (T30_fit_ratio - T25_fit_ratio)^2),
         dist10 = sqrt((T35_stabil_potential - T25_stabil_potential)^2 + (T35_fit_ratio - T25_fit_ratio)^2),
         dist20 = sqrt((T45_stabil_potential - T25_stabil_potential)^2 + (T45_fit_ratio - T25_fit_ratio)^2)) %>% 
  dplyr::select(iteration, dist5, dist10, dist20)

#get set of parameters for each iteration
iter_parms <- rrc_n %>% 
  dplyr::select(iteration, ref_temp:m2_b) %>% 
  unique()

rrcn_dist <- left_join(iter_parms, rrc_n_e)

#euclidean distance of end point as a function of different parameters ######
rrcn_dist %>% 
  dplyr::select(-ref_temp) %>% 
  pivot_longer(cols = c(dist5, dist10, dist20), names_to = "warm_dist", values_to = "eucl_dist") %>%  
  pivot_longer(cols = c(r_EaN:m2_b), names_to = "param", values_to = "value") %>% 
  filter(!str_detect(param, "_b")) %>% 
  ggplot(aes(x = value, y = eucl_dist, colour = warm_dist)) +
  geom_point() + 
  facet_wrap(~param, scales = "free") + 
  labs(x = "Parameter value", y = "Euclidean distance from start point")


#trying to see if there is are certain combinations of parameters than generate large euclidean distances
rrcn_dist %>% 
  dplyr::select(-ref_temp) %>% 
  pivot_longer(cols = c(dist5, dist10, dist20), names_to = "warm_dist", values_to = "eucl_dist") %>%  
  filter(warm_dist == "dist20") %>% 
  # ggplot(aes(x = c_Ea1P, y = c_Ea2P, colour = eucl_dist, size = eucl_dist)) +
  ggplot(aes(x = r_EaN, y = r_EaP, colour = eucl_dist, size = eucl_dist)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey") + 
  geom_point() 
#big differences in r_EaN vs. r_EaP can be associated with large euclidean distances, though not always

#### explore thermal asymmetries as a source of large euclidean distances ######
##### first with reciprocal resource use, N grows faster #####
#### Note: in rrc_n, sp 1 prefers P & species 2 prefers N

#the base tassel plot again, just for orientation
#base pompom for comparison
ggplot() +
  geom_path(data = rrc_n, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0.01, 0.6, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(rrc_n, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 6, shape = "triangle") + #N grows faster
  geom_point(data = filter(rrc_n, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "turquoise", size = 5, shape = "triangle") + #N grows faster
  geom_point(data = rrc_n_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 5) +
  geom_point(data = rrc_n_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
  scale_y_log10() + 
  scale_x_log10()

# # just look at euclidean distance at 20C warming 
# rrcn_dist1 <- rrcn_dist %>% 
#   dplyr::select(-c(dist5, dist10))

hist(rrcn_dist1$dist20) #these are not normally distributed!!! most dists are small, few are large. I guess we knew that from the plots w mean position overlaid on the tassel

#look at diffs between:
# -- EA c of two species
# -- EA c of preferred resources for each species
# -- EA of the two resource growth rates
# -- EA c and EA r at the same time?
# -- EA K of two growth rates vs EA r of two growth rates

rrcn_ta <- rrcn_dist %>% 
  mutate(c_N_ta = abs(c_Ea1N - c_Ea2N),
         c_P_ta = abs(c_Ea1P - c_Ea2P),
         c_pref_ta = abs(c_Ea1P - c_Ea2N),
         R_r_ta = abs(r_EaN - r_EaP),
         R_K_ta = abs(K_EaN - K_EaP))

# effect of TA in c_EAs on euclidean distance. No real relationship for TA between two species on a single resource. Seems to be an effect of the TA between preferred resources
rrcn_ta %>% 
  pivot_longer(cols = c(dist5, dist10, dist20), names_to = "warm_dist", values_to = "eucl_dist") %>% 
  ggplot(aes(x = c_pref_ta, y = eucl_dist)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~warm_dist) 
#biggest euclidean distances happen at 20C warming -- not surprising because of the exponentials
hist(rrcn_ta$dist20)
summary(rrcn_ta$dist20)

#plot just eucl dist after 20C warming and plot the 3rd quantile of distance over the graph so we can see which points are on the longer end
rrcn_ta %>% 
  ggplot(aes(x = c_pref_ta, y = dist20)) + 
  geom_point(aes(colour = dist20 > summary(dist20)[4])) + 
  geom_smooth(method = "lm") +
  labs(colour = "Euclidean distance > \nmean distance, T/F")

lm(dist20 ~ c_pref_ta, data = rrcn_ta) #m = 0.8414

## -- EA of the two resource growth rates
rrcn_ta %>% 
  ggplot(aes(x = R_r_ta, y = dist20)) + 
  geom_point() + 
  # geom_point(aes(colour = dist20 > summary(dist20)[4])) +
  geom_smooth(method = "lm")
  # labs(colour = "Euclidean distance > \nmean distance, T/F")

lm(dist20 ~ R_r_ta, data = rrcn_ta) #m = 1.33

#make the whole df long and look at all the asymmetries
rrcn_ta %>% 
  pivot_longer(cols = c(c_N_ta, c_P_ta, c_pref_ta, R_r_ta, R_K_ta), names_to = "param_therm_asym", values_to = "abs_diff") %>% 
  mutate(param_therm_asym = str_replace(param_therm_asym, "_ta", "")) %>% 
  ggplot(aes(x = abs_diff, y = dist20)) + 
  geom_point() + 
  geom_smooth(method = "lm", colour = "red") +
  facet_wrap(~param_therm_asym) + 
  geom_hline(yintercept = summary(rrcn_ta$dist20)[4], linetype = "dashed") + #hline at mean eucl dist
  labs(y = "Euclidean distance from 25C to 45C\n (i.e. Start point to 20C warming)", x = "|EA1 - EA2|")
#amazing! thermal asymmetries in consumption rate of preferred resource, and of resource growth rates have the clearest relationships with euclidean distance to 20C

hist(rrcn_ta$c_pref_ta) #kind of skewed
hist(rrcn_ta$c_N_ta) #very skewed #why would distributions of c_pref_ta and c_N_ta have different shapes?
hist(rrcn_ta$R_K_ta) #very skewed

##### now with generalist-specialist trade-off #####
#### Note: in gs_ext, sp 1 specializes on N, and species 2 is a generalist

#Calculate Euclidean distance at 5C, 10C, 20C #####
gs_ext_e <- gs_ext %>% 
  filter(T %in% c(25, 30, 35, 45)) %>% 
  dplyr::select(-c(new_stabil_potential, a11:g2, new_fit_ratio:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(stabil_potential, fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist5 = sqrt((T30_stabil_potential - T25_stabil_potential)^2 + (T30_fit_ratio - T25_fit_ratio)^2),
         dist10 = sqrt((T35_stabil_potential - T25_stabil_potential)^2 + (T35_fit_ratio - T25_fit_ratio)^2),
         dist20 = sqrt((T45_stabil_potential - T25_stabil_potential)^2 + (T45_fit_ratio - T25_fit_ratio)^2)) %>% 
  dplyr::select(iteration, dist5, dist10, dist20)

#get set of parameters for each iteration
iter_parms <- gs_ext %>% 
  dplyr::select(iteration, ref_temp:m2_b) %>% 
  unique()

gs_ext_dist <- left_join(iter_parms, gs_ext_e)

#the base tassel plot again, just for orientation
#base pompom for comparison
ggplot() +
  geom_path(data = gs_ext, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0.01, 0.6, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  
  geom_point(data = filter(gs_ext, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 6, shape = "square") + #N grows faster
  geom_point(data = filter(gs_ext, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "salmon", size = 5, shape = "square") + #N grows faster
  geom_point(data = gs_ext_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 5) +
  geom_point(data = gs_ext_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 4) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
  scale_y_log10() + 
  scale_x_log10()


hist(gs_ext_dist$dist20) #these are not normally distributed!!! most dists are small, few are large. I guess we knew that from the plots w mean position overlaid on the tassel

#look at diffs between:
# -- EA c of two species
# -- EA c of preferred resources for each species
# -- EA of the two resource growth rates
# -- EA c and EA r at the same time?
# -- EA K of two growth rates vs EA r of two growth rates

gs_ext_ta <- gs_ext_dist %>% 
  mutate(c_N_ta = abs(c_Ea1N - c_Ea2N),
         c_P_ta = abs(c_Ea1P - c_Ea2P),
         c_pref_ta = abs(c_Ea1P - c_Ea2N),
         R_r_ta = abs(r_EaN - r_EaP),
         R_K_ta = abs(K_EaN - K_EaP))

# effect of TA in c_EAs on euclidean distance. No real relationship for TA between two species on a single resource, or preffered resource
gs_ext_ta %>% 
  pivot_longer(cols = c(dist5, dist10, dist20), names_to = "warm_dist", values_to = "eucl_dist") %>% 
  ggplot(aes(x = c_P_ta, y = eucl_dist)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~warm_dist) 
#biggest euclidean distances happen at 20C warming -- not surprising because of the exponentials
hist(gs_ext_ta$dist20)
summary(gs_ext_ta$dist20)

#plot just eucl dist after 20C warming and plot the 3rd quantile of distance over the graph so we can see which points are on the longer end
gs_ext_ta %>% 
  ggplot(aes(x = c_pref_ta, y = dist20)) + 
  geom_point(aes(colour = dist20 > summary(dist20)[4])) + 
  geom_smooth(method = "lm") +
  labs(colour = "Euclidean distance > \nmean distance, T/F")

lm(dist20 ~ c_pref_ta, data = gs_ext_ta) #m = -0.04

## -- EA of the two resource growth rates
gs_ext_ta %>% 
  ggplot(aes(x = R_r_ta, y = dist20)) + 
  geom_point() + 
  # geom_point(aes(colour = dist20 > summary(dist20)[4])) +
  geom_smooth(method = "lm")
# labs(colour = "Euclidean distance > \nmean distance, T/F")

lm(dist20 ~ R_r_ta, data = gs_ext_ta) #m = 0.59

#make the whole df long and look at all the asymmetries
gs_ext_ta %>% 
  pivot_longer(cols = c(c_N_ta, c_P_ta, c_pref_ta, R_r_ta, R_K_ta), names_to = "param_therm_asym", values_to = "abs_diff") %>% 
  mutate(param_therm_asym = str_replace(param_therm_asym, "_ta", "")) %>% 
  ggplot(aes(x = abs_diff, y = dist20)) + 
  geom_point() + 
  geom_smooth(method = "lm", colour = "red") +
  facet_wrap(~param_therm_asym) + 
  geom_hline(yintercept = summary(gs_ext_ta$dist20)[4], linetype = "dashed") + #hline at mean eucl dist
  labs(y = "Euclidean distance from 25C to 45C\n (i.e. Start point to 20C warming)", x = "|EA1 - EA2|")
#amazing! thermal asymmetries in consumption rate of resource growth rates (and maybe P consumption rate?) have the clearest relationships with euclidean distance to 20C

hist(gs_ext_ta$c_pref_ta) #very skewed
hist(gs_ext_ta$c_N_ta) #very skewed #why would distributions of c_pref_ta and c_N_ta have different shapes?
hist(gs_ext_ta$R_K_ta) #very skewed

### plot euclidean distance of start point against euclidean distance of end points #####
#### ISSUE: 3/20 plane -- this is not working properly for the reciprocal resource scenarios -- not sure why. Re-run and see if it behaves? 
big <- bind_rows(rrc_e, rrc_n, rrc_p, gs_mod, gs_ext) %>% #combine all simulation outputs
  filter(T %in% c(50, 35, 30)) %>% #just temps of interest
  mutate(dist_origin = sqrt((stabil_potential - 0)^2 + (fit_ratio - 1)^2)) %>%  #get distance from origin
  left_join(., starts_e, by = c("scenario", "sub_scenario")) #bind starting positions onto the end

#plot start position against end position
big %>% 
  ggplot(aes(x = start_dist_origin, y = dist_origin, colour = scenario, shape = sub_scenario)) +
  geom_jitter(width = 0.01) +
  labs(x = "Start position distance from origin", y = "End posiiton distance from origin") + 
  facet_wrap(~T)

### plot shifts in stab potential based on start point ####
big1 <- big %>% 
  mutate(fit_shift = fit_ratio - start_fit_ratio,
         stab_shift = stabil_potential - start_stabil_potential,
         std_fit_shift = as.vector(scale(fit_shift)),
         std_stab_shift = as.vector(scale(stab_shift)))

big1 %>% 
  ggplot(aes(x = start_dist_origin, y = std_fit_shift)) + 
  geom_point(aes(colour = T, shape = scenario), size = 4) + 
  geom_smooth(method = "lm") +
  labs(x = "Start position distance from origin", y = "Standardized shift in fitness ratio")

big1 %>% 
  ggplot(aes(x = start_dist_origin, y = std_stab_shift)) + 
  geom_point(aes(colour = T, shape = scenario), size = 4) + 
  geom_smooth(method = "lm") +
  labs(x = "Start position distance from origin", y = "Standardized shift in stabilization potential")

big1 %>% 
  ggplot(aes(x = std_stab_shift, y = std_fit_shift, colour = T, shape = sub_scenario)) + 
  geom_point(size = 4) + 
  labs(x = "Standardized shift in stabilization potential", y = "Standardized shift in fitness ratio") + 
  facet_wrap(~scenario)


