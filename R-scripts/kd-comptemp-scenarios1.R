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
  mutate(scenario = "borderline",
         sub_scenario = "N-faster")

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
  mutate(scenario = "borderline",
         sub_scenario = "P-faster")

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
  mutate(scenario = "borderline",
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
  ggplot(aes(x = stabil_potential, y = fit_ratio)) + 
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)), #can log these in order to make the coexistence zone more symmetrical
              fill = "darkgrey", color = "black", alpha = 0.2) +
  geom_point(aes(shape = scenario, colour = sub_scenario), size = 4) +
  geom_hline(yintercept = 1, linetype=5) +
  geom_point(x = 0.39, y = 0.35, size = 4, shape = 21) + 
  coord_cartesian(ylim = c(0.2, 2)) +
  # guides(colour = "none") +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
  scale_x_log10() + 
  scale_y_log10()

#get euclidean distance of start points to "origin"
starts_e <- starts %>% 
  dplyr::select(c1N_b:m2_b, stabil_potential, fit_ratio) %>%
  mutate(dist_origin = sqrt((stabil_potential - 0)^2 + (fit_ratio - 1)^2))
#need to bind these onto the right hand side of the final dataset? Or can cram all the simulations together into a massive df and make a column for scenario and then calculate start euclidean distance there.



############################################################################################
################## RECIPROCAL RESOURCE CONSUMPTION AT T_REF SCENARIOS ######################
############################################################################################

##### draw all params at random ##############
jb <- data.frame()
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
  jb <- bind_rows(jb, hold) 
}

jb %>% filter(rho == 1)

#get average change in position after 5, 10, 20C warming
# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
jb_avg <- jb %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#base pompom for comparison
ggplot() +
  geom_path(data = jb, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = filter(jb, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = filter(jb, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4, shape = "triangle") + #N grows faster
  geom_point(data = filter(jb, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "turquoise", size = 3, shape = "triangle") + #N grows faster
  geom_point(data = jb_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = jb_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 3) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) + 
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
jb %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(5, 10, 15, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally()

#get quantiles for all parameters
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
  geom_point(data = filter(jb, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
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
gs <- data.frame()
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
  gs <- bind_rows(gs, hold) 
}

#get average change in position after 5, 10, 20C warming
# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
gs_avg <- gs %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#pompom plot
ggplot() +
  geom_path(data = gs, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(0.0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = filter(gs, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = filter(gs, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4, shape = "square") + #N grows faster
  geom_point(data = filter(gs, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "salmon", size = 3, shape = "square") + #N grows faster
  geom_point(data = gs_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = gs_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T),  size = 3) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  labs(colour = "Degrees warming") +
  # coord_cartesian(ylim=c(0.2, 2), xlim = c(0, 0.5)) +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
  scale_x_log10() + 
  scale_y_log10()

# QUESTION: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
gs %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(5, 10, 15, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally()

#Calculate Euclidean distance at 5C, 10C, 20C
gs_e <- gs %>% 
  filter(T %in% c(25, 30, 35, 45)) %>% 
  dplyr::select(-c(log_stabil_potential, a11:g2, log_fit_ratio:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(stabil_potential, fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  mutate(dist5 = sqrt((T30_stabil_potential - T25_stabil_potential)^2 + (T30_fit_ratio - T25_fit_ratio)^2),
         dist10 = sqrt((T35_stabil_potential - T25_stabil_potential)^2 + (T35_fit_ratio - T25_fit_ratio)^2),
         dist20 = sqrt((T45_stabil_potential - T25_stabil_potential)^2 + (T45_fit_ratio - T25_fit_ratio)^2)) %>% 
  dplyr::select(iteration, dist5, dist10, dist20)

#get set of parameters for each iteration
iter_parms <- gs %>% 
  dplyr::select(iteration, ref_temp:m2_b) %>% 
  unique()

gs_dist <- left_join(iter_parms, gs_e)

#distance as a function of different parameters
gs_dist %>% 
  dplyr::select(-ref_temp) %>% 
  pivot_longer(cols = c(dist5, dist10, dist20), names_to = "warm_dist", values_to = "eucl_dist") %>%  
  pivot_longer(cols = c(r_EaN:m2_b), names_to = "param", values_to = "value") %>% 
  filter(!str_detect(param, "_b")) %>% 
  ggplot(aes(x = value, y = eucl_dist, colour = warm_dist)) +
  geom_point() + 
  facet_wrap(~param, scales = "free")

#trying to see if there is are certain combinations of parameters than generate large euclidean distances
gs_dist %>% 
  dplyr::select(-ref_temp) %>% 
  pivot_longer(cols = c(dist5, dist10, dist20), names_to = "warm_dist", values_to = "eucl_dist") %>%  
  filter(warm_dist == "dist20") %>% 
  ggplot(aes(x = v_EaP, y = c_Ea2P, colour = eucl_dist)) +
  geom_point(size = 4) 
#big differences in r_EaN vs. r_EaP can be associated with large euclidean distances, though not always

hist(gs$rho)
hist(jb$rho)




