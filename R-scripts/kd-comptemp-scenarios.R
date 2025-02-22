#this script is to play around with some specific biological scenarios, and see how coex/compex might play out with warming in those scenarios

#Questions I'd like to answer for each scenario:
## 1. How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming
####### In order to answer this I need to calculate the coordinates of these lines? Or can I just set up the inequality and see who satisfies is at these different time points?
## 2. Do niche or fitness differences change more with warming?
####### Set up a scaled position axis and see what values are after 25 deg warming or cooling along each axis
## 3. Are there more lines that move net W (decrease in stab potential) or net E (increase stab potential) after 20 degrees of warming; net N vs net S
######## Count number of species pairs that end up more than and less than the starting stab potential and fitness ratio values

#script DOB: Feb 19, 2025
#author: Kaleigh Davis, University of Guelph postdoc

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
source("R-scripts/temp-dep-macarthur-KD.R") #this contains the macarthur translation function, exactly as Joey's functions did, but has all parameters flexibly defined in the function for assigning at time of use. Nothing is hard coded in.

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

###### SCENARIO 1: all params are temp dependent, drawn from empirical distribution ##########
##### run model with Joey's starting param values ##############
jb <- data.frame()
for(f in 1:100){ #was 200
  hold = temp_dep_mac(T = seq(0, 50, by = 0.1), #was by 0.1
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
  jb <- bind_rows(jb, hold) 
}

##### change magnitude of difference between consumption vectors to see if I can move the starting point along the threshold ####
#there should be an analytical way of identifying pairs of c_b values that will fall along the barrier line -- I need to figure this out

#distribution of points at T == 25 (== ref temp)
jb %>% 
  filter(T == 25) %>% 
  ggplot(aes(x = stabil_potential, y = fit_ratio, colour = coexist)) + 
  geom_ribbon(data = data.frame(x = seq(0.1, 0.5, 0.001)),
  aes(x = x,
      y = NULL,
      ymin = 1-x,
      ymax = 1/(1-x)),
fill = "grey", color = "black", alpha = 0.2) +
  # geom_ribbon(data = data.frame(x = seq(min(jb$stabil_potential)*0.99, max(jb$stabil_potential)*1.01, 0.001)),
              # aes(x = x,
              #     y = NULL,
              #     ymin = 1-x,
              #     ymax = 1/(1-x)),
              # fill = "grey", color = "black", alpha = 0.2) +
  geom_point() + 
  geom_hline(yintercept = 1, linetype=5) + 
  guides(colour = "none") + 
  ggtitle("Coexistence trait pairs for all iterations, at T = 25")

#base pompom for comparison
ggplot() +
  geom_path(data = jb, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(min(jb$stabil_potential)*0.99, max(jb$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION 1: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
jb_qualres1 <- jb %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(-20, -10, -5, 5, 10, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally() #these numbers are suspiciously symmetrical... they change # of true/false per temp if I switch the scenario, but it's still the same exact number of Ts and Fs for each temperature, and switches when it goes from cooling to warming

# QUESTION 2. Do niche or fitness differences change more with warming? ####
#get starting position
jb_start_stab <- jb %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(stabil_potential)) %>% #only one value in dataset
  unlist()

jb_start_fitrat <- jb %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(fit_ratio)) %>% 
  unlist()

#get dist from starting positions to end positions (i.e. 25C of cooling or warming)
jb_qualres2_25 <- jb %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -25 | rel_T == 25) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - jb_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - jb_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -25, "cooling", "warming")),
         shift = "25C")
  
jb_qualres2_25 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

jb_qualres2_25 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#get dist from starting positions to 5 degrees warmer or cooler
jb_qualres2_5 <- jb %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -5 | rel_T == 5) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - jb_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - jb_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -5, "cooling", "warming")),
         shift = "5C")

jb_qualres2_5 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

jb_qualres2_5 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#bind results together from 25 and 5 degree shifts
jb_qualres2_allshifts <- bind_rows(jb_qualres2_5, jb_qualres2_25)

#just to make sure this worked, there should be 4 rows per iteration now -- cooling 25, cooling 5, warming 5, warming 25. 
jb_qualres2_allshifts %>% 
  group_by(iteration) %>% 
  tally() %>% 
  summarise(mean_n = mean(n), sd_n = sd(n)) #good
  
# QUESTION 3. Are there more lines that move net W (decrease in stab potential) or net E (increase stab potential) after 20 degrees of warming; net N vs net S #####
jb_shift_counts <- jb_qualres2_allshifts %>% 
  ungroup() %>% 
  group_by(rel_T) %>% 
  summarise(stab_shift_pos = sum(stab_shift > 0),
            stab_shift_none = sum(stab_shift == 0),
            stab_shift_neg = sum(stab_shift < 0),
            fr_shift_pos = sum(fr_shift > 0),
            fr_shift_none = sum(fr_shift == 0),
            fr_shift_neg = sum(fr_shift < 0))
#under both cooling scenarios: most species pairs shift down (negative) in stab potential and fitness ratios
#under both warming scenarios: most species pairs shift down in stab potential. At 25 degrees warming, most fitness ratios also decrease, but at 5C warming, most fitness ratios increase.

###### N supply is temperature dependent but P supply is not ##########
# hard code in EAs #####
neanop <- data.frame()
for(f in 1:200){
  hold = temp_dep_mac(T = seq(0, 50, by = 0.5),
                      ref_temp = 25,
                      r_EaN = 1.1, #draw all EAs from empirical distributions above
                      r_EaP = 0.1, 
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
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  neanop <- bind_rows(neanop, hold) 
}

#quick pompom plot
ggplot() +
  geom_path(data = neanop, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), size = 2) +
  geom_ribbon(data = data.frame(x = seq(min(neanop$stabil_potential)*0.99, max(neanop$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# draw N EA from dist #######
neanop1 <- data.frame()
for(f in 1:200){
  hold = temp_dep_mac(T = seq(0, 50, by = 0.5),
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = 0, 
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
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  neanop1 <- bind_rows(neanop1, hold) 
}

#quick pompom plot
ggplot() +
  geom_path(data = neanop1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), size = 2) +
  geom_ribbon(data = data.frame(x = seq(min(neanop1$stabil_potential)*0.99, max(neanop1$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 


###### SCENARIO 2: Generalist-specialist trade-off #######
gs1 <- data.frame()
for(f in 1:200){
  hold = temp_dep_mac(T = seq(0, 50, by = 0.5),
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
                      c1N_b = 0.4, c1P_b = 0.6, #spec 1 consumes much more P than N
                      c2N_b = 0.5, c2P_b = 0.5, #spec 2 consumes N and P equally
                      r_N_b = 0.1, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.4, v1P_b = 0.6, #sp 1 converts P much more efficiently than N
                      v2N_b = 0.5, v2P_b = 0.5, #sp 2 converts N & P equally well
                      m1_b = 0.1, m2_b = 0.01) #same for both species
  hold$iteration <- f
  gs1 <- bind_rows(gs1, hold) 
}

#distribution of points at T == 25 (== ref temp)
gs1 %>% 
  filter(T == 25) %>% 
  ggplot(aes(x = stabil_potential, y = fit_ratio, colour = coexist)) + 
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_ribbon(data = data.frame(x = seq(min(jb$stabil_potential)*0.99, max(jb$stabil_potential)*1.01, 0.001)),
  # aes(x = x,
  #     y = NULL,
  #     ymin = 1-x,
  #     ymax = 1/(1-x)),
  # fill = "grey", color = "black", alpha = 0.2) +
  geom_point() + 
  geom_hline(yintercept = 1, linetype=5) + 
  guides(colour = "none") + 
  ggtitle("Coexistence trait pairs for all iterations, at T = 25")

#quick pompom plot
ggplot() +
  geom_path(data = gs1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), size = 2) +
  geom_ribbon(data = data.frame(x = seq(min(gs1$stabil_potential)*0.99, max(gs1$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

#make consumption rates temp dep ####
gs2 <- data.frame()
for(f in 1:200){
  hold = temp_dep_mac(T = seq(0, 50, by = 0.5),
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = 0.586, #bottom quartile
                      c_Ea1P = 0.544, #top quartile
                      c_Ea2N = 0.565, #mean
                      c_Ea2P = 0.565, #mean
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.4, c1P_b = 0.6, #spec 1 consumes much more P than N
                      c2N_b = 0.5, c2P_b = 0.5, #spec 2 consumes N and P equally
                      r_N_b = 0.1, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.4, v1P_b = 0.6, #sp 1 converts P much more efficiently than N
                      v2N_b = 0.5, v2P_b = 0.5, #sp 2 converts N & P equally well
                      m1_b = 0.1, m2_b = 0.01) #same for both species
  hold$iteration <- f
  gs2 <- bind_rows(gs2, hold) 
}

#distribution of points at T == 25 (== ref temp)
gs2 %>% 
  filter(T == 25) %>% 
  ggplot(aes(x = stabil_potential, y = fit_ratio, colour = coexist)) + 
  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_ribbon(data = data.frame(x = seq(min(jb$stabil_potential)*0.99, max(jb$stabil_potential)*1.01, 0.001)),
  # aes(x = x,
  #     y = NULL,
  #     ymin = 1-x,
  #     ymax = 1/(1-x)),
  # fill = "grey", color = "black", alpha = 0.2) +
  geom_point() + 
  geom_hline(yintercept = 1, linetype=5) + 
  guides(colour = "none") + 
  ggtitle("Coexistence trait pairs for all iterations, at T = 25")

#quick pompom plot
ggplot() +
  geom_path(data = gs2, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), size = 2) +
  geom_ribbon(data = data.frame(x = seq(min(gs2$stabil_potential)*0.99, max(gs2$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  # geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

