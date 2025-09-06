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
source("R-scripts/kd-analysis-scripts/02-temp-dep-macarthur-KD.R") #this contains the macarthur translation function, exactly as Joey's functions did, but has all parameters flexibly defined in the function for assigning at time of use. Nothing is hard coded in.

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
for(f in 1:200){ #was 200
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
  # geom_ribbon(data = data.frame(x = seq(-1, 1, 0.001)),
  #             aes(x = x,
  #                 y = NULL,
  #                 ymin = 1-x,
  #                 ymax = 1/(1-x)),
  #             fill = "darkgrey", color = "black", alpha = 0.2) +
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
  geom_point(colour = "black", size = 4) +
  geom_hline(yintercept = 1, linetype=5) +
  coord_cartesian(ylim = c(0, 2.5)) +
  guides(colour = "none") +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 
  # ggtitle("Coexistence trait pairs for all iterations, at T = 25") #exporting at 

#base pompom for comparison
ggplot() +
  geom_path(data = jb, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2) +
  geom_ribbon(data = data.frame(x = seq(min(jb$stabil_potential)*0.99, max(jb$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(jb, T==25), aes(x = stabil_potential, y = fit_ratio), fill = "pink", size = 3) +
  geom_point(data = filter(jb, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION 1: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
# jb_qualres1 <- 
  jb %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(-20, -10, -5, 5, 10, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally() #these numbers are suspiciously symmetrical... they change # of true/false per temp if I switch the scenario, but it's still the same exact number of Ts and Fs for each temperature, and switches when it goes from cooling to warming

# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
jb_avg <- jb %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))


#plot em over the pompom
ggplot() +
  geom_path(data = jb, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2, alpha = 0.5) +
  geom_ribbon(data = data.frame(x = seq(min(jb$stabil_potential)*0.99, max(jb$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(jb, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = jb_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = jb_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T), size = 3) + 
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

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
jb_qualres2_20 <- jb %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -20 | rel_T == 20) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - jb_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - jb_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -20, "20C cooling", "20C warming")),
         shift = "20C")
  
jb_qualres2_20 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("20C cooling" = "blue", "20C warming" = "red"))

jb_qualres2_20 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#get dist from starting positions to 10 degrees warmer or cooler
jb_qualres2_10 <- jb %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -10 | rel_T == 10) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - jb_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - jb_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -10, "10C cooling", "10C warming")),
         shift = "10C")

jb_qualres2_10 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("10C cooling" = "blue", "10C warming" = "red"))

jb_qualres2_10 %>% 
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
  mutate(t_scen = as.factor(ifelse(rel_T == -5, "5C cooling", "5C warming")),
         shift = "5C")

jb_qualres2_5 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("5C cooling" = "blue", "5C warming" = "red"))

jb_qualres2_5 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#bind results together from 25 and 5 degree shifts
jb_qualres2_allshifts <- bind_rows(jb_qualres2_5, jb_qualres2_10, jb_qualres2_20)

#just to make sure this worked, there should be 6 rows per iteration now -- cooling 25, cooling 5, warming 5, warming 25. 
jb_qualres2_allshifts %>% 
  group_by(iteration) %>% 
  tally() %>% 
  summarise(mean_n = mean(n), sd_n = sd(n)) #good

#now plot the different distributions over each other, warming scenarios only
jb_qualres2_allshifts %>%
  filter(str_detect(.$t_scen, "warming")) %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  viridis::scale_fill_viridis(discrete = TRUE) +
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") 
  # scale_fill_manual(values = ifelse(grepl("cooling", jb_qualres2_allshifts$t_scen), "blue", ifelse(grepl("warming", jb_qualres2_allshifts$t_scen), "red", "gray")))

  
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

###### SCENARIO 2: N supply is temperature dependent but P supply is not ##########
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
                      r_N_b = 0.05, r_P_b = 0.1, #growth rate for each resource at ref temp. Note: setting these equal to each other, rather than with r_N = 1 and r_P = 0.05 moves the starting point from the upper line down to perfectly on the dotted line (fitness diffs == 1)
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
  geom_point(data = filter(neanop, T == 25), aes(x = stabil_potential, y = fit_ratio), color = "black", size = 4) +
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

# QUESTION 1: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
neanop1_1 <- neanop1 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(-20, -10, -5, 5, 10, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally() #these numbers are suspiciously symmetrical... they change # of true/false per temp if I switch the scenario, but it's still the same exact number of Ts and Fs for each temperature, and switches when it goes from cooling to warming

# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
neanop1_avg <- neanop1 %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#plot em over the pompom
ggplot() +
  geom_path(data = neanop1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2, alpha = 0.5) +
  geom_ribbon(data = data.frame(x = seq(min(neanop1$stabil_potential)*0.99, max(neanop1$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(neanop1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = neanop1_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = neanop1_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T), size = 3) + 
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION 2. Do niche or fitness differences change more with warming? ####
#get starting position
nea_start_stab <- neanop1 %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(stabil_potential)) %>% #only one value in dataset
  unlist()

nea_start_fitrat <- neanop1 %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(fit_ratio)) %>% 
  unlist()

#get dist from starting positions to end positions (i.e. 25C of cooling or warming)
nea_qualres2_20 <- neanop1 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -20 | rel_T == 20) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - nea_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - nea_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -20, "20C cooling", "20C warming")),
         shift = "20C")

nea_qualres2_20 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("20C cooling" = "blue", "20C warming" = "red"))

nea_qualres2_20 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#get dist from starting positions to 10 degrees warmer or cooler
nea_qualres2_10 <- neanop1 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -10 | rel_T == 10) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - nea_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - nea_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -10, "10C cooling", "10C warming")),
         shift = "10C")

nea_qualres2_10 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("10C cooling" = "blue", "10C warming" = "red"))

nea_qualres2_10 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)


#get dist from starting positions to 5 degrees warmer or cooler
nea_qualres2_5 <- neanop1 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -5 | rel_T == 5) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - nea_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - nea_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -5, "5C cooling", "5C warming")),
         shift = "5C")

nea_qualres2_5 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("5C cooling" = "blue", "5C warming" = "red"))

nea_qualres2_5 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#bind results together from 25 and 5 degree shifts
nea_qualres2_allshifts <- bind_rows(nea_qualres2_5, nea_qualres2_10, nea_qualres2_20)

#just to make sure this worked, there should be 6 rows per iteration now -- cooling 25, cooling 5, warming 5, warming 25. 
nea_qualres2_allshifts %>% 
  group_by(iteration) %>% 
  tally() %>% 
  summarise(mean_n = mean(n), sd_n = sd(n)) #good

#now plot the different distributions over each other, warming scenarios only
nea_qualres2_allshifts %>%
  # filter(str_detect(.$t_scen, "warming")) %>% 
  filter(str_detect(.$t_scen, "cooling")) %>%
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  viridis::scale_fill_viridis(discrete = TRUE) +
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") 
# scale_fill_manual(values = ifelse(grepl("cooling", nea_qualres2_allshifts$t_scen), "blue", ifelse(grepl("warming", nea_qualres2_allshifts$t_scen), "red", "gray")))


###### SCENARIO 3: Generalist-specialist trade-off #######
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
                      m1_b = 0.01, m2_b = 0.01) #same for both species
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

# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
gs1_avg <- gs1 %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#plot em over the pompom
ggplot() +
  geom_path(data = gs1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2, alpha = 0.5) +
  geom_ribbon(data = data.frame(x = seq(min(gs1$stabil_potential)*0.99, max(gs1$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(gs1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = gs1_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = gs1_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T), size = 3) + 
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
                      c_Ea1N = 0.544, #bottom quartile ## I HAD THESE SWITCHED BEFORE! 9/5/25
                      c_Ea1P = 0.586, #top quartile
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
                      m1_b = 0.01, m2_b = 0.01) #same for both species
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
  geom_point(data = filter(gs2, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) + 
  # geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION 1: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
gs2_1 <- gs2 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(-20, -10, -5, 5, 10, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally() #these numbers are suspiciously symmetrical... they change # of true/false per temp if I switch the scenario, but it's still the same exact number of Ts and Fs for each temperature, and switches when it goes from cooling to warming

# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
gs2_avg <- gs2 %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#plot em over the pompom
ggplot() +
  geom_path(data = gs2, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2, alpha = 0.5) +
  geom_ribbon(data = data.frame(x = seq(min(gs2$stabil_potential)*0.99, max(gs2$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(gs2, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = gs2_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = gs2_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T), size = 3) + 
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION 2. Do niche or fitness differences change more with warming? ####
#get starting position
gs_start_stab <- gs2 %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(stabil_potential)) %>% #only one value in dataset
  unlist()

gs_start_fitrat <- gs2 %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(fit_ratio)) %>% 
  unlist()

#get dist from starting positions to end positions (i.e. 25C of cooling or warming)
gs_qualres2_20 <- gs2 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -20 | rel_T == 20) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - gs_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - gs_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -20, "20C cooling", "20C warming")),
         shift = "20C")

gs_qualres2_20 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("20C cooling" = "blue", "20C warming" = "red"))

gs_qualres2_20 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#get dist from starting positions to 10 degrees warmer or cooler
gs_qualres2_10 <- gs2 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -10 | rel_T == 10) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - gs_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - gs_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -10, "10C cooling", "10C warming")),
         shift = "10C")

gs_qualres2_10 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("10C cooling" = "blue", "10C warming" = "red"))

gs_qualres2_10 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)


#get dist from starting positions to 5 degrees warmer or cooler
gs_qualres2_5 <- gs2 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -5 | rel_T == 5) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - gs_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - gs_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -5, "5C cooling", "5C warming")),
         shift = "5C")

gs_qualres2_5 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("5C cooling" = "blue", "5C warming" = "red"))

gs_qualres2_5 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#bind results together from 25 and 5 degree shifts
gs_qualres2_allshifts <- bind_rows(gs_qualres2_5, gs_qualres2_10, gs_qualres2_20)

#just to make sure this worked, there should be 6 rows per iteration now -- cooling 25, cooling 5, warming 5, warming 25. 
gs_qualres2_allshifts %>% 
  group_by(iteration) %>% 
  tally() %>% 
  summarise(mean_n = mean(n), sd_n = sd(n)) #good

#now plot the different distributions over each other, warming scenarios only
gs_qualres2_allshifts %>%
  filter(str_detect(.$t_scen, "warming")) %>%
  # filter(str_detect(.$t_scen, "cooling")) %>%
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  viridis::scale_fill_viridis(discrete = TRUE) +
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario")

###### SCENARIO 4: Resident that does well on dominant resource, invader that is very temperature sensitive in its consumption #######
ri <- data.frame()
for(f in 1:200){
  hold = temp_dep_mac(T = seq(0, 50, by = 0.5),
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, #draw all EAs from empirical distributions above
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2N = 0.6,
                      c_Ea2P = 0.6,
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.8, c1P_b = 0.2, #spec 1 consumes much more P than N
                      c2N_b = 0.5, c2P_b = 0.5, #spec 2 consumes N and P equally
                      r_N_b = 0.1, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.8, v1P_b = 0.2, #sp 1 converts P much more efficiently than N
                      v2N_b = 0.5, v2P_b = 0.5, #sp 2 converts N & P equally well
                      m1_b = 0.1, m2_b = 0.01) #same for both species
  hold$iteration <- f
  ri <- bind_rows(ri, hold) 
}

#distribution of points at T == 25 (== ref temp)
ri %>% 
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
  geom_path(data = ri, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), size = 2) +
  geom_ribbon(data = data.frame(x = seq(min(ri$stabil_potential)*0.99, max(ri$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(ri, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  # geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 


# QUESTION 1: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
ri_1 <- ri %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(-20, -10, -5, 5, 10, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally() #these numbers are suspiciously symmetrical... they change # of true/false per temp if I switch the scenario, but it's still the same exact number of Ts and Fs for each temperature, and switches when it goes from cooling to warming

# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
ri_avg <- ri %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#plot em over the pompom
ggplot() +
  geom_path(data = ri, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2, alpha = 0.5) +
  geom_ribbon(data = data.frame(x = seq(min(ri$stabil_potential)*0.99, max(ri$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(ri, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = ri_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = ri_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T), size = 3) + 
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION 2. Do niche or fitness differences change more with warming? ####
#get starting position
ri_start_stab <- ri %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(stabil_potential)) %>% #only one value in dataset
  unlist()

ri_start_fitrat <- ri %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(fit_ratio)) %>% 
  unlist()

#get dist from starting positions to end positions (i.e. 25C of cooling or warming)
ri_qualres2_20 <- ri %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -20 | rel_T == 20) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - ri_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - ri_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -20, "20C cooling", "20C warming")),
         shift = "20C")

ri_qualres2_20 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("20C cooling" = "blue", "20C warming" = "red"))

ri_qualres2_20 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#get dist from starting positions to 10 degrees warmer or cooler
ri_qualres2_10 <- ri %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -10 | rel_T == 10) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - ri_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - ri_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -10, "10C cooling", "10C warming")),
         shift = "10C")

ri_qualres2_10 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("10C cooling" = "blue", "10C warming" = "red"))

ri_qualres2_10 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)


#get dist from starting positions to 5 degrees warmer or cooler
ri_qualres2_5 <- ri %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -5 | rel_T == 5) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - ri_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - ri_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -5, "5C cooling", "5C warming")),
         shift = "5C")

ri_qualres2_5 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario") +
  scale_fill_manual(values = c("5C cooling" = "blue", "5C warming" = "red"))

ri_qualres2_5 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#bind results together from 25 and 5 degree shifts
ri_qualres2_allshifts <- bind_rows(ri_qualres2_5, ri_qualres2_10, ri_qualres2_20)

#just to make sure this worked, there should be 6 rows per iteration now -- cooling 25, cooling 5, warming 5, warming 25. 
ri_qualres2_allshifts %>% 
  group_by(iteration) %>% 
  tally() %>% 
  summarise(mean_n = mean(n), sd_n = sd(n)) #good

#now plot the different distributions over each other, warming scenarios only
ri_qualres2_allshifts %>%
  filter(str_detect(.$t_scen, "warming")) %>%
  # filter(str_detect(.$t_scen, "cooling")) %>%
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  viridis::scale_fill_viridis(discrete = TRUE) +
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenario")


#Scenario 4, alternate parameterization -- make dominant resource have higher growth rate #####
ri1 <- data.frame()
for(f in 1:200){
  hold = temp_dep_mac(T = seq(0, 50, by = 0.5),
                      ref_temp = 25,
                      r_EaN = sample_n(rgr_post_dist, size = 1)$intercept, 
                      r_EaP = sample_n(rgr_post_dist, size = 1)$intercept, 
                      c_Ea1N = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea1P = sample_n(c_post_dist, size = 1)$intercept,
                      c_Ea2N = 0.6,
                      c_Ea2P = 0.6,
                      K_EaN = sample_n(k_post_dist, size = 1)$intercept, 
                      K_EaP = sample_n(k_post_dist, size = 1)$intercept, 
                      v_EaN = sample_n(v_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(v_post_dist, size = 1)$intercept, 
                      m_Ea1 = sample_n(m_post_dist, size = 1)$intercept, 
                      m_Ea2 = sample_n(m_post_dist, size = 1)$intercept,
                      c1N_b = 0.8, c1P_b = 0.2, #spec 1 consumes much more P than N
                      c2N_b = 0.5, c2P_b = 0.5, #spec 2 consumes N and P equally
                      r_N_b = 0.3, r_P_b = 0.05, #growth rate for each resource at ref temp
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp
                      v1N_b = 0.8, v1P_b = 0.2, #sp 1 converts P much more efficiently than N
                      v2N_b = 0.5, v2P_b = 0.5, #sp 2 converts N & P equally well
                      m1_b = 0.1, m2_b = 0.01) #same for both species
  hold$iteration <- f
  ri1 <- bind_rows(ri1, hold) 
}

#distribution of points at T == 25 (== ref temp)
ri1 %>% 
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
  geom_path(data = ri1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), size = 2) +
  geom_ribbon(data = data.frame(x = seq(min(ri1$stabil_potential)*0.99, max(ri1$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(ri1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  # geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION 1: How many species pairs end up in coexistence vs competitive exclusion after a) 20 degrees cooling, b) 10 degrees cooling, c) 5 degrees warming, d) 10 degrees warming, e) 20 degrees warming ####
ri1_1 <- ri1 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T %in% c(-20, -10, -5, 5, 10, 20)) %>% 
  group_by(rel_T, coexist) %>% 
  tally() #these numbers are suspiciously symmetri1cal... they change # of true/false per temp if I switch the scenari1o, but it's still the same exact number of Ts and Fs for each temperature, and switches when it goes from cooling to warming

# QUESTION. What is the average position of the dot after 5C, 10C, 20C warming? #####
ri1_avg <- ri1 %>% 
  mutate(rel_T = T-25) %>% 
  filter(rel_T %in% c(5, 10, 20)) %>% 
  group_by(rel_T) %>% 
  summarise(mean_stab_pot = mean(stabil_potential),
            mean_fit_rat = mean(fit_ratio))

#plot em over the pompom
ggplot() +
  geom_path(data = ri1, aes(x = stabil_potential, y = fit_ratio, color = T-25, group = iteration), linewidth = 2, alpha = 0.5) +
  geom_ribbon(data = data.frame(x = seq(min(ri1$stabil_potential)*0.99, max(ri1$stabil_potential)*1.01, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", color = "black", alpha = 0.2) +
  geom_point(data = filter(ri1, T==25), aes(x = stabil_potential, y = fit_ratio), colour = "black", size = 4) +
  geom_point(data = ri1_avg, aes(x = mean_stab_pot, y = mean_fit_rat), colour = "black",  size = 4) +
  geom_point(data = ri1_avg, aes(x = mean_stab_pot, y = mean_fit_rat, colour = rel_T), size = 3) + 
  geom_hline(yintercept = 1, linetype=5) + 
  scale_colour_continuous_diverging() +
  xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
  ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

# QUESTION 2. Do niche or fitness differences change more with warming? ####
#get starting position
ri1_start_stab <- ri1 %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(stabil_potential)) %>% #only one value in dataset
  unlist()

ri1_start_fitrat <- ri1 %>% 
  filter(T == 25) %>% 
  summarise(mean = mean(fit_ratio)) %>% 
  unlist()

#get dist from starting positions to end positions (i.e. 25C of cooling or warming)
ri1_qualres2_20 <- ri1 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -20 | rel_T == 20) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - ri1_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - ri1_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -20, "20C cooling", "20C warming")),
         shift = "20C")

ri1_qualres2_20 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenari1o") +
  scale_fill_manual(values = c("20C cooling" = "blue", "20C warming" = "red"))

ri1_qualres2_20 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#get dist from starting positions to 10 degrees warmer or cooler
ri1_qualres2_10 <- ri1 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -10 | rel_T == 10) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - ri1_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - ri1_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -10, "10C cooling", "10C warming")),
         shift = "10C")

ri1_qualres2_10 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenari1o") +
  scale_fill_manual(values = c("10C cooling" = "blue", "10C warming" = "red"))

ri1_qualres2_10 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)


#get dist from starting positions to 5 degrees warmer or cooler
ri1_qualres2_5 <- ri1 %>% 
  mutate(rel_T = T - 25) %>% 
  filter(rel_T == -5 | rel_T == 5) %>% 
  mutate(scaled_stab = as.vector(scale(stabil_potential)),
         scaled_fitratio = as.vector(scale(fit_ratio)),
         stab_shift = stabil_potential - ri1_start_stab, #this gives a different answer than the scaled value; need to figure out what the scaled value thing is actually doing...
         fr_shift = fit_ratio - ri1_start_fitrat) %>% 
  mutate(t_scen = as.factor(ifelse(rel_T == -5, "5C cooling", "5C warming")),
         shift = "5C")

ri1_qualres2_5 %>% 
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenari1o") +
  scale_fill_manual(values = c("5C cooling" = "blue", "5C warming" = "red"))

ri1_qualres2_5 %>% 
  dplyr::select(-scaled_stab, -scaled_fitratio) %>% 
  pivot_longer(cols = c(stab_shift, fr_shift), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var)

#bind results together from 25 and 5 degree shifts
ri1_qualres2_allshifts <- bind_rows(ri1_qualres2_5, ri1_qualres2_10, ri1_qualres2_20)

#just to make sure this worked, there should be 6 rows per iteration now -- cooling 25, cooling 5, warming 5, warming 25. 
ri1_qualres2_allshifts %>% 
  group_by(iteration) %>% 
  tally() %>% 
  summarise(mean_n = mean(n), sd_n = sd(n)) #good

#now plot the different distri1butions over each other, warming scenari1os only
ri1_qualres2_allshifts %>%
  filter(str_detect(.$t_scen, "warming")) %>%
  # filter(str_detect(.$t_scen, "cooling")) %>%
  dplyr::select(-stab_shift, -fr_shift) %>% 
  pivot_longer(cols = c(scaled_stab, scaled_fitratio), names_to = "var", values_to = "value") %>% 
  ggplot() + 
  geom_histogram(aes(x = value, fill = t_scen)) + 
  facet_wrap(~ var, labeller = labeller(var = c("scaled_stab" = "Stability", "scaled_fitratio" = "Fitness Ratio"))) + 
  viridis::scale_fill_viridis(discrete = TRUE) +
  xlab("Shift from value at ref temp") + 
  labs(fill = "Temperature \nscenari1o") #this seems wrong -- it seems like these should be moving fitness ratio up and stability to the left
