#this script is to visualize the explore the behaviour of each parameter as temperatures increase, to clarify the behaviour where species pairs tend toward neutrality

# Author: Kaleigh Davis, PDF University of Guelph
# Script DOB: 13 January 2025

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
##### draw all param EAs at random ##############
rrc <- data.frame()
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
  rrc <- bind_rows(rrc, hold) 
}

# develop relevant composite variables
rrc1 <-
  rrc %>% 
  mutate(KN_over_rN = KN/rN,
         KP_over_rP = KP/rP,
         intra_inter1 = a11/a12,
         intra_inter2 = a22/a21,
         inter_intra1 = a12/a11,
         inter_intra2 = a21/a22) %>% 
  select(-ref_temp, -coexist) %>% 
  select(iteration, T, everything()) 

#plot KN/rN as a function of KN
rrc1 %>% 
  ggplot() + 
  geom_line(aes(x = KN, y = KN_over_rN, group = iteration, colour = K_EaN))

#because K has a negative temperature dependence, K approaches 0 exponentially as systems heat up. That means that K/r grows rapidly as the system heats up

#KN temp response relative to KN/rN
rrc1 %>% 
  ggplot() + 
  geom_line(aes(x = T, y = KN, group = iteration), colour = "blue") + #negative exponential
  geom_line(aes(x = T, y = KN_over_rN, group = iteration), colour = "red") #negative exponential divided by a positive exponential

#get summary trajectories for each parameter
rrc1_add_sum <- rrc1_add %>% 
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

rrc_sums <- rrc1_add_sum %>% 
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
rrc_sums_to25 <- rrc_sums %>% 
  filter(T > 9, T < 26)

# trajectories of each param over 200C
warming_tendencies200 <-
  rrc_sums %>% 
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

# ggsave(plot = warming_tendencies200, filename = "figures/200C_param_trajectories.pdf", width = 12, height = 8)

#summaries of each param over the first 15C
warming_tendencies15 <- 
  rrc_sums_to25 %>% 
  filter(statistic == "median") %>% 
  filter(! parameter %in% c("inter_intra1", "inter_intra2")) %>% 
  ggplot() + 
  # geom_hline(yintercept = 1, colour = "green4") + 
  geom_point(aes(x = T, y = value)) +
  facet_wrap(~parameter, 
             labeller = as_labeller(param_labels, label_parsed), 
             scales = "free")

# ggsave(plot = warming_tendencies15, filename = "figures/15C_param_trajectories.pdf", width = 12, height = 8)