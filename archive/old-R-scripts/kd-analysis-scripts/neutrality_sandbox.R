#this script is to explore the behaviour of each parameter as temperatures increase, in an attempt to clarify the behaviour where species pairs tend toward neutrality

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

#make v only positive to see effects on the pompom
vpos_post_dist <- param_vals %>% 
  filter(parameter == "conversion_efficiency" & intercept > 0) %>% 
  mutate(parameter1 = "vpos") %>% 
  dplyr::select(-parameter) %>% 
  rename(parameter = parameter1)
  

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


#again with vpos
v_pos %>%
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
  hold = temp_dep_mac(T = seq(10, 210, by = 0.5), #was by 0.1
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

rrc1 <- rrc %>% 
  mutate(amtN = KN/rN,
         amtP = KP/rP)

# absolute competition coefficients
rrc1_long <- 
  rrc1 %>% 
  select(-ref_temp, -coexist) %>% 
  select(iteration, T, everything()) %>% 
  pivot_longer(cols = c(r_EaN:beta12), names_to = "parameter", values_to = "value")

  
# 53 variables *75500temp*iteration pairs = 3850500 #check

rrc_long1 <- 
  rrc1 %>% 
  mutate(KN_over_rN = KN/rN,
         KP_over_rP = KP/rP) %>% 
  select(-ref_temp, -coexist) %>% 
  select(iteration, T, everything()) %>% 
  pivot_longer(cols = c(r_EaN:KP_over_rP), names_to = "parameter", values_to = "value")

rrc_long1 %>% 
  filter(str_detect(parameter, "amt")) %>% 
  ggplot() +
  geom_line(aes(x = T, y = value, group = iteration)) +
  facet_wrap(~parameter)

rrc1 %>% 
  ggplot() + 
  geom_line(aes(x = KN, y = amtN, group = iteration, colour = K_EaN))

#because K has a negative temperature dependence, K approaches 0 exponentially as systems heatup. That means that K/r grows rapidly as the system heats up

rrc1 %>% 
  ggplot() + 
  geom_line(aes(x = T, y = KN, group = iteration), colour = "blue") + #negative exponential
  geom_line(aes(x = T, y = amtN, group = iteration), colour = "red") #negative exponential divided by a positive exponential

rrc_long3 <- 
  rrc1 %>% 
  mutate(KN_over_rN = KN/rN,
         KP_over_rP = KP/rP,
         intra_inter1 = a11/a12,
         intra_inter2 = a22/a21,
         inter_intra1 = a12/a11,
         inter_intra2 = a21/a22) %>% 
  select(-ref_temp, -coexist) %>% 
  select(iteration, T, everything()) %>% 
  pivot_longer(cols = c(a11:intra_inter2), names_to = "parameter", values_to = "value")

#Species 1 prefers P, species 2 prefers N, N grows faster at ref temp
rrc_long3 %>% 
  filter(parameter %in% c("KN", "rN", "KN_over_rN", "a11", "a12", "a22", "a21", "intra_inter1", "intra_inter2", "inter_intra1", "inter_intra2")) %>% 
  ggplot() + 
  geom_line(aes(x = T, y = value, group = iteration, colour = r_EaN)) + 
  facet_wrap(~parameter, scales = "free_y")

rrc_long3 %>% 
  filter(parameter %in% c("a11", "a12", "a22", "a21")) %>%
  ggplot() + 
  geom_line(aes(x = T, y = value, group = iteration)) + 
  facet_wrap(~parameter) + 
  coord_cartesian(ylim = c(-10, 10))

rrc_long3 %>% 
  filter(parameter %in% c("KN", "rN", "KN_over_rN", "a11", "a12", "a22", "a21", "new_fit_ratio", "new_stabil_potential")) %>% 
  ggplot() + 
  geom_line(aes(x = T, y = value, group = iteration, colour = K_EaN/r_EaN)) + 
  facet_wrap(~parameter, scales = "free_y") 
  # coord_cartesian(ylim = c(-5, 5), xlim = c(0, 100))

rrc_long3 %>% 
  filter(parameter %in% c("KN", "KN_over_rN", "K_EaN", "r_EaN")) %>% 
  pivot_wider(names_from = parameter, values_from = value) %>% 
  ggplot() + 
  geom_line(aes(x = KN, y = KN_over_rN, colour = K_EaN/r_EaN, group = iteration)) +
  geom_abline(slope = 1)

rrc_wide3 <- rrc_long3 %>% 
  pivot_wider(names_from = parameter, values_from = value)

rrc_wide3_sum <- rrc_wide3 %>% 
  group_by(T) %>% 
  summarise(across(
    .cols = c(a11, a12, a21, a22, KN, KN_over_rN, new_fit_ratio, new_stabil_potential, rN, intra_inter1, intra_inter2, inter_intra1, inter_intra2),
    .fns = list(
      mean   = ~mean(.x, na.rm = TRUE),
      median = ~median(.x, na.rm = TRUE),
      sd     = ~sd(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ),
  .groups = "drop"
  )

rrc_sums <- rrc_wide3_sum %>% 
  pivot_longer(cols = c(a11_mean:inter_intra2_sd), names_to = "parameter", values_to = "value") %>% 
  separate(parameter, sep = "_(?=[^_]+$)", into = c("parameter", "statistic")) %>%  
  filter(parameter %in% c("a11", "a12", "a21", "a22", "KN", "KN_over_rN", "new_fit_ratio", "new_stabil_potential", "rN", "intra_inter1", "intra_inter2", "inter_intra1", "inter_intra2") & statistic %in% c("median", "mean")) %>% 
  mutate(parameter = forcats::fct_relevel(parameter, "KN", "rN", "KN_over_rN", "a11", "a12", "a21", "a22", "intra_inter1", "intra_inter2", "inter_intra1", "inter_intra2", "new_fit_ratio", "new_stabil_potential")) 
  
param_labels <- c(
  a11 = "alpha[11]",
  a12 = "alpha[12]",
  a21 = "alpha[21]",
  a22 = "alpha[22]",
  KN  = "K[N]",
  KN_over_rN = "K[N]/r[N]",
  new_fit_ratio = "Fitness~ratio",
  new_stabil_potential = "Stabilization~potential",
  rN = "r[N]",
  intra_inter1 = "alpha[11]/alpha[12]",
  intra_inter2 = "alpha[22]/alpha[21]",
  inter_intra1 = "alpha[12]/alpha[11]",
  inter_intra2 = "alpha[21]/alpha[22]"
)

rrc_sums_to25 <- rrc_sums %>% 
filter(T > 9, T < 26)
  
warming_tendencies200 <-
rrc_sums %>% 
ggplot() + 
  geom_point(aes(x = T, y = value, colour = statistic)) +
  facet_wrap(~parameter, 
             labeller = as_labeller(param_labels, label_parsed), 
             scales = "free")

# ggsave(plot = warming_tendencies200, filename = "figures/200C_param_trajectories.pdf", width = 12, height = 8)

warming_tendencies15 <- rrc_sums_to25 %>% 
  ggplot() + 
  geom_point(aes(x = T, y = value, colour = statistic)) +
  facet_wrap(~parameter, scales = "free")

# ggsave(plot = warming_tendencies15, filename = "figures/15C_param_trajectories.pdf", width = 12, height = 8)

# restrict transfer efficiencies to all positive values and see how that affects trajectories
# set K to mean value, to remove effects of thermal asymmetries to observe those effects. We need to understand why K and v turn back non-linearly to understand the neutrality finding

##### repeat analysis with only positive side of the v distribution ##############
rrcv <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1), #was by 0.1
                      ref_temp = 10,
                      r_EaN = 0, #draw all EAs from empirical distributions above
                      r_EaP = 0, 
                      c_Ea1N = 0,
                      c_Ea1P = 0, 
                      c_Ea2N = 0,
                      c_Ea2P = 0, 
                      K_EaN = 0, 
                      K_EaP = 0, 
                      v_EaN = sample_n(vpos_post_dist, size = 1)$intercept,
                      v_EaP = sample_n(vpos_post_dist, size = 1)$intercept, 
                      m_Ea1 = 0, 
                      m_Ea2 = 0,
                      c1N_b = 0.5, c1P_b = 1, #spec 1 consumes more P 
                      c2N_b = 1, c2P_b = 0.5, #spec 2 consumes more N 
                      r_N_b = 1, r_P_b = 0.5, #growth rate for each resource at ref temp 
                      K_N_b= 2000, K_P_b = 2000, #carrying capacity for each resource at ref temp 
                      v1N_b = 0.5, v1P_b = 1, #sp 1 converts P more efficiently 
                      v2N_b = 1, v2P_b = 0.5, #sp 2 converts N more efficiently 
                      m1_b = 0.01, m2_b = 0.01) #same for both species
  hold$iteration <- f
  rrcv <- bind_rows(rrcv, hold) 
}

#pompom
log_pomv <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 0.75, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = rrcv, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(rrcv, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(rrcv, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  # geom_point(data = rrcv_avg_new, x = 0, y = 0, colour = "black", size = 6) +
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

  # get euclidean distances for each species pair
rrcv_e <- rrcv %>% 
  filter(T %in% c(10, 25)) %>%
  dplyr::select(-c(a11:g2, m1:beta12)) %>% 
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>% 
  #thermal asymmetries
  mutate(abs_v_ta = abs(v_EaN - v_EaP)) %>% 
  #shifts in competition
  mutate(dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(rrcv_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#which thermal asymmetries drive large shifts in competition?
rrcv_ed <- rrcv_e %>% 
  filter(response_var == "dist15")

#unscaled TA - euclidean distance plot - v
rrcv_plot_e2_v <-
  rrcv_ed %>% 
  ggplot(aes(x = abs_v_ta, y = value)) +
  geom_point(aes(colour = value), size = 3) +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", colour = "red") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.45)) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[va] * "- E"[vb]*"|"))
