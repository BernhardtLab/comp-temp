# This is a daughter script of script 06 to explore the reviewers' analysis suggestion of exploring more starting points.

# This script is to reproduce the main results figures of the paper -- figs 3 and 5 -- with different start points. We achieve these different start points by manipulating the resource preferences (analogous to niche widths) of each species at the models starting conditions, i.e. ambient temperature. Scenario 1 is given in the main text (script 04) and features two consumers that each specialize on one of two resources and they have equally strong preference for this resource. The two resources have uneven growth rates under ambient conditions, which places the species pair on the boundary of coexistence under ambient conditions. In scenario 2, consumer preferences are the same as in Scenario 1, but the resources grow at the same rate under ambient conditions, which moves the start point from the coexistence boundary to the middle of the coexistence region. In scenario 3, species have very different resource preferences, where one consumer has a strong preference for its preferred resource and the other consumer has a weak preference for its preferred resources. Both resources grow at the same rate under ambient conditions. Parameters defining each of these starting conditions are given in each simulation below, and in summary in Table S1. In each simulation, each MacArthur consumer-resource parameter is given by an Arrhenius function, with a temperature sensitivity (activation energy, slope) term and an intercept term, which determines the value of the function at ambient temperatures (Tref, ref temp). In each simulation, temperature sensitivities are defined as "{parameter_EAik}", where ik captures the relevant consumer, resource, or both, and intercepts are defined as "{parameter-ik_b}". Consumers are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively, throughout the script. The script simulates the effects of warming when each parameter is given a temperature sensitivity, randomly drawn from the parameter's empirical distribution (generated in 01-param-dists), simultaneously.

#script DOB: 12/15/2025
#author: Kaleigh Davis, Postdoc UoG with Joey Bernhardt

#### packages and referencing #####
#load necessary pkgs
library(tidyverse)
library(janitor)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)
library(purrr)
library(viridis)
library(beepr)

# get referencing set up for MacArthur temp dependence function
source("R-scripts/02-temp-dep-macarthur.R") #this contains the MacArthur translation function, with all parameters flexibly defined in the function for assigning at time of use
source("R-scripts/03-arrhenius.R") #arrhenius function

#load in distributions for parameter values.
# these are continuous distributions generated from empirical data using MCMC regression, in 01-param-dists.R
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

#### POMPOM; uneven reciprocal preference, equal growth rates - scenario 3 from first submission ####
urrce_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
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
                      c1N_b = 0.1, c1P_b = 0.9, 
                      c2N_b = 0.6, c2P_b = 0.4, 
                      r_N_b = 0.5, r_P_b = 0.5, 
                      K_N_b= 2000, K_P_b = 2000, 
                      v1N_b = 0.1, v1P_b = 0.9, 
                      v2N_b = 0.6, v2P_b = 0.4, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce_all <- bind_rows(urrce_all, hold) 
}

#get average change in position after 15C warming
urrce_all_avg_new <- urrce_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.35, y = -0.08, label = "Coexistence", size = 9, fontface = 2) +
  annotate("text", x = 0.05, y = -0.2, label = "Species 1 wins", size = 9, fontface = 2) +
  annotate("text", x = 0.05, y = 0.7, label = "Species 2 wins", size = 9, fontface = 2) +
  annotate("text", x = -0.015, y = 0.05, label = "Neutrality", size = 9, fontface = 2) +
  theme_cowplot(font_size = 26)

#repeat pompom, but without annotations
#pompom
log_pom_urrce_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  coord_cartesian(ylim = c(-0.27, 0.8), xlim = c(-0.022, 1.02)) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  theme_cowplot(font_size = 26) +
  ggtitle("Full model")

#for how many simulations does ND decrease and does FD decrease

urrce_start <- urrce_all %>% 
  filter(T == 10) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_sd_stab_pot = sd(new_stabil_potential),
            new_sd_fit_rat = sd(new_fit_ratio)) %>% 
  select(new_mean_stab_pot, new_mean_fit_rat)

urrce_shifts <- urrce_all %>% 
  select(iteration, T, new_stabil_potential, new_fit_ratio) %>% 
  mutate(stab_pot_start = urrce_start$new_mean_stab_pot,
         fit_rat_start = urrce_start$new_mean_fit_rat) %>% 
  mutate(fd_shift = ifelse(new_fit_ratio < fit_rat_start, "decrease", 
                           ifelse(new_fit_ratio > fit_rat_start, "increase",
                                  ifelse(new_fit_ratio == fit_rat_start, "no change", "potato"))),
         nd_shift = ifelse(new_stabil_potential < stab_pot_start, "decrease",
                           ifelse(new_stabil_potential > stab_pot_start, "increase",
                                  ifelse(new_stabil_potential == stab_pot_start, "no change","potato")))) %>% 
  filter(T == 10 | T == 25)

urrce_shifts %>% 
  group_by(T, fd_shift, nd_shift) %>% 
  tally()

#FD decrease in 253/500 cases; 2nd sim: 254/500; 3rd:279/500
#ND decrease in 250/500 cases; 2nd sim: 255/500; 3rd: 263/500

# calculate euclidean distances at 25C for each iteration
urrce_all_e <- urrce_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end? No, great.
nrow(urrce_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential))

#which thermal asymmetries drive large shifts in competition?
urrce_all_ed <- urrce_all_e %>% 
  filter(response_var == "dist15")

#unscaled TA - euclidean distance plot - r
urrce_all_plot_e2_r <-
  urrce_all_ed %>% 
  ggplot(aes(x = abs_r_ta, y = value)) +
  geom_point(aes(colour = r_EaP > r_EaN), size = 3) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) + 
  annotate("text", x = 0.8, y = 0.35, label = "Resource \ngrowth rate, r", size = 5.5) + 
  annotate("text", x = 0.8, y = 0.15, label = "m = 0.62, \np = <0.001, adj. r2 = 0.13", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[ra] * "- E"[rb]*"|"))

summary(lm(value ~ abs_r_ta, data = urrce_all_ed)) #S

#unscaled TA - euclidean distance plot 
urrce_all_plot_e2_c <-
  urrce_all_ed %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_c1_ta, y = value)) +
  geom_point(aes(colour = c_Ea1P > c_Ea2N), size = 3) + 
  geom_smooth(method = "lm", linetype = "dashed", se = F, colour = "black") +
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.8, y = 0.35, label = "Consumption \nrate, c", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[c1a] * "- E"[c2b]*"|"))

#none significant
summary(lm(value~abs_c1_ta, data = urrce_all_ed))
summary(lm(value~abs_c2_ta, data = urrce_all_ed))
summary(lm(value~abs_c3_ta, data = urrce_all_ed))
summary(lm(value~abs_c4_ta, data = urrce_all_ed))
summary(lm(value~abs_c5_ta, data = urrce_all_ed))
summary(lm(value~abs_c6_ta, data = urrce_all_ed)) #none significant

urrce_all_plot_e2_k <-
  urrce_all_ed %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_k_ta, y = value)) +
  geom_point(aes(colour = K_EaN > K_EaP), size = 3) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.97, y = 0.35, label = "Carrying \ncapacity, K", size = 5.5) + 
  annotate("text", x = 0.97, y = 0.15, label = "m = 0.24, \np = <0.001, \nadj. r2 = 0.07", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[Ka] * "- E"[Kb]*"|"))

summary(lm(value~abs_k_ta, data = urrce_all_ed)) #S

urrce_all_plot_e2_v <-
  urrce_all_ed %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_v_ta, y = value)) +
  geom_point(aes(colour = v_EaN > v_EaP), size = 3) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.92, y = 0.25, label = "Conversion \nefficiency, v", size = 5.5) + 
  annotate("text", x = 0.92, y = 0.1, label = "m = 0.32, \np = <0.001,\nadj. r2 = 0.25", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[va] * "- E"[vb]*"|"))

summary(lm(value~abs_v_ta, data = urrce_all_ed)) #S

urrce_all_plot_e2_m <-
  urrce_all_ed %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_m_ta, y = value)) +
  geom_point(aes(colour = m_Ea1 > m_Ea2), size = 3) + 
  geom_smooth(method = "lm", linetype = "dashed", se = F, colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.8, y = 0.5, label = "Mortality \nrate, m", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[m1] * "- E"[m2]*"|"))

summary(lm(value~abs_m_ta, data = urrce_all_ed)) #NS
hist(urrce_all_ed$value)

urrce_tas <- urrce_all_plot_e2_c +urrce_all_plot_e2_r + urrce_all_plot_e2_k + urrce_all_plot_e2_v +  urrce_all_plot_e2_m + log_pom_urrce_all_noanno + plot_annotation(tag_levels = "A", title = "Thermal asymmetry links to shifts in competition: Example 1")

# ggsave(plot = urrce_tas, filename = "figures/TA-ED-full1.pdf", width = 18, height = 12)


#### POMPOM 1; uneven reciprocal preference, unequal growth rates ####
urrce1_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
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
                      c1N_b = 0.6, c1P_b = 0.9, 
                      c2N_b = 0.5, c2P_b = 0.2, 
                      r_N_b = 0.7, r_P_b = 0.2, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.3, v1P_b = 0.9, 
                      v2N_b = 0.3, v2P_b = 0.1, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce1_all <- bind_rows(urrce1_all, hold) 
}

#get average change in position after 15C warming
urrce1_all_avg_new <- urrce1_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce1_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce1_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce1_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce1_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce1_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce1_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce1_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = -0.8, 
           label = "c1N_b = 0.6, c1P_b = 0.9, 
c2N_b = 0.5, c2P_b = 0.2, 
r_N_b = 0.7, r_P_b = 0.2, 
K_N_b= 1000, K_P_b = 1000, 
v1N_b = 0.3, v1P_b = 0.9, 
v2N_b = 0.3, v2P_b = 0.1, 
m1_b = 0.01, m2_b = 0.01", 
           size = 9, fontface = 1, hjust = 0) 

#make the pompom again, but without annotations
#pompom
log_pom_urrce1_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce1_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce1_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce1_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce1_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce1_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce1_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  ggtitle("Full model")

# calculate euclidean distances at 25C for each iteration
urrce1_all_e <- urrce1_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end?
nrow(urrce1_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential)) #No, great.

urrce1_all_ed <- urrce1_all_e %>% 
  filter(response_var == "dist15")

#which thermal asymmetries drive large shifts in competition?
#unscaled TA - euclidean distance plot - r
urrce1_all_plot_e2_r <-
  urrce1_all_ed %>% 
  ggplot(aes(x = abs_r_ta, y = value)) +
  geom_point(aes(colour = r_EaP > r_EaN), size = 3) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) + 
  annotate("text", x = 0.8, y = 0.35, label = "Resource \ngrowth rate, r", size = 5.5) + 
  annotate("text", x = 0.8, y = 0.15, label = "m = 0.87, \np = <0.001, r2 = 0.41", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[ra] * "- E"[rb]*"|"))

summary(lm(value ~ abs_r_ta, data = urrce1_all_ed)) #S

#unscaled TA - euclidean distance plot 
urrce1_all_plot_e2_c <-
  urrce1_all_ed %>%
  ggplot(aes(x = abs_c1_ta, y = value)) +
  geom_point(aes(colour = c_Ea1P > c_Ea1N), size = 3) + 
  geom_smooth(method = "lm", colour = "black") +
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)", title = expression("|E"[c1a] * "- E"[c1b]*"|*"), subtitle = "*Only this pair of consumption rate thermal \nasymmetries had a significant effect on \nEuclidean distance. All other \npairs not significant.") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.8, y = 0.35, label = "Consumption \nrate, c", size = 5.5) + 
  annotate("text", x = 0.8, y = 0.15, label = "m = 0.47, \np = 0.008, r2 = 0.011", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") 

summary(lm(value~abs_c1_ta, data = urrce1_all_ed)) #NS
summary(lm(value~abs_c2_ta, data = urrce1_all_ed)) #S
summary(lm(value~abs_c3_ta, data = urrce1_all_ed)) #NS
summary(lm(value~abs_c4_ta, data = urrce1_all_ed)) #NS
summary(lm(value~abs_c5_ta, data = urrce1_all_ed)) #NS
summary(lm(value~abs_c6_ta, data = urrce1_all_ed)) #NS

urrce1_all_plot_e2_k <-
  urrce1_all_ed %>% 
  ggplot(aes(x = abs_k_ta, y = value)) +
  geom_point(aes(colour = K_EaN > K_EaP), size = 3) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.97, y = 0.78, label = "Carrying \ncapacity, K", size = 5.5) + 
  annotate("text", x = 0.97, y = 0.55, label = "m = 0.14, \np = < 0.001, \nr2 = 0.037", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[Ka] * "- E"[Kb]*"|"))

summary(lm(value~abs_k_ta, data = urrce1_all_ed)) #S

urrce1_all_plot_e2_v <-
  urrce1_all_ed %>% 
  ggplot(aes(x = abs_v_ta, y = value)) +
  geom_point(aes(colour = v_EaN > v_EaP), size = 3) +
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.92, y = 0.78, label = "Conversion \nefficiency, v", size = 5.5) + 
  annotate("text", x = 0.92, y = 0.6, label = "m = 0.18, \np = <0.001, r2 = 0.11", size = 5.5) +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[va] * "- E"[vb]*"|"))

summary(lm(value~abs_v_ta, data = urrce1_all_ed)) #S
hist(urrce1_all_ed$value)

urrce1_all_plot_e2_m <-
  urrce1_all_ed %>% 
  ggplot(aes(x = abs_m_ta, y = value)) +
  geom_point(aes(colour = m_Ea1 > m_Ea2), size = 3) + 
  geom_smooth(method = "lm", linetype = "dashed", se = F, colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.8, y = 0.5, label = "Mortality \nrate, m", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[m1] * "- E"[m2]*"|"))

summary(lm(value~abs_m_ta, data = urrce1_all_ed)) #NS
hist(urrce1_all_ed$value)

urrce1_tas <- urrce1_all_plot_e2_c + urrce1_all_plot_e2_r + urrce1_all_plot_e2_k +  urrce1_all_plot_e2_v + urrce1_all_plot_e2_m + log_pom_urrce1_all_noanno + plot_annotation(tag_levels = "A", title = "Thermal asymmetry links to shifts in competition: Example 2")

ggsave(plot = urrce1_tas, filename = "figures/TA-ED-full2.pdf", width = 18, height = 12)

#### POMPOM 2; uneven reciprocal preference, unequal growth rates ####
urrce2_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
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
                      c1N_b = 0.5, c1P_b = 0.2, 
                      c2N_b = 0.9, c2P_b = 0.9, 
                      r_N_b = 0.2, r_P_b = 0.7, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.1, v1P_b = 0.6, 
                      v2N_b = 0.4, v2P_b = 0.7, 
                      m1_b = 0.01, m2_b = 0.01) 
  hold$iteration <- f
  urrce2_all <- bind_rows(urrce2_all, hold) 
}

#get average change in position after 15C warming
urrce2_all_avg_new <- urrce2_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce2_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce2_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce2_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce2_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce2_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce2_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce2_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = 0.8, 
           label = "c1N_b = 0.5, c1P_b = 0.2, 
c2N_b = 0.9, c2P_b = 0.9, 
r_N_b = 0.2, r_P_b = 0.7, 
K_N_b= 1000, K_P_b = 1000, 
v1N_b = 0.1, v1P_b = 0.6, 
v2N_b = 0.4, v2P_b = 0.7, 
m1_b = 0.01, m2_b = 0.01", 
           size = 9, fontface = 1, hjust = 0) 

#### POMPOM 3; uneven reciprocal preference, unequal growth rates ####
urrce3_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
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
                      c1N_b = 1.3, c1P_b = 0.9, 
                      c2N_b = 0.9, c2P_b = 0.2, 
                      r_N_b = 0.7, r_P_b = 2, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.9, 
                      v2N_b = 0.9, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce3_all <- bind_rows(urrce3_all, hold) 
}

#get average change in position after 15C warming
urrce3_all_avg_new <- urrce3_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce3_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce3_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce3_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce3_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce3_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce3_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce3_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = 0.8, 
           label = "c1N_b = 1.3, c1P_b = 0.9, 
c2N_b = 0.9, c2P_b = 0.2, 
r_N_b = 0.7, r_P_b = 2, 
K_N_b= 1000, K_P_b = 1000, 
v1N_b = 0.9, v1P_b = 0.9, 
v2N_b = 0.9, v2P_b = 0.9, 
m1_b = 0.1, m2_b = 0.1", 
           size = 9, fontface = 1, hjust = 0) 

#make the pompom again, but without annotations
#pompom
log_pom_urrce3_all_noanno <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce3_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce3_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce3_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce3_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce3_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce3_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  ggtitle("Full model")

# calculate euclidean distances at 25C for each iteration
urrce3_all_e <- urrce3_all %>% 
  filter(T %in% c(10, 25)) %>% 
  dplyr::select(-c(a11:g2, coexist:beta12)) %>%
  pivot_wider(id_cols = c(ref_temp:m2_b, iteration),
              names_from = T,
              values_from = c(new_stabil_potential, new_fit_ratio),
              names_glue = "T{T}_{.value}") %>%
  mutate(abs_r_ta = abs(r_EaN - r_EaP),
         r_ta = r_EaN - r_EaP,
         abs_c1_ta = abs(c_Ea1P - c_Ea2N),
         abs_c2_ta = abs(c_Ea1P - c_Ea1N),
         abs_c3_ta = abs(c_Ea1P - c_Ea2P), 
         abs_c4_ta = abs(c_Ea2P - c_Ea2N), 
         abs_c5_ta = abs(c_Ea2P - c_Ea1N), 
         abs_c6_ta = abs(c_Ea1N - c_Ea2N), 
         abs_k_ta = abs(K_EaN - K_EaP),
         abs_v_ta = abs(v_EaN - v_EaP),
         abs_m_ta = abs(m_Ea1 - m_Ea2),
         dist15 = sqrt((T25_new_stabil_potential - T10_new_stabil_potential)^2 + (T25_new_fit_ratio - T10_new_fit_ratio)^2),
         shift_fitrat = T25_new_fit_ratio - T10_new_fit_ratio,
         shift_nichediffs = T25_new_stabil_potential - T10_new_stabil_potential) %>% 
  pivot_longer(cols = c(dist15, shift_fitrat, shift_nichediffs), names_to = "response_var", values_to = "value") 

#anywhere where ND and FD are exactly the same at end?
nrow(urrce3_all_e %>% filter(T == 25 & T25_new_fit_ratio == T25_new_stabil_potential)) #No, great.

#which thermal asymmetries drive large shifts in competition?
urrce3_all_ed <- urrce3_all_e %>% 
  filter(response_var == "dist15")

#unscaled TA - euclidean distance plot - r
urrce3_all_plot_e2_r <-
  urrce3_all_ed %>% 
  ggplot(aes(x = abs_r_ta, y = value)) +
  geom_point(aes(colour = r_EaP > r_EaN), size = 3) + 
  geom_smooth(method = "lm", colour = "black", linetype = "dashed", se = F) + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) + 
  annotate("text", x = 0.8, y = 0.35, label = "Resource \ngrowth rate, r", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[ra] * "- E"[rb]*"|"))

summary(lm(value ~ abs_r_ta, data = urrce3_all_ed)) #NS

#unscaled TA - euclidean distance plot 
urrce3_all_plot_e2_c <-
  urrce3_all_ed %>% 
  ggplot(aes(x = abs_c1_ta, y = value)) +
  geom_point(aes(colour = c_Ea1P > c_Ea2N), size = 3) + 
  geom_smooth(method = "lm", linetype = "dashed", se = F, colour = "black") +
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.8, y = 0.35, label = "Consumption \nrate, c", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[c1a] * "- E"[c2b]*"|"))

#none significant
summary(lm(value~abs_c1_ta, data = urrce3_all_ed))
summary(lm(value~abs_c2_ta, data = urrce3_all_ed))
summary(lm(value~abs_c3_ta, data = urrce3_all_ed))
summary(lm(value~abs_c4_ta, data = urrce3_all_ed))
summary(lm(value~abs_c5_ta, data = urrce3_all_ed))
summary(lm(value~abs_c6_ta, data = urrce3_all_ed)) #none significant

urrce3_all_plot_e2_k <-
  urrce3_all_ed %>% 
  ggplot(aes(x = abs_k_ta, y = value)) +
  geom_point(aes(colour = K_EaN > K_EaP), size = 3) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.97, y = 0.75, label = "Carrying \ncapacity, K", size = 5.5) + 
  annotate("text", x = 0.92, y = 0.55, label = "m = 0.12, \np = <0.001, r2 = 0.064", size = 5.5) +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[Ka] * "- E"[Kb]*"|"))

summary(lm(value~abs_k_ta, data = urrce3_all_ed)) #S

urrce3_all_plot_e2_v <-
  urrce3_all_ed %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_v_ta, y = value)) +
  geom_point(aes(colour = v_EaN > v_EaP), size = 3) + 
  geom_smooth(method = "lm", colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.92, y = 0.75, label = "Conversion \nefficiency, v", size = 5.5) + 
  annotate("text", x = 0.92, y = 0.55, label = "m = 0.15, \np = <0.001, r2 = 0.22", size = 5.5) +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[va] * "- E"[vb]*"|"))

summary(lm(value~abs_v_ta, data = urrce3_all_ed)) #S


urrce3_all_plot_e2_m <-
  urrce3_all_ed %>% 
  filter(response_var == "dist15") %>% 
  ggplot(aes(x = abs_m_ta, y = value)) +
  geom_point(aes(colour = m_Ea1 > m_Ea2), size = 3) + 
  geom_smooth(method = "lm", linetype = "dashed", se = F, colour = "black") + 
  labs(x = "Magnitude \nof thermal asymmetry", y = "Displacement of species pair with \nwarming (Euclidean distance)") + 
  coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 0.8)) +
  annotate("text", x = 0.8, y = 0.5, label = "Mortality \nrate, m", size = 5.5) + 
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  ggtitle(expression("|E"[m1] * "- E"[m2]*"|"))

summary(lm(value~abs_m_ta, data = urrce3_all_ed)) #NS

urrce3_tas <- urrce3_all_plot_e2_c + urrce3_all_plot_e2_r +  urrce3_all_plot_e2_k +  urrce3_all_plot_e2_v + urrce3_all_plot_e2_m + log_pom_urrce3_all_noanno + plot_annotation(tag_levels = "A", title = "Thermal asymmetry links to shifts in competition: Example 3")

# ggsave(plot = urrce3_tas, filename = "figures/TA-ED-full3.pdf", width = 18, height = 12)

#### POMPOM 4; uneven reciprocal preference, unequal growth rates ####
urrce4_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
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
                      c1N_b = 1.3, c1P_b = 0.4, 
                      c2N_b = 0.9, c2P_b = 0.4, 
                      r_N_b = 0.7, r_P_b = 1, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.4, 
                      v2N_b = 0.9, v2P_b = 0.4, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce4_all <- bind_rows(urrce4_all, hold) 
}

#get average change in position after 15C warming
urrce4_all_avg_new <- urrce4_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce4_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce4_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce4_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce4_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce4_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce4_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce4_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = -0.8, 
           label = "c1N_b = 1.3, c1P_b = 0.4, 
c2N_b = 0.9, c2P_b = 0.4, 
r_N_b = 0.7, r_P_b = 1, 
K_N_b= 1000, K_P_b = 1000, 
v1N_b = 0.9, v1P_b = 0.4, 
v2N_b = 0.9, v2P_b = 0.4, 
m1_b = 0.1, m2_b = 0.1", 
           size = 9, fontface = 1, hjust = 0) 

#### POMPOM 5; uneven reciprocal preference, unequal growth rates ####
urrce5_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
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
                      c1N_b = 1.3, c1P_b = 0.8, 
                      c2N_b = 0.3, c2P_b = 2, 
                      r_N_b = 0.8, r_P_b = 0.8, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.8, 
                      v2N_b = 0.1, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce5_all <- bind_rows(urrce5_all, hold) 
}

#get average change in position after 15C warming
urrce5_all_avg_new <- urrce5_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce5_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce5_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce5_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce5_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce5_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce5_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce5_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  theme_cowplot(font_size = 28) + 
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = 0.8, 
           label = "c1N_b = 1.3, c1P_b = 0.8, 
c2N_b = 0.3, c2P_b = 2, 
r_N_b = 0.8, r_P_b = 0.8, 
K_N_b= 1000, K_P_b = 1000, 
v1N_b = 0.9, v1P_b = 0.8, 
v2N_b = 0.1, v2P_b = 0.9, 
m1_b = 0.1, m2_b = 0.1", 
           size = 9, fontface = 1, hjust = 0)


#### POMPOM 6; uneven reciprocal preference, unequal growth rates ####
urrce6_all <- data.frame()
for(f in 1:500){ 
  hold = temp_dep_mac(T = seq(10, 25, by = 0.1),
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
                      c1N_b = 1.3, c1P_b = 0.6, 
                      c2N_b = 0.3, c2P_b = 2, 
                      r_N_b = 0.8, r_P_b = 1.5, 
                      K_N_b= 1000, K_P_b = 1000, 
                      v1N_b = 0.9, v1P_b = 0.8, 
                      v2N_b = 0.1, v2P_b = 0.9, 
                      m1_b = 0.1, m2_b = 0.1) 
  hold$iteration <- f
  urrce6_all <- bind_rows(urrce6_all, hold) 
}

#get average change in position after 15C warming
urrce6_all_avg_new <- urrce6_all %>% 
  mutate(rel_T = T-10) %>% 
  filter(rel_T == 15) %>% 
  group_by(rel_T) %>% 
  summarise(new_mean_stab_pot = mean(new_stabil_potential),
            new_mean_fit_rat = mean(new_fit_ratio),
            new_med_stab_pot = median(new_stabil_potential),
            new_med_fit_rat = median(new_fit_ratio))

#pompom
log_pom_urrce6_all <-
  ggplot() +
  # coexist area
  geom_ribbon(data = data.frame(x = seq(0, 1.25, 0.001)),
              aes(x = x,
                  y = NULL,
                  ymin = -x,
                  ymax = x),
              fill = "grey", color = "black", alpha = 0.2) +
  # sim paths
  geom_path(data = urrce6_all, aes(x = new_stabil_potential, y = new_fit_ratio, color = T-10, group = iteration), linewidth = 3) +
  # position before warming
  geom_point(data = filter(urrce6_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio), colour = "black", size = 7.5) +
  geom_point(data = filter(urrce6_all, T==10), aes(x = new_stabil_potential, y = new_fit_ratio, colour = T-10), size = 6) +
  # position after 15C warming
  geom_point(data = urrce6_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat), colour = "black",  size = 7.5) +
  geom_point(data = urrce6_all_avg_new, aes(x = new_med_stab_pot, y = new_med_fit_rat, colour = rel_T),  size = 6) +
  geom_hline(yintercept = 0, linetype=5) +
  geom_point(data = urrce6_all, x = 0, y = 0, colour = "black", size = 6) +
  #aesthetic customization
  scale_colour_viridis_c(option = "magma", begin = 0.53, end = 1, direction = -1) +
  xlab(expression(paste("Niche differences (-log(", rho, "))"))) +
  ylab(expression(paste("Fitness differences (log(", f[2], "/", f[1], "))"))) +
  labs(colour = "°C Warming") +
  theme_cowplot(font_size = 26) + 
  annotate("text", x = 0, y = -0.8, 
           label = "c1N_b = 1.3, c1P_b = 0.6,
c2N_b = 0.3, c2P_b = 2,
r_N_b = 0.8, r_P_b = 1.5,
K_N_b= 1000, K_P_b = 1000,
v1N_b = 0.9, v1P_b = 0.8,
v2N_b = 0.1, v2P_b = 0.9,
m1_b = 0.1, m2_b = 0.1", 
           size = 9, fontface = 1, hjust = 0) 

diff_start_points <- (log_pom_urrce2_all + log_pom_urrce3_all +  log_pom_urrce5_all) / (log_pom_urrce4_all + log_pom_urrce1_all + log_pom_urrce6_all) + plot_annotation(tag_levels = "A")

# ggsave(plot = diff_start_points, filename = "figures/extra-start-points.pdf", height = 24, width = 30)
