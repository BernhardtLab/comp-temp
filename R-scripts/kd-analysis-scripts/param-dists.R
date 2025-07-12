# This script is to generate distributions of activation energy for each MacArthur parameter based on empirical data collected by Grace Wells 2022(?)

# author: Kaleigh Davis, Postdoc with JB at U of Guelph
# script DOB: 2/12/2025

#load in necessary pkgs
library(tidyverse)
library(janitor)
library(MCMCpack)
library(bayesplot)
library(MCMCvis)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)
library(colorspace)

#get data for published parameter estimates of relevant traits
mac_means <- read_csv("data/mac-means.csv") %>% 
  clean_names()

## Get EA distributions for each parameter ######
### mortality rates -------------------------------

lm_mort <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "mortality rate"), burnin = 1000) %>% 
  as.data.frame() %>% 
  clean_names() %>% 
  mutate(parameter = "mortality_rate")

mortality_rates <- mac_means %>% 
  filter(simple_parameter == "mortality rate") %>% 
  dplyr::select(activation_energy)


# resource growth rate ----------------------------------------------------

lm_rgr <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "resource growth rate"), burnin = 1000) %>% 
  as.data.frame() %>% 
  clean_names() %>% 
  mutate(parameter = "resource_growth_rate")

growth_rates <- mac_means %>% 
  filter(simple_parameter == "resource growth rate") %>% 
  dplyr::select(activation_energy)


# conversion efficiency ---------------------------------------------------
lm_conv_eff <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "conversion efficiency"), burnin = 1000) %>% 
  as.data.frame() %>% 
  clean_names() %>% 
  mutate(parameter = "conversion_efficiency")

conv_rates <- mac_means %>% 
  filter(simple_parameter == "conversion efficiency") %>% 
  dplyr::select(activation_energy)

# resource carrying capacity  ---------------------------------------------------
lm_carrying_capacity <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "resource carrying capacity"), burnin = 1000) %>% 
  as.data.frame() %>% 
  clean_names() %>% 
  mutate(parameter = "carrying_capacity")

carrying_capacity <- mac_means %>% 
  filter(simple_parameter == "resource carrying capacity") %>% 
  dplyr::select(activation_energy)

# consumption rate  ---------------------------------------------------
lm_consumption_rate <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "consumption rate"), burnin = 1000) %>% 
  as.data.frame() %>% 
  clean_names() %>% 
  mutate(parameter = "consumption rate")

consumption_rate <- mac_means %>% 
  filter(simple_parameter == "consumption rate") %>% 
  dplyr::select(activation_energy)

#stitch all these dfs together for use in other scripts
# bind_rows(lm_mort, lm_rgr, lm_conv_eff, lm_carrying_capacity, lm_consumption_rate) %>% 
# write_csv(., "data/processed-data/param_post_dists.csv")

### plot all distributions from saved distributions #####
data <- read_csv("data/processed-data/param_post_dists.csv")

#plot distribution over original data 
mort_ea_plot <- data %>% 
  filter(parameter == "mortality_rate") %>% 
  ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
  geom_point(aes(x = activation_energy, y = 0), data = mortality_rates, color = "orange", size = 3) +
  geom_point(aes(x = activation_energy, y = 0), data = mortality_rates, color = "black", size = 3, shape = 1) +
  # geom_vline(aes(xintercept = mean(intercept))) +
  geom_vline(aes(xintercept = median(intercept)), color = "darkred") +
  coord_cartesian(ylim = c(0, 5000)) + 
  labs(y = "Count", x = "Temperature Sensitivity (eV)") +
  # coord_flip() + 
  theme_cowplot(font_size = 20) + 
  # label = expression("Conversion \nefficiency," ~ italic(v))
  annotate("text", x = 0.75, y = 4500, label = expression("Consumer \nmortality rate," ~italic(m)), size = 6)

rgr_plot <- data %>% 
  filter(parameter == "resource_growth_rate") %>% 
  ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
  geom_point(aes(x = activation_energy, y = 0), data = growth_rates, color = "orange", size = 3) +
  geom_point(aes(x = activation_energy, y = 0), data = growth_rates, color = "black", size = 3, shape = 1) +
  # geom_vline(aes(xintercept = mean(intercept))) +
  geom_vline(aes(xintercept = median(intercept)), color = "darkred") +
  coord_cartesian(ylim = c(0, 5000)) + 
  labs(y = "Count", x = "Temperature Sensitivity (eV)") +
  # coord_flip() + 
  theme_cowplot(font_size = 20) + 
  # label = expression("Conversion \nefficiency," ~ italic(v))
  annotate("text", x = 1.5, y = 4500, label = expression("Resource \ngrowth rate," ~italic(r)), size = 6)

conv_eff_plot <- data %>% 
  filter(parameter == "conversion_efficiency") %>% 
  ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
  geom_point(aes(x = activation_energy, y = 0), data = conv_rates, color = "orange", size = 3) +
  geom_point(aes(x = activation_energy, y = 0), data = conv_rates, color = "black", size = 3, shape = 1) +
  # geom_vline(aes(xintercept = mean(intercept))) +
  geom_vline(aes(xintercept = median(intercept)), color = "darkred") +
  # coord_flip() + 
  coord_cartesian(ylim = c(0, 5000)) + 
  labs(y = "Count", x = "Temperature Sensitivity (eV)") +
  theme_cowplot(font_size = 20) + 
  annotate("text", x = 1.2, y = 4500, label = expression("Conversion \nefficiency," ~italic(v)), size = 6)

carrying_capacity_plot <- data %>% 
  filter(parameter == "carrying_capacity") %>% 
  ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
  geom_point(aes(x = activation_energy, y = 0), data = carrying_capacity, color = "orange", size = 3) +
  geom_point(aes(x = activation_energy, y = 0), data = carrying_capacity, color = "black", size = 3, shape = 1) +
  # geom_vline(aes(xintercept = mean(intercept))) +
  geom_vline(aes(xintercept = median(intercept)), color = "darkred") +
  # coord_flip() + 
  coord_cartesian(ylim = c(0, 5000)) + 
  labs(y = "Count", x = "Temperature Sensitivity (eV)") +
  theme_cowplot(font_size = 20) + 
  annotate("text", x = 0.05, y = 4000, label = expression("Resource \ncarrying \ncapacity," ~italic(K)), size = 6)

consumption_rate_plot <- data %>% 
  filter(parameter == "consumption rate") %>% 
  ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
  geom_point(aes(x = activation_energy, y = 0), data = consumption_rate, color = "orange", size = 3) +
  geom_point(aes(x = activation_energy, y = 0), data = consumption_rate, color = "black", size = 3, shape = 1) +
  # geom_vline(aes(xintercept = mean(intercept))) +
  geom_vline(aes(xintercept = median(intercept)), color = "darkred") +
  # coord_flip() + 
  coord_cartesian(ylim = c(0, 5000)) + 
  labs(y = "Count", x = "Temperature Sensitivity (eV)") +
  theme_cowplot(font_size = 20) + 
  annotate("text", x = 1.3, y = 4500, label = expression("Consumption \nrate," ~italic(c)), size = 6)


ea_plots <-
  consumption_rate_plot + rgr_plot + carrying_capacity_plot + conv_eff_plot + mort_ea_plot +
  plot_annotation(tag_levels = "A")
ggsave(filename = "figures/kd-figs/ea-plots1.pdf", ea_plots, width = 15, height = 12)
