#this script is a daughter script to script 01, and creates separate datasets for unicells and multicells in order to identify systematic differences in activation energies, and test for their effects

# This script is to generate posterior distributions of temperature sensitivity for each MacArthur parameter from synthesized empirical data. This scripts output is used in each subsequent analysis (scripts 03, 04, 05).

# author: Kaleigh Davis, Postdoc with Joey Bernhardt at U of Guelph
# script DOB: 1/7/2025

# inputs:
#   data/raw-data/param-eas.csv — synthesized published estimates of temperature
#     sensitivities (activation energies) for MacArthur model parameters

# outputs:
#   data/processed-data/unimulti_param_post_dists.csv — posterior distributions
#     of temperature sensitivity split by cellularity (unicellular vs. multicellular)
#   figures/unimulti-ea-plots1.pdf — supplementary figure S7

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

# get data for published parameter estimates of relevant traits
mac_means <- read_csv("data/raw-data/param-eas.csv") %>%
  clean_names()
#here, temperature sensitivity is stored in a column called "activation_energy" in accordance with the Arrhenius model

# create uni/multi columns
#isolate the ones with not enough info
mac_means_nas <- mac_means %>% 
  filter(is.na(general_taxon_grouping_1)) %>% 
  dplyr::select(general_taxon_grouping_1:species)

mac_cell <- mac_means %>% 
  mutate(cellularity = ifelse(general_taxon_grouping_1 == "Bacteria", "uni",
                              ifelse(general_taxon_grouping_1 %in% c("Invertebrates", "Vertebrates"), "multi", "potato")))

mac_cell1 <- mac_cell %>% 
  mutate(cellularity1 = ifelse(!is.na(cellularity), cellularity,
                               ifelse(general_taxon_grouping_3 %in% c("Algae", "Cilliate", "Phytoplankton", "Protist"), "uni",
                                      ifelse(str_detect(general_taxon_grouping_3, "Duckweed") | str_detect(general_taxon_grouping_3, "zooplankton"), "multi", "potato"))))

#double check
mac_cell1 %>% 
   distinct(species, general_taxon_grouping_1, general_taxon_grouping_2, general_taxon_grouping_3, cellularity1) %>% view()

mac_cell2 <- mac_cell1 %>% 
  dplyr::select(-cellularity) %>% 
  rename(cellularity = cellularity1)

mac_cell2 %>% 
  ungroup() %>% 
  group_by(cellularity, simple_parameter) %>% 
  tally()

## Generate EA distributions for each parameter ######
### mortality rates -------------------------------
lm_mort <- mac_cell2 %>%
  filter(simple_parameter == "mortality rate") %>% 
  group_by(cellularity) %>% #no unicellular estimates
  group_modify(~ {
    MCMCregress(
      activation_energy ~ 1,
      data = .x,
      burnin = 1000
    ) %>%
      as.data.frame() %>%
      clean_names()
  }) %>%
  ungroup() %>%
  mutate(parameter = "mortality_rate")

mortality_rates <- mac_cell2 %>% 
  filter(simple_parameter == "mortality rate") %>% 
  dplyr::select(activation_energy, cellularity)

mortality_rates %>% 
  group_by(cellularity) %>% 
  summarise(min = min(activation_energy),
            max = max(activation_energy))

# resource growth rate ----------------------------------------------------
lm_rgr <- mac_cell2 %>%
  filter(simple_parameter == "resource growth rate") %>% 
  group_by(cellularity) %>%
  group_modify(~ {
    MCMCregress(
      activation_energy ~ 1,
      data = .x,
      burnin = 1000
    ) %>%
      as.data.frame() %>%
      clean_names()
  }) %>%
  slice_sample(n = 5000) %>% #sample 5000 per cellularity group so total distribution is the same size as other parameters
  ungroup() %>%
  mutate(parameter = "resource_growth_rate")

growth_rates <- mac_cell2 %>% 
  filter(simple_parameter == "resource growth rate") %>% 
  dplyr::select(cellularity, activation_energy)

growth_rates %>% 
  group_by(cellularity) %>% 
  summarise(min = min(activation_energy),
            max = max(activation_energy))

# conversion efficiency ---------------------------------------------------
lm_conv_eff <- mac_cell2 %>%
  filter(simple_parameter == "conversion efficiency") %>% 
  # group_by(cellularity) %>% #this may not be appropriate -- come back to this after I run the sim
  group_modify(~ {
    MCMCregress(
      activation_energy ~ 1,
      data = .x,
      burnin = 1000
    ) %>%
      as.data.frame() %>%
      clean_names()
  }) %>%
  ungroup() %>%
  mutate(parameter = "conversion_efficiency")

conv_rates <- mac_cell2 %>% 
  filter(simple_parameter == "conversion efficiency") %>% 
  dplyr::select(activation_energy, cellularity)

conv_rates %>% 
  group_by(cellularity) %>% 
  summarise(min = min(activation_energy),
            max = max(activation_energy))

# resource carrying capacity  ---------------------------------------------------
lm_carrying_capacity <- mac_cell2 %>%
  filter(simple_parameter == "resource carrying capacity") %>% 
  # group_by(cellularity) %>% #only one estimate for multicellular
  group_modify(~ {
    MCMCregress(
      activation_energy ~ 1,
      data = .x,
      burnin = 1000
    ) %>%
      as.data.frame() %>%
      clean_names()
  }) %>%
  ungroup() %>%
  mutate(parameter = "carrying_capacity")

carrying_capacity <- mac_cell2 %>% 
  filter(simple_parameter == "resource carrying capacity") %>% 
  dplyr::select(activation_energy, cellularity)

carrying_capacity %>%
  group_by(cellularity) %>% 
  summarise(min = min(activation_energy),
            max = max(activation_energy))

# consumption rate  ---------------------------------------------------
lm_consumption_rate <- mac_cell2 %>%
  filter(simple_parameter == "consumption rate") %>% 
  # group_by(cellularity) %>% #only 1 unicellular observation
  group_modify(~ {
    MCMCregress(
      activation_energy ~ 1,
      data = .x,
      burnin = 1000
    ) %>%
      as.data.frame() %>%
      clean_names()
  }) %>%
  ungroup() %>%
  mutate(parameter = "consumption rate")

consumption_rate <- mac_cell2 %>% 
  filter(simple_parameter == "consumption rate") %>% 
  dplyr::select(activation_energy, cellularity)

consumption_rate %>% 
  group_by(cellularity) %>% 
  summarise(min = min(activation_energy),
            max = max(activation_energy))

# stitch all these dfs together for use in other scripts
bind_rows(lm_mort, lm_rgr, lm_conv_eff, lm_carrying_capacity, lm_consumption_rate) %>% 
write_csv(., "data/processed-data/unimulti_param_post_dists.csv")

### get summary stats for each parameter #####
#in order to prevent small changes to parameter distribution values with each simulation, I do not re-run the regressions each time.
data <- read_csv("data/processed-data/unimulti_param_post_dists.csv") 
  # filter(!(parameter == "carrying_capacity" & cellularity == "multi")) %>% 
  # filter(!(parameter == "consumption rate" & cellularity == "uni")) #this is hacky, but the best way forward right now. One alternative is to use the single value for the other rate as a stand in. This is possible for the two params with only one cellularity level present, as well.
data_all <- read_csv("data/processed-data/param_post_dists.csv")

param_sum <- data %>%
  group_by(parameter, cellularity) %>% 
  summarize(
    across(
      intercept,
      list(
        Mean = mean,
        ci_low = ~quantile(., 0.025),
        ci_up = ~quantile(., 0.975),
        Q1 = ~quantile(., 0.25),
        Median = median,
        Q3 = ~quantile(., 0.75),
        Min = min,
        Max = max,
        sd = sd
      ),
      .names = "{.fn}" )
  ) 

### plot all distributions  #####
# plot distribution over original data 
mort_ea_plot <-
  ggplot() + 
  # geom_density(aes(x = intercept, y = after_stat(density * n), group = cellularity), fill = "#F8766D", alpha = 0.5, data = filter(data, parameter == "mortality_rate")) + 
  geom_density(aes(x = intercept, group = cellularity), fill = "#F8766D", alpha = 0.5, data = filter(data, parameter == "mortality_rate")) + 
  geom_point(aes(x = activation_energy, y = 0), data = mortality_rates, color = "#F8766D", size = 5, alpha = 0.5) +
  geom_point(aes(x = activation_energy, y = 0), data = mortality_rates, color = "#252525", size = 5, shape = 1) +
  geom_vline(aes(xintercept = mean(intercept)), color = "#F8766D", data = filter(data, parameter == "mortality_rate")) +
  # coord_cartesian(xlim = c(-1.5, 2)) +
  coord_cartesian(ylim = c(0, 12), xlim = c(-1.5, 2)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12)) +
  labs(y = "Density", x = "Temperature Sensitivity (eV)") +
  theme_cowplot(font_size = 20) +
  annotate("text", x = -0.3, y = 9, label = expression("Consumer \nmortality rate," ~ italic(m)[italic(i)]), size = 6)

rgr_plot <-
  ggplot() + 
  # geom_density(aes(x = intercept, y = after_stat(density * n), group = cellularity, fill = cellularity), alpha = 0.5, data = filter(data, parameter == "resource_growth_rate")) +
  # geom_density(aes(x = intercept, y = after_stat(density * n)), fill = "grey", alpha = 0.5, data = filter(data_all, parameter == "resource_growth_rate")) +
  geom_density(aes(x = intercept, group = cellularity, fill = cellularity), alpha = 0.5, data = filter(data, parameter == "resource_growth_rate")) +
  geom_density(aes(x = intercept), fill = "grey", alpha = 0.5, data = filter(data_all, parameter == "resource_growth_rate")) +
  geom_point(aes(x = activation_energy, y = 0), data = filter(growth_rates, cellularity == "multi"), color = "#F8766D", size = 5, alpha = 0.5) +
    geom_point(aes(x = activation_energy, y = 0), data = filter(growth_rates, cellularity == "uni"), color = "#00BFC4", size = 5, alpha = 0.5) +
  geom_point(aes(x = activation_energy, y = 0), data = growth_rates, color = "#252525", size = 5, shape = 1) +
  geom_vline(aes(xintercept = mean(intercept)), color = "red3", data = filter(data, parameter == "resource_growth_rate" & cellularity == "multi")) +
  geom_vline(aes(xintercept = mean(intercept)), color = "darkblue", data = filter(data, parameter == "resource_growth_rate" & cellularity == "uni")) +
  geom_vline(aes(xintercept = mean(intercept)), color = "black", data = filter(data_all, parameter == "resource_growth_rate")) +
  # coord_cartesian(xlim = c(-1.5, 2)) +
  coord_cartesian(ylim = c(0, 12), xlim = c(-1.5, 2)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12)) +
  labs(y = "Density", x = "Temperature Sensitivity (eV)") +
  theme_cowplot(font_size = 20) + 
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 9, label = expression("Resource \ngrowth rate,"~italic(r)[italic(k)]), size = 6)

conv_eff_plot <- 
  ggplot() + 
  # geom_density(aes(x = intercept, y = after_stat(density * n)), fill = "grey", alpha = 0.5, data = filter(data, parameter == "conversion_efficiency")) + 
  geom_density(aes(x = intercept), fill = "grey", alpha = 0.5, data = filter(data, parameter == "conversion_efficiency")) + 
  geom_point(aes(x = activation_energy, y = 0), data = filter(conv_rates, cellularity == "multi"), color = "#F8766D", size = 5, alpha = 0.5) +
  geom_point(aes(x = activation_energy, y = 0), data = filter(conv_rates, cellularity == "uni"), color = "#00BFC4", size = 5, alpha = 0.5) +
  geom_point(aes(x = activation_energy, y = 0), data = conv_rates, color = "#252525", size = 5, shape = 1) +
  geom_vline(aes(xintercept = mean(intercept)), data = filter(data, parameter == "conversion_efficiency")) +
  # coord_cartesian(xlim = c(-1.5, 2)) +
  coord_cartesian(ylim = c(0, 12), xlim = c(-1.5, 2)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12)) +
  labs(y = "Density", x = "Temperature Sensitivity (eV)") +
  theme_cowplot(font_size = 20) +
  annotate("text", x = 1.2, y = 9, label = expression("Conversion \nefficiency,"), size = 6) +
  # Second line: italic c with subscript
  annotate(
    "text",
    x = 1.2, y = 8.5,  # adjust vertical spacing
    parse = TRUE,
    label = "italic(v)[i*k]",
    size = 6,
    hjust = 0
  )

#problem with unicellular estimates ruining the plot
carrying_capacity_plot <-
  ggplot() + 
  # geom_density(aes(x = intercept, y = after_stat(density * n)), fill = "grey", alpha = 0.5, data = filter(data, parameter == "carrying_capacity")) + 
  geom_density(aes(x = intercept), fill = "grey", alpha = 0.5, data = filter(data, parameter == "carrying_capacity")) + 
  geom_point(aes(x = activation_energy, y = 0), data = filter(carrying_capacity, cellularity == "uni"), color = "#00BFC4", size = 5, alpha = 0.5) +
  geom_point(aes(x = activation_energy, y = 0), data = filter(carrying_capacity, cellularity == "multi"), color = "#F8766D", size = 5, alpha = 0.5) +
  geom_point(aes(x = activation_energy, y = 0), data = carrying_capacity, color = "#252525", size = 5, shape = 1) +
  geom_vline(aes(xintercept = mean(intercept)), color = "black", data = filter(data, parameter == "carrying_capacity")) +
  # coord_cartesian(xlim = c(-1.5, 2)) +
  coord_cartesian(ylim = c(0, 12), xlim = c(-1.5, 2)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12)) +
  labs(y = "Density", x = "Temperature Sensitivity (eV)") +
  theme_cowplot(font_size = 20) +
    # First line: plain text
    annotate(
      "text",
      x = -1, y = 9,
      label = "Resource carrying",
      size = 6,
      hjust = 0
    ) +
    # Second line: "carrying capacity, K_k" with italic K and subscript
    annotate(
      "text",
      x = -1, y = 8,  # adjust spacing
      parse = TRUE,
      label = '"capacity, "*italic(K)[k]',
      size = 6,
      hjust = 0
    )

consumption_rate_plot <-
  ggplot() + 
  # geom_density(aes(x = intercept, y = after_stat(density * n)), fill = "grey", alpha = 0.5, data = filter(data, parameter == "consumption rate")) + 
  geom_density(aes(x = intercept), fill = "grey", alpha = 0.5, data = filter(data, parameter == "consumption rate")) + 
  geom_point(aes(x = activation_energy, y = 0), data = filter(consumption_rate, cellularity == "multi"), color = "#F8766D", size = 5, alpha = 0.5) +
  geom_point(aes(x = activation_energy, y = 0), data = filter(consumption_rate, cellularity == "uni"), color = "#00BFC4", size = 5) +
  geom_point(aes(x = activation_energy, y = 0), data = consumption_rate, color = "#252525", size = 5, shape = 1) +
  geom_vline(aes(xintercept = mean(intercept)), color = "black", data = filter(data, parameter == "consumption rate")) +
  # coord_cartesian(xlim = c(-1.5, 2)) +
  coord_cartesian(ylim = c(0, 12), xlim = c(-1.5, 2)) +
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12)) +
  labs(y = "Density", x = "Temperature Sensitivity (eV)") +
  theme_cowplot(font_size = 20) +
    annotate(
      "text",
      x = -0.2, y = 9,
      parse = TRUE,
      label = '"Consumption "',
      size = 6
    ) +
    annotate(
      "text",
      x = -0.2, y = 8,   # adjust y to move it below
      parse = TRUE,
      label = '"rate, "*italic(c)[i*k]',
      size = 6
    )

# plot inter-process TAs
mac_cell2 %>% 
  group_by(simple_parameter, cellularity) %>% 
  summarise(n = n()) %>% view()

single_obs <- mac_cell2 %>% 
  filter(simple_parameter == "resource carrying capacity" & cellularity == "multi" | simple_parameter == "consumption rate" & cellularity == "uni") %>% 
  dplyr::select(simple_parameter, activation_energy)
 
#plot each estimate and CI bars - mortality rate and growth rate?
interTAs <-
  param_sum %>% 
  clean_names() %>% 
  mutate(parameter1 = str_replace_all(parameter, "_", " "),
         parameter1 = str_to_title(parameter1, local = "en"),
         parameter1 = str_replace(parameter1, "Mortality Rate", "Consumer Mortality Rate"),
         parameter1 = str_replace(parameter1, "Carrying", "Resource Carrying"),
         parameter1 = fct_reorder(parameter1, mean),
         Cellularity = ifelse(cellularity == "multi", "multicellular",
                               ifelse(cellularity == "uni", "unicellular", "potato"))) %>% 
  ggplot() + 
  geom_point(aes(x = parameter1, y = mean, group = Cellularity, colour = Cellularity, fill = Cellularity), size = 3) + 
  geom_errorbar(aes(x = parameter1, ymin = ci_low, ymax = ci_up, colour = Cellularity), width = 0.2) + 
  theme_cowplot(font_size = 20) +
  scale_x_discrete(
      labels = c(expression(italic(K)[italic(k)]),
                 expression(italic(c)[italic(ik)]),
                 expression(italic(v)[italic(ik)]),
                 expression(italic(m)[italic(i)]),
                 expression(italic(r)[italic(k)])
                 )) +
    labs(x = "Parameter", y = "Mean Temperature Senstivitity \n(eV)")

#### save multipanel - Fig S7 #######
ea_plots <-
  consumption_rate_plot + rgr_plot + carrying_capacity_plot + conv_eff_plot + mort_ea_plot + interTAs +
  plot_annotation(tag_levels = "A") #there are no unicellular mortality rate estimates, so I'm not sure why there are two dots there

# ggsave(filename = "figures/unimulti-ea-plots1.pdf", ea_plots, width = 16, height = 12)