# tpcs for conceptual figure

library(tidyverse)
library(cowplot)
library(rTPC)

theme <- theme_set(theme_cowplot(font_size = 30))

#make some data
model_data <- seq(from = 0, to = 40, length.out = 200) %>%
    tibble(temp = .) %>%
    mutate(y1 = sharpeschoolhigh_1981(temp = temp, tref = 20, r_tref = 5, e = 0.25, eh = 5, th = 35),
           y2 = sharpeschoolhigh_1981(temp = temp, tref = 20, r_tref = 5, e = 0.5, eh = 5, th = 35),
           y3 = sharpeschoolhigh_1981(temp = temp, tref = 20, r_tref = 5, e = 0.35, eh = 5, th = 35))

#### plotting thermal asymmetries with TPCs #######
#plot a TPC
model_data %>% 
  ggplot() + 
  # geom_line(aes(x = temp, y = y1), colour = "dodgerblue2", linewidth = 1.5) +
  # geom_line(aes(x = temp, y = y2), colour = "orchid", linewidth = 1.5) +
  geom_line(aes(x = temp, y = log(y2)), linewidth = 1.5) +
  coord_cartesian(xlim = c(0,40), ylim = c(0, 3)) + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) + 
  labs(x = "Temperature (°C)",
         y = "ln(Rate)") 
  
# ggsave(plot = last_plot(), file = "figures/kd-figs/tpc_generic.png", bg = "transparent")
# ggsave(plot = last_plot(), file = "figures/kd-figs/inter-process.png", bg = "transparent")

#plot inter-process TA 
model_data %>% 
  ggplot() + 
  annotate("rect",
           xmin = 0, xmax = 18, ymin = -1, ymax = 3, 
           fill = "lightgrey", alpha = 0.6) +
  geom_line(aes(x = temp, y = log(y1)), colour = "dodgerblue4", linewidth = 1.5) + 
  geom_line(aes(x = temp, y = log(y2)), colour = "goldenrod", linewidth = 1.5) + 
  coord_cartesian(xlim = c(2,40), ylim = c(0, 3)) +
  theme_cowplot(font_size = 20) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) + 
  labs(x = "Temperature (°C)",
       y = "ln(Rate)") 

# ggsave(plot = last_plot(), file = "figures/kd-figs/inter-process.png", bg = "transparent", height = 5, width = 5, units = "in")

#plot intra-process TA
model_data %>% 
  ggplot() + 
  annotate("rect",
           xmin = 0, xmax = 18, ymin = -1, ymax = 3, 
           fill = "lightgrey", alpha = 0.6) +
  geom_line(aes(x = temp, y = log(y1)), colour = "dodgerblue4", linewidth = 1.5) + 
  geom_line(aes(x = temp, y = log(y3)), colour = "dodgerblue2", linewidth = 1.5) + 
  coord_cartesian(xlim = c(2,40), ylim = c(0, 3)) + 
  theme_cowplot(font_size = 20) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) + 
  labs(x = "Temperature (°C)",
       y = "ln(Rate)") 

# ggsave(plot = last_plot(), file = "figures/kd-figs/intra-process.png", bg = "transparent", height = 5, width = 5, units = "in")

### plotting thermal asymmetries in distributions ####
#data
# simulate some data
data <- data.frame(x = rnorm(100000, mean = 0.01, sd = 0.22),
                     x1 = rnorm(100000, mean = 0.5, sd = 0.12))

ggplot(data) + 
  geom_density(aes(x = x1), bw = 0.1, alpha = 0.5, colour = "goldenrod" , linewidth = 1.5) +
  geom_density(aes(x = x), bw = 0.1, alpha = 0.5, colour = "dodgerblue4", linewidth = 1.5) + 
  theme_cowplot() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  labs(x = "Temperature dependence", y = "Frequency")
# ggsave(plot = last_plot(), file = "figures/kd-figs/inter_gaussian.png", bg = "transparent", height = 5, width = 5, units = "in")

ggplot(data) + 
  geom_density(aes(x = x1), bw = 0.1, fill = NA, colour = alpha("goldenrod", 0) , linewidth = 1.5) +
  geom_density(aes(x = x), bw = 0.1, alpha = 0.5, colour = "dodgerblue4", linewidth = 1.5) +
  theme_cowplot() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  labs(x = "Temperature dependence", y = "Frequency")

# ggsave(plot = last_plot(), file = "figures/kd-figs/intra_gaussian.png", bg = "transparent", height = 5, width = 5, units = "in")

### resource use preference conceptual #####
# simulate some data
pref_data <- data.frame(x = rnorm(100000, mean = 0.5, sd = 0.7),
                   x1 = rnorm(100000, mean = 1.5, sd = 0.7))

ggplot(pref_data) + 
  geom_density(aes(x = x1), bw = 0.1, alpha = 0.5, colour = "black" , linewidth = 2.5) +
  geom_density(aes(x = x), bw = 0.1, alpha = 0.5, colour = "darkgrey", linewidth = 2.5) + 
  theme_cowplot() + 
  coord_cartesian(ylim = c(0, 1), xlim = c(-0.75, 2.75)) +
  labs(x = "Resource use", y = "Consumption rate") + 
  theme_cowplot(font_size = 30) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) 

### uneven preference conceptual ####
pref_data <- pref_data %>% 
  mutate(x2 = rnorm(100000, mean = 0.7, sd = 0.7), #weak preference species - N2
         x3 = rnorm(100000, mean = 1.9, sd = 0.3)) #strong preference species - N1

ggplot(pref_data) + 
  geom_density(aes(x = x2), bw = 0.1, alpha = 0.5, colour = "darkgrey" , linewidth = 2.5) +
  geom_density(aes(x = x3), bw = 0.1, alpha = 0.5, colour = "black", linewidth = 2.5) + 
  theme_cowplot() + 
  coord_cartesian(ylim = c(0, 1.5), xlim = c(-0.75, 2.75)) +
  labs(x = "Resource use", y = "Consumption rate") + 
  theme_cowplot(font_size = 30) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# resource preference axes
prefs1 <- data.frame(resource = c("A", "B")) %>% 
  mutate(N1 = ifelse(resource == "A", 0.5, 1),
         N2 = ifelse(resource == "B", 0.5, 1),
         N3 = ifelse(resource == "A", 0.1, 0.9),
         N4 = ifelse(resource == "A", 0.6, 0.4)) %>% 
  pivot_longer(cols = N1:N4, names_to = "species", values_to = "consumption_rate") %>% 
  mutate(scenario = ifelse(species %in% c("N1", "N2"), "even", "uneven")) %>% 
  mutate(species = ifelse(scenario == "uneven" & species == "N3", "N1",
                          ifelse(scenario == "uneven" & species == "N4", "N2", species)))


prefs1 %>% 
  filter(scenario == "even") %>% 
  ggplot() + 
  geom_point(aes(x = resource, y = consumption_rate, colour = species), size = 6) + 
  geom_line(aes(x = resource, y = consumption_rate, colour = species, group = species), linewidth = 1.5) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_cowplot(font_size = 20) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(y = "Consumption rate", x = "Resource", colour = "Species")
  
ggsave(plot = last_plot(), file = "figures/conceptual/errce.png", width = 4.5, height = 3)

prefs1 %>% 
  filter(scenario == "uneven") %>% 
  ggplot() + 
  geom_point(aes(x = resource, y = consumption_rate, colour = species), size = 6) + 
  geom_line(aes(x = resource, y = consumption_rate, colour = species, group = species), linewidth = 1.5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_cowplot(font_size = 20) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(y = "Consumption rate", x = "Resource", colour = "Species")

ggsave(plot = last_plot(), file = "figures/conceptual/urrce.png", width = 4.5, height = 3)
