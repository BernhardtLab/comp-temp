# tpcs for conceptual figure

library(tidyverse)
library(cowplot)
library(rTPC)

theme <- theme_set(theme_cowplot(font_size = 30))

#make some data
model_data <- seq(from = 0, to = 40, length.out = 200) %>%
    tibble(temp = .) %>%
    mutate(y1 = sharpeschoolhigh_1981(temp = temp, tref = 20, r_tref = 5, e = 0.15, eh = 5, th = 35),
           y2 = sharpeschoolhigh_1981(temp = temp, tref = 20, r_tref = 5, e = 0.5, eh = 5, th = 35),
           y3 = sharpeschoolhigh_1981(temp = temp, tref = 20, r_tref = 5, e = 0.25, eh = 5, th = 35))

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
  geom_line(aes(x = temp, y = log(y1)), colour = "dodgerblue2", linewidth = 1.5) + 
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
  geom_line(aes(x = temp, y = log(y1)), colour = "dodgerblue4", linewidth = 1.5) + 
  geom_line(aes(x = temp, y = log(y3)), colour = "dodgerblue2", linewidth = 1.5) + 
  coord_cartesian(xlim = c(2,40), ylim = c(0, 3)) + 
  theme_cowplot(font_size = 20) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) + 
  labs(x = "Temperature (°C)",
       y = "ln(Rate)") 

ggsave(plot = last_plot(), file = "figures/kd-figs/intra-process.png", bg = "transparent", height = 5, width = 5, units = "in")




### plotting thermal asymmetries in distibutions ####
#data
# simulate some data
mydata <- data.frame(x = rnorm(100000, mean = 0.01, sd = 0.22),
                     x1 = rnorm(100000, mean = 0.5, sd = 0.12))

ggplot(mydata) + 
  geom_density(aes(x = x1), bw = 0.1, alpha = 0.5, colour = "goldenrod" , linewidth = 1.5) +
  geom_density(aes(x = x), bw = 0.1, alpha = 0.5, colour = "dodgerblue2", linewidth = 1.5) + 
  theme_cowplot() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  labs(x = "Temperature dependence", y = "Frequency")
ggsave(plot = last_plot(), file = "figures/kd-figs/dummy.png", bg = "transparent")
# ggsave(plot = last_plot(), file = "figures/kd-figs/inter_gaussian.png", bg = "transparent")

ggplot(mydata) + 
  geom_density(aes(x = x1), bw = 0.1, fill = NA, colour = alpha("goldenrod", 0.25) , linewidth = 1.5) +
  geom_density(aes(x = x), bw = 0.1, alpha = 0.5, colour = "dodgerblue2", linewidth = 1.5) +
  theme_cowplot() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  labs(x = "Temperature dependence", y = "Frequency")

# ggsave(plot = last_plot(), file = "figures/kd-figs/intra_gaussian.png", bg = "transparent")
