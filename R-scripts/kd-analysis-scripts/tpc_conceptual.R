# tpcs for conceptual figure

library(tidyverse)
library(cowplot)
library(rTPC)

theme <- theme_set(theme_cowplot(font_size = 30))

#make some data
model_data <- seq(from = 0, to = 40, length.out = 200) %>%
    tibble(temp = .) %>%
    mutate(y1 = sharpeschoolhigh_1981(temp = temp, tref = 20, r_tref = 5, e = 0.3, eh = 5, th = 35),
           y2 = sharpeschoolhigh_1981(temp = temp, tref = 20, r_tref = 5, e = 0.6, eh = 5, th = 35))

#plot it
model_data %>% 
  ggplot() + 
  geom_line(aes(x = temp, y = y1), colour = "dodgerblue2", linewidth = 1.5) + 
  geom_line(aes(x = temp, y = y2), colour = "orchid", linewidth = 1.5) + 
  coord_cartesian(xlim = c(0,40), ylim = c(0, 12)) + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) + 
  labs(x = "Temperature (°C)",
         y = "Rate") 
  

ggsave(plot = last_plot(), file = "figures/kd-figs/inter-process.png", bg = "transparent")

#plot it
model_data %>% 
  ggplot() + 
  geom_line(aes(x = temp, y = y1), colour = "goldenrod2", linewidth = 1.5) + 
  geom_line(aes(x = temp, y = y2), colour = "tomato3", linewidth = 1.5) + 
  coord_cartesian(xlim = c(0,40), ylim = c(0, 12)) + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) + 
  labs(x = "Temperature (°C)",
       y = "Rate")

ggsave(plot = last_plot(), file = "figures/kd-figs/intra-process.png", bg = "transparent")
