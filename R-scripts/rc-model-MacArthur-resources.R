### PJK Feb 16 2020
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())



### We use the Arrhenius function to model the temperature dependence
arrhenius_function <- function(Temp, E, b1, ref_temp = 1) {
  k <- 8.62e-05 #Boltzmann's constant
  E <- E # 0.6 # activation energy (eV)
  T <- Temp+273.15 #range of temp in K
  Tc <- ref_temp+273.15 #reference temperature
  
  metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
  return(metabolism)
}

# ### to visualize what the arrhenius function looks like
# b <- arrhenius_function(Temp = 0:40, E = 0.6, b1 = 0.5, ref_temp = 1)
# d <- 0.25 + arrhenius_function(Temp = 0:40, E = 0.75, b1 = 0.25, ref_temp = 1)
# plot(b-d, ylim = c(0,max(b-d)), type = "l")
# 
# b2 <- arrhenius_function(Temp = 0:40, E = 0.6, b1 = 0.5, ref_temp = -2)
# d2 <- 0.25 + arrhenius_function(Temp = 0:40, E = 0.75, b1 = 0.25, ref_temp = -2)
# lines(b2-d2, col = 2)



### Here we make the species differt in the intercepts of the their thermal performance curves
temp_dependences_MacArthur <- function(T = 25, r_Ea = 0.6, c_Ea = 0.6, K_Ea = 0, v_Ea = 0.0, m_Ea = 0.2, ref_temp2 = 0){
  
  # resource growth rates
  r1 = arrhenius_function(Temp = T, E = r_Ea, b1 = 0.05); 
  r2 = arrhenius_function(Temp = T, E = r_Ea, b1 = 0.05, ref_temp = ref_temp2) 
  
  # resource carrying capacity
  K1 = arrhenius_function(Temp = T, E = K_Ea, b1 = 2) 
  K2 = arrhenius_function(Temp = T, E = K_Ea, b1 = 2) 
  
  # cij = per capita consumption of comsumer i on resource j
  c11 = arrhenius_function(Temp = T, E = c_Ea, b1 = 0.2)
  c12 = arrhenius_function(Temp = T, E = c_Ea, b1 = 0.4)
  c21 = arrhenius_function(Temp = T, E = c_Ea, b1 = 0.4)
  c22 = arrhenius_function(Temp = T, E = c_Ea, b1 = 0.2) 
  
  # vij = conversion factor that converts resource j into biomass of consumer i
  v11 = arrhenius_function(Temp = T, E = v_Ea, b1 = 0.2) 
  v12 = arrhenius_function(Temp = T, E = v_Ea, b1 = 0.4) 
  v21 = arrhenius_function(Temp = T, E = v_Ea, b1 = 0.4) 
  v22 = arrhenius_function(Temp = T, E = v_Ea, b1 = 0.2)
  
  # mortality rates
  m1 = arrhenius_function(Temp = T, E = m_Ea, b1 = 0.01)
  m2 = arrhenius_function(Temp = T, E = m_Ea, b1 = 0.01)
  
  # Absolute competition coefficients
  beta11 = v11 * c11 * (K1/r1) * c11 + v12 * c12 * (K2/r2) * c12
  beta12 = v11 * c11 * (K1/r1) * c21 + v12 * c12 * (K2/r2) * c22
  beta21 = v21 * c21 * (K1/r1) * c11 + v22 * c22 * (K2/r2) * c12
  beta22 = v21 * c21 * (K1/r1) * c21 + v22 * c22 * (K2/r2) * c22
  g1 = v11 * c11 * K1 + v12 * c12 * K2 - m1
  g2 = v21 * c21 * K1 + v22 * c22 * K2 - m2
  
  # Relative competition coefficients 
  a11 = beta11 / g1
  a12 = beta12 / g1
  a22 = beta22 / g2
  a21 = beta21 / g2
  
  # MCT components
  rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
  stabil_potential <- 1 - rho #stabilizing potential
  fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio -- ask Patrick why he scaled his by the r's
  coexist <- rho < fit_ratio &  fit_ratio < 1/rho
  
  # report results
  data.frame(a11 = a11, a12 = a12, a22 = a22, a21 = a21, g1 = g1, g2 = g2, T = T, m1 = m1, m2 = m2, r1 = r1, r2 = r2, c11 = c11,  c12 = c12,  c21 = c21, c22 = c22, stabil_potential = stabil_potential, fit_ratio = fit_ratio, rho = rho, coexist = coexist)
}



### Explore different temperature dependences
results_MacArthur <- data.frame()
for(i in 1:6){
  for(j in 1:6){
    r_Ea <- 0.3
    c_Ea <- 0.1*(i-1)
    m_Ea <- 0.1*(j-1)
    hold <- temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), r_Ea = r_Ea , c_Ea = c_Ea, m_Ea = m_Ea, ref_temp2 = 1)
    hold$r_Ea <- r_Ea
    hold$c_Ea <- c_Ea
    hold$m_Ea <- m_Ea
    results_MacArthur <- bind_rows(results_MacArthur, hold)
  }
}



### Plot
# log(Fitness ratio) ~ Stabilization potential
# Fitness ratio ~ Stabilization potential
results_MacArthur %>%
  filter((g1 > 0) & (g2 > 0)) %>%
  ggplot(aes(x = stabil_potential, y = fit_ratio, col=T))+
  geom_point()+
  # geom_ribbon(data = data.frame(x = seq(min(results_MacArthur$stabil_potential)*0.8, max(results_MacArthur$stabil_potential)*1.2, 0.001)),
  #             aes(x = x,
  #                 ymin = 1-x,
  #                 ymax = 1/(1-x)),
  #             fill = "grey", color = "black", alpha = 0.3) +
  geom_hline(yintercept = 1, linetype=5) +
  scale_linetype(guide = F) +
  scale_color_viridis_c(end = 0.8)+
  facet_grid(rows=vars(c_Ea), col=vars(m_Ea)) +
  ylab(label = "Fitness ratio") + xlab("Stabilization potential")

# results_MacArthur %>% 
#   filter((g1 > 0) & (g2 > 0)) %>%
#   ggplot(aes(x = stabil_potential, y = log(fit_ratio), col=T))+
#   geom_point()+
#   geom_ribbon(data = data.frame(x = seq(min(results_MacArthur$stabil_potential)*0.8, max(results_MacArthur$stabil_potential)*1.2, 0.001)),
#               aes(x = x, 
#                   ymin = log(1-x), 
#                   ymax = log(1/(1-x))),
#               fill = "grey", color = "black", alpha = 0.3) +
#   geom_hline(yintercept = 0, linetype=5) +
#   scale_linetype(guide = F) +
#   scale_color_viridis_c(end = 0.8)+
#   facet_grid(c_Ea ~ m_Ea) +
#   ylab(label = "log(Fitness ratio)") + xlab("Stabilization potential")

# # Stabilization potential ~ Temperature
# results_MacArthur %>% 
#   filter((g1 > 0) & (g2 > 0)) %>%
#   ggplot(aes(x = T, y = stabil_potential, group = m_Ea, color = m_Ea))+
#   geom_line()+
#   scale_linetype(guide = F)+
#   scale_color_viridis_c(end = 0.8)+
#   facet_grid(c_Ea ~ m_Ea) +
#   ylab(label = "Stabilization potential") + xlab("Temperature (°C)")

# # Fitness ratio ~ Temperature
# results_MacArthur %>% 
#   filter((g1 > 0) & (g2 > 0)) %>%
#   ggplot(aes(x = T, y = fit_ratio, group = m_Ea, color = m_Ea))+
#   geom_line()+
#   scale_linetype(guide = F)+
#   scale_color_viridis_c(end = 0.8)+
#   facet_grid(c_Ea ~ m_Ea) +
#   ylab(label = "Fitness ratio") + xlab("Temperature (°C)")



### Code for Pretty plots
### Grab data from results_MacArthur
data = results_MacArthur[((results_MacArthur$c_Ea==0.0) & (results_MacArthur$m_Ea==0.5)), ]
components.all = data.frame(scenario = c("base", "warm"), 
                            ND = c(head(data)[1, ]$stabil_potential, tail(data)[6, ]$stabil_potential), 
                            FR = c(head(data)[1, ]$fit_ratio, tail(data)[6, ]$fit_ratio), 
                            Coexist = c(TRUE, FALSE))
components.base = components.all[components.all$scenario == "base", ]
components = components.all[components.all$scenario == "warm", ]



### Prepare plotting objects
line.avoid.radius = 0.015   
arrow.avoid.radius = 0.015 

### Calculate arrow coordinates (slightly offset from start/end coordinates)
scenario.arrows = components %>%
  
  # All scenarios start at the sterile point
  mutate(ND0 = components.base$ND, FR0 = components.base$FR) %>%
  
  # Calculate change in ND and FR, as well as length of line
  mutate(dND = ND - ND0, dFR = FR - FR0, length = sqrt(dND^2+dFR^2)) %>%
  
  # Offset arrow start/end by the required length
  mutate(ND1 = ND0 + (line.avoid.radius / length) * dND,
         FR1 = FR0 + (line.avoid.radius / length) * dFR,
         ND2 = ND - (arrow.avoid.radius / length) * dND,
         FR2 = FR - (arrow.avoid.radius / length) * dFR)

### Put start and end points in tidy format
scenario.points <- components %>%
  mutate(ND = components.base$ND, FR = components.base$FR, point.type = 'start') %>%
  bind_rows(mutate(components, point.type = 'end'))



### Plot
ggplot() +
  
  # Plot Chesson parameter space
  geom_ribbon(data = data.frame(x = seq(min(components.all$ND)*0.8, max(components.all$ND)*1.2, 0.001)),
              aes(x = x, 
                  ymin = log(1-x), 
                  ymax = log(1/(1-x))),
              fill = "grey", alpha = 0.3) +
  
  # Add starting and end points
  geom_point(data = rbind(components.base, components),
             aes(x = ND, 
                 y = log(FR), 
                 color = scenario),
             size = 2, stroke = 1) +
  
  # Add arrows
  geom_segment(data = scenario.arrows, 
               aes(x = ND1, 
                   y = log(FR1), 
                   xend = ND2, 
                   yend = log(FR2), 
                   color = scenario),
               size = 1, linejoin = 'mitre',
               arrow = arrow(type = 'closed', length = unit(7, 'pt')),
               show.legend = F) +
  
  # Color and text setting
  scale_color_manual('Temperature', 
                     values = c("base" = "black", 
                                "warm" = "#fe9929",
                                "warmer" = "#8c2d04"), 
                     breaks = c('base', 'warm', 'warmer'),
                     guide = F) +
  
  # Coordinate and axes settings
  coord_cartesian(xlim = c(0.15, 0.20), ylim = log(c(0.5, 1.6))) +
  ylab(expression(paste("fitness ratio, ", italic(f[2]), "/", italic(f[1])))) +
  xlab(expression(paste("niche difference, ", 1-rho))) +
  theme_bw()






