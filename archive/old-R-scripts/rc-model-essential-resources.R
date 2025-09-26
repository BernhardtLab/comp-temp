### Joey April 16 2019
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# The essential resources model from soup to nuts -------------------------

## this is the case where Each species is limited by a different resource

### mortality rates
m1 <- 0.2; m2 <- 0.2

## max population growth rates
r1 <- 1.5; r2 <- 1.5

## half saturation constants
k12 <- 0.1; k21 <- 0.1
k11 <- 0.1; k22 <- 0.1

## consumption vectors

c11 <- 0.2
c12 <- 0.4
c21 <- 0.2
c22 <- 0.4

## R star
R12 <- (m1 * k12)/ (r1 - m1)
R11 <- (m1 * k11)/ (r1 - m1)
R21 <- (m2 * k21)/ (r2 - m2)
R22 <- (m2 * k22)/ (r2 - m2)
D <- 0.25

# N_1_star <- (D / c_12) * (S_2 - R_12) - (c_22 / c_12) * N_2_star


S1 <- 10 #0.248
S2 <- 10  #0.245


a11 <- c12 / (D * (S2 - R12))
a12 <- c22 / (D * S2 - R12)
a21 <- c11 / (D * S1 - R21)
a22 <- c21 / (D * S1 - R21)


rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
stabil_potential <- 1 - rho #stabilizing potential
fit_ratio <- sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D) #fitness ratio -- ask Patrick why he scaled his by the r's
coexist <- rho < fit_ratio &  fit_ratio < 1/rho


### We use the Arrhenius function to model the temperature dependence

joeys_function <- function(temperature, activation_energy = 2) {
	metabolism <- temperature*activation_energy
	return(metabolism)
}

temperatures <- seq(from = 0, to = 35)


metabolic_rates <- sapply(temperatures, joeys_function)
metabolic_rates

arrhenius_function <- function(Temp, E = 0.65, b1 = 10, ref_temp = 1) {
  k<-8.62e-05 #Boltzmann's constant
  E <- E # 0.6 # activation energy (eV)
  T<-Temp+273.15 #range of temp in K
  Tc <- ref_temp+273.15 #reference temperature
  
  metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
  return(metabolism)
}

growth_rate <- function(birth_rate, death_rate) {
	growth <- birth_rate - death_rate
}

metabolic_rates <- data.frame(temperature = temperatures, birth_rate = sapply(temperatures, arrhenius_function), death_rate = sapply(temperatures, arrhenius_function))
View(metabolic_rates)

metabolic_rates %>% 
	ggplot(aes(x = temperature, y = metabolic_rate)) + geom_point()

### to visualize what the arrhenius function looks like
b <- arrhenius_function(Temp = 0:40, E = 0.6, b1 = 0.5, ref_temp = 1)
d <- 0.25 + arrhenius_function(Temp = 0:40, E = 0.75, b1 = 0.25, ref_temp = 1)
plot(b-d, ylim = c(0,max(b-d)), type = "l")

b2 <- arrhenius_function(Temp = 0:40, E = 0.6, b1 = 0.5, ref_temp = -2)
d2 <- 0.25 + arrhenius_function(Temp = 0:40, E = 0.75, b1 = 0.25, ref_temp = -2)
lines(b2-d2, col = 2)

library(tidyverse)
library(cowplot)

#' Setup Tilman's consumer resource model parameters 

### Here we make the species differt in the intercepts of the their thermal performance curves
temp_dependences <- function(T = 25, r_Ea = 0.6, c_Ea = 0.6, k_Ea = 0, m_Ea = 0.2, ref_temp2 = 0){
	c11 = arrhenius_function(Temp = T, E = c_Ea, b1 = 0.2); c12 = arrhenius_function(Temp = T, E = c_Ea, b1 = 0.4); c21 = arrhenius_function(Temp = T, E = c_Ea, b1 = 0.5); c22 = arrhenius_function(Temp = T, E = c_Ea, b1 = 0.2) ### cij = per capita consumption of comsumer i on resource j
  k12 = arrhenius_function(Temp = T, E = k_Ea, b1 = 0.2); k21 = arrhenius_function(Temp = T, E = k_Ea, b1 = 0.1) #half saturation constant for N resource consumption
  k11 = arrhenius_function(Temp = T, E = k_Ea, b1 = 0.1); k22 = arrhenius_function(Temp = T, E = k_Ea, b1 = 0.1)
  r1 = arrhenius_function(Temp = T, E = r_Ea, b1 = 0.05); r2 = arrhenius_function(Temp = T, E = r_Ea, b1 = 0.05, ref_temp = ref_temp2) #population growth rates
  m1 = arrhenius_function(Temp = T, E = m_Ea, b1 = 0.01); m2 = arrhenius_function(Temp = T, E = m_Ea, b1 = 0.01) # mortality rates
  
  
  #' Supply concentrations of resources
  S1 <- 0.3 #0.248
  S2 <- 0.3  #0.245
  
  
  R12 <- (m1 * k12)/ (r1 - m1)
  R11 <- (m1 * k11)/ (r1 - m1)
  R21 <- (m2 * k21)/ (r2 - m2)
  R22 <- (m2 * k22)/ (r2 - m2)
  
  D <- 0.25
  
  a11 <- c12 / (D * (S2 - R12))
  a12 <- c22 / (D * (S2 - R12))
  a21 <- c11 / (D * (S1 - R21))
  a22 <- c21 / (D * (S1 - R21))
  
  
  rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
  stabil_potential <- 1 - rho #stabilizing potential
  # fit_ratio <- sqrt((a11*a12)/(a22*a21)) * (r2-D)/(r1-D) #fitness ratio -- ask Patrick why he scaled his by the r's
  fit_ratio <- sqrt((a11*a12)/(a22*a21))
  coexist <- rho < fit_ratio &  fit_ratio < 1/rho
  
  
  data.frame(a11 = a11, a12 = a12, a22 = a22, a21 = a21, R12 = R12, R11 = R11, R22 = R22, R21= R21, K1 = (r1)/a11, K2 = (r2)/a22, T = T, m1 = m1, m2 = m2, r1 = r1, r2 = r2, c11 = c11,  c12 = c12,  c21 = c21, c22 = c22, stabil_potential = stabil_potential, fit_ratio = fit_ratio, rho = rho, coexist = coexist)
}

### explore different temperature dependences
results <- data.frame()
for(i in 1:6){
  for(j in 1:6){
  r_Ea <- 0.3
  c_Ea <- 0.1*(i-1)
  m_Ea <- 0.1*(j-1)
  hold <- temp_dependences(T = seq(0, 35, by = 0.1), r_Ea = r_Ea , c_Ea = c_Ea, m_Ea = m_Ea, ref_temp2 = 1)
  hold$r_Ea <- r_Ea
  hold$c_Ea <- c_Ea
  hold$m_Ea <- m_Ea
  results <- bind_rows(results, hold)
  }}


## save first round of results


ggplot(results, aes(x = T, y = r1))+
	geom_line()+
	scale_linetype("w_Ea")+
	scale_color_viridis_c(end = 0.8, "log ratio\nc/r Ea")+
	coord_cartesian(ylim = c(0, max(results$r1)))+
	theme(legend.justification=c(0,1), legend.position=c(0,1))

ggplot(results, aes(x = T, y = m1, color = m_Ea, group = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	theme(legend.justification=c(0,1), legend.position=c(0,1))


results %>% 
	filter(r1 > 0) %>%
	ggplot(aes(x = T, y = K1, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "K_1") + xlab("Temperature (°C)")

results %>% 
	filter(r1 > 0) %>%
	ggplot(aes(x = T, y = K2, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "K_2") + xlab("Temperature (°C)")
ggsave("figures/temperature-dependent-K-c-m.png", width = 12, height = 12)

results %>% 
	filter(r1 > 0) %>%
	ggplot(aes(x = T, y = K2, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "K_2") + xlab("Temperature (°C)")

results %>% 
	filter(r1 > 0) %>%
	ggplot(aes(x = T, y = a11, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "a11") + xlab("Temperature (°C)")
ggsave("figures/essential-resources-a11-temperature-consx-my.png", width = 7, height = 5)

results %>% 
	filter(r1 > 0) %>%
	ggplot(aes(x = T, y = R21, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "R*_21") + xlab("Temperature (°C)")
ggsave("figures/essential-resources-r-star-temperature-cx-my.png", width = 7, height = 5)

results %>% 
	filter(r1 > 0) %>%
	ggplot(aes(x = T, y = R12, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "R*_12") + xlab("Temperature (°C)")

results %>% 
	filter(r1 > 0) %>%
	ggplot(aes(x = T, y = K1, group = m_Ea, color = m_Ea)) +
	geom_line() +
	scale_color_viridis_c(end = 0.8) +
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "K") + xlab("Temperature (°C)")
ggsave("figures/essential-resources-K-temperature.png", width = 7, height = 5)

results %>% 
	filter(r1 > 0) %>%
	ggplot(aes(x = T, y = a12, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_linetype(guide = F)+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea) +
	ylab(label = "alpha_21") + xlab("Temperature (°C)")


### can we make the square zngi plots from these results?

## how do we define the regions? Well the lines that delineate the sections start at the ZNGIs 
### are run parallel to the consumption vectors. 


### now need to write a script to first figure out what region we are in, then calculate the alphas etc.

### maybe first figure out how much it is possible for the consumption vectors and R*s to change with temperature?

## species i resource j

results2 <- results %>% 
	rename(n_star = R12,
		   p_star = R21)


results %>% 
	filter(T %in% c(0, 10, 20, 30, 40)) %>% 
	ggplot(aes(x = R12, y = R11, color = T)) + geom_point() +
	geom_segment(aes(x = R12, y = R11, xend = R12 + c12/20, yend = R11 + c11/20, color = T),
				 size = 0.5, linetype = 1, alpha = 1) +
	geom_segment(aes(x = R12, y = R11, xend = R12, yend = 0.11, color = T),
				 size = 0.5, linetype = 1, alpha = 1) +
	geom_segment(aes(x = R12, y = R11, xend = 0.11, yend = R11, color = T),
				 size = 0.5, linetype = 1, alpha = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	scale_colour_viridis_c() +
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	coord_cartesian() +
	xlim(0, 0.15) + ylim(0, 0.15) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 

ggsave("figures/zingis-temperature.png", width = 12, height = 10)


