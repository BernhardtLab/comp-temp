### PJK Apr 30 2020
library(tidyverse)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())
library(plotrix)



### Arrhenius function to model the temperature dependence
arrhenius_function <- function(Temp, E, b1, ref_temp = 1) {
  k <- 8.62e-05 #Boltzmann's constant
  E <- E # 0.6 # activation energy (eV)
  T <- Temp+273.15 #range of temp in K
  Tc <- ref_temp+273.15 #reference temperature
  
  metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
  return(metabolism)
}


### JB modifying to have two consumer species, 1 and 2, and two resource species, N and P
### Function to calculate MCT components under different temperatures
temp_dependences_MacArthur <- function(T = 25, ref_temp2 = 1,
                                       r_EaN = 0.5, r_EaP = 0.5, 
                                       c_Ea1N = 0.8, c_Ea1P = 0.8, 
                                       c_Ea2N = 0.8, c_Ea2P = 0.8, 
                                       K_EaN = -0.3, K_EaP = -0.3, 
                                       v_EaN = 0.0, v_EaP = 0.0, 
                                       m_Ea1 = 0.6, m_Ea2 = 0.6,
									   C1N_b = 0.2, c2P_b = 0.2,
									   c1P_b = 0.4, c2N_b = 0.4
									   ){
  
  # resource growth rates
  rN = arrhenius_function(Temp = T, E = r_EaN, b1 = 0.05) 
  rP = arrhenius_function(Temp = T, E = r_EaP, b1 = 0.05) 
  
  # resource carrying capacity
  KN = arrhenius_function(Temp = T, E = K_EaN, b1 = 2) 
  KP = arrhenius_function(Temp = T, E = K_EaP, b1 = 2) 
  
  # cij = per capita consumption of comsumer i on resource j
  c1N = arrhenius_function(Temp = T, E = c_Ea1N, b1 = C1N_b)
  c1P = arrhenius_function(Temp = T, E = c_Ea1P, b1 = c1P_b) ## species 1 consumes more P than N
  c2N = arrhenius_function(Temp = T, E = c_Ea2N, b1 = c2N_b) ## species 2 consumes more N than P
  c2P = arrhenius_function(Temp = T, E = c_Ea2P, b1 = c2P_b) 
  
  # vij = conversion factor that converts resource j into biomass of consumer i
  v1N = arrhenius_function(Temp = T, E = v_EaN, b1 = 0.2)
  v2N = arrhenius_function(Temp = T, E = v_EaN, b1 = 0.4) ## species 2 converts N more efficiently 
  v1P = arrhenius_function(Temp = T, E = v_EaP, b1 = 0.4) ## species 1 converts P more efficiently 
  v2P = arrhenius_function(Temp = T, E = v_EaP, b1 = 0.2)
  
  # mortality rates
  m1 = arrhenius_function(Temp = T, E = m_Ea1, b1 = 0.01)
  m2 = arrhenius_function(Temp = T, E = m_Ea2, b1 = 0.01)
  
  # Absolute competition coefficients
  beta11 = v1N * c1N * (KN/rN) * c1N + v1P * c1P * (KP/rP) * c1P ### intra
  beta12 = v1N * c1N * (KN/rN) * c2N + v1P * c1P * (KP/rP) * c2P ### inter
  beta22 = v2N * c2N * (KN/rN) * c2N + v2P * c2P * (KP/rP) * c2P ### intra
  beta21 = v2N * c2N * (KN/rN) * c1N + v2P * c2P * (KP/rP) * c1P ### inter
  
  
  g1 = v1N * c1N * KN + v1P * c1P * KP - m1 ### growth rate of the consumer 1
  g2 = v2N * c2N * KN + v2P * c2P * KP - m2 ### growth rate of consumer 2
  
  # Relative competition coefficients 
  a11 = beta11 / g1
  a21 = beta21 / g2
  a22 = beta22 / g2
  a12 = beta12 / g1
  
  # MCT components
  rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
  stabil_potential <- 1 - rho #stabilizing potential
  fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio
  coexist <- rho < fit_ratio &  fit_ratio < 1/rho
  
  # report results
  data.frame(T = T, 
             a11 = a11, a12 = a12, a22 = a22, a21 = a21, g1 = g1, g2 = g2,
             stabil_potential = stabil_potential, fit_ratio = fit_ratio, rho = rho, coexist = coexist, 
             m1 = m1, m2 = m2, rN = rN, rP = rP, KN = KN, KP = KP,
             c1N = c1N,  c1P = c1P,  c2N = c2N, c2P = c2P)
}



### Function to plot: parameter values, consumption vector slopes, ND & FD
plot_MacArthur <- function(Data.parameter, Data.temperature){
  
  # panel a: parameter value ~ temperature | species
  panel.a = 
    ggplot(Data.parameter, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    xlab("Temperature") +
    ylab(label = "Parameter value") + 
    ggtitle("(a)") +
    theme(legend.position = "bottom")
  
  # panel b: consumption vector slope ~ temperature | species
  Temp = data.frame(T = Data.temperature$T, 
                    slope1 = Data.temperature$c1P / Data.temperature$c1N, 
                    slope2 = Data.temperature$c2P / Data.temperature$c2N)
  Temp = Temp %>% gather(value=value, key=parameter, -T)
  panel.b = 
    ggplot(Temp, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype=5) +
    xlab("Temperature") +
    ylab("Parameter value") + 
    ggtitle("(b)") +
    theme(legend.position = "bottom")
  
  # panel c: ND ~ temperature
  panel.c = 
    Data.temperature %>%
    filter((g1 > 0) & (g2 > 0)) %>%
    ggplot(aes(x = T, y = stabil_potential, col = T)) +
    geom_point() +
    scale_color_viridis_c("Temperature", end = 0.8) +
    xlab("Temperature") +
    ylab(expression(paste("Stabilization potential (1-", rho, ")"))) + 
    ggtitle("(c)") +
    theme(legend.position = "bottom")
  
  # panel d: FD ~ temperature
  panel.d = 
    Data.temperature %>%
    filter((g1 > 0) & (g2 > 0)) %>%
    ggplot(aes(x = T, y = fit_ratio, col = T)) +
    geom_point() +
    scale_color_viridis_c("Temperature", end = 0.8) +
    xlab("Temperature") +
    ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
    ggtitle("(d)") +
    theme(legend.position = "bottom")
  
  # panel e: FD ~ ND
  panel.e = 
    Data.temperature %>%
    filter((g1 > 0) & (g2 > 0)) %>%
    ggplot(aes(x = stabil_potential, y = fit_ratio, col = T)) +
    geom_point() +
    geom_ribbon(data = data.frame(x = seq(min(Data.temperature$stabil_potential)*0.99, max(Data.temperature$stabil_potential)*1.01, 0.001)),
                aes(x = x,
                    y = NULL, 
                    ymin = 1-x,
                    ymax = 1/(1-x)),
                fill = "grey", color = "black", alpha = 0.3) +
    geom_hline(yintercept = 1, linetype=5) +
    scale_x_continuous(expand=c(0, 0)) + 
    scale_color_viridis_c("Temperature", end = 0.8) +
    xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
    ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
    ggtitle("(e)") +
    theme(legend.position = "bottom")
  
  # Combine panels
  return(panel.a + panel.b + panel.c + panel.d + panel.e + plot_layout(nrow=1))
  
}



### Function to plot: alpha, competitive impact, sensitivity to competition
plot_MacArthur_alpha <- function(Data.temperature){
  
  # panel a: alpha ~ temperature | species
  Data.alpha = Data.temperature[, c("T", "a11", "a12", "a22", "a21")] %>% gather(value=value, key=parameter, -T)
  panel.a = 
    ggplot(Data.alpha, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    scale_color_manual(values=c("a11"="#08519c", "a21"="#6baed6", "a22"="#a50f15", "a12"="#fb6a4a")) + 
    xlab("Temperature") +
    ylab(expression(alpha[ij])) + 
    ggtitle("(a)") +
    theme(legend.position = "bottom")
  
  # panel b: relative impact ~ temperature | species
  Temp = data.frame(T = Data.temperature$T, 
                    rel.impact1 = Data.temperature$a21 / Data.temperature$a11, 
                    rel.impact2 = Data.temperature$a12 / Data.temperature$a22)
  Temp = Temp %>% gather(value=value, key=parameter, -T)
  panel.b = 
    ggplot(Temp, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    scale_color_manual(values=c("rel.impact1"="#08519c", "rel.impact2"="#a50f15")) + 
    xlab("Temperature") +
    ylab(expression(paste("Relative impact (", alpha[ji], "/", alpha[ii], ")"))) + 
    ggtitle("(b)") +
    theme(legend.position = "bottom")
  
  # panel c: sensitivity to competition ~ temperature | species
  Temp2 = data.frame(T = Data.temperature$T, 
                     sen1 = Data.temperature$a12 * Data.temperature$a11, 
                     sen2 = Data.temperature$a21 * Data.temperature$a22)
  Temp2 = Temp2 %>% gather(value=value, key=parameter, -T)
  panel.c = 
    ggplot(Temp2, 
           aes(x = T, y = value, col = parameter)) +
    geom_point() +
    scale_color_manual(values=c("sen1"="#08519c", "sen2"="#a50f15")) + 
    xlab("Temperature") +
    ylab(expression(paste("Sensitivity (", alpha[ii], "*", alpha[ij], ")"))) + 
    ggtitle("(c)") +
    theme(legend.position = "bottom")
  
  # Combine panels
  return(panel.a + panel.b + panel.c + plot_layout(nrow=1))
  
}



### SIMULATE !
### r varies with temperature & species have different activation energy for r (benefitting Sp1: r_Ea2 > r_Ea1)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.r = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
#                                                 r_Ea1 = 0.3, r_Ea2 = 0.5, 
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.3, 
#                                                 c_Ea11 = 0.8, c_Ea12 = 0.8, 
#                                                 c_Ea21 = 0.8, c_Ea22 = 0.8, 
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0, 
#                                                 m_Ea1 = 0.6, m_Ea2 = 0.6)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism 
### Patterns are similar with previous simulations
Data.temperature.r = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_EaN = 0.3, r_EaP = 0.5, 
                                                K_EaN = -0.0, K_EaP = -0.0, 
                                                c_Ea1N = 0.0, c_Ea1P = 0.0, 
                                                c_Ea2N = 0.0, c_Ea2P = 0.0, 
                                                v_EaN = 0.0, v_EaP = 0.0, 
                                                m_Ea1 = 0.0, m_Ea2 = 0.0)
Data.parameter.r = Data.temperature.r[, c("T", "rN", "rP")] %>% gather(value=value, key=parameter, -T)
Plot.r = plot_MacArthur(Data.parameter.r, Data.temperature.r)
Plot.alpha.r = plot_MacArthur_alpha(Data.temperature.r)
ggsave(file="figures/joey-exploring/MacArthur-scenario-r-unequal-vs.png", plot = Plot.r, device = "png", width=17.5, height=4)
ggsave(file="figures/joey-exploring/MacArthur-scenario-r-alpha.png", plot = Plot.alpha.r, device = "png", width=10.5, height=4)



### K varies with temperature & species have different activation energy for K (benefitting Sp1: K_Ea2 > K_Ea1)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.K = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
#                                                 r_Ea1 = 0.5, r_Ea2 = 0.5, 
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.1, 
#                                                 c_Ea11 = 0.8, c_Ea12 = 0.8, 
#                                                 c_Ea21 = 0.8, c_Ea22 = 0.8, 
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0, 
#                                                 m_Ea1 = 0.6, m_Ea2 = 0.6)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism
### Patterns are similar with previous simulations
Data.temperature.K = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_EaN = 0.0, r_EaP = 0.0, 
                                                K_EaN = -0.3, K_EaP = -0.1, 
                                                c_Ea1N = 0.0, c_Ea1P = 0.0, 
                                                c_Ea2N = 0.0, c_Ea2P = 0.0, 
                                                v_EaN = 0.0, v_EaP = 0.0, 
                                                m_Ea1 = 0.0, m_Ea2 = 0.0)
Data.parameter.K = Data.temperature.K[, c("T", "KN", "KP")] %>% gather(value=value, key=parameter, -T)
Plot.K = plot_MacArthur(Data.parameter.K, Data.temperature.K)
Plot.alpha.K = plot_MacArthur_alpha(Data.temperature.K)
ggsave(file="figures/joey-exploring/MacArthur-scenario-K.png", plot=Plot.K, device = "png", width=17.5, height=4)
ggsave(file="figures/joey-exploring/MacArthur-scenario-K-alpha.png", plot=Plot.alpha.K, device = "png", width=10.5, height=4)



### m varies with temperature & species have different activation energy for m (benefitting Sp1: m_Ea1 < m_Ea2)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.m = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
#                                                 r_Ea1 = 0.5, r_Ea2 = 0.5, 
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.3, 
#                                                 c_Ea11 = 0.8, c_Ea12 = 0.8, 
#                                                 c_Ea21 = 0.8, c_Ea22 = 0.8, 
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0, 
#                                                 m_Ea1 = 0.4, m_Ea2 = 0.6)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism
### Patterns are similar with previous simulations
Data.temperature.m = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_EaN = 0.0, r_EaP = 0.0, 
                                                K_EaN = -0.0, K_EaP = -0.0, 
                                                c_Ea1N = 0.0, c_Ea1P = 0.0, 
                                                c_Ea2N = 0.0, c_Ea2P = 0.0, 
                                                v_EaN = 0.0, v_EaP = 0.0, 
                                                m_Ea1 = 0.4, m_Ea2 = 0.6)
Data.parameter.m = Data.temperature.m[, c("T", "m1", "m2")] %>% gather(value=value, key=parameter, -T)
Plot.m = plot_MacArthur(Data.parameter.m, Data.temperature.m)
Plot.alpha.m = plot_MacArthur_alpha(Data.temperature.m)
ggsave(file="figures/joey-exploring/MacArthur-scenario-m.png", plot=Plot.m, device = "png", width=17.5, height=4)
ggsave(file="figures/joey-exploring/MacArthur-scenario-m-alpha.png", plot=Plot.alpha.m, device = "png", width=10.5, height=4)



### c varies with temperature & species have different activation energy for c (benefitting Sp1: c_Ea1i > c_Ea2i)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.c = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1),
#                                                 r_Ea1 = 0.5, r_Ea2 = 0.5,
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.3,
#                                                 c_Ea11 = 0.8, c_Ea12 = 0.8,
#                                                 c_Ea21 = 0.6, c_Ea22 = 0.6,
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0,
#                                                 m_Ea1 = 0.6, m_Ea2 = 0.6)
# Data.parameter.c = Data.temperature.c[, c("T", "c11", "c12", "c21", "c22")] %>% gather(value=value, key=parameter, -T)
# Plot.c = plot_MacArthur(Data.parameter.c, Data.temperature.c)
# Plot.alpha.c = plot_MacArthur_alpha(Data.temperature.c)
# ggsave(file="../figures/MacArthur-scenario-c-2.png", plot=Plot.c, device = "png", width=17.5, height=4)
# ggsave(file="../figures/MacArthur-scenario-c-alpha-2.png", plot=Plot.alpha.c, device = "png", width=10.5, height=4)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism
### Patterns are different with previous simulations, indicating interactive effect among parameters
rnorm(100,mean=0.2,sd=0.05)

Data.temperature.c = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_EaN = 0.0, r_EaP = 0.0, 
                                                K_EaN = -0.0, K_EaP = -0.0, 
                                                c_Ea1N = 0.8, c_Ea1P = 0.8, 
                                                c_Ea2N = 0.2, c_Ea2P = 0.6, 
                                                v_EaN = 0.0, v_EaP = 0.0, 
                                                m_Ea1 = 0.0, m_Ea2 = 0.0)
Data.parameter.c = Data.temperature.c[, c("T", "c1N", "c1P", "c2N", "c2P")] %>% gather(value=value, key=parameter, -T)
Plot.c = plot_MacArthur(Data.parameter.c, Data.temperature.c)
Plot.alpha.c = plot_MacArthur_alpha(Data.temperature.c)
ggsave(file="figures/joey-exploring/MacArthur-scenario-c.png", plot=Plot.c, device = "png", width=17.5, height=4)
ggsave(file="figures/joey-exploring/MacArthur-scenario-c-alpha.png", plot=Plot.alpha.c, device = "png", width=10.5, height=4)



### preference varies with temperature & species have different activation energy for c (c_Ea11 > c_Ea12 & c_Ea22 > c_Ea21)
# ### Other parameters vary with temperature (same AE between species for other parameters) 
# Data.temperature.p = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1),
#                                                 r_Ea1 = 0.5, r_Ea2 = 0.5,
#                                                 K_Ea1 = -0.3, K_Ea2 = -0.3,
#                                                 c_Ea11 = 0.95, c_Ea12 = 0.65,
#                                                 c_Ea21 = 0.6, c_Ea22 = 0.8,
#                                                 v_Ea1 = 0.0, v_Ea2 = 0.0,
#                                                 m_Ea1 = 0.6, m_Ea2 = 0.6)

### Other parameters do not vary with temperature (AE=0) to disentangle underlying mechanism
### Patterns are similar with previous simulations
Data.temperature.p = temp_dependences_MacArthur(T = seq(0, 35, by = 0.1), 
                                                r_EaN = 0.0, r_EaP = 0.0, 
                                                K_EaN = -0.0, K_EaP = -0.0, 
                                                c_Ea1N = 0.95, c_Ea1P = 0.65, 
                                                c_Ea2N = 0.6, c_Ea2P = 0.8, 
                                                v_EaN = 0.0, v_EaP = 0.0, 
                                                m_Ea1 = 0.0, m_Ea2 = 0.0)
Data.parameter.p = Data.temperature.p[, c("T", "c1N", "c1P", "c2N", "c2P")] %>% gather(value=value, key=parameter, -T)
Plot.p = plot_MacArthur(Data.parameter.p, Data.temperature.p)
Plot.alpha.p = plot_MacArthur_alpha(Data.temperature.p)
ggsave(file="figures/joey-exploring/MacArthur-scenario-p.png", plot=Plot.p, device = "png", width=14, height=4)
ggsave(file="figures/joey-exploring/MacArthur-scenario-p-alpha.png", plot=Plot.alpha.p, device = "png", width=10.5, height=4)





results_range <- data.frame()
for(f in seq(0, 1, by = 0.5)){
	for(b in seq(0, 1, by = 0.5)){			
	for(j in seq(0, 1, by = 0.5)){
				for(i in seq(0, 1, by = 0.5)){
					for(q in seq(0, 1, by = 0.5)){
						for(r in seq(0, 1, by = 0.5)){
							for(s in seq(-0.6, 0, by = 0.3)){
								for(l in seq(-0.6, 0, by = 0.3)){
									for(r1 in seq(0, 1, by = 0.5)){
										for(r2 in seq(0, 1, by = 0.5)){
											
hold = temp_dependences_MacArthur(T = seq(0, 34, by = 0.1), 
												r_EaN = r1, r_EaP = r2, 
												K_EaN = s, K_EaP = l, 
												c_Ea1N = q, c_Ea1P = r, 
												c_Ea2N = j, c_Ea2P = i, 
												v_EaN = 0.0, v_EaP = 0.0, 
												m_Ea1 = f, m_Ea2 = b) %>% 
	mutate(c_EA2N = j) %>% 
	mutate(c_EA2P = i) %>% 
	mutate(c_EA1N = q) %>% 
	mutate(c_EA1P = r) %>% 
	mutate(K_EAN = s) %>% 
	mutate(K_EAP = l) %>% 
	mutate(r_EAN = r1) %>% 
	mutate(r_EAP = r2)  %>% 
	mutate(m_EAN = f) %>% 
	mutate(m_EAP = b) 
	

	
results_range <- bind_rows(results_range, hold)

				}
					}
				}}}}}}}}
View(results_range)
write_csv(results_range, "data-processed/results-range.csv")
write_csv(results_range, "data-processed/results-range-0.1.csv")
library(plotrix)



# results simulations -----------------------------------------------------
#### this is where I'm coming back, Feb 2 2022

results_range <- read_csv("data-processed/results-range.csv")
results_range_by_0.1 <- read_csv("data-processed/results-range-0.1.csv")
results_range1 <- results_range
results_range2 <- results_range
results_range3 <- results_range

names(results_range)

results_sum <- results_range %>%
	# mutate(temp_deg = round(T, digits = 0)) %>% 
	# filter(g1 > 0) %>% 
	# filter(g2 > 0) %>% 
	# group_by(temp_deg) %>% 
	group_by(T) %>% 
	summarise_all(funs(median, mean, std.error)) 

rs2 <- results_sum %>% 
	mutate(krN_mean = KN_mean / rN_mean) %>% 
	mutate(krP_mean = KP_mean / rP_mean) %>% 
	mutate(v1N_mean =  0.2) %>% 
	mutate(v2N_mean = 0.4) %>% 
	mutate(v1P_mean = 0.4 ) %>% 
	mutate(v2P_mean =  0.2) %>% 
	mutate(beta_11_first_mean = v1N_mean * c1N_mean * (KN_mean/rN_mean) * c1N_mean) %>% 
	mutate(beta_11_second_mean = v1P_mean * c1P_mean * (KP_mean/rP_mean) * c1P_mean) %>% 
	mutate(beta_11_mean = beta_11_first_mean + beta_11_second_mean) %>% 
	mutate(alpha11_cal_mean = (beta_11_first_mean + beta_11_second_mean) / g1_mean) 

names(results_range)

# results_sum <- results_range %>% 
# 	filter(g1 > 0) %>% 
# 	filter(g2 > 0) %>% 
# 	group_by(T) %>% 
# 	summarise_each(funs(mean, std.error, median), stabil_potential, fit_ratio)

results_range %>% 
	ggplot(aes(x = T, y = m1)) + geom_point() +
	geom_point(aes(x = T, y = m2), color = "pink")

neg_alphas <- results_range %>% 
	filter(a21 <0)


results_sum %>% 
	ggplot(aes(x = T, y = stabil_potential_mean)) + geom_line() +
	geom_ribbon(aes(x = T, ymin = stabil_potential_mean - stabil_potential_std.error, ymax = stabil_potential_mean + stabil_potential_std.error), alpha = 0.2) +
	ylab("Stabilization potential") + xlab("Temperature (°C)")
ggsave("figures/stabil-potent-param-space.png", width = 6, height = 4)


results_sum %>% 
	ggplot(aes(x = T, y = fit_ratio_mean)) + geom_line() +
	geom_ribbon(aes(x = T, ymin = fit_ratio_mean - fit_ratio_std.error, ymax = fit_ratio_mean + fit_ratio_std.error), alpha = 0.2) +
	ylab("Fitness ratio") + xlab("Temperature (°C)")
ggsave("figures/fitness-ratio-param-space.png", width = 6, height = 4)

rs2 %>% 
	dplyr::select(T, contains("mean")) %>% 
	gather(contains("mean"), value = value, key = parameter) %>% 
	ggplot(aes(x = T, y = value, color = parameter)) + geom_line() +
	xlab("Temperature (°C)") + ylab("Parameter value") +
	facet_wrap( ~ parameter, scales = "free")
ggsave("figures/alphas-mean.png", width = 6, height = 4)
ggsave("figures/params-mean.png", width = 20, height = 15)



results_sum %>% 
	ggplot(aes(x = T, y = a11_mean)) + geom_line() +
	geom_line(aes(x = T, y = a12_mean), color = "pink") +
	geom_line(aes(x = T, y = a22_mean), color = "blue") +
	geom_line(aes(x = T, y = a21_mean), color = "green")
	
names(results_sum)	
results_sum %>% 
	ggplot(aes(x = T, y = a21_mean)) + geom_line() +
	geom_line(aes(x = T, y = a12_mean), color = "pink") 
	

results_range %>% 
	ggplot(aes(x = c_EA2N)) + geom_histogram()

results_range %>% 
	ggplot(aes(x = T, y = a11)) + geom_point()


p1 <- results_range %>% 
	ggplot(aes(x = T, y = fit_ratio, color = factor(r_EAP))) + geom_point(alpha = 0.5)

ggsave("figures/p1.pdf", width = 8, height = 6)	


### OK I honestly don't know what's going on here with the weird spiky values of alpha. I can't figure it out.
### in the mean time, I'm going to try to take the Bayesian approach.
### How would we do this?
### I think we need to fit a model to the data from the dataset, i.e. estimate the mean of the activation energies (this would be an intercept only model?)
### Then we get the posterior distribution
### Ok how do I do this??
### step 1 find the data.





