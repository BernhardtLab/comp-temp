#This script is to implement a unimodal model in the temperature dependence function, in order to explore how that affects competitive outcomes. It is a daughter script to 02-temp-dep-macarthur.R

#from parent script: This script defines the MacArthur consumer-resource function, where all parameters are assignable at the function deploying stage. Temperature sensitivities are referred to here as activation energy are parameters containing "EA." Intercept terms, given with the "_b" notation in parameter names, determine the value of the each function at the ambient temperature (Tref, ref_temp). Consumers are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively. These functions are called in all subsequent analysis scripts that simulate warming: Scripts 03, 04, and 05.

# note: this script defines its own copy of arrhenius_function() (see bottom
#   of script), which is functionally identical to the one in
#   R-scripts/03-arrhenius.R. This script is therefore self-contained and
#   does not need to source script 03.

#script author: Kaleigh Davis, UoG postdoc with Joey Bernhardt
#script DOB: [date]


#johnson-Lewin model for temperature dependence 
# #41 in Kontopoulos et al 2024 supplemental
jl_function2 <- function(Temp, E, b1, ED, Topt_C, ref_temp) {
  T <- Temp + 273.15 #Temp, Kelvin
  k <- 8.62e-05 #Boltzmann's constant
  E <- E
  ED <- ED
  Topt <- Topt_C + 273.15 #Topt in Kelvin
  Tc = ref_temp + 273.15
  
  increase <- b1*((exp(E*(1/(k*Tc) - 1/(k*T))))) 
  decrease <- 1 + (E/(ED - E))*exp((ED/k)*(((1/Topt) - (1/Tc))-((1/T) - (1/Tc))))
  
  metabolism <- increase / decrease
  return(metabolism)
}


# testing that wrote the JL model properly
model1 <- tibble(temp = seq(from = 5, to = 25, by = 0.01)) %>% 
  mutate(metabolism = map_dbl(temp, 
                              ~jl_function2(Temp = .x, 
                                           E = -0.65, 
                                           b1 = 1,
                                           ED = 4.5,
                                           Topt_C = 19,
                                           ref_temp = 10)))

ggplot(model1, aes(x = temp, y = metabolism)) + 
  geom_point() + 
  labs(x = "Temperature", y = "Metabolic rate")

# ggsave(filename = "figures/jl-tpc.pdf", plot = last_plot(), device = "pdf")

#temperature dependent CR function with diff Topts for two consumer species, Topt at max temp for resource; mortality increases exponentially, not unimodally with warming
uni_temp_dep_mac_spec_diffs <- function(T, ED, Topt_C1, Topt_C2, Topt_Cr, ref_temp,
                             r_EaN, r_EaP, #activation energy for resource growth rate N and P
                             c_Ea1N, c_Ea1P, #activation energy for the consumption rate N and P, species 1
                             c_Ea2N, c_Ea2P, #activation energy consumption rate N and P, species 2
                             K_EaN, K_EaP, #activation energy carrying capacity N and P
                             v_EaN, v_EaP, #activation energy conversion efficiency N & P (same for both consumer species)
                             m_Ea1, m_Ea2, #activation energy mortality rate, species 1 and 2
                             c1N_b, c1P_b, #consumption rate of N and P at ref temp for species 1
                             c2N_b, c2P_b, #consumption rate of N and P at ref temp for species 2
                             r_N_b, r_P_b, #growth rate for each resource at ref temp
                             K_N_b, K_P_b, #carrying capacity for each resource at ref temp
                             v1N_b, v1P_b, #conversion efficiency for each resource at ref temp for species 1
                             v2N_b, v2P_b, #conversion efficiency for each resource at ref temp for species 2
                             m1_b, m2_b){ #mortality rate at ref temp for each species
  
  # resource growth rates
  rN = jl_function2(Temp = T, E = r_EaN, b1 = r_N_b, ED = ED, Topt_C = Topt_Cr, ref_temp = ref_temp)
  rP = jl_function2(Temp = T, E = r_EaP, b1 = r_P_b, ED = ED, Topt_C = Topt_Cr, ref_temp = ref_temp)
  
  # resource carrying capacity
  KN = arrhenius_function(Temp = T, E = K_EaN, b1 = K_N_b, ref_temp = ref_temp)
  KP = arrhenius_function(Temp = T, E = K_EaP, b1 = K_P_b, ref_temp = ref_temp)
  
  # cij = per capita consumption of consumer i on resource j
  c1N = jl_function2(Temp = T, E = c_Ea1N, b1 = c1N_b, ED = ED, Topt_C = Topt_C1, ref_temp = ref_temp)
  c1P = jl_function2(Temp = T, E = c_Ea1P, b1 = c1P_b, ED = ED, Topt_C = Topt_C1, ref_temp = ref_temp)
  c2N = jl_function2(Temp = T, E = c_Ea2N, b1 = c2N_b, ED = ED, Topt_C = Topt_C2, ref_temp = ref_temp) 
  c2P = jl_function2(Temp = T, E = c_Ea2P, b1 = c2P_b, ED = ED, Topt_C = Topt_C2, ref_temp = ref_temp)
  
  # vij = conversion factor that converts resource j into biomass of consumer i
  v1N = jl_function2(Temp = T, E = v_EaN, b1 = v1N_b, ED = ED, Topt_C = Topt_C1, ref_temp = ref_temp)
  v2N = jl_function2(Temp = T, E = v_EaN, b1 = v2N_b, ED = ED, Topt_C = Topt_C1, ref_temp = ref_temp) 
  v1P = jl_function2(Temp = T, E = v_EaP, b1 = v1P_b, ED = ED, Topt_C = Topt_C2, ref_temp = ref_temp)
  v2P = jl_function2(Temp = T, E = v_EaP, b1 = v2P_b, ED = ED, Topt_C = Topt_C2, ref_temp = ref_temp)
  
  # mortality rates
  m1 = arrhenius_function(Temp = T, E = m_Ea1, b1 = m1_b, ref_temp = ref_temp)
  m2 = arrhenius_function(Temp = T, E = m_Ea2, b1 = m2_b, ref_temp = ref_temp)
  
  # Absolute competition coefficients
  beta11 = v1N * c1N * (KN/rN) * c1N + v1P * c1P * (KP/rP) * c1P ### intra
  beta12 = v1N * c1N * (KN/rN) * c2N + v1P * c1P * (KP/rP) * c2P ### inter
  beta22 = v2N * c2N * (KN/rN) * c2N + v2P * c2P * (KP/rP) * c2P ### intra
  beta21 = v2N * c2N * (KN/rN) * c1N + v2P * c2P * (KP/rP) * c1P ### inter
  
  #In Song et al 2019, this is r_i
  g1 = v1N * c1N * KN + v1P * c1P * KP - m1 ### growth rate of consumer 1
  g2 = v2N * c2N * KN + v2P * c2P * KP - m2 ### growth rate of consumer 2
  
  # Relative competition coefficients
  a11 = beta11 / g1 #increased growth rate --> decreased alpha
  a21 = beta21 / g2
  a22 = beta22 / g2
  a12 = beta12 / g1
  
  # MCT components
  rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
  stabil_potential <- 1 - rho #stabilizing potential
  new_stabil_potential <- -log(rho)
  fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio = k2/k1
  new_fit_ratio <- log(fit_ratio)
  coexist <- rho < fit_ratio &  fit_ratio < 1/rho
  
  # report results
  data.frame(T = T, ED = ED, Topt_Cr = Topt_Cr, Topt_C1 = Topt_C1, Topt_C2 = Topt_C2, ref_temp = ref_temp, 
             r_EaN = r_EaN, r_EaP = r_EaP,
             c_Ea1N = c_Ea1N, c_Ea1P = c_Ea1P,
             c_Ea2N = c_Ea2N, c_Ea2P = c_Ea2P,
             K_EaN = K_EaN, K_EaP = K_EaP,
             v_EaN = v_EaN, v_EaP = v_EaP,
             m_Ea1 = m_Ea1, m_Ea2 = m_Ea2,
             c1N_b = c1N_b, c2P_b = c2P_b,
             c1P_b = c1P_b, c2N_b = c2N_b,
             r_N_b = r_N_b, r_P_b = r_P_b,
             K_N_b = K_N_b, K_P_b = K_P_b,
             v1N_b = v1N_b, v1P_b = v1P_b,
             v2N_b = v2N_b, v2P_b = v2P_b,
             m1_b = m1_b, m2_b = m2_b,
             a11 = a11, a12 = a12, a22 = a22, a21 = a21, g1 = g1, g2 = g2,
             stabil_potential = stabil_potential, new_stabil_potential = new_stabil_potential, fit_ratio = fit_ratio, new_fit_ratio = new_fit_ratio, rho = rho, coexist = coexist,
             m1 = m1, m2 = m2, rN = rN, rP = rP, KN = KN, KP = KP,
             c1N = c1N,  c1P = c1P,  c2N = c2N, c2P = c2P, beta11 = beta11, beta21 = beta21, beta22 = beta22, beta12 = beta12)}


#arrhenius -- mortality
arrhenius_function <- function(Temp, E, b1, ref_temp) {
	k <- 8.62e-05 #Boltzmann's constant
	E <- E # 0.6 # activation energy (eV)
	T <- Temp+273.15 #range of temp in K
	Tc <- ref_temp+273.15 #reference temperature
	
	metabolism <- (b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
	return(metabolism)
}




