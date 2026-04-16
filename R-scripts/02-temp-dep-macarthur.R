#This script defines the MacArthur consumer-resource function, where all parameters are assignable at the function deploying stage, which feeds directly into the MacArthur consumer resource function. Temperature sensitivities are referred to here as activation energy, and are parameters containing "EA." Intercept terms, given with the "_b" notation in parameter names, determine the value of each function at the ambient temperature (Tref, ref_temp). Consumers are given by the numbers 1 and 2 and substitutable resources (referred to as a and b in the manuscript text) are referred to here as N and P, respectively. This temperature-dependent consumer resource function is used in all subsequent analysis scripts that simulate warming: Scripts 03, 04, and 05.

#script author: Kaleigh Davis, UoG postdoc with Joey Bernhardt

temp_dep_mac <- function(T, ref_temp, 
                         r_EaN, r_EaP, #activation energies for resource growth rate  for resources N and P
                         c_Ea1N, c_Ea1P, #activation energies for the consumption rate of resources N and P, consumer species 1
                         c_Ea2N, c_Ea2P, #activation energy for the consumption rate of resources N and P, consumer species 2
                         K_EaN, K_EaP, #activation energies carrying capacity of resources N and P
                         v_EaN, v_EaP, #activation energies for the conversion efficiency of resources N & P (same for both consumer species)
                         m_Ea1, m_Ea2, #activation energies for mortality rate, species 1 and 2
                         c1N_b, c1P_b, #consumption rate of resources N and P at ref temp for species 1
                         c2N_b, c2P_b, #consumption rate of resources N and P at ref temp for species 2
                         r_N_b, r_P_b, #growth rate for each resource, N and P, at the ref temp
                         K_N_b, K_P_b, #carrying capacity for each resource, N and P, at ref temp
                         v1N_b, v1P_b, #conversion efficiency for each resource, N and P, at ref temp for species 1
                         v2N_b, v2P_b, #conversion efficiency for each resource, N and P, at ref temp for species 2
                         m1_b, m2_b){ #mortality rate at ref temp for each consumer species
  
  # resource growth rates
  rN = arrhenius_function(Temp = T, E = r_EaN, b1 = r_N_b, ref_temp = ref_temp)
  rP = arrhenius_function(Temp = T, E = r_EaP, b1 = r_P_b, ref_temp = ref_temp)
  
  # resource carrying capacity
  KN = arrhenius_function(Temp = T, E = K_EaN, b1 = K_N_b, ref_temp = ref_temp)
  KP = arrhenius_function(Temp = T, E = K_EaP, b1 = K_P_b, ref_temp = ref_temp)
  
  # cij = per capita consumption of consumer i on resource j
  c1N = arrhenius_function(Temp = T, E = c_Ea1N, b1 = c1N_b, ref_temp = ref_temp)
  c1P = arrhenius_function(Temp = T, E = c_Ea1P, b1 = c1P_b, ref_temp = ref_temp)
  c2N = arrhenius_function(Temp = T, E = c_Ea2N, b1 = c2N_b, ref_temp = ref_temp) 
  c2P = arrhenius_function(Temp = T, E = c_Ea2P, b1 = c2P_b, ref_temp = ref_temp)
  
  # vij = conversion factor that converts resource j into biomass of consumer i
  v1N = arrhenius_function(Temp = T, E = v_EaN, b1 = v1N_b, ref_temp = ref_temp)
  v2N = arrhenius_function(Temp = T, E = v_EaN, b1 = v2N_b, ref_temp = ref_temp) 
  v1P = arrhenius_function(Temp = T, E = v_EaP, b1 = v1P_b, ref_temp = ref_temp)
  v2P = arrhenius_function(Temp = T, E = v_EaP, b1 = v2P_b, ref_temp = ref_temp)
  
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
  data.frame(T = T, ref_temp = ref_temp,
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
             m1 = m1, m2 = m2, rN = rN, rP = rP, KN = KN, KP = KP, v1N = v1N, v2N = v2N, v1P = v1P, v2P = v2P,
             c1N = c1N,  c1P = c1P,  c2N = c2N, c2P = c2P, beta11 = beta11, beta21 = beta21, beta22 = beta22, beta12 = beta12)}



