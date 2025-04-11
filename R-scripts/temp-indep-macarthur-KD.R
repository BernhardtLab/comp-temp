#KD macarthur function where all EAs are == 0
temp_indep_mac <- function(r_EaN = 0, r_EaP = 0, #activation energy growth rate N and P
                         c_Ea1N = 0, c_Ea1P = 0, #activation energy consumption rate N and P, species 1
                         c_Ea2N = 0, c_Ea2P = 0, #activation energy consumption rate N and P, species 2
                         K_EaN = 0, K_EaP = 0, #activation energy carrying capacity N and P
                         v_EaN = 0, v_EaP = 0, #activation energy conversion efficiency N & P (same for both species)
                         m_Ea1 = 0, m_Ea2 = 0, #activation energy mortality rate, species 1 and 2
                         c1N_b, c1P_b, #consumption rate of N and P at ref temp for species 1
                         c2N_b, c2P_b, #consumption rate of N and P at ref temp for species 2
                         r_N_b, r_P_b, #growth rate for each resource at ref temp
                         K_N_b, K_P_b, #carrying capacity for each resource at ref temp
                         v1N_b, v1P_b, #conversion efficiency for each resource at ref temp for species 1
                         v2N_b, v2P_b, #conversion efficiency for each resource at ref temp for species 2
                         m1_b, m2_b){ #mortality rate at ref temp for each species
  
  # resource growth rates
  rN = arrhenius_function(E = r_EaN, b1 = r_N_b)
  rP = arrhenius_function(E = r_EaP, b1 = r_P_b)
  
  # resource carrying capacity
  KN = arrhenius_function(E = K_EaN, b1 = K_N_b)
  KP = arrhenius_function(E = K_EaP, b1 = K_P_b)
  
  # cij = per capita consumption of comsumer i on resource j
  c1N = arrhenius_function(E = c_Ea1N, b1 = c1N_b)
  c1P = arrhenius_function(E = c_Ea1P, b1 = c1P_b) ## species 1 consumes more P than N
  c2N = arrhenius_function(E = c_Ea2N, b1 = c2N_b) ## species 2 consumes more N than P
  c2P = arrhenius_function(E = c_Ea2P, b1 = c2P_b)
  
  # vij = conversion factor that converts resource j into biomass of consumer i
  v1N = arrhenius_function(E = v_EaN, b1 = v1N_b)
  v2N = arrhenius_function(E = v_EaN, b1 = v2N_b) ## species 2 converts N more efficiently
  v1P = arrhenius_function(E = v_EaP, b1 = v1P_b) ## species 1 converts P more efficiently
  v2P = arrhenius_function(E = v_EaP, b1 = v2P_b)
  
  # mortality rates
  m1 = arrhenius_function(E = m_Ea1, b1 = m1_b)
  m2 = arrhenius_function(E = m_Ea2, b1 = m2_b)
  
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
  new_stabil_potential <- -log(rho)
  fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio
  new_fit_ratio <- log(fit_ratio)
  coexist <- rho < fit_ratio &  fit_ratio < 1/rho
  
  # report results
  data.frame(r_EaN = r_EaN, r_EaP = r_EaP,
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

#make this temperature independent -- not expecting a T
arrhenius_function <- function(E, b1) {
	k <- 8.62e-05 #Boltzmann's constant
	E <- E # 0.6 # activation energy (eV)
	T <- 298.15 #range of temp in K
	
	metabolism <- (b1*exp(1)^(E*(1/(k*T))))
	return(metabolism)
}




