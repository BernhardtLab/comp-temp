#This script defines the MacArthur consumer-resource function, where all parameters are assignable at the function deploying stage, and the Arrhenius temperature-depenence function, which feeds directly into the MacArthur C-R function. Temperature senstivities are referred to here as activation energy are parameters containing "EA." Intercept terms, given with the "_b" notation in parameter names, determine the value of the each function at the ambient temperature (Tref, ref_temp). Consumers are given by the numbers 1 and 2 and substitutable resources a and b are referred to as N and P, respectively. These functions are called in all subsequent analysis scripts that simulate warming: Scripts 03, 04, and 05.

#script author: Kaleigh Davis, UoG postdoc with Joey Bernhardt

temp_dep_mac_8r <- function(T, ref_temp, 
                        # temperature sensitivities
                         r_EaN, r_EaP, r_EaC, r_EaD, r_EaE, r_EaF, r_EaG, r_EaH, 
                         c_Ea1N, c_Ea1P, c_Ea1C, c_Ea1D, c_Ea1E, c_Ea1F, c_Ea1G, c_Ea1H,  
                         c_Ea2N, c_Ea2P, c_Ea2C, c_Ea2D, c_Ea2E, c_Ea2F, c_Ea2G, c_Ea2H, 
                         K_EaN, K_EaP, K_EaC, K_EaD, K_EaE, K_EaF, K_EaG, K_EaH,
                         v_EaN, v_EaP, v_EaC, v_EaD, v_EaE, v_EaF, v_EaG, v_EaH,
                         m_Ea1, m_Ea2, 
                        # base rates at reference temperature
                         c1N_b, c1P_b, c1C_b, c1D_b,  c1E_b, c1F_b, c1G_b, c1H_b, 
                         c2N_b, c2P_b, c2C_b, c2D_b,  c2E_b, c2F_b, c2G_b, c2H_b,
                         r_N_b, r_P_b, r_C_b, r_D_b,  r_E_b, r_F_b, r_G_b, r_H_b, 
                         K_N_b, K_P_b, K_C_b, K_D_b,  K_E_b, K_F_b, K_G_b, K_H_b, 
                         v1N_b, v1P_b, v1C_b, v1D_b,  v1E_b, v1F_b, v1G_b, v1H_b, 
                         v2N_b, v2P_b, v2C_b, v2D_b,  v2E_b, v2F_b, v2G_b, v2H_b, 
                         m1_b, m2_b) { 
  
  # resource growth rates
  rN = arrhenius_function(Temp = T, E = r_EaN, b1 = r_N_b, ref_temp = ref_temp)
  rP = arrhenius_function(Temp = T, E = r_EaP, b1 = r_P_b, ref_temp = ref_temp)
  rC = arrhenius_function(Temp = T, E = r_EaC, b1 = r_C_b, ref_temp = ref_temp)
  rD = arrhenius_function(Temp = T, E = r_EaD, b1 = r_D_b, ref_temp = ref_temp)
  rE = arrhenius_function(Temp = T, E = r_EaE, b1 = r_E_b, ref_temp = ref_temp)
  rF = arrhenius_function(Temp = T, E = r_EaF, b1 = r_F_b, ref_temp = ref_temp)
  rG = arrhenius_function(Temp = T, E = r_EaG, b1 = r_G_b, ref_temp = ref_temp)
  rH = arrhenius_function(Temp = T, E = r_EaH, b1 = r_H_b, ref_temp = ref_temp)
  
  # resource carrying capacity
  KN = arrhenius_function(Temp = T, E = K_EaN, b1 = K_N_b, ref_temp = ref_temp)
  KP = arrhenius_function(Temp = T, E = K_EaP, b1 = K_P_b, ref_temp = ref_temp)
  KC = arrhenius_function(Temp = T, E = K_EaC, b1 = K_C_b, ref_temp = ref_temp)
  KD = arrhenius_function(Temp = T, E = K_EaD, b1 = K_D_b, ref_temp = ref_temp)
  KE = arrhenius_function(Temp = T, E = K_EaE, b1 = K_E_b, ref_temp = ref_temp)
  KF = arrhenius_function(Temp = T, E = K_EaF, b1 = K_F_b, ref_temp = ref_temp)
  KG = arrhenius_function(Temp = T, E = K_EaG, b1 = K_G_b, ref_temp = ref_temp)
  KH = arrhenius_function(Temp = T, E = K_EaH, b1 = K_H_b, ref_temp = ref_temp)
  
  # cij = per capita consumption of consumer i on resource j
  c1N = arrhenius_function(Temp = T, E = c_Ea1N, b1 = c1N_b, ref_temp = ref_temp)
  c1P = arrhenius_function(Temp = T, E = c_Ea1P, b1 = c1P_b, ref_temp = ref_temp)
  c1C = arrhenius_function(Temp = T, E = c_Ea1C, b1 = c1C_b, ref_temp = ref_temp)
  c1D = arrhenius_function(Temp = T, E = c_Ea1D, b1 = c1D_b, ref_temp = ref_temp)
  c1E = arrhenius_function(Temp = T, E = c_Ea1E, b1 = c1E_b, ref_temp = ref_temp)
  c1F = arrhenius_function(Temp = T, E = c_Ea1F, b1 = c1F_b, ref_temp = ref_temp)
  c1G = arrhenius_function(Temp = T, E = c_Ea1G, b1 = c1G_b, ref_temp = ref_temp)
  c1H = arrhenius_function(Temp = T, E = c_Ea1H, b1 = c1H_b, ref_temp = ref_temp)
  
  c2N = arrhenius_function(Temp = T, E = c_Ea2N, b1 = c2N_b, ref_temp = ref_temp) 
  c2P = arrhenius_function(Temp = T, E = c_Ea2P, b1 = c2P_b, ref_temp = ref_temp)
  c2C = arrhenius_function(Temp = T, E = c_Ea2C, b1 = c2C_b, ref_temp = ref_temp)
  c2D = arrhenius_function(Temp = T, E = c_Ea2D, b1 = c2D_b, ref_temp = ref_temp)
  c2E = arrhenius_function(Temp = T, E = c_Ea2E, b1 = c2E_b, ref_temp = ref_temp)
  c2F = arrhenius_function(Temp = T, E = c_Ea2F, b1 = c2F_b, ref_temp = ref_temp)
  c2G = arrhenius_function(Temp = T, E = c_Ea2G, b1 = c2G_b, ref_temp = ref_temp)
  c2H = arrhenius_function(Temp = T, E = c_Ea2H, b1 = c2H_b, ref_temp = ref_temp)
  
  # vij = conversion factor that converts resource j into biomass of consumer i
  v1N = arrhenius_function(Temp = T, E = v_EaN, b1 = v1N_b, ref_temp = ref_temp)
  v1P = arrhenius_function(Temp = T, E = v_EaP, b1 = v1P_b, ref_temp = ref_temp)
  v1C = arrhenius_function(Temp = T, E = v_EaC, b1 = v1C_b, ref_temp = ref_temp)
  v1D = arrhenius_function(Temp = T, E = v_EaD, b1 = v1D_b, ref_temp = ref_temp)
  v1E = arrhenius_function(Temp = T, E = v_EaE, b1 = v1E_b, ref_temp = ref_temp)
  v1F = arrhenius_function(Temp = T, E = v_EaF, b1 = v1F_b, ref_temp = ref_temp)
  v1G = arrhenius_function(Temp = T, E = v_EaG, b1 = v1G_b, ref_temp = ref_temp)
  v1H = arrhenius_function(Temp = T, E = v_EaH, b1 = v1H_b, ref_temp = ref_temp)
  
  v2N = arrhenius_function(Temp = T, E = v_EaN, b1 = v2N_b, ref_temp = ref_temp) 
  v2P = arrhenius_function(Temp = T, E = v_EaP, b1 = v2P_b, ref_temp = ref_temp)
  v2C = arrhenius_function(Temp = T, E = v_EaC, b1 = v2C_b, ref_temp = ref_temp)
  v2D = arrhenius_function(Temp = T, E = v_EaD, b1 = v2D_b, ref_temp = ref_temp)
  v2E = arrhenius_function(Temp = T, E = v_EaE, b1 = v2E_b, ref_temp = ref_temp)
  v2F = arrhenius_function(Temp = T, E = v_EaF, b1 = v2F_b, ref_temp = ref_temp)
  v2G = arrhenius_function(Temp = T, E = v_EaG, b1 = v2G_b, ref_temp = ref_temp)
  v2H = arrhenius_function(Temp = T, E = v_EaH, b1 = v2H_b, ref_temp = ref_temp)
  
  # mortality rates
  m1 = arrhenius_function(Temp = T, E = m_Ea1, b1 = m1_b, ref_temp = ref_temp)
  m2 = arrhenius_function(Temp = T, E = m_Ea2, b1 = m2_b, ref_temp = ref_temp)
  
  # Absolute competition coefficients
  ### intra 11
  beta11 = v1N * c1N * (KN/rN) * c1N + v1P * c1P * (KP/rP) * c1P + v1C * c1C * (KC/rC) * c1C + v1D * c1D * (KD/rD) * c1D + v1E * c1E * (KE/rE) * c1E + v1F * c1F * (KF/rF) * c1F + v1G * c1G * (KG/rG) * c1G + v1H * c1H * (KH/rH) * c1H  
  
  ### inter 12
  beta12 = v1N * c1N * (KN/rN) * c2N + v1P * c1P * (KP/rP) * c2P + v1C * c1C * (KC/rC) * c2C + v1D * c1D * (KD/rD) * c2D + v1E * c1E * (KE/rE) * c2E + v1F * c1F * (KF/rF) * c2F + v1G * c1G * (KG/rG) * c2G + v1H * c1H * (KH/rH) * c2H
  
  ### intra 22
  beta22 = v2N * c2N * (KN/rN) * c2N + v2P * c2P * (KP/rP) * c2P + v2C * c2C * (KC/rC) * c2C + v2D * c2D * (KD/rD) * c2D + v2E * c2E * (KE/rE) * c2E + v2F * c2F * (KF/rF) * c2F + v2G * c2G * (KG/rG) * c2G + v2H * c2H * (KH/rH) * c2H  
  
  ### inter 21
  beta21 = v2N * c2N * (KN/rN) * c1N + v2P * c2P * (KP/rP) * c1P + v2C * c2C * (KC/rC) * c1C + v2D * c2D * (KD/rD) * c1D + v2E * c2E * (KE/rE) * c1E + v2F * c2F * (KF/rF) * c1F + v2G * c2G * (KG/rG) * c1G + v2H * c2H * (KH/rH) * c1H
  
  #In Song et al 2019, this is r_i
  g1 = v1N * c1N * KN + v1P * c1P * KP + v1C * c1C * KC + v1D * c1D * KD + v1E * c1E * KE + v1F * c1F * KF + v1G * c1G * KG + v1H * c1H * KH - m1 
  
  g2 = v2N * c2N * KN + v2P * c2P * KP + v2C * c2C * KC + v2D * c2D * KD + v2E * c2E * KE + v2F * c2F * KF + v2G * c2G * KG + v2H * c2H * KH - m2 
  
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
             r_EaN = r_EaN, r_EaP = r_EaP, r_EaC = r_EaC, r_EaD = r_EaD, r_EaE = r_EaE, r_EaF = r_EaF, r_EaG = r_EaG, r_EaH = r_EaH, 
             c_Ea1N = c_Ea1N, c_Ea1P = c_Ea1P, c_Ea1C = c_Ea1C, c_Ea1D = c_Ea1D, c_Ea1E = c_Ea1E, c_Ea1F = c_Ea1F, c_Ea1G = c_Ea1G, c_Ea1H = c_Ea1H,
             c_Ea2N = c_Ea2N, c_Ea2P = c_Ea2P, c_Ea2C = c_Ea2C, c_Ea2D = c_Ea2D, c_Ea2E = c_Ea2E, c_Ea2F = c_Ea2F, c_Ea2G = c_Ea2G, c_Ea2H = c_Ea2H,
             K_EaN = K_EaN, K_EaP = K_EaP, K_EaC = K_EaC, K_EaD = K_EaD, K_EaE = K_EaE, K_EaF = K_EaF, K_EaG = K_EaG, K_EaH = K_EaH,
             v_EaN = v_EaN, v_EaP = v_EaP, v_EaC = v_EaC, v_EaD = v_EaD, v_EaE = v_EaE, v_EaF = v_EaF, v_EaG = v_EaG, v_EaH = v_EaH, 
             m_Ea1 = m_Ea1, m_Ea2 = m_Ea2,
             c1N_b = c1N_b, c2P_b = c2P_b, c1C_b = c1C_b, c1D_b = c1D_b, c1E_b = c1E_b, c1F_b = c1F_b, c1G_b = c1G_b, c1H_b = c1H_b,
             c1P_b = c1P_b, c2N_b = c2N_b, c2C_b = c2C_b, c2D_b = c2D_b, c2E_b = c2E_b, c2F_b = c2F_b, c2G_b = c2G_b, c2H_b = c2H_b,
             r_N_b = r_N_b, r_P_b = r_P_b, r_C_b = r_C_b, r_D_b = r_D_b, r_E_b = r_E_b, r_F_b = r_F_b, r_G_b = r_G_b, r_H_b = r_H_b,
             K_N_b = K_N_b, K_P_b = K_P_b, K_C_b = K_C_b, K_D_b = K_D_b, K_E_b = K_E_b, K_F_b = K_F_b, K_G_b = K_G_b, K_H_b = K_H_b,
             v1N_b = v1N_b, v1P_b = v1P_b, v1C_b = v1C_b, v1D_b = v1D_b, v1E_b = v1E_b, v1F_b = v1F_b, v1G_b = v1G_b, v1H_b = v1H_b,
             v2N_b = v2N_b, v2P_b = v2P_b, v2C_b = v2C_b, v2D_b = v2D_b, v2E_b = v2E_b, v2F_b = v2F_b, v2G_b = v2G_b, v2H_b = v2H_b,
             m1_b = m1_b, m2_b = m2_b,
             a11 = a11, a12 = a12, a22 = a22, a21 = a21, g1 = g1, g2 = g2,
             stabil_potential = stabil_potential, new_stabil_potential = new_stabil_potential, fit_ratio = fit_ratio, new_fit_ratio = new_fit_ratio, rho = rho, coexist = coexist,
             m1 = m1, m2 = m2, rN = rN, rP = rP, rC = rC, rD = rD, rE = rE, rF = rF, rG = rG, rH = rH, KN = KN, KP = KP, KC = KC, KD = KD, KE = KE, KF = KF, KG = KG, KH = KH,
             c1N = c1N,  c1P = c1P,  c1C = c1C, c1D = c1D, c1E = c1E, c1F = c1F, c1G = c1G, c1H = c1H, c2N = c2N, c2P = c2P, c2C = c2C, c2D = c2D, c2E = c2E, c2F = c2F, c2G = c2G, c2H = c2H, beta11 = beta11, beta21 = beta21, beta22 = beta22, beta12 = beta12)}



