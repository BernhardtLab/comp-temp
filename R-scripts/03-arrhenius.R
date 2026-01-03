#this script defines the Arrhenius function, which describes the rate of increase of parameters over a temperature gradient. This function feeds into the various different temp-dep-macarthur functions

#author: KD
#script DOB - 1/2/2026

arrhenius_function <- function(Temp, E, b1, ref_temp) {
  k <- 8.62e-05 #Boltzmann's constant
  E <- E # activation energy (eV)
  T <- Temp+273.15 #range of temp in K
  Tc <- ref_temp+273.15 #reference temperature
  
  metabolism <- (b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
  return(metabolism)
}
