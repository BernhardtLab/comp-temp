

arrhenius_function <- function(Temp, E, b1, ref_temp = 1) {
	k <- 8.62e-05 #Boltzmann's constant
	E <- E # 0.6 # activation energy (eV)
	T <- Temp+273.15 #range of temp in K
	Tc <- ref_temp+273.15 #reference temperature
	
	metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
	return(metabolism)
}