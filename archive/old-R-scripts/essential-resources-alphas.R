
library(tidyverse)
all_rstars <- read_csv("data/all-rstars.csv")

#### goal here is to find the alphas depending on what zone we are in

### We use Arrhenius function to model the temperature dependence
arrhenius_function <- function(Temp, E, b1, ref_temp = 20) {
	k<-8.62e-05 #Boltzmann's constant
	E <- E # 0.6 # activation energy (eV)
	temperature <- Temp+273.15 #range of temp in K
	Tc <- ref_temp+273.15 #reference temperature
	
	metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*temperature))))
	return(metabolism)
}


find_alphas <- function(temperature, r_Ea = 0.6, c_Ea = 0.6, k_Ea = 0, m_Ea = 0.2, ref_temp = 20){
	snippet <- all_rstars %>% 
		filter(ancestor_id %in% c("anc4", "anc5")) %>% 
		filter(treatment == "none")
	
	pop1 <- snippet$population[[1]]
	pop2 <- snippet$population[[2]]
	
	c1P_b <- snippet$pc[snippet$population == pop1]
	c2P_b <- snippet$pc[snippet$population == pop2]
	c1N_b <- snippet$nc[snippet$population == pop1]
	c2N_b <- snippet$nc[snippet$population == pop2]
	R1P_b <- snippet$p_star[snippet$population == pop1]
	R2N_b <- snippet$n_star[snippet$population == pop2]
	R2P_b <- snippet$p_star[snippet$population == pop2]
	R1N_b <- snippet$n_star[snippet$population == pop1]
	k1N_b <- snippet$n_ks[snippet$population == pop1]
	k2N_b <- snippet$n_ks[snippet$population == pop2]
	k1P_b <- snippet$p_ks[snippet$population == pop1]
	k2P_b <- snippet$p_ks[snippet$population == pop2]
	D <- 0.5
	
	## how can we define the consumption vector lines?
	
	SN <- 1000/300
	SP <- 140/350
	
	r1_b <- max(snippet$p_umax[snippet$population == pop1], snippet$n_umax[snippet$population == pop1])
	r2_b <- max(snippet$p_umax[snippet$population == pop2], snippet$n_umax[snippet$population == pop2])
	
	c_Ea = 0.3
	k_Ea = 0.5
	r_Ea = 0.3
	m_Ea = 0.65
	m1_b <- 0.1
	m2_b <- 0.1
	
	c1P = arrhenius_function(Temp = temperature, E = c_Ea, b1 = c1P_b)
	c1N = arrhenius_function(Temp = temperature, E = c_Ea, b1 = c1N_b)
	c2P = arrhenius_function(Temp = temperature, E = c_Ea, b1 = c2P_b)
	c2N = arrhenius_function(Temp = temperature, E = c_Ea, b1 = c2N_b) ### cij = per capita consumption of comsumer i on resource j
	k1N = arrhenius_function(Temp = temperature, E = k_Ea, b1 = k1N_b)
	k2P = arrhenius_function(Temp = temperature, E = k_Ea, b1= k2P_b) #half saturation constant for N resource consumption
	k1P = arrhenius_function(Temp = temperature, E = k_Ea, b1 = k1P_b)
	k2N = arrhenius_function(Temp = temperature, E = k_Ea, b1= k2N_b)
	r1 = arrhenius_function(Temp = temperature, E = r_Ea, b1 = r1_b)
	r2 = arrhenius_function(Temp = temperature, E = r_Ea, b1= r2_b, ref_temp = ref_temp2) #population growth rates
	m1 = arrhenius_function(Temp = temperature, E = m_Ea, b1 = m1_b)
	m2 = arrhenius_function(Temp = temperature, E = m_Ea, b1 = m2_b) # mortality rates
	
	R1N <- (m1 * k1N)/ (r1 - m1)
	R1P <- (m1 * k1P)/ (r1 - m1)
	R2P <- (m2 * k2P)/ (r2 - m2)
	R2N <- (m2 * k2N)/ (r2 - m2)
	
	cons_vec1_intercept <- max(R1P, R2P) -(c1P/c1N)*max(R1N, R2N)
	cons_vec2_intercept <- max(R1P, R2P) -(c2P/c2N)*max(R1N, R2N)
	supply_vec <- SP/SN
	
	cons_vec1_fun <- function(x){
		y <- (c1P/c1N)*x + cons_vec1_intercept
		return(y)
	}
	
	cons_vec2_fun <- function(x){
		y <- (c2P/c2N)*x + cons_vec2_intercept
		return(y)
	}
	
	supply_vec_fun <- function(x){
		y <- (SP/SN)*x
		return(y)
	}
	
	zone_middle <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN) | cons_vec1_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec2_fun(SN)
	zone_bottom <- cons_vec2_fun(SN) >= supply_vec_fun(SN) & supply_vec_fun(SN) <= cons_vec1_fun(SN)
	zone_top <- cons_vec2_fun(SN) <= supply_vec_fun(SN) & supply_vec_fun(SN) >= cons_vec1_fun(SN)
	
	zones <- c("zone_middle", "zone_bottom", "zone_top")
	zone <- zones[c(zone_middle, zone_bottom, zone_top)]
	
	alphas <- function(zone) {
		if (zone == "zone_middle") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_top") {
			a11 <- c1N / (D * (SN - R1N))
			a12 <- c2N / (D * (SN - R1N))
			a21 <- c1N / (D * (SN - R2N))
			a22 <- c2N / (D * (SN - R2N))
		} else if (zone == "zone_bottom") {
			a11 <- c1P / (D * (SP - R1P))
			a12 <- c2P / (D * (SP - R1P))
			a21 <- c1P / (D * (SP - R2P))
			a22 <- c2P / (D * (SP - R2P))
		}
		alphas1 <- data.frame(a11 = a11, a12 = a12, a21 = a21, a22 = a22)
		return(alphas1)
	}
	
	alphas_calc <- alphas(zone)
	r_s <- data.frame(r1 = r1, r2 = r2)
	comps <- data.frame(pop1 = pop1, pop2 = pop2, D = D, treatment = snippet$treatment[1], zone = zone)
	
	
	
	rho <- sqrt((alphas_calc$a12*alphas_calc$a21)/(alphas_calc$a11*alphas_calc$a22)) #niche overlap
	stabil_potential <- 1 - rho #stabilizing potential
	fit_ratio <- sqrt((alphas_calc$a11*alphas_calc$a12)/(alphas_calc$a22*alphas_calc$a21))  #fitness ratio -- ask Patrick why he scaled his by the r's
	coexist <- rho < fit_ratio &  fit_ratio < 1/rho
	
	
	rs <- data.frame(R1N = R1N, R1P = R1P, R2N = R2N, R2P= R2P,
			   K1 = (r1)/alphas_calc$a11, K2 = (r2)/alphas_calc$a22, temperature = temperature, m1 = m1, m2 = m2, c1P = c1P,  c1N = c1N,
			   c2P = c2P, c2N = c2N, stabil_potential = stabil_potential,
			   fit_ratio = fit_ratio, rho = rho, coexist = coexist,
					 a11 = alphas_calc$a11, a12 = alphas_calc$a12, a21 = alphas_calc$a21, a22 = alphas_calc$a22)
	
	alphas_calc2 <- bind_cols(r_s, comps, rs)
	return(alphas_calc2)
}



results <- data.frame()
for(i in 1:350){
	for(j in 1:6){
		for(z in 1:6){
		r_Ea <- 0.3
	 	c_Ea <- 0.1*(z-1)
		# c_Ea <- 0.5
		m_Ea <- 0.1*(j-1)
		temperature <- 0.1*(i-1)
		hold <- find_alphas(temperature = temperature, r_Ea = r_Ea , c_Ea = c_Ea, m_Ea = m_Ea, ref_temp = 20)
		hold$r_Ea <- r_Ea
		hold$c_Ea <- c_Ea
		hold$m_Ea <- m_Ea
		results <- bind_rows(results, hold)
	}}}


unique(results$rho)
length(unique(results$a11))
max(results$fit_ratio)
min(results$fit_ratio)

results %>% 
	# filter(fit_ratio < 2) %>% 
	# filter(T %in% c(0, 10, 20, 30, 40)) %>% 
	ggplot(aes(x = R1N, y = R1P, color = temperature)) + geom_point() +
	geom_segment(aes(x = R1N, y = R1P, xend = R1N + c1N/20, yend = R1P + c1P/20, color = temperature),
				 size = 0.5, linetype = 1, alpha = 1) +
	geom_segment(aes(x = R1N, y = R1P, xend = R1N, yend = 0.2, color = temperature),
				 size = 0.5, linetype = 1, alpha = 1) +
	geom_segment(aes(x = R1N, y = R1P, xend = 5, yend = R1P, color = temperature),
				 size = 0.5, linetype = 1, alpha = 1) +
	ylab("P (uM)") + xlab("N (uM)") +
	scale_colour_viridis_c() +
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	coord_cartesian() +
	# xlim(0, 0.15) + ylim(0, 0.15) +
	theme( 
		plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
		axis.text = element_text(size=13),
		axis.title=element_text(size=14)) +
	panel_border(colour = "black") 
ggsave("figures/essential-resources-tilman-plots.png", width = 14, height = 8)


x <- seq(0.5,1.5,length = 100)
lines.df <- data.frame(x = x, y = 1/(x), y2 = x)

results %>% 
	filter(fit_ratio < 2) %>% 
ggplot(aes(x = stabil_potential, y = fit_ratio, color = temperature, shape = coexist))+
	geom_hline(yintercept = 1, color = "grey", linetype = 2)+
	geom_vline(xintercept = 0, color = "grey", linetype = 2)+
	geom_line(data = lines.df, aes(x = 1-x, y = y, color = NULL, shape = NULL))+
	geom_line(data = lines.df, aes(x = 1-x, y = y2, color = NULL, shape = NULL))+
	geom_point()+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	scale_color_viridis_c(option = "magma", end = 0.8)+
	scale_shape_manual(values = c(0,19)) +
	ylab("Fitness ratio") + xlab("Niche differences") +
	xlim(-0.1, 0.1)
ggsave("figures/essential-resources-chesson-plots.png", width = 14, height = 8)

results %>% 
	# filter(r1 > 0) %>%
	ggplot(aes(x = temperature, y = R2N, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "N* (species 2)") + xlab("Temperature (°C)")

results %>% 
	# filter(r1 > 0) %>%
	ggplot(aes(x = temperature, y = R2P, group = m_Ea, color = m_Ea))+
	geom_line()+
	scale_color_viridis_c(end = 0.8)+
	facet_grid(c_Ea ~ m_Ea, scales = "free") +
	ylab(label = "P* (species 1)") + xlab("Temperature (°C)")
ggsave("figures/essential-resources-P-star-plots.png", width = 14, height = 8)
