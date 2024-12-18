

library(tidyverse)
library(janitor)
library(MCMCpack)
library(bayesplot)
library(MCMCvis)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)




### Arrhenius function to model the temperature dependence
arrhenius_function <- function(Temp, E, b1, ref_temp = 25) {
	k <- 8.62e-05 #Boltzmann's constant
	E <- E # 0.6 # activation energy (eV)
	T <- Temp+273.15 #range of temp in K
	Tc <- ref_temp+273.15 #reference temperature
	
	metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
	return(metabolism)
}


mac_means <- read_csv("data/mac-means.csv") %>% 
	clean_names()




### mortality rates

lm_mort <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "mortality rate"), burnin = 1000) %>% 
	as.data.frame() %>% 
	clean_names() %>% 
	mutate(parameter = "mortality_rate")

mortality_rates <- mac_means %>% 
	filter(simple_parameter == "mortality rate") %>% 
	dplyr::select(activation_energy)

mort_ea_plot <- lm_mort %>% 
	ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
	xlab("Activation energy for mortality rate") +
	geom_point(aes(x = activation_energy, y = 0), data = mortality_rates, color = "orange", size = 3) +
	geom_point(aes(x = activation_energy, y = 0), data = mortality_rates, color = "black", size = 3, shape = 1) +
	geom_vline(aes(xintercept = mean(lm_mort$intercept))) +
	geom_vline(aes(xintercept = median(lm_mort$intercept)), color = "pink") +
	coord_flip()
# resource growth rate ----------------------------------------------------

lm_rgr <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "resource growth rate"), burnin = 1000) %>% 
	as.data.frame() %>% 
	clean_names() %>% 
	mutate(parameter = "resource_growth_rate")

growth_rates <- mac_means %>% 
	filter(simple_parameter == "resource growth rate") %>% 
	dplyr::select(activation_energy)

rgr_plot <- lm_rgr %>% 
	ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
	xlab("Activation energy for resource growth rate") +
	geom_point(aes(x = activation_energy, y = 0), data = growth_rates, color = "orange", size = 3) +
	geom_point(aes(x = activation_energy, y = 0), data = growth_rates, color = "black", size = 3, shape = 1) +
	geom_vline(aes(xintercept = mean(lm_rgr$intercept))) +
	geom_vline(aes(xintercept = median(lm_rgr$intercept)), color = "pink") +
	coord_flip()

# conversion efficiency ---------------------------------------------------
lm_conv_eff <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "conversion efficiency"), burnin = 1000) %>% 
	as.data.frame() %>% 
	clean_names() %>% 
	mutate(parameter = "conversion_efficiency")

conv_rates <- mac_means %>% 
	filter(simple_parameter == "conversion efficiency") %>% 
	dplyr::select(activation_energy)

conv_eff_plot <- lm_conv_eff %>% 
	ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
	xlab("Activation energy for conversion efficiency") +
	geom_point(aes(x = activation_energy, y = 0), data = conv_rates, color = "orange", size = 3) +
	geom_point(aes(x = activation_energy, y = 0), data = conv_rates, color = "black", size = 3, shape = 1) +
	coord_flip()

# resource carrying capacity  ---------------------------------------------------
lm_carrying_capacity <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "resource carrying capacity"), burnin = 1000) %>% 
	as.data.frame() %>% 
	clean_names() %>% 
	mutate(parameter = "carrying_capacity")

carrying_capacity <- mac_means %>% 
	filter(simple_parameter == "resource carrying capacity") %>% 
	dplyr::select(activation_energy)

carrying_capacity_plot <- lm_carrying_capacity %>% 
	ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
	xlab("Activation energy for carrying capacity") +
	geom_point(aes(x = activation_energy, y = 0), data = carrying_capacity, color = "orange", size = 3) +
	geom_point(aes(x = activation_energy, y = 0), data = carrying_capacity, color = "black", size = 3, shape = 1) +
	geom_vline(aes(xintercept = mean(lm_carrying_capacity$intercept))) +
	geom_vline(aes(xintercept = median(lm_carrying_capacity$intercept)), color = "pink") +
	coord_flip()



# consumption rate  ---------------------------------------------------
lm_consumption_rate <- MCMCregress(activation_energy ~ 1, data = filter(mac_means, simple_parameter == "consumption rate"), burnin = 1000) %>% 
	as.data.frame() %>% 
	clean_names() %>% 
	mutate(parameter = "consumption rate")

consumption_rate <- mac_means %>% 
	filter(simple_parameter == "consumption rate") %>% 
	dplyr::select(activation_energy)

consumption_rate_plot <- lm_consumption_rate %>% 
	ggplot(aes(x = intercept)) + geom_histogram(fill = "lightgrey") +
	xlab("Activation energy for consumption rate") +
	geom_point(aes(x = activation_energy, y = 0), data = consumption_rate, color = "orange", size = 3) +
	geom_point(aes(x = activation_energy, y = 0), data = consumption_rate, color = "black", size = 3, shape = 1) +
	geom_vline(aes(xintercept = mean(lm_consumption_rate$intercept))) +
	geom_vline(aes(xintercept = median(lm_consumption_rate$intercept)), color = "pink") +
	
	coord_flip()


ea_plots <- mort_ea_plot +  rgr_plot + conv_eff_plot + carrying_capacity_plot + consumption_rate_plot
# ggsave(filename = "figures/ea-plots.pdf", ea_plots, width = 15, height = 12)


# simulate fitness and niche params with 0.1 interval ---------------------

results_b25 <- data.frame()
for(f in 1:200){
	hold = temp_dependences_MacArthur(T = seq(0, 50, by = 0.1),
									  r_EaN = sample_n(lm_rgr, size = 1)$intercept,
									  r_EaP = sample_n(lm_rgr, size = 1)$intercept, 
									  K_EaN = sample_n(lm_carrying_capacity, size = 1)$intercept, 
									  K_EaP = sample_n(lm_carrying_capacity, size = 1)$intercept, 
									  c_Ea1N = sample_n(lm_consumption_rate, size = 1)$intercept,
									  c_Ea1P = sample_n(lm_consumption_rate, size = 1)$intercept, 
									  c_Ea2N = sample_n(lm_consumption_rate, size = 1)$intercept,
									  c_Ea2P = sample_n(lm_consumption_rate, size = 1)$intercept, 
									  v_EaN = sample_n(lm_conv_eff, size = 1)$intercept,
									  v_EaP = sample_n(lm_conv_eff, size = 1)$intercept, 
									  m_Ea1 = sample_n(lm_mort, size = 1)$intercept, 
									  m_Ea2 = sample_n(lm_mort, size = 1)$intercept)
	hold$iteration <- f
	results_b25 <- bind_rows(results_b25, hold)
}
# results05 <- results_b
# results10 <- results_b
# results15 <- results_b
# results20 <- results_b
# results25 <- results_b	
# results30 <- results_b	
# results35 <- results_b


write_csv(results_b25, "data-processed/results_b25.csv")

### now save activation energies and b's
results_b25b <- data.frame()
for(f in 1:200){
	hold = temp_dependences_MacArthurb(T = seq(0, 50, by = 0.1),
									  r_EaN = sample_n(lm_rgr, size = 1)$intercept,
									  r_EaP = sample_n(lm_rgr, size = 1)$intercept, 
									  K_EaN = sample_n(lm_carrying_capacity, size = 1)$intercept, 
									  K_EaP = sample_n(lm_carrying_capacity, size = 1)$intercept, 
									  c_Ea1N = sample_n(lm_consumption_rate, size = 1)$intercept,
									  c_Ea1P = sample_n(lm_consumption_rate, size = 1)$intercept, 
									  c_Ea2N = sample_n(lm_consumption_rate, size = 1)$intercept,
									  c_Ea2P = sample_n(lm_consumption_rate, size = 1)$intercept, 
									  v_EaN = sample_n(lm_conv_eff, size = 1)$intercept,
									  v_EaP = sample_n(lm_conv_eff, size = 1)$intercept, 
									  m_Ea1 = sample_n(lm_mort, size = 1)$intercept, 
									  m_Ea2 = sample_n(lm_mort, size = 1)$intercept)
	hold$iteration <- f
	# hold$r_EaN <- r_EaN
	results_b25b <- bind_rows(results_b25b, hold)
}

View(results_b25b)
head(results_b25b)
write_csv(results_b25b, "data-processed/results_b25b.csv")

# results05 <- results_b
# results10 <- results_b
# results15 <- results_b
# results20 <- results_b
# results25 <- results_b	
# results30 <- results_b	
# results35 <- results_b


write_csv(results_b25, "data-processed/results_b25.csv")




results_b %>% 
	ggplot(aes(x = T, y = abs(log(fit_ratio)), group = iteration)) + geom_line(alpha = 0.1) +
	ylab("Absolute value of Log Fitness ratio") + xlab("Temperature")
# ggsave("figures/fitness-ratio-lines-log-v-low-mortality.jpeg", width = 8, height = 6)
# ggsave("figures/fitness-ratio-lines-log-high-K.jpeg", width = 8, height = 6)

results_b %>% 
	ggplot(aes(x =(T-25), y = fit_ratio, group = iteration)) + geom_path(alpha = 0.1) +
	ylab("Fitness ratio") + xlab("(temperature - ref temp)")
# ggsave("figures/fitness-ratio-reftemp-20-relative-to-ref-temp.jpeg", width = 8, height = 6)


#### working graphs
res_sum <- results_b25b %>% ### update dec 18 2024
	group_by(T) %>% 
	summarise(mean_fit = mean(fit_ratio),
			  mean_stabil = mean(stabil_potential))

results_b25b %>% 
	ggplot() +
	geom_path(aes(x =(T-25), y = log(fit_ratio), group = iteration), alpha = 0.1) +
	geom_point(aes(x = T-25, y = log(mean_fit)), data = res_sum, color = "red", size = 0.5) +
	ylab("Fitness ratio") + xlab("(temperature - ref temp)")
ggsave("figures/averages.png", width = 8, height = 6)

results_b %>% 
	ggplot() +
	geom_path(aes(x =(T-25), y = (stabil_potential), group = iteration), alpha = 0.1) +
	geom_point(aes(x = T-25, y = (mean_stabil)), data = res_sum, color = "red", size = 0.5) +
	ylab("Stabilization potential") + xlab("(temperature - ref temp)")
ggsave("figures/averages-stabil.png", width = 8, height = 6)


results_b %>% 
	ggplot(aes(x = T, y = stabil_potential, group = iteration)) + geom_line(alpha = 0.5) +
	ylab("Stabilization potential") + xlab("Temperature")
# ggsave("figures/stabil-reftemp-20.jpeg", width = 8, height = 6)
ggsave("figures/stabil-reftemp-25.jpeg", width = 8, height = 6)


results_b %>% 
	ggplot(aes(x = T, y = fit_ratio, group = iteration, color = coexist)) + geom_path(alpha = 0.6) +
	ylab("Fitness ratio") + xlab("Temperature")
ggsave("figures/fitness-reftemp-25.jpeg", width = 8, height = 6)

results_b %>% 
	filter(iteration %in% c(1:45)) %>% 
	ggplot(aes(x = T, y = fit_ratio, group = iteration, color = coexist)) + geom_path(alpha = 0.6) +
	ylab("Fitness ratio") + xlab("Temperature") +
	facet_wrap( ~ iteration)
# ggsave("figures/fitness-reftemp-20.jpeg", width = 8, height = 6)
# ggsave("figures/fitness-reftemp-20.jpeg", width = 8, height = 6)


results_b %>% 
	# filter(iteration %in% c(1:50)) %>% 
	ggplot(aes(x = stabil_potential, y = fit_ratio, group = iteration, color = coexist)) + geom_path(size = 2) +
	ylab("Fitness ratio") + xlab("stabilization potential") +
	facet_wrap( ~ iteration)
ggsave("figures/fitness-stabil-reftemp-20.jpeg", width = 20, height = 16)
ggsave("figures/fitness-stabil-reftemp-20-1.2.jpeg", width = 20, height = 16)


results05$ref_temp <- 5
results10$ref_temp <- 10
results15$ref_temp <- 15
results20$ref_temp <- 20
results25$ref_temp <- 25
results35$ref_temp <- 35

all_results <- bind_rows(results05, results15, results25, results35, results10, results20)

results20 %>% 
	group_by(coexist) %>% 
	tally()

all_results %>% 
	group_by(iteration, coexist, ref_temp) %>% 
	tally() %>% 
	ggplot(aes(x = n, fill = coexist)) + geom_density(alpha = 0.5) + 
	facet_wrap( ~ ref_temp, scales = "free")
ggsave("figures/coexist-reftemp-25.jpeg", width = 20, height = 16)
ggsave("figures/coexist-reftemp-range.jpeg", width = 20, height = 16)
ggsave("figures/coexist-reftemp-range-2.jpeg", width = 20, height = 16)


all_results %>% 
	group_by(coexist, ref_temp) %>% 
	count() %>% 
	ggplot(aes(x = ref_temp, y = n, group = coexist, fill = coexist)) + geom_bar(stat = "identity")

results25 <- results_b	

#### calculate distances
res20 <- results_b %>% 
	filter(T == 0)

results_b <- results_b25


resb2 <- results_b %>% 
	mutate(stand_temp = T - 25) 

resb2b <- results_b %>% 
	mutate(stand_temp = T - 25) 

library(colorspace)

ggplot() +
	geom_path(data = resb2, aes(x = stabil_potential, y = fit_ratio, color = stand_temp, group = iteration), size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(resb2$stabil_potential)*0.99, max(resb2$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	scale_colour_continuous_diverging() +
		xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
		ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 
ggsave("figures/stabil-fitness-ref-temp-25.jpeg", width = 8, height = 6)
ggsave("figures/stabil-fitness-ref-temp-25-diverging.jpeg", width = 8, height = 6)
ggsave("figures/stabil-fitness-ref-temp-1.jpeg", width = 8, height = 6)

res3_warming <- resb2 %>% 
	filter(stand_temp > 0)

ggplot() +
	geom_path(data = res3_warming, aes(x = stabil_potential, y = fit_ratio, color = stand_temp, group = iteration), size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(res3_warming$stabil_potential)*0.99, max(res3_warming$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	scale_colour_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 
ggsave("figures/stabil-fitness-ref-temp-25-warming.jpeg", width = 8, height = 6)

res3_cooling <- resb2 %>% 
	filter(stand_temp < 0)

ggplot() +
	geom_path(data = res3_cooling, aes(x = stabil_potential, y = fit_ratio, color = stand_temp, group = iteration), size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(res3_cooling$stabil_potential)*0.99, max(res3_cooling$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	scale_colour_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 
ggsave("figures/stabil-fitness-ref-temp-25-cooling.jpeg", width = 8, height = 6)






res20b <- res20 %>% 
	mutate(std_temp = T-0)

results_bb <- results_b %>% 
	mutate(std_temp = T-0)

str(results_bb)


resmat <- results_bb %>% 
	dplyr::select(fit_ratio, stabil_potential) %>% 
	as.matrix()

res20b_mat <- res20b %>% 
	dplyr::select(fit_ratio, stabil_potential) %>% 
	distinct() %>% 
	as.matrix()
library(fields)

distances <-  rdist(res20b_mat, resmat) %>% 
	as.data.frame() %>% 
	gather()

all_distances <- bind_cols(results_bb, distances)

all_distances %>% 
	ggplot(aes(x = std_temp, y = value)) + geom_point(alpha = 0.1)
ggsave("figures/distances-ref0-all.png", width = 6, height = 4)

all_distances %>% 
	group_by(std_temp) %>% 
	summarise(mean_distance = mean(value)) %>% 
	ggplot(aes(x = std_temp, y = mean_distance)) + geom_point()
ggsave("figures/distances-ref25.png", width = 6, height = 4)
ggsave("figures/distances-ref1.png", width = 6, height = 4)
ggsave("figures/distances-ref0.png", width = 6, height = 4)


all_distances %>% 
	ggplot(aes(x = std_temp, y = value, color = coexist)) + geom_point(alpha = 0.1)+
	facet_wrap( ~ iteration, scales = "free")
ggsave("figures/distances-ref0-all-coexist-facs-free.png", width = 20, height = 20)

### come back here to do with ref temp = 1
	

ggplot() +
	geom_path(data = results_bb, aes(x = stabil_potential, y = fit_ratio, color = std_temp, group = iteration), size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(results_b$stabil_potential)*0.99, max(results_b$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	geom_point(data = res20b, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + scale_color_viridis_c() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 
ggsave("figures/stabil-fitness-ref-temp-25.jpeg", width = 8, height = 6)





#### Now find the number of cases where they coexist at 20, but then with warming they no longer co-exist

results_b %>%
	filter(T != 20) %>% 
	group_by(coexist) %>% 
	tally()

results_b

# rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
# stabil_potential <- 1 - rho #stabilizing potential
# fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio
# coexist <- rho < fit_ratio &  fit_ratio < 1/rho


# now do the thing where I manipulate each param by a percent -------------
## manipulate each by 25% of the max, and then manipulate each by  25% of its value
## start by finding out: which one has the biggest range?



results_b %>% 
	# filter(iteration == "415") %>% 
	gather(key = parameter, value = value, 2:26) %>% 
	ggplot(aes(x = T, y = value, color = parameter)) + geom_line() +
	facet_wrap(~ parameter, scales = "free")
ggsave("figures/params-vlow-mort.jpeg", width = 12, height = 12)


results_b %>% 
	filter(iteration == "415") %>% 
	dplyr::select(T, a11, beta11, g1) %>% 
	mutate(new_alpha = beta11/g1) %>% 
	gather(key = parameter, value = value, a11, beta11, g1, new_alpha) %>% 
	ggplot(aes(x = T, y = value, color = parameter)) + geom_line() +
	facet_wrap(~ parameter, scales = "free")
ggsave("figures/params-415-a11.jpeg", width = 8, height = 4)

results_b %>% 
	filter(iteration == "415") %>% 
	dplyr::select(T, a11, beta11, g1) %>% 
	mutate(new_alpha = beta11/g1) %>% 
	gather(key = parameter, value = value, beta11, g1) %>% 
	ggplot(aes(x = T, y = value, color = parameter)) + geom_line(size = 2) +
	geom_hline(yintercept = 0)
ggsave("figures/params-415-a11.jpeg", width = 8, height = 4)


# resultsb <- results_b
res_sum <- results_b %>% 
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(fit_ratio = log(fit_ratio)) %>% ### log transforming the fitness ratio
	mutate(fit_ratio = abs(fit_ratio)) %>%  ### taking the absolute value of the fitness ratio
	mutate(mean_fit = mean(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = mean(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 
# write_csv(res_sum, "data-processed/sum_bayes_param40_0.1.csv")
write_csv(res_sum, "data-processed/sum_bayes_param40_0.1_vlow_mortality.csv")
write_csv(res_sum, "data-processed/sum_bayes_param-high-k.csv")




res_sum <- read_csv("data-processed/sum_bayes_param40.csv")

### 50 percentiles

plot1 <- res_sum %>%
	ggplot() + 
	# geom_point(alpha = 0.1) +
	# geom_line(aes(x = T, y = lower_q50), color = "green", alpha = 1) +
	# geom_line(aes(x = T, y = upper_q50), color = "green", alpha = 1) +
	# geom_line(aes(x = T, y = lower_q95), color = "blue", alpha = 1) +
	# geom_line(aes(x = T, y = upper_q95), color = "blue", alpha = 1) +
	
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Log fitness ratio, absolute value") + xlab("Temperature")

plot2 <- res_sum %>% 
	ggplot() + 
	# geom_point(alpha = 0.1) +
	# geom_line(aes(x = T, y = lower_q50), color = "green", alpha = 1) +
	# geom_line(aes(x = T, y = upper_q50), color = "green", alpha = 1) +
	# geom_line(aes(x = T, y = lower_q95), color = "blue", alpha = 1) +
	# geom_line(aes(x = T, y = upper_q95), color = "blue", alpha = 1) +
	
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")
	

plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-0.1-low-mortality.jpeg", plot = plots, width = 10, height = 5)
ggsave(filename = "figures/bayes-parms-high-k.jpeg", plot = plots, width = 10, height = 5)
ggsave(filename = "figures/bayes-parms-high-k-log.jpeg", plot = plots, width = 10, height = 5)
ggsave(filename = "figures/bayes-parms-high-k-log-abs.jpeg", plot = plots, width = 10, height = 5)

# ggsave(filename = "figures/bayes-parms-0.1-50per.jpeg", plot = plots, width = 10, height = 5)


res_sum %>% 
	ggplot() + 
	# geom_point(alpha = 0.1) +
	geom_ribbon(aes(x = T, ymin = lower_q, ymax  = upper_q), fill = "grey", alpha = 0.1) +
	# geom_line(aes(x = T, y = upper_q), color = "blue", alpha = 1) +
	geom_line(aes(x = T, y = mean_fit), color = "pink", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")


### come back after dinner to try the thing where we only vary one of the parameters at a time, keep the others constant.


### update May 7 2024 -- coming back here because this seems promising
### Ok i need to look at this function 'temp_dependences_MacArthur' to figure out how I set up the scenario


# REAP --------------------------------------------------------------------


results_r_EaP <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  r_EaP  = sample_n(lm_rgr, size = 1)$intercept)
	results_r_EaP   <- bind_rows(results_r_EaP, hold)
}

res_r_EaP  <- results_r_EaP  %>% 
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 
write_csv(res_r_EaP, "data-processed/res_r_EaP_zero.csv")


plot1 <- res_r_EaP  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_r_EaP  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-r_EaP_zero.jpeg", plot = plots, width = 10, height = 5)


ggplot() +
	geom_point(data = res_r_EaP, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_r_EaP, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_r_EaP$stabil_potential)*0.99, max(res_r_EaP$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_r_EaP-biplot-afternoon-2.jpeg", width = 8, height= 6)


### ok what if I try Po-Ju's new definition of niche and fitness differences
#niche difference as -log(\rho) and fitness difference as log(f1/f2).
# fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio

res_r_EaP2  <- results_r_EaP  %>% 
	mutate(iteration = rownames(.)) %>% 
	mutate(stabil_potential = -log(rho)) %>% 
	mutate(fit_ratio = log(fit_ratio)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 

plot1 <- res_r_EaP2  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_r_EaP2  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")

plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-r_EaP_new_zero.jpeg", plot = plots, width = 10, height = 5)



res_r_EaP2 %>%
	filter((g1 > 0) & (g2 > 0)) %>%
	ggplot(aes(x = stabil_potential, y = fit_ratio, col = T)) +
	geom_point() +
	geom_point(aes(x = stabil_potential, y = fit_ratio), data = filter(res_r_EaP2, T == 25), color = "red", size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(res_r_EaP2$stabil_potential)*0.99, max(res_r_EaP2$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL, 
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.3) +
	geom_hline(yintercept = 1, linetype=5) +
	scale_x_continuous(expand=c(0, 0)) + 
	scale_color_viridis_c("Temperature", end = 0.8) +
	xlab(expression(paste("Niche difference, -log(", rho, ")"))) +
	ylab(expression(paste("Fitness difference, log(", f[2], "/", f[1], ")"))) + 
	ggtitle("(e)") +
	theme(legend.position = "bottom")
ggsave(filename = "figures/bayes-parms-r_EaP_new_coexist_zero.jpeg", width = 6, height = 5)

res_r_EaP2 %>%
	filter((g1 > 0) & (g2 > 0)) %>%
	ggplot(aes(x = stabil_potential, y = fit_ratio, col = T)) +
	geom_point() +
	geom_point(aes(x = stabil_potential, y = fit_ratio), data = filter(res_r_EaP2, T == 25), color = "red", size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(res_r_EaP2$stabil_potential)*0.99, max(res_r_EaP2$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.3) +
	geom_hline(yintercept = 1, linetype=5) +
	geom_abline(slope =1, intercept = 0) +
	scale_x_continuous(expand=c(0, 0)) + 
	scale_color_viridis_c("Temperature", end = 0.8) +
	xlab(expression(paste("Niche difference, -log(", rho, ")"))) +
	ylab(expression(paste("Fitness difference, log(", f[2], "/", f[1], ")"))) + 
	ggtitle("(e)") +
	theme(legend.position = "bottom")
ggsave(filename = "figures/bayes-parms-r_EaP_new_coexist-slope1_zero.jpeg", width = 6, height = 5)





# REAN --------------------------------------------------------------------


results_r_EaN <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  r_EaN  = sample_n(lm_rgr, size = 1)$intercept)
	results_r_EaN   <- bind_rows(results_r_EaN, hold)
}

res_r_EaN  <- results_r_EaN  %>% 
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975))
write_csv(res_r_EaN, "data-processed/res_r_EaN_zero.csv")


plot1 <- res_r_EaN  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_r_EaN  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-r_EaN_zero.jpeg", plot = plots, width = 10, height = 5)


ggplot() +
	geom_point(data = res_r_EaN, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_r_EaN, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_r_EaN$stabil_potential)*0.99, max(res_r_EaN$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_r_EaN-biplot-afternoon.jpeg", width = 8, height= 6)




# kEAN --------------------------------------------------------------------


results_K_EaN <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  K_EaN = sample_n(lm_carrying_capacity, size = 1)$intercept) 
	results_K_EaN  <- bind_rows(results_K_EaN , hold)
}

res_K_EaN <- results_K_EaN %>% 
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 


plot1 <- res_K_EaN %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")


plot2 <- res_K_EaN %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-K_EaN_zero.jpeg", plot = plots, width = 10, height = 5)

ggplot() +
	geom_point(data = res_K_EaN, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_K_EaN, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_K_EaN$stabil_potential)*0.99, max(res_K_EaN$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_K_EaN-biplot-afternoon.jpeg", width = 8, height= 6)

# K_EaP -------------------------------------------------------------------

results_K_EaP <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  K_EaP = sample_n(lm_carrying_capacity, size = 1)$intercept)
	hold$iteration <- f
	results_K_EaP  <- bind_rows(results_K_EaP , hold)
}

res_K_EaP <- results_K_EaP %>% 
	# mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 


plot1 <- res_K_EaP %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_K_EaP %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-K_EaP_zero.jpeg", plot = plots, width = 10, height = 5)


ggplot() +
	geom_point(data = res_K_EaP, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_K_EaP, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_K_EaP$stabil_potential)*0.99, max(res_K_EaP$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_K_EaP-biplot-afternoon.jpeg", width = 8, height= 6)



res_K_EaP2 <- results_K_EaP %>% 
	mutate(iteration = rownames(.)) 

ggplot() +
	geom_point(data = results_K_EaP, aes(x = stabil_potential, y = fit_ratio, color = T), size = 2) +
	# geom_point(data = filter(res_K_EaP, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	# geom_ribbon(data = data.frame(x = seq(min(results_K_EaP$stabil_potential)*0.99, max(results_K_EaP$stabil_potential)*1.01, 0.001)),
	# 			aes(x = x,
	# 				y = NULL,
	# 				ymin = 1-x,
	# 				ymax = 1/(1-x)),
	# 			fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))



results_K_EaP %>% 
	ggplot(aes(x = T, y = fit_ratio)) + geom_point(size = 0.25)


results_K_EaP %>% 
	ggplot(aes(x = T, y = fit_ratio, color = T, group = iteration)) + geom_path(size = 0.25, alpha = 0.25)
results_K_EaP %>% 
	ggplot(aes(x = T, y = stabil_potential, color = T, group = iteration)) + geom_line(size = 0.25, alpha = 0.25)

results_K_EaP %>% 
	ggplot(aes(x = stabil_potential, y = fit_ratio, color = T)) + geom_point(size = 0.25, alpha = 0.25)



#### now we need to do the same thing for the other metrics



# c_Ea1P -------------------------------------------------------------------

results_c_Ea1P <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  c_Ea1P = sample_n(lm_consumption_rate, size = 1)$intercept)
	results_c_Ea1P  <- bind_rows(results_c_Ea1P , hold)
}

res_c_Ea1P <- results_c_Ea1P %>%
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 


plot1 <- res_c_Ea1P %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_c_Ea1P %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-c_Ea1P_zero.jpeg", plot = plots, width = 10, height = 5)


res_c_Ea1P %>% 
	ggplot() + 
	# geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	# geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")

ggplot() +
	geom_point(data = res_c_Ea1P, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_c_Ea1P, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_c_Ea1P$stabil_potential)*0.99, max(res_c_Ea1P$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_c_Ea1P-biplot.jpeg", width = 8, height= 6)




# c_Ea1N -------------------------------------------------------------------

results_c_Ea1N <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  c_Ea1N = sample_n(lm_consumption_rate, size = 1)$intercept)
	results_c_Ea1N  <- bind_rows(results_c_Ea1N, hold)
}

res_c_Ea1N <- results_c_Ea1N %>%
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 


plot1 <- res_c_Ea1N %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_c_Ea1N %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-c_Ea1N_zero.jpeg", plot = plots, width = 10, height = 5)



# v_EaN -------------------------------------------------------------------

results_v_EaN <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  v_EaN = sample_n(lm_conv_eff, size = 1)$intercept)
	results_v_EaN  <- bind_rows(results_v_EaN, hold)
}

res_v_EaN <- results_v_EaN %>%
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 


plot1 <- res_v_EaN %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_v_EaN %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-v_EaN_zero.jpeg", plot = plots, width = 10, height = 5)



# v_EaP -------------------------------------------------------------------

results_v_EaP <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  v_EaP = sample_n(lm_conv_eff, size = 1)$intercept)
	results_v_EaP  <- bind_rows(results_v_EaP, hold)
}

res_v_EaP <- results_v_EaP %>%
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 


plot1 <- res_v_EaP %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_v_EaP %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-v_EaP_zero.jpeg", plot = plots, width = 10, height = 5)






# m_Ea1 -------------------------------------------------------------------

results_m_Ea1 <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  m_Ea1 = sample_n(lm_mort, size = 1)$intercept)
	results_m_Ea1  <- bind_rows(results_m_Ea1, hold)
}

res_m_Ea1 <- results_m_Ea1 %>%
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 


plot1 <- res_m_Ea1 %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_m_Ea1 %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-m_Ea1_zero.jpeg", plot = plots, width = 10, height = 5)



ggplot() +
	geom_point(data = res_m_Ea1, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_m_Ea1, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_m_Ea1$stabil_potential)*0.99, max(res_m_Ea1$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_m_Ea1-biplot.jpeg", width = 8, height= 6)


# m_Ea2 -------------------------------------------------------------------

results_m_Ea2 <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_zero(T = seq(0, 40, by = 0.1), 
									  m_Ea2 = sample_n(lm_mort, size = 1)$intercept)
	results_m_Ea2  <- bind_rows(results_m_Ea2, hold)
}

res_m_Ea2 <- results_m_Ea2 %>%
	mutate(iteration = rownames(.)) %>% 
	group_by(T) %>% 
	mutate(mean_fit = median(fit_ratio)) %>% 
	mutate(lower_q50 = quantile(fit_ratio, prob = 0.25)) %>% 
	mutate(upper_q50 = quantile(fit_ratio, prob = 0.75))  %>% 
	mutate(lower_q95 = quantile(fit_ratio, prob = 0.025)) %>% 
	mutate(upper_q95 = quantile(fit_ratio, prob = 0.975)) %>% 
	mutate(mean_fit_s = median(stabil_potential)) %>% 
	mutate(lower_q50s = quantile(stabil_potential, prob = 0.25)) %>% 
	mutate(upper_q50s = quantile(stabil_potential, prob = 0.75))  %>% 
	mutate(lower_q95s = quantile(stabil_potential, prob = 0.025)) %>% 
	mutate(upper_q95s = quantile(stabil_potential, prob = 0.975)) 


plot1 <- res_m_Ea2 %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_m_Ea2 %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-m_Ea2_zero.jpeg", plot = plots, width = 10, height = 5)


ggplot() +
	geom_point(data = res_m_Ea2, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_m_Ea2, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_m_Ea2$stabil_potential)*0.99, max(res_m_Ea2$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_m_Ea2-biplot.jpeg", width = 8, height= 6)




#### now plot them all on a common scale to see how the effects compare!

res_c_Ea1N <- res_c_Ea1N %>% 
	mutate(param = "res_c_Ea1N")
res_c_Ea1P <- res_c_Ea1P %>% 
	mutate(param = "res_c_Ea1P")
res_K_EaN <- res_K_EaN %>% 
	mutate(param = "res_K_EaN")
res_K_EaP <- res_K_EaP %>% 
	mutate(param = "res_K_EaP")
res_m_Ea1 <- res_m_Ea1 %>% 
	mutate(param = "res_m_Ea1")
res_m_Ea2 <- res_m_Ea2 %>% 
	mutate(param = "res_m_Ea2")
res_r_EaN <- res_r_EaN %>% 
	mutate(param = "res_r_EaN")
res_r_EaP <- res_r_EaP %>% 
	mutate(param = "res_r_EaP")
res_v_EaN <- res_v_EaN %>% 
	mutate(param = "res_v_EaN")
res_v_EaP <- res_v_EaP %>% 
	mutate(param = "res_v_EaP")


all_results <- bind_rows(res_c_Ea1N,
						 res_c_Ea1P,
						 res_K_EaN,
						 res_K_EaP,
						 res_m_Ea1,
						 res_m_Ea2,
						 res_r_EaN,
						 res_r_EaP,
						 res_v_EaN,
						 res_v_EaP)



all_results %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature") + facet_wrap( ~ param)
ggsave(filename = "figures/fitness_diffs_zero.jpeg", width = 20, height = 25)

all_results %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential")  + xlab("Temperature") + facet_wrap( ~ param)
ggsave(filename = "figures/stabilization_diffs_zero.jpeg", width = 20, height = 25)




### note from Po-Ju
### In this past couple years, more and more people are starting to define niche difference as -log(\rho) and fitness difference as log(f1/f2).



# now make pom-pom plot with the empirical data ---------------------------
### Update May 20 2024
ggplot() +
	geom_point(data = res_r_EaP, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_r_EaP, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_r_EaP$stabil_potential)*0.99, max(res_r_EaP$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_r_EaP-biplot.jpeg", width = 8, height= 6)
### something is weird here!



ggplot() +
	geom_point(data = res_K_EaP, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_K_EaP, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_K_EaP$stabil_potential)*0.99, max(res_K_EaP$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_K_EaP-biplot.jpeg", width = 8, height= 6)



ggplot() +
	geom_point(data = res_K_EaN, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_K_EaN, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_K_EaN$stabil_potential)*0.99, max(res_K_EaN$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_K_EaN-biplot.jpeg", width = 8, height= 6)



ggplot() +
	geom_point(data = res_r_EaN, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_r_EaN, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_r_EaN$stabil_potential)*0.99, max(res_r_EaN$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
	geom_hline(yintercept = 1.5, linetype=5, color = "pink")
ggsave("figures/res_r_EaN-biplot.jpeg", width = 8, height= 6)



ggplot() +
	geom_point(data = res_c_Ea1N, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_c_Ea1N, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_c_Ea1N$stabil_potential)*0.99, max(res_c_Ea1N$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_c_Ea1N-biplot.jpeg", width = 8, height= 6)

ggplot() +
	geom_point(data = res_c_Ea1P, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_c_Ea1P, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_c_Ea1P$stabil_potential)*0.99, max(res_c_Ea1P$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_c_Ea1P-biplot.jpeg", width = 8, height= 6)

ggplot() +
	geom_point(data = res_v_EaN, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_v_EaN, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_v_EaN$stabil_potential)*0.99, max(res_v_EaN$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")")))
ggsave("figures/res_v_EaN-biplot.jpeg", width = 8, height= 6)

ggplot() +
	geom_point(data = res_v_EaP, aes(x = mean_fit_s, y = mean_fit, color = T), size = 2) +
	geom_point(data = filter(res_v_EaP, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(res_v_EaP$stabil_potential)*0.99, max(res_v_EaP$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
	xlim(min(res25c4$stabil_potential)*0.99, max(res25c4$stabil_potential)*1.01)
ggsave("figures/res_v_EaP-biplot-2.jpeg", width = 8, height= 6)



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




# let's do a scenario that favours one species and make all the ot --------
## May 24
#let's favour species 2, which prefers N

### Come back here December 18 2024 (figure out where the temp_dependences_MacArthur_zero comes from)

results_favourite <- data.frame()
for(f in 1:100){
	hold = temp_dependences_MacArthur_ref0(T = seq(0, 50, by = 0.1), ### Dec 2024; I can't find this function --- (temp_dependences_MacArthur_zero); so I'm swapping it for (temp_dependences_MacArthur_ref0)
										   c_Ea2N = sample_n(lm_consumption_rate, size = 1)$intercept,
										v_EaN = sample_n(lm_conv_eff, size = 1)$intercept,
										r_EaN = sample_n(lm_rgr, size = 1)$intercept)
	hold$iteration <- f
	results_favourite   <- bind_rows(results_favourite, hold)
}


ggplot() +
	geom_path(data = results_favourite, aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 1) +
	# geom_point(data = filter(res_summary, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(results_favourite$stabil_potential)*0.99, max(results_favourite$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
	scale_color_viridis_c() +
	theme(axis.text=element_text(size=24),
		  axis.title=element_text(size=24,face="bold"))
ggsave("figures/pompom-favourite-1000.png",  width = 12, height = 10) ### December 13 2024 this looks like latest script


results_favourite_warm <- results_favourite %>% 
	filter(T > 25)
results_favourite_warm2 <- results_favourite %>% 
	filter(T > 25) %>% 
	group_by(T) %>% 
	summarise(mean_fit = mean(fit_ratio),
			  mean_stabil = mean(stabil_potential))


end <- results_favourite_warm %>% 
	filter(T == 50) %>%
	summarise(mean_fit = mean(fit_ratio),
			  mean_stabil = mean(stabil_potential))

start <- results_favourite_warm %>% 
	filter(T == 25.1) %>%
	summarise(mean_fit = mean(fit_ratio),
			  mean_stabil = mean(stabil_potential))
	
# 
# group_by(coexist) %>% 
# 	tally()
# (6/1000)*100
# (994/1000)*100

#niche difference as -log(\rho) and fitness difference as log(f1/f2).


ggplot() +
	geom_path(data = results_favourite_warm, aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 1) +
	geom_point(data = filter(results_favourite_warm2, T == 35), aes(x = mean_stabil, y = mean_fit), size = 5, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(results_favourite_warm$stabil_potential)*0.99, max(results_favourite_warm$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
	scale_color_viridis_c() +
	theme(axis.text=element_text(size=24),
		  axis.title=element_text(size=24,face="bold"))
ggsave("figures/pompom-favourite-warm-1000-35.png",  width = 12, height = 10)

ggplot() +
	geom_path(data = results_favourite_warm, aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 1) +
	geom_point(data = filter(results_favourite_warm2, T == 25.1), aes(x = mean_stabil, y = mean_fit), size = 5, color = "pink") +
	geom_point(data = end, aes(x = mean_stabil, y = mean_fit), size = 5, color = "red") +
	geom_ribbon(data = data.frame(x = seq(min(results_favourite_warm$stabil_potential)*0.99, max(results_favourite_warm$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
	scale_color_viridis_c() +
	theme(axis.text=element_text(size=24),
		  axis.title=element_text(size=24,face="bold"))



results_loser <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur_ref0(T = seq(0, 50, by = 0.1), ## swapped this function -- don't think it's the same one I used before
										   c_Ea2P = sample_n(lm_consumption_rate, size = 1)$intercept,
										   v_EaP = sample_n(lm_conv_eff, size = 1)$intercept,
										   r_EaP = sample_n(lm_rgr, size = 1)$intercept)
	hold$iteration <- f
	results_loser   <- bind_rows(results_loser, hold)
}


ggplot() +
	geom_path(data = results_loser, aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 1) +
	# geom_point(data = filter(res_summary, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(results_loser$stabil_potential)*0.99, max(results_loser$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
	scale_color_viridis_c() +
	theme(axis.text=element_text(size=24),
		  axis.title=element_text(size=24,face="bold"))
ggsave("figures/pompom-loser-1000.png",  width = 12, height = 10)


results_loser_warm <- results_loser %>% 
	filter(T > 25)

results_loser_warm2 <- results_loser %>% 
	filter(T > 25) %>% 
	group_by(T) %>% 
	summarise(mean_fit = mean(fit_ratio),
			  mean_stabil = mean(stabil_potential))

ggplot() +
	geom_path(data = results_loser_warm, aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 1) +
	geom_point(data = filter(results_loser_warm2, T == 50), aes(x = mean_stabil, y = mean_fit), size = 5, color = "pink") +
	geom_ribbon(data = data.frame(x = seq(min(results_loser_warm$stabil_potential)*0.99, max(results_loser_warm$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) +
	scale_color_viridis_c() +
	theme(axis.text=element_text(size=24),
		  axis.title=element_text(size=24,face="bold"))
ggsave("figures/pompom-loser-warm-1000-50.png",  width = 12, height = 10)

