


results_c_Ea2N <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur(T = seq(0, 40, by = 0.1), 
									  c_Ea2N  = sample_n(lm_consumption_rate, size = 1)$intercept)
	results_c_Ea2N   <- bind_rows(results_c_Ea2N, hold)
}

res_c_Ea2N  <- results_c_Ea2N  %>% 
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


plot1 <- res_c_Ea2N  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Log fitness ratio, absolute value") + xlab("Temperature")

plot2 <- res_c_Ea2N  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-c_Ea2N.jpeg", plot = plots, width = 10, height = 5)
