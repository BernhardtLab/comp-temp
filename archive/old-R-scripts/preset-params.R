


## all


results_all <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur(T = seq(0, 40, by = 0.1))
	results_all   <- bind_rows(results_all, hold)
}

res_all  <- results_all  %>% 
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


plot1 <- res_all  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_all  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-all.jpeg", plot = plots, width = 10, height = 5)



results_all %>%
	dplyr::select(T, a11, a12, a21, a22, rho) %>% 
	gather(key = parameter, value = value, a11, a12, a21, a22, rho) %>%
	ggplot(aes(x = T, y = value, color = parameter)) + geom_line()
