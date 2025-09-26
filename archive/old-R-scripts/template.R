


results_parm <- data.frame()
for(f in 1:1000){
	hold = temp_dependences_MacArthur(T = seq(0, 40, by = 0.1), 
									  parm  = sample_n(lm_rgr, size = 1)$intercept)
	results_parm   <- bind_rows(results_parm, hold)
}

res_parm  <- results_parm  %>% 
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


plot1 <- res_parm  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Fitness difference") + xlab("Temperature")

plot2 <- res_parm  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature")


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parm.jpeg", plot = plots, width = 10, height = 5)