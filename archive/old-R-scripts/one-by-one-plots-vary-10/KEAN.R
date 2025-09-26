# KEAN --------------------------------------------------------------------


results_k_EaN <- data.frame()
for(f in 1:200){
	hold = temp_dependences_MacArthur(T = seq(0, 40, by = 0.1), 
									  K_EaN  = -0.3 * runif(1, 0.5, 1.5))
	results_k_EaN   <- bind_rows(results_k_EaN, hold)
}

res_k_EaN  <- results_k_EaN  %>% 
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


plot1 <- res_k_EaN  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50, ymax  = upper_q50), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95, ymax  = upper_q95), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit), color = "red", alpha = 1) +
	ylab("Log fitness ratio, absolute value") + xlab("Temperature")

plot2 <- res_k_EaN  %>% 
	ggplot() + 
	geom_ribbon(aes(x = T, ymin = lower_q50s, ymax  = upper_q50s), fill = "pink", alpha = 0.5) +
	geom_ribbon(aes(x = T, ymin = lower_q95s, ymax  = upper_q95s), fill = "lightblue", alpha = 0.5) +
	geom_line(aes(x = T, y = mean_fit_s), color = "red", alpha = 1) +
	ylab("Stabilization potential") + xlab("Temperature") 


plots <- plot1 + plot2
ggsave(filename = "figures/bayes-parms-k_EaN-10-ds.jpeg", plot = plots, width = 10, height = 5)


results_k_EaN %>% 
	ggplot(aes(x = T, y = fit_ratio)) + geom_path(alpha = 0.1)

results_k_EaN %>% 
	ggplot(aes(x = T, y = stabil_potential)) + geom_path(alpha = 0.1)


ggplot() +
	geom_path(data = results_k_EaN, aes(x = stabil_potential, y = fit_ratio, color = T), size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(results_k_EaN$stabil_potential)*0.99, max(results_k_EaN$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	geom_point(data = results_k_EaN, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + scale_color_viridis_c() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 

ggsave("figures/k-scenario-0.5.png", width = 8, height = 6)
