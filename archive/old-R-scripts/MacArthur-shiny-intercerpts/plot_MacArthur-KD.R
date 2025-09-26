plot_MacArthur <- function(Data.parameter, Data.temperature){
	
	# panel a: parameter value ~ temperature | species
	# panel.a = 
	# 	ggplot(Data.parameter, 
	# 		   aes(x = T, y = value, col = parameter)) +
	# 	geom_point() +
	# 	xlab("Temperature") +
	# 	ylab(label = "Parameter value") + 
	# 	ggtitle("(a)") +
	# 	theme(legend.position = "bottom")
	# 
	# # panel b: consumption vector slope ~ temperature | species
	# Temp = data.frame(T = Data.temperature$T, 
	# 				  slope1 = Data.temperature$c1P / Data.temperature$c1N, 
	# 				  slope2 = Data.temperature$c2P / Data.temperature$c2N)
	# Temp = Temp %>% gather(value=value, key=parameter, -T)
	# panel.b = 
	# 	ggplot(Temp, 
	# 		   aes(x = T, y = value, col = parameter)) +
	# 	geom_point() +
	# 	geom_hline(yintercept = 1, linetype=5) +
	# 	xlab("Temperature") +
	# 	ylab("Consumption of P relative to N") + 
	# 	ggtitle("(b)") +
	# 	theme(legend.position = "bottom")
	# 
	# # panel c: ND ~ temperature
	# panel.c = 
	# 	Data.temperature %>%
	# 	filter((g1 > 0) & (g2 > 0)) %>%
	# 	ggplot(aes(x = T, y = stabil_potential, col = T)) +
	# 	geom_point() +
	# 	scale_color_viridis_c("Temperature", end = 0.8) +
	# 	xlab("Temperature") +
	# 	ylab(expression(paste("Stabilization potential (1-", rho, ")"))) + 
	# 	ggtitle("(c)") +
	# 	theme(legend.position = "bottom")
	# 
	# # panel d: FD ~ temperature
	# panel.d = 
	# 	Data.temperature %>%
	# 	filter((g1 > 0) & (g2 > 0)) %>%
	# 	ggplot(aes(x = T, y = fit_ratio, col = T)) +
	# 	geom_point() +
	# 	scale_color_viridis_c("Temperature", end = 0.8) +
	# 	xlab("Temperature") +
	# 	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
	# 	ggtitle("(d)") +
	# 	theme(legend.position = "bottom")
	# 
	# panel e: FD ~ ND
	panel.e = 
		Data.temperature %>%
		filter((g1 > 0) & (g2 > 0)) %>%
		ggplot(aes(x = stabil_potential, y = fit_ratio, col = T)) +
		geom_point(size = 4) +
	  geom_ribbon(data = data.frame(x = seq(0, 0.5, 0.001)),
	                            aes(x = x,
	                                y = NULL,
	                                ymin = exp(-x),
	                                ymax = 1/(exp(-x))),
	              fill = "grey", color = "black", alpha = 0.3) +
	#   geom_ribbon(data = data.frame(x = seq(min(Data.temperature$stabil_potential)*0.99, max(Data.temperature$stabil_potential)*1.01, 0.001)),
	# 				aes(x = x,
	# 					y = NULL, 
	# 					ymin = 1-x,
	# 					ymax = 1/(1-x)),
	# 				fill = "grey", color = "black", alpha = 0.3) +
		geom_hline(yintercept = 1, linetype=5) +
		# scale_x_continuous(expand=c(0, 0)) + 
		# scale_color_viridis_c("Temperature", end = 0.8) +
	    # coord_cartesian(ylim=c(-0.5, 2.5), xlim = c(0, 0.5)) +
	    xlab(expression(paste("Stabilization potential (-log(", rho, "))"))) +
	    ylab(expression(paste("Fitness difference (log(", f[2], "/", f[1], "))"))) +
		# xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
		# ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
		ggtitle("(e)") +
		theme(legend.position = "bottom")
	
	# Combine panels
	return(panel.e)
	
}

