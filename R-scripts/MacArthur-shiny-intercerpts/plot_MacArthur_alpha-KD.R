plot_MacArthur_alpha <- function(Data.temperature){
	
	# panel a: alpha ~ temperature | species
	Data.alpha = Data.temperature[, c("T", "a11", "a12", "a22", "a21")] %>% gather(value=value, key=parameter, -T)
	panel.a = 
		ggplot(Data.alpha, 
			   aes(x = T, y = value, col = parameter)) +
		geom_point() +
		scale_color_manual(values=c("a11"="#08519c", "a21"="#6baed6", "a22"="#a50f15", "a12"="#fb6a4a")) + 
		xlab("Temperature") +
		ylab(expression(alpha[ij])) + 
		ggtitle("(a)") +
		theme(legend.position = "bottom")
	
	# panel b: relative impact ~ temperature | species
	Temp = data.frame(T = Data.temperature$T, 
					  rel.impact1 = Data.temperature$a21 / Data.temperature$a11, 
					  rel.impact2 = Data.temperature$a12 / Data.temperature$a22)
	Temp = Temp %>% gather(value=value, key=parameter, -T)
	panel.b = 
		ggplot(Temp, 
			   aes(x = T, y = value, col = parameter)) +
		geom_point() +
		scale_color_manual(values=c("rel.impact1"="#08519c", "rel.impact2"="#a50f15")) + 
		xlab("Temperature") +
		ylab(expression(paste("Relative impact (", alpha[ji], "/", alpha[ii], ")"))) + 
		ggtitle("(b)") +
		theme(legend.position = "bottom")
	
	# panel c: sensitivity to competition ~ temperature | species
	Temp2 = data.frame(T = Data.temperature$T, 
					   sen1 = Data.temperature$a12 * Data.temperature$a11, 
					   sen2 = Data.temperature$a21 * Data.temperature$a22)
	Temp2 = Temp2 %>% gather(value=value, key=parameter, -T)
	panel.c = 
		ggplot(Temp2, 
			   aes(x = T, y = value, col = parameter)) +
		geom_point() +
		scale_color_manual(values=c("sen1"="#08519c", "sen2"="#a50f15")) + 
		xlab("Temperature") +
		ylab(expression(paste("Sensitivity (", alpha[ii], "*", alpha[ij], ")"))) + 
		ggtitle("(c)") +
		theme(legend.position = "bottom")
	
	# Combine panels
	return(panel.a + panel.b + panel.c + plot_layout(nrow=1))
	
}

