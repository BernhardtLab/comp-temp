
### update March 13 trying to understand what's going on
#### Plotting
library(colorspace)
library(tidyverse)

res25 <- read_csv("data-processed/results_b25.csv")
res25 <- read_csv("data-processed/results_b25b.csv")
res25 <- data.frame()
res25 <- read_csv("data-processed/results_b25b-May.csv")
View(res25)
res0 <- read_csv("data-processed/results_ref0.csv")


res0b <- res0 %>% 
	mutate(diff_from25 = T - 25) %>% 
	mutate(stand_temp = T - 0)

res25b <- res25 %>% 
	mutate(diff_from25 = T - 25) %>% ## I am guessing this is assuming a reference temperature of 25C
	mutate(stand_temp = T - 25) ### this is getting us standardized temperature -- I guess a distance away from 25C


# other stuff -------------------------------------------------------------


View(res0b)

str(res25b)
class(res25b$diff_from25)


res0b_warming <- res0b %>% 
	filter(stand_temp > 0)

res0b_cooling <- res0b %>% 
	filter(stand_temp < 0)

ggplot() +
	geom_path(data = res0b_warming, aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(res0b_warming$stabil_potential)*0.99, max(res0b_warming$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 
ggsave("figures/ref0-graph-from25-warming.png", width = 8, height = 6)



ggplot() +
	geom_path(data = res0b_warming, aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(res0b_warming$stabil_potential)*0.99, max(res0b_warming$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + 
	scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 



### ok now try the fitness difference * rho as the distance


res25c <- res25b %>%
	mutate(distance = rho * fit_ratio) %>% 
	group_by(T) %>% 
	mutate(mean_distance= mean(distance))

res25c_sum <- res25c %>% 
	group_by(T) %>% 
	summarise(distance = mean(distance),
			  mean_stabil = mean(stabil_potential),
			  mean_fit = mean(fit_ratio))


res25c %>% 
	ggplot() + geom_line(aes(x = T, y = distance, group = iteration), size = 1) +
	geom_line(aes(x = T, y = mean_distance), color = "red") +ylab("rho * fitness ratio")
ggsave("figures/distances.png", width = 8, height = 6)


res25c %>% 
	ggplot(aes(x = T, y = distance)) + geom_point() +
		   	ylab("rho * fitness ratio") + geom_smooth(method = "lm")
ggsave("figures/distances-slope.png", width = 8, height = 6)


res25c %>% 
	ggplot() + geom_line(aes(x = T, y = stabil_potential, group = iteration), size = 1, alpha = 0.7) +
	geom_line(aes(x = T, y = mean_stabil), color = "red", data = res25c_sum) +ylab("Stabilization potential")
ggsave("figures/stabils-ref25.png", width = 8, height = 6)
	


View(res25c)


plot1 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = stabil_potential, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")
ggsave("figures/stabils-ref25-facet10.png", width = 14, height = 2)


res25c %>% 
	filter(iteration %in% c(2)) %>% 
	ggplot(aes( x= T, y = stabil_potential, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")


plot2 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = a11, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

plot2b <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = a12, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")
plot2c <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = a21, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

plot2d <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = a22, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

plot3 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = g1, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

plot4 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = c1P, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

plot4b <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = c1N, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")


plot5 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = beta11, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

plot5b <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = beta12, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")
plot5c <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = beta21, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")
plot5d <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = beta22, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

plot6 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = KN, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")
plot7 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = KP, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

plot8 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = a12*a21, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")
plot9 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = a11*a22, group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")
plot10 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = (a12*a21)/(a11*a22), group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")
plot10b <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot() +
	geom_line(aes( x= T, y = (a12*a21), group = iteration)) +
	geom_line(aes( x= T, y = (a11*a22), group = iteration)) +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

res25c %>% 
	filter(iteration %in% c(2)) %>% 
	ggplot() +
	geom_line(aes( x= T, y = (a12*a21), group = iteration)) +
	geom_line(aes( x= T, y = (a11*a22), group = iteration)) +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")

### ok let's report out the Eas and b's that caused that

res25c %>% 
	filter(iteration %in% c(2)) %>% 
	dplyr::select(1:17, iteration) %>% 
	dplyr::select(iteration, 1:17) %>% 
	ungroup() %>% 
	dplyr::select(-T) %>% 
	gather(key = parameter, value = value, 2:17) %>%
	distinct() %>% 
	View


plot11 <- res25c %>% 
	filter(iteration %in% c(1:10)) %>% 
	ggplot(aes( x= T, y = 1- sqrt((a12*a21)/(a11*a22)), group = iteration)) + geom_line() +
	facet_grid(~ iteration)  + geom_hline(yintercept = 0, color = "red")



library(patchwork)

plots <- plot1 / plot2 /plot2b/plot2c/plot2d/ plot3 /plot4/plot4b /plot5/plot5b/ plot5c/ plot5d/plot6/plot7/plot8/plot9/plot10/plot10b/plot11
ggsave('figures/stabils-ref25-facet10-params.png', plots, width = 14, height = 36)
ggsave('figures/stabils-ref25-facet10-params-b.png', plots, width = 14, height = 36)
ggsave('figures/stabils-ref25-facet10-params-c.png', plots, width = 14, height = 36)



### Make a figure that’s like a phase space diagram that outlines when we should expect to see non-monotonic relationships in stabilization potential, add to this figure the toy scenario 

### first pull out all the non-monotonic iterations
res25d <- res25c %>% 
	group_by(iteration) %>%
	filter(stabil_potential == max(stabil_potential)) %>% 
	mutate(monotonic = ifelse(T %in% c(0, 50), "monotonic", "nonmonotonic")) ### ok this makes sense

res25d %>% View

plot_mon <- res25d %>% 
	dplyr::select(2:17, monotonic) %>%
	gather(2:17, key = parameter, value = value) %>%
	filter(!str_detect(parameter, "b")) %>%
	group_by(monotonic, parameter) %>% 
	mutate(mean_param_value = mean(value)) %>% 
	ungroup() %>%
	ggplot(aes(x = value, fill = monotonic)) + geom_density() +
	geom_vline(aes(xintercept = mean_param_value, color = monotonic)) +
	facet_wrap( ~ parameter,scales = "free")
ggsave(plot = plot_mon, "figures/params-monotonic.png", width = 8, height = 8)

res25d2 <- res25d %>% 
	dplyr::select(iteration, monotonic)
res25c2 <- full_join(res25c, res25d2)
res25c3 <- res25c2 %>% 
	group_by(T) %>% 
	summarise(mean_stabil = mean(stabil_potential))
View(res25c2)
View(res25d)


	ggplot() + 
	geom_line(aes(x = T, y = stabil_potential, group = iteration, color = monotonic), data = res25c2) +
	geom_line(aes(x = T, y = mean_stabil), data = res25c3, color = "black")
	ggsave("figures/stabil-pot-summary.png", width = 10, height = 8)
	
	res25c2 %>% 
		ungroup() %>% 
		filter(iteration %in% c(1: 10)) %>% 
		dplyr::select(monotonic, iteration, 2:17) %>% 
		gather(3:18, key = parameter, value = value) %>% 
		distinct() %>% 
		ggplot(aes(x = parameter, y = value, fill = monotonic)) + geom_violin()
	ggsave("figures/stabil-pot-summary-params-b.png", width = 12, height = 4)
	ggsave("figures/stabil-pot-summary-params-bc-violin.png", width = 12, height = 4)
	
	res25c2 %>% 
		ungroup() %>% 
		filter(iteration %in% c(1: 10)) %>% 
		ggplot() +
		geom_line(aes(x = T, y = stabil_potential, group = iteration, color = monotonic), data = filter(res25c2, iteration %in% c(1: 10)))
	ggsave("figures/stabil-pot-summary-params-bc-lines.png", width = 12, height = 4)
	
	res25c2 %>% 
		ungroup() %>% 
		filter(iteration %in% c(1:10)) %>% 
		dplyr::select(monotonic, iteration, 2:17) %>% 
		gather(3:18, key = parameter, value = value) %>% 
		distinct() %>% 
		ggplot(aes(x = parameter, y = value, color = monotonic)) + geom_point()
	ggsave("figures/stabil-pot-summary-params-bc.png", width = 12, height = 4)
	
	
	#### ok let's try with focusing on the really high r Ea values
	
	res25c2 %>%
		dplyr::select(r_EaN, r_EaP, monotonic, stabil_potential) %>% 
		gather(r_EaN, r_EaP, stabil_potential, key = parameter, value = value) %>%
		ggplot() +
		geom_line(aes(x = T, y = value, color = monotonic), alpha = 0.2) +
		facet_grid( ~ parameter)
	ggsave("figures/stabil-pot-summary-params-facet.png", width = 12, height = 4)
	
	res25c2 %>%
		dplyr::select(r_EaN, r_EaP, monotonic, iteration) %>% 
		gather(r_EaN, r_EaP, key = parameter, value = value) %>%
		ggplot() +
		geom_line(aes(x = parameter, y = value, color = monotonic, group = iteration)) 
	ggsave("figures/stabil-pot-summary-params-rparams.png", width = 8, height = 8)
	
	res25c2 %>%
		dplyr::select(c_Ea1N, c_Ea2N, monotonic, iteration) %>% 
		gather(c_Ea1N, c_Ea2N, key = parameter, value = value) %>%
		ggplot() +
		geom_line(aes(x = parameter, y = value, color = monotonic, group = iteration)) 
	ggsave("figures/stabil-pot-summary-params-cparams.png", width = 8, height = 8)
	
	
	### ok now let's take the monotonic ones, and split them out by increasing or decreasing and see if we can figure out what explains the difference between increasing and decreasing
	
	res25c3 <- res25c2 %>%
		filter(monotonic == "monotonic") %>% 
		group_by(iteration) %>% 
		filter(stabil_potential == max(stabil_potential)) %>% 
		mutate(increasing = ifelse(T == 50, "increasing", "decreasing")) %>% 
		dplyr::select(increasing, monotonic, iteration) %>% 
		ungroup() %>% 
		distinct()
	
	
	monotonic <- res25c2 %>%
		filter(monotonic == "monotonic") 
	res25c4 <- left_join(monotonic, res25c3)
	
	View(res25c3)
	
	res25c4 %>% 
		ggplot(aes(x = T, y = stabil_potential, color = increasing, group = iteration)) + geom_line()
	ggsave("figures/decreasing-increasing-stabil.png", width = 8, height = 6)
	
	library(plotrix)
	mean_res25c4 <- res25c4 %>% 
		dplyr::select(iteration, increasing, 2:17) %>% 
		gather(4:19, key = parameter, value = value) %>% 
		group_by(increasing, parameter) %>% 
		summarise(mean_parameter_value = mean(value),
					   std_err = std.error(value)) %>% 
		filter(!str_detect(parameter, "b"))
	
	
	mean_res25c4 %>% 
		ggplot() +
		geom_pointrange(aes(color = increasing, x = parameter, y = mean_parameter_value, ymin = mean_parameter_value - std_err, ymax = mean_parameter_value + std_err)) +
		facet_wrap( ~ parameter, scales = "free") 
	ggsave("figures/params-increasing-decreasing-facet.png", width = 10, height = 8)
	
	res25c4 %>% 
		dplyr::select(iteration, increasing, 2:17) %>% 
		gather(4:19, key = parameter, value = value) %>% 
		filter(!str_detect(parameter, "b")) %>%
		ggplot(aes(x = value, fill = increasing)) + geom_density(alpha = 0.5) +
		geom_vline(aes(xintercept = mean_parameter_value, color = increasing), data = mean_res25c4) +
		facet_wrap( ~ parameter, scales = "free")
	ggsave("figures/decreasing-increasing-stabil-density.png", width = 12, height = 10)

	
#### OK now let's make the pom pom but highlighting the increasing and decreasing cases
View(res25c4)	
#	reformulate as increasing or decreasing relative to the reference temp. This would allow you to have all the red lines going to the left and all the blue to the right.

write_csv(res25c4, "data-processed/res25c4.csv")
	

ggplot() +
		geom_path(data = filter(res25c4, T > 25), aes(x = stabil_potential, y = fit_ratio, color = increasing, group = iteration), size = 2) +
		geom_ribbon(data = data.frame(x = seq(min(res25c4$stabil_potential)*0.99, max(res25c4$stabil_potential)*1.01, 0.001)),
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
	 ggsave("figures/decreasing-increasing-pompom.png", width = 12, height = 10)
	 ggsave("figures/decreasing-increasing-pompom-higher-than-ref.png", width = 12, height = 10)
	 
### ok this is with the empirical parameters all at the same time	
	  ggplot() +
	 	geom_path(data = filter(res25c4), aes(x = stabil_potential, y = fit_ratio, color = increasing, group = iteration), size = 2) +
	 	geom_ribbon(data = data.frame(x = seq(min(res25c4$stabil_potential)*0.99, max(res25c4$stabil_potential)*1.01, 0.001)),
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
	 	ylim(-0, 3) + xlim(min(res25c4$stabil_potential)*0.99, max(res25c4$stabil_potential)*1.01)
	 ggsave("figures/decreasing-increasing-pompom-all.png", width = 12, height = 10)
	 

# start again here --------------------------------------------------------

	 
	res25b <- res25 %>% 
		mutate(diff_from_25 = T-25) 
	 
	 ### ok this is with the empirical parameters all at the same time	(May 20 2024 afternoon)
	 ggplot() +
	 	geom_path(data = res25b, aes(x = stabil_potential, y = fit_ratio, group = iteration, color = diff_from_25), size = 2) +
	 	geom_point(data = filter(res25b, diff_from_25 == 0), aes(x = stabil_potential, y = fit_ratio), size = 4, color = "green") +
	 	# geom_point(data = filter(res25, T == 0), aes(x = stabil_potential, y = fit_ratio), size = 2, color = "red") +
	 	geom_ribbon(data = data.frame(x = seq(min(res25b$stabil_potential)*0.99, max(res25b$stabil_potential)*1.01, 0.001)),
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
	 	scale_color_distiller(palette = 'RdBu')
	 ggsave("figures/decreasing-increasing-pompom-all-afternoon.png", width = 12, height = 10)
	 

## find the average end point for the hot and cool	 
	 
View(res25b)	 

### ok what if I try Po-Ju's new definition of niche and fitness differences
#niche difference as -log(\rho) and fitness difference as log(f1/f2).
# fit_ratio <- sqrt((a11*a12)/(a22*a21))  #fitness ratio
	 
r2 <- res25b %>% 
	group_by(T) %>% 
	# mutate(fit_ratio = log(fit_ratio)) %>% 
	# mutate(stabil_potential = -log(rho)) %>% 
	summarise(mean_fit = mean(fit_ratio),
			  mean_stabil = mean(stabil_potential),
			  median_fit = median(fit_ratio),
			  median_stabil = median(stabil_potential),
			  max_fit = max(fit_ratio),
			  max_stabil = max(stabil_potential),
			  min_fit = min(fit_ratio),
			  min_stabil = min(stabil_potential))

r3 <- res25b %>% 
	# group_by(T) %>% 
	# mutate(fit_ratio = log(fit_ratio)) %>% 
	# mutate(stabil_potential = -log(rho)) %>% 
	summarise(mean_fit = mean(fit_ratio),
			  mean_stabil = mean(stabil_potential),
			  median_fit = median(fit_ratio),
			  median_stabil = median(stabil_potential),
			  max_fit = max(fit_ratio),
			  max_stabil = max(stabil_potential),
			  min_fit = min(fit_ratio),
			  min_stabil = min(stabil_potential))

	 

r2 %>% 
	ggplot(aes(x = T, y = mean_fit)) + geom_point()

r2 %>% 
	ggplot(aes(x = T, y = mean_stabil)) + geom_point()


ggplot() +
	geom_path(data = res25b, aes(x = stabil_potential, y = fit_ratio, group = iteration, color = diff_from_25), size = 2) +
	geom_point(data = filter(res25b, diff_from_25 == 0), aes(x = stabil_potential, y = fit_ratio), size = 4, color = "green") +
	# geom_point(data = filter(res25, T == 0), aes(x = stabil_potential, y = fit_ratio), size = 2, color = "red") +
	geom_ribbon(data = data.frame(x = seq(min(res25b$stabil_potential)*0.99, max(res25b$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = filter(r2, T == 0), aes(x = mean_stabil, y = mean_fit), color = "blue", size = 5) +
	# geom_point(data = filter(r2, T == 50), aes(x = mean_stabil, y = mean_fit), color = "red", size = 5) +
	# geom_point(data = filter(r2, T == 25), aes(x = mean_stabil, y = mean_fit), color = "black", size = 5) +
	geom_point(data = filter(r2, T == 0), aes(x = median_stabil, y = median_fit), color = "blue", size = 5) +
	geom_point(data = filter(r2, T == 50), aes(x = median_stabil, y = median_fit), color = "red", size = 5) +
	geom_point(data = filter(r2, T == 25), aes(x = median_stabil, y = median_fit), color = "black", size = 5) +
	geom_point(data = filter(r3), aes(x = min_stabil, y = min_fit), color = "lightblue", size = 5) +
	geom_point(data = filter(r3), aes(x = min_stabil, y = min_fit), color = "pink", size = 5) +
	# geom_point(data = filter(r2, T == 25), aes(x = median_stabil, y = median_fit), color = "black", size = 5) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
	scale_color_distiller(palette = 'RdBu') +
	theme(axis.text=element_text(size=16),
		  axis.title=element_text(size=16,face="bold"))
ggsave("figures/decreasing-increasing-pompom-all-afternoon-averages.png", width = 12, height = 10)
ggsave("figures/decreasing-increasing-pompom-all-afternoon-median.png", width = 12, height = 10)
ggsave("figures/decreasing-increasing-pompom-all-afternoon-median-min.png", width = 12, height = 10)





ggplot() +
	geom_path(data = r2, aes(x = mean_stabil, y = mean_fit), size = 1) +
	# geom_point(data = filter(res25b, diff_from_25 == 0), aes(x = stabil_potential, y = fit_ratio), size = 4, color = "green") +
	# geom_point(data = filter(res25, T == 0), aes(x = stabil_potential, y = fit_ratio), size = 2, color = "red") +
	geom_ribbon(data = data.frame(x = seq(min(r2$mean_stabil)*0.99, max(r2$mean_stabil)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	geom_point(data = filter(r2, T == 0), aes(x = mean_stabil, y = mean_fit), color = "blue", size = 5) +
	geom_point(data = filter(r2, T == 50), aes(x = mean_stabil, y = mean_fit), color = "red", size = 5) +
	geom_point(data = filter(r2, T == 25), aes(x = mean_stabil, y = mean_fit), color = "black", size = 5) +
	geom_hline(yintercept = 1, linetype=5) + 
	# scale_color_continuous_diverging() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) + 
	scale_color_distiller(palette = 'RdBu')



	 library(viridis)
	 ggplot() +
	 	geom_path(data = filter(res25b, diff_from_25 < 0), aes(x = stabil_potential, y = fit_ratio, group = iteration, color = diff_from_25), size = 2) +
	 	geom_point(data = filter(res25b, diff_from_25 == 0), aes(x = stabil_potential, y = fit_ratio), size = 4, color = "green") +
	 	# geom_point(data = filter(res25, T == 0), aes(x = stabil_potential, y = fit_ratio), size = 2, color = "red") +
	 	geom_ribbon(data = data.frame(x = seq(min(res25b$stabil_potential)*0.99, max(res25b$stabil_potential)*1.01, 0.001)),
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
	 	scale_color_viridis(begin = 0, end = 1)
	  ggsave("figures/decreasing-increasing-pompom-all-afternoon-lessthan25.png", width = 12, height = 10)	
	  
	  ggplot() +
	  	geom_path(data = filter(res25b, diff_from_25 > 0), aes(x = stabil_potential, y = fit_ratio, group = iteration, color = diff_from_25), size = 2) +
	  	geom_point(data = filter(res25b, diff_from_25 == 0), aes(x = stabil_potential, y = fit_ratio), size = 4, color = "green") +
	  	# geom_point(data = filter(res25, T == 0), aes(x = stabil_potential, y = fit_ratio), size = 2, color = "red") +
	  	geom_ribbon(data = data.frame(x = seq(min(res25b$stabil_potential)*0.99, max(res25b$stabil_potential)*1.01, 0.001)),
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
	  	scale_color_viridis(begin = 0, end = 1)
	  ggsave("figures/decreasing-increasing-pompom-all-afternoon-morethan25.png", width = 12, height = 10)
	  ggsave("figures/decreasing-increasing-pompom-all-afternoon-morethan25-all.png", width = 12, height = 10)
	 
	 
	  warming <- res25b %>% 
	  	filter(diff_from_25 > 0) 
	  
	  w2 <- warming %>% 
	  	summarise(median_fit = median(fit_ratio),
	  			  median_stabl = median(stabil_potential))
	  w3 <- warming %>% 
	  	filter(diff_from_25 == 25) %>% 
	  	summarise(median_fit = mean(fit_ratio),
	  			  median_stabl = mean(stabil_potential))
	  w4 <- warming %>% 
	  	filter(diff_from_25 == min(diff_from_25)) %>% 
	  	summarise(median_fit = mean(fit_ratio),
	  			  median_stabl = mean(stabil_potential))
	  
	  
	  ggplot() +
	  	geom_path(data = filter(res25b, diff_from_25 > 0), aes(x = stabil_potential, y = fit_ratio, group = iteration, color = diff_from_25), size = 2) +
	  	geom_point(data = filter(res25b, diff_from_25 == 0), aes(x = stabil_potential, y = fit_ratio), size = 4, color = "green") +
	  	geom_point(data = w4, aes(x = median_stabl, y = median_fit), size = 4, color = "blue") +
	  	geom_point(data = w3, aes(x = median_stabl, y = median_fit), size = 4, color = "red") +
	  	# geom_point(data = filter(res25, T == 0), aes(x = stabil_potential, y = fit_ratio), size = 2, color = "red") +
	  	geom_ribbon(data = data.frame(x = seq(min(res25b$stabil_potential)*0.99, max(res25b$stabil_potential)*1.01, 0.001)),
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
	  	scale_color_viridis(begin = 0, end = 1) +
	  	theme(axis.text=element_text(size=24),
	  		  axis.title=element_text(size=24,face="bold"))
	  ggsave("figures/decreasing-increasing-pompom-all-afternoon-morethan25-all-averages.png", width = 12, height = 10)
	  
	  
	w5 <- warming %>% 
		group_by(coexist) %>% 
		tally()
	
	w5 %>% 
		ggplot(aes(x = coexist, y = n)) + geom_bar(stat = "identity")
	
	((125500 - 124500) / (125500 + 124500))*100
	
	
	
	coex <- warming %>% 
		filter(coexist == "TRUE")
	
	exclude <- warming %>% 
		filter(coexist == "FALSE")
	
	
	exclude %>% 
		ggplot(aes(x = c_Ea1N, y = stabil_potential)) + geom_point() +
		geom_smooth(method = "lm")
	
	exclude %>% 
		ggplot(aes(x = c_Ea1P, y = stabil_potential)) + geom_point() +
		geom_smooth(method = "lm")
	
	exclude %>% 
		ggplot(aes(x = K_EaN, y = stabil_potential)) + geom_point() +
		geom_smooth(method = "lm")
	
	exclude %>% 
		ggplot(aes(x = K_EaP, y = stabil_potential)) + geom_point() +
		geom_smooth(method = "lm")
	
	 
	 res_summary <- res25  %>% 
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
	 
	 
	 
	 
	 ggplot() +
	 	geom_path(data = res_summary, aes(x = mean_fit_s, y = mean_fit, group = iteration, color = T), size = 2) +
	 	geom_point(data = filter(res_summary, T == 25), aes(x = mean_fit_s, y = mean_fit), size = 2, color = "pink") +
	 	geom_ribbon(data = data.frame(x = seq(min(res_summary$mean_fit_s)*0.99, max(res_summary$mean_fit_s)*1.01, 0.001)),
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
	 ggsave("figures/decreasing-increasing-pompom-all.png", width = 12, height = 10)	 
	 
	 
	 
	 ggplot() +
	 	geom_path(data = filter(res25c4, T <= 40), aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 2) +
	 	geom_point(data = filter(res25c4, T == 25), aes(x = stabil_potential, y = fit_ratio), size = 2, color = "pink") +
	 	geom_ribbon(data = data.frame(x = seq(min(res25c4$stabil_potential)*0.99, max(res25c4$stabil_potential)*1.01, 0.001)),
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
	 	ylim(-0, 3) + xlim(min(res25c4$stabil_potential)*0.99, max(res25c4$stabil_potential)*1.01)
	 ggsave("figures/decreasing-increasing-pompom-temperatures.png", width = 12, height = 10)
	 
	 
### Now we want to do this with the empirical data	 
	 
	
	 
### Now we want to know what combinations of parameters cause the increase or decrease
	 
	### start at the highest level equations and go down from there
	 res25c4 %>% 
	 	filter(T > 25) %>%
	 	mutate(numerator = a12*a21) %>% 
	 	mutate(denominator = a11*a22) %>% 
	 	ggplot(aes(x = T, y = stabil_potential, color = increasing)) + geom_point(size = .75)
	 
	 res25c4 %>% 
	 	filter(T > 25) %>%
	 	mutate(numerator = a12*a21) %>% 
	 mutate(denominator = a11*a22) %>% 
	 	ggplot(aes(x = T, y = denominator, color = increasing)) + geom_point(size = .75) ### all denominators are decreasing
	
	 res25c4 %>% 
		filter(T > 25) %>%
		mutate(numerator = a12*a21) %>% 
		mutate(denominator = a11*a22) %>% 
		ggplot(aes(x = T, y = numerator, color = increasing)) + geom_point(size = .75)### all numerators are decreasing
	
	res25c4 %>% 
		filter(T > 25) %>%
		mutate(numerator = a12*a21) %>% 
		mutate(denominator = a11*a22) %>% 
		ggplot(aes(x = T, y = a22, color = increasing)) + geom_point(size = .75) ## one is going up then down
	
	res25c4 %>% 
		filter(T > 25) %>%
		mutate(numerator = a12*a21) %>% 
		mutate(denominator = a11*a22) %>% 
		ggplot(aes(x = T, y = a12, color = increasing)) + geom_point(size = .75) ## all decreasing
	
	res25c4 %>% 
		filter(T > 25) %>%
		mutate(numerator = a12*a21) %>% 
		mutate(denominator = a11*a22) %>% 
		ggplot(aes(x = T, y = a11, color = increasing)) + geom_point(size = .75) ## all decreasing
	res25c4 %>% 
		filter(T > 25) %>%
		mutate(numerator = a12*a21) %>% 
		mutate(denominator = a11*a22) %>% 
		ggplot(aes(x = T, y = a21, color = increasing)) + geom_point(size = .75) ## one increasing
	res25c4 %>% 
		filter(T > 25) %>%
		mutate(numerator = a12*a21) %>% 
		mutate(denominator = a11*a22) %>% 
		ggplot(aes(x = T, y = a22, color = increasing)) + geom_point(size = .75) ## one increasing
	
	res25c4 %>% 
		filter(T > 25) %>%
		ggplot(aes(x = T, y = g1, color = increasing)) + geom_point(size = .75)  ## some increaseing, some decreasing
	res25c4 %>% 
		filter(T > 25) %>%
		ggplot(aes(x = T, y = g2, color = increasing)) + geom_point(size = .75) ## some increasing, some decreasing
	res25c4 %>% 
		filter(T > 25) %>%
		ggplot(aes(x = T, y = beta11, color = increasing)) + geom_point(size = .75) 
	res25c4 %>% 
		filter(T > 25) %>%
		ggplot(aes(x = T, y = beta21, color = increasing)) + geom_point(size = .75) 
	res25c4 %>% 
		filter(T > 25) %>%
		ggplot(aes(x = T, y = beta22, color = increasing)) + geom_point(size = .75) 
	res25c4 %>% 
		filter(T > 25) %>%
		ggplot(aes(x = T, y = beta12, color = increasing)) + geom_point(size = .75) 
	
	 		   	
	 
# rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
	 
	 
	 
	 

## Note from Patrick: 
#	reformulate as increasing or decreasing relative to the reference temp. This would allow you to have all the red lines going to the left and all the blue to the right.
	#Then attempt to distinguish the minimum set of parameters that allows us to determine whether the lines decrease or increase. 
	#I recommend starting at the highest level equations and then moving down level by level as far as it still makes sense.
	
	
# Suggestions for next time:
# 	Make a figure that’s like a phase space diagram that outlines when we should expect to see non-monotonic relationships in stabilization potential, add to this figure the toy scenario 
# 
# ​​what does it take to break the flat relationship? (i.e. where are we likely to see non-random combinations of species)
# 
# Based on the scenario that we outline, subset the random trait space for parameter values that align with that scenario and plot on the pom pom (this is likely something like warm adapted species has higher temp dependence than cooler adapted species - subset all values in pom pom that fit that and hopefully those all show a similar pattern in the pom pom). 
# 
# Start with density independent parameters, like growth rate – see if that makes sense on the pom-pom then if not… point out the need for looking at density independent.
# 
# Think about other possible scenario – inspired by Amy Angert’s mimulus work on the relationship between thermal niche and range size

### ok come back to the target vs. reference species stuff -- think that's going to be promising
	
	
### read in the empirical data
library(janitor)
library(cowplot)
theme_set(theme_cowplot())	
	
data <- read_csv("data/mac-means.csv") %>% 
	clean_names()

View(data)

data %>% 
	ggplot(aes(x = activation_energy)) + geom_histogram() + 
	facet_grid(~simple_parameter)
