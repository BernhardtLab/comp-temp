
library(fields)
library(tidyverse)
# res_ref0 <- results_b
res_ref25 <- results_b25

# write_csv(res_ref0, "data-processed/sim-results-reftemp0.csv")
write_csv(res_ref25, "data-processed/sim-results-reftemp25.csv")


res_ref0 <- read_csv("data-processed/sim-results-reftemp0.csv")
res_ref25 <- read_csv("data-processed/sim-results-reftemp25.csv")

ggplot() +
	geom_path(data = res_ref25, aes(x = stabil_potential, y = fit_ratio, color = T, group = iteration), size = 2) +
	geom_ribbon(data = data.frame(x = seq(min(results_b$stabil_potential)*0.99, max(results_b$stabil_potential)*1.01, 0.001)),
				aes(x = x,
					y = NULL,
					ymin = 1-x,
					ymax = 1/(1-x)),
				fill = "grey", color = "black", alpha = 0.2) +
	# geom_point(data = res20, aes(x = stabil_potential, y = fit_ratio), color = "pink", size = 4, shape = 1) +
	geom_hline(yintercept = 1, linetype=5) + scale_color_viridis_c() +
	xlab(expression(paste("Stabilization potential (1-", rho, ")"))) +
	ylab(expression(paste("Fitness difference (", f[2], "/", f[1], ")"))) 



res25_ref <- res_ref25 %>% 
	filter(T == 25)

res0_ref <- res_ref0 %>% 
	filter(T == 0)


#### for the ref 25

resmat25 <- res_ref25 %>% 
	dplyr::select(fit_ratio, stabil_potential) %>% 
	as.matrix()

res25_mat_ref <- res_ref25 %>% 
	filter(T == 25) %>% 
	dplyr::select(fit_ratio, stabil_potential) %>% 
	distinct() %>% 
	as.matrix()


distances25 <-  rdist(res25_mat_ref, resmat25) %>% 
	as.data.frame() %>% 
	gather()

all_distances25 <- bind_cols(res_ref25, distances25) %>% 
	mutate(ref_temp = 25) %>% 
	mutate(std_temp = T - 25)

all_distances25 %>% 
	ggplot(aes(x = T, y = value)) + geom_point(alpha = 0.1)
all_distances25 %>% 
	ggplot(aes(x = std_temp, y = value)) + geom_point(alpha = 0.1)

all_distances25 %>% 
	ggplot(aes(x = std_temp, y = value, color = coexist)) + geom_point(alpha = 0.1)+
	facet_wrap( ~ iteration, scales = "free")
ggsave("figures/distances-ref25-all-coexist-facs-free.png", width = 20, height = 20)


all_distances25 %>% 
	ggplot(aes(x = T, y = value, color = coexist)) + geom_point(alpha = 0.1)+
	facet_wrap( ~ iteration, scales = "free")
ggsave("figures/distances-ref25-all-coexist-facs-free-from25.png", width = 20, height = 20)


all_distances25 %>% 
	mutate(direction = ifelse(std_temp > 0, "warming", "cooling")) %>% 
	ggplot(aes(x = value)) + geom_density(fill = "blue") +
	facet_wrap( ~ direction, ncol = 1)
ggsave("figures/distances-density.png", width = 6, height = 4)


# for the ref 0 -----------------------------------------------------------


resmat0 <- res_ref0 %>% 
	dplyr::select(fit_ratio, stabil_potential) %>% 
	as.matrix()

res0_mat_ref <- res_ref0 %>% 
	filter(T == 0) %>% 
	dplyr::select(fit_ratio, stabil_potential) %>% 
	distinct() %>% 
	as.matrix()


distances0 <-  rdist(res0_mat_ref, resmat0) %>% 
	as.data.frame() %>% 
	gather()

all_distances0 <- bind_cols(res_ref0, distances0) %>% 
	mutate(ref_temp = 0) %>% 
	mutate(std_temp = T-0)

all_distances0 %>% 
	ggplot(aes(x = T, y = value)) + geom_point(alpha = 0.1)

all_distances0 %>% 
	ggplot(aes(x = std_temp, y = value, color = coexist)) + geom_point(alpha = 0.1)+
	facet_wrap( ~ iteration, scales = "free")
ggsave("figures/distances-reft0-all-coexist-facs-free.png", width = 20, height = 20)







