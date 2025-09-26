


### how does the supply point impact the fitness differences?

a11 <- function(SP) c1P / (D * (SP - R1P))
a12 <- c2P / (D * (SP - R1P))
a21 <- c1N / (D * (SN - R2N))
a22 <- c2N / (D * (SN - R2N))


rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
stabil_potential <- 1 - rho #stabilizing potential
fit_ratio <- sqrt((a11*a12)/(a22*a21)) #fitness ratio 
coexist <- rho < fit_ratio &  fit_ratio < 1/rho


a11 <- sapply(10:1000, function(SP) c1P / (D * (SP - R1P)))
a12 <- sapply(10:1000, function(SP) c2P / (D * (SP - R1P)))
a21_consN <- sapply(rep(20, times = length(10:1000)), function(SN) c1N / (D * (SN - R2N)))
a22_consN <- sapply(rep(20, times = length(10:1000)), function(SN) c2N / (D * (SN - R2N)))
a21_varN <- sapply(10:1000, function(SN) c1N / (D * (SN - R2N)))
a22_varN <- sapply(10:1000, function(SN) c2N / (D * (SN - R2N)))



df <- data.frame(a11 = a11, a12 = a12, a21 = a21, a22 = a22,
				 SN = rep(20, times = length(10:1000)), SP = 10:1000, SN = 10:1000, a21_varN = a21_varN, a22_varN, a22_varN) %>% 
	mutate(fitness_ratio = sqrt((a11*a12)/(a22*a21))) %>% 
	mutate(fitness_ratio_varN = sqrt((a11*a12)/(a22_varN*a21_varN)))

df %>% 
	ggplot(aes(x = SP, y = fitness_ratio_varN)) + geom_point() +
	geom_point(aes(x = SP, y = fitness_ratio), color = "pink")
