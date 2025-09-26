

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())


arrhenius<-function(Temp, E, b1, ref_temp = 0) {
  k<-8.62e-05 #Boltzmann's constant
  E <- E # 0.6 # activation energy (eV)
  T<-Temp+273.15 #range of temp in K
  Tc <- ref_temp+273.15 #reference temperature
  
  metabolism<-(b1*exp(1)^(E*(1/(k*Tc)-1/(k*T))))
  return(metabolism)
}

b <- arrhenius(Temp = 0:40, E = 0.6, b1 = 0.5, ref_temp = 1)
d <- 0.25 + arrhenius(Temp = 0:40, E = 0.75, b1 = 0.25, ref_temp = )
plot(b-d, ylim = c(0,max(b-d)), type = "l")

b2 <- arrhenius(Temp = 0:40, E = 0.6, b1 = 0.5, ref_temp = -2)
d2 <- 0.25 + arrhenius(Temp = 0:40, E = 0.75, b1 = 0.25, ref_temp = -2)
lines(b2-d2, col = 2)

library(tidyverse)
library(cowplot)

#' Setup Tilman's consumer resource model parameters 

### Here species 1 is RED species 2 is BLUE
temp_RC <- function(T = 25, r_eA = 0.6, c_eA = 0.6, w_eA = 0, k_eA = 0, ref_temp2 = 0){
  D = 0.25 #0.7 #mortality rate
  c11 = arrhenius(Temp = T, E = c_eA, b1 = 0.2); c12 = arrhenius(Temp = T, E = c_eA, b1 = 0.4); c21 = arrhenius(Temp = T, E = c_eA, b1 = 0.5); c22 = arrhenius(Temp = T, E = c_eA, b1 = 0.2) ### cij = per capita consumption of comsumer i on resource j
  w11 = arrhenius(Temp = T, E = w_eA, b1 = 0.2); w12 = arrhenius(Temp = T, E = w_eA, b1 = 0.4); w21 = arrhenius(Temp = T, E = w_eA, b1 = 0.4); w22 = arrhenius(Temp = T, E = w_eA, b1 = 0.2) ## wij = weighting factor that converts availability of resource j into consumer Ni
  k1 = arrhenius(Temp = T, E = k_eA, b1 = 0.1); k2 = arrhenius(Temp = T, E = k_eA, b1 = 0.1) #half saturation constant for N resource consumption
  r1 = arrhenius(Temp = T, E = r_eA, b1 = 0.5) - arrhenius(Temp = T, E = 0.74, b1 = 0.25); 
  r2 = arrhenius(Temp = T, E = r_eA, b1 = 0.5, ref_temp = ref_temp2) - arrhenius(Temp = T, E = 0.74, b1 = 0.25, ref_temp = ref_temp2) #population growth rates
  T1 = 1; T2 = 1 #minimum amount of resource required for growth
  
  species1_consumption <- c12/c11
  species2_consumption <- c22/c21
  
  #' Calculate slope and intercept for the ZNGI's (get this from equation S3 in Ke and Letten)
  y.inter.1 <-  (D * (k1 - T1) + r1 * T1) / (w12 * (r1 - D))
  y.inter.2 <-  (w21 / w22) * (D * (k2 - T2) + r2 * T2) / (w21 * (r2 - D))
  
  
  #' Calculate equilibrium resource values (when the ZNGIs cross each other, i.e. such that Ke and Letten equations S3.1 = S3.2)
  B1 <-  y.inter.1
  B2 <-  y.inter.2
  Lamda.1 <-  (w11 / w12)
  Lamda.2 <- (w21 / w22)
  R1.star <-  (B1 - B2) / (Lamda.1 - Lamda.2)
  R2.star <-  (B2 * Lamda.1 - B1 * Lamda.2) / (Lamda.1 - Lamda.2)  
  
  
  #' Supply concentrations of resources
  S1 <- 10 #0.248
  S2 <- 12  #0.245
  
  
  #Calculate alphas
  a11 <- (c12 + c11 * (w11/w12)) / (D * (S2 + (w11/w12) * S1 - B1))
  a12 <- (c22 + c21 * (w11/w12)) / (D * (S2 + (w11/w12) * S1 - B1))
  
  a22 <- (c22 + c21 * (w21/w22)) / (D *(S2 + (w21/w22) * S1 - B2))
  a21 <- (c12 + c11 * (w21/w22)) / (D *(S2 + (w21/w22) * S1 - B2))
  
  p <- sqrt((a12*a21)/(a11*a22)) #niche overlap
  stabil_potential <- 1 - p #stabilizing potential
  fit_ratio <- sqrt((a11*a12)/(a22*a21)) * (r2)/(r1) #fitness ratio
  coexist <- p < fit_ratio &  fit_ratio < 1/p
  
  data.frame(a11 = a11, a12 = a12, a22 = a22, a21 = a21, R1.star = R1.star, R2.star= R2.star, K1 = (r1)/a11, K2 = (r2)/a22, T = T, r1 = r1, r2 = r2, c = c11, stabil_potential = stabil_potential, fit_ratio = fit_ratio, p = p , coexist = coexist)
}

results<-data.frame()
for(i in 1:20){
  for(j in 1:5){
  r_eA <- 0.6
  c_eA <- 0.1*(i-1)
  w_eA <- 0.1*(j-1)
  hold<-temp_RC(T = seq(0, 35,length = 100), r_eA = r_eA , c_eA = c_eA, w_eA = w_eA, ref_temp2 = 1)
  hold$r_eA <- r_eA
  hold$c_eA <- c_eA
  hold$w_eA <- w_eA
  hold$r_c_ratio <- log(c_eA/r_eA)
  results <- bind_rows(results, hold)
  }}


plot_grid(ggplot(results, aes(x = T, y = r1, group = interaction(r_c_ratio, w_eA), color = r_c_ratio, linetype = factor(w_eA)))+
            geom_line()+
            scale_linetype("w_eA")+
            scale_color_viridis_c(end = 0.8, "log ratio\nc/r eA")+
            coord_cartesian(ylim = c(0, max(results$r1)))+
            theme(legend.justification=c(0,1), legend.position=c(0,1)),
          
          ggplot(filter(results,r1>0), aes(x = T, y = R1.star, group = interaction(r_c_ratio, w_eA), color = r_c_ratio, linetype = factor(w_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            coord_cartesian()+
            ylab(label = "R*"),
          
          ggplot(results, aes(x = T, y = a22, group = interaction(r_c_ratio, w_eA), color = r_c_ratio, linetype = factor(w_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            coord_cartesian()+
            ylab("alpha"),
          
          ggplot(results, aes(x = T, y = a22, group = interaction(r_c_ratio, w_eA), color = r_c_ratio, linetype = factor(w_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            scale_y_log10()+
            ylab("alpha"),
          
          ggplot(results, aes(x = T, y = K1, group = interaction(r_c_ratio, w_eA), color = r_c_ratio, linetype = factor(w_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            coord_cartesian(ylim = c(0, max(results$K1,na.rm = T)))+
            ylab(label = "K"),
          
          ggplot(results, aes(x = T, y = K1, group = interaction(r_c_ratio, w_eA), color = r_c_ratio, linetype = factor(w_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            #coord_cartesian(ylim = c(0, max(results$K1)))+
            ylab(label = "K")+
            scale_y_log10(), nrow = 3)
ggsave("figures/metabolic_explore.png", height = 11, width = 10)

ggplot(results, aes(x = T, y = a22, group = interaction(r_c_ratio, w_eA), color = r_c_ratio, linetype = factor(w_eA)))+
  geom_line()+
  facet_wrap(~r_c_ratio, scales = "free_y")+
  scale_linetype(guide = F)+
  scale_color_viridis_c(end = 0.8, guide = F)+
  ylab("alpha")

results %>% ggplot(aes(x = T, y = fit_ratio, color = w_eA, shape = coexist))+
  geom_hline(yintercept = 1)+
  geom_point()+
  scale_color_viridis_c(option = "B", end = 0.8)+
  scale_y_log10()+
  scale_shape_manual(values = c(1,19))

results2<-data.frame()
for(i in 1:10){
  for(j in 1:5){
    r_eA <- 0.6
    c_eA <- 0.1*(i-1)
    k_eA <- 0.1*(j-1)
    hold<-temp_RC(T = seq(0, 35,length = 100), r_eA = r_eA , c_eA = c_eA, w_eA = 0, k_eA = k_eA, ref_temp2 = 1)
    hold$r_eA <- r_eA
    hold$c_eA <- c_eA
    hold$k_eA <- k_eA
    results2 <- bind_rows(results2, hold)
  }}


#results2 <- results2 %>% filter(r >0)

plot_grid(ggplot(results2, aes(x = T, y = r1, group = interaction(c_eA, k_eA), color = c_eA, linetype = factor(k_eA)))+
            geom_line()+
            scale_linetype("k_eA")+
            scale_color_viridis_c(end = 0.8, "c eA")+
            coord_cartesian(ylim = c(0, max(results$r1)))+
            theme(legend.justification=c(0,1), legend.position=c(0,1)),
          
          ggplot(filter(results2,r1>0), aes(x = T, y = R1.star, group = interaction(c_eA, k_eA), color = c_eA, linetype = factor(k_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            ylab(label = "R*"),
          
          ggplot(results2, aes(x = T, y = a22, group = interaction(c_eA, k_eA), color = c_eA, linetype = factor(k_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            coord_cartesian()+
            ylab("alpha"),
          
          ggplot(results2, aes(x = T, y = a22, group = interaction(c_eA, k_eA), color = c_eA, linetype = factor(k_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            scale_y_log10()+
            ylab("alpha"),
          
          ggplot(results2, aes(x = T, y = K1, group = interaction(c_eA, k_eA), color = c_eA, linetype = factor(k_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            coord_cartesian(ylim = c(0, max(results$K1,na.rm = T)))+
            ylab(label = "K"),
          
          ggplot(results2, aes(x = T, y = K1, group = interaction(c_eA, k_eA), color = c_eA, linetype = factor(k_eA)))+
            geom_line()+
            scale_linetype(guide = F)+
            scale_color_viridis_c(end = 0.8, guide = F)+
            #coord_cartesian(ylim = c(0, max(results$K1)))+
            ylab(label = "K")+
            scale_y_log10(), nrow = 3)
ggsave("figures/metabolic_explore2.png", height = 11, width = 10)


results3<-temp_RC(T = seq(0, 35,length = 100), ref_temp2 = -2)

plot_grid(
results3 %>% 
  select(r1, r2, T) %>% 
  gather(key = species, value = r, r1:r2) %>% 
  ggplot(aes(x = T, y = r, color = species))+
  geom_line()+
  coord_cartesian(ylim = c(0,1)),

results3 %>% ggplot(aes(x = T, y = fit_ratio, shape = coexist))+
  geom_hline(yintercept = 1)+
  geom_point()+
  scale_shape_manual(values = c(1,19))+
  scale_y_log10(),nrow = 2)
ggsave("figures/metabolic_explore3.png", height = 6, width = 5)

#now vary Cs and Ws separately with T#####
temp_RC_2 <- function(T = 25, r_eA = 0.6, c_eA11 = 0.6, c_eA12 = 0.6, c_eA22 = 0.6, c_eA21 = 0.6, w_eA = 0, k_eA = 0, ref_temp2 = 0){
  D = 0.25 
  c11 = arrhenius(Temp = T, E = c_eA11, b1 = 0.2); c12 = arrhenius(Temp = T, E = c_eA12, b1 = 0.4); c21 = arrhenius(Temp = T, E = c_eA21, b1 = 0.5); c22 = arrhenius(Temp = T, E = c_eA22, b1 = 0.2) ### cij = per capita consumption of comsumer i on resource j
  w11 = arrhenius(Temp = T, E = w_eA, b1 = 0.2); w12 = arrhenius(Temp = T, E = w_eA, b1 = 0.4); w21 = arrhenius(Temp = T, E = w_eA, b1 = 0.4); w22 = arrhenius(Temp = T, E = w_eA, b1 = 0.2) ## wij = weighting factor that converts availability of resource j into consumer Ni
  k1 = arrhenius(Temp = T, E = k_eA, b1 = 0.1); k2 = arrhenius(Temp = T, E = k_eA, b1 = 0.1) #half saturation constant for N resource consumption
  r1 = arrhenius(Temp = T, E = r_eA, b1 = 0.5) - arrhenius(Temp = T, E = 0.74, b1 = 0.25); 
  r2 = arrhenius(Temp = T, E = r_eA, b1 = 0.5, ref_temp = ref_temp2) - arrhenius(Temp = T, E = 0.74, b1 = 0.25, ref_temp = ref_temp2) #population growth rates
  T1 = 1; T2 = 1 #minimum amount of resource required for growth
  
  species1_consumption <- c12/c11
  species2_consumption <- c22/c21
  
  #' Calculate slope and intercept for the ZNGI's (get this from equation S3 in Ke and Letten)
  y.inter.1 <-  (D * (k1 - T1) + r1 * T1) / (w12 * (r1 - D))
  y.inter.2 <-  (w21 / w22) * (D * (k2 - T2) + r2 * T2) / (w21 * (r2 - D))
  
  
  #' Calculate equilibrium resource values (when the ZNGIs cross each other, i.e. such that Ke and Letten equations S3.1 = S3.2)
  B1 <-  y.inter.1
  B2 <-  y.inter.2
  Lamda.1 <-  (w11 / w12)
  Lamda.2 <- (w21 / w22)
  R1.star <-  (B1 - B2) / (Lamda.1 - Lamda.2)
  R2.star <-  (B2 * Lamda.1 - B1 * Lamda.2) / (Lamda.1 - Lamda.2)  
  
  
  #' Supply concentrations of resources
  S1 <- 10 #0.248
  S2 <- 12  #0.245
  
  
  #Calculate alphas
  a11 <- (c12 + c11 * (w11/w12)) / (D * (S2 + (w11/w12) * S1 - B1))
  a12 <- (c22 + c21 * (w11/w12)) / (D * (S2 + (w11/w12) * S1 - B1))
  
  a22 <- (c22 + c21 * (w21/w22)) / (D *(S2 + (w21/w22) * S1 - B2))
  a21 <- (c12 + c11 * (w21/w22)) / (D *(S2 + (w21/w22) * S1 - B2))
  
  p <- sqrt((a12*a21)/(a11*a22)) #niche overlap
  stabil_potential <- 1 - p #stabilizing potential
  fit_ratio <- sqrt((a11*a12)/(a22*a21)) * (r2)/(r1) #fitness ratio
  coexist <- p < fit_ratio &  fit_ratio < 1/p
  
  data.frame(a11 = a11, a12 = a12, a22 = a22, a21 = a21, R1.star = R1.star, R2.star= R2.star, K1 = (r1)/a11, K2 = (r2)/a22, T = T, r1 = r1, r2 = r2, c = c11, stabil_potential = stabil_potential, fit_ratio = fit_ratio, p = p , coexist = coexist, species1_consumption = species1_consumption, species2_consumption = species2_consumption)
}

results4<-data.frame()
for(i in 1:5){
  for(j in 1:5){
    for(k in 1:5){
      for(l in 1:5){
        r_eA <- 0.6
        c_eA11 <- 0.1*(i-1)
        c_eA22 <- 0.1*(j-1)
        c_eA21 <- 0.1*(k-1)
        c_eA12 <- 0.1*(l-1)
        hold<-temp_RC_2(T = seq(0, 35,by = 2), r_eA = r_eA , c_eA11 = c_eA11, c_eA12 = c_eA12, c_eA22 = c_eA22, c_eA21 = c_eA21)
        hold$r_eA <- r_eA
        hold$c_eA11 <- c_eA11
        hold$c_eA12 <- c_eA12
        hold$c_eA22 <- c_eA22
        hold$c_eA21 <- c_eA21
        results4 <- bind_rows(results4, hold)
      }}}}

x <- seq(0.5,1.5,length = 100)

lines.df <- data.frame(x = x, y = 1/(x), y2 = x)

ggplot(results4, aes(x = stabil_potential, y = fit_ratio, color = T, shape = coexist))+
  geom_hline(yintercept = 1, color = "grey", linetype = 2)+
  geom_vline(xintercept = 0, color = "grey", linetype = 2)+
  geom_line(data = lines.df, aes(x = 1-x, y = y, color = NULL, shape = NULL))+
  geom_line(data = lines.df, aes(x = 1-x, y = y2, color = NULL, shape = NULL))+
  geom_point()+
  facet_grid(c_eA11~c_eA12)+
  scale_y_log10()+
  scale_color_viridis_c(option = "B", end = 0.8)+
  scale_shape_manual(values = c(1,19))
ggsave("figures/vary_cii.png", height = 8, width = 10)

ggplot(results4, aes(x = stabil_potential, y = fit_ratio, color = T, shape = coexist))+
  geom_hline(yintercept = 1, color = "grey", linetype = 2)+
  geom_vline(xintercept = 0, color = "grey", linetype = 2)+
  geom_line(data = lines.df, aes(x = 1-x, y = y, color = NULL, shape = NULL))+
  geom_line(data = lines.df, aes(x = 1-x, y = y2, color = NULL, shape = NULL))+
  geom_point()+
  facet_grid(~c_eA11)+
  scale_y_log10()+
  scale_color_viridis_c(option = "B", end = 0.8)+
  scale_shape_manual(values = c(1,19))
ggsave("figures/coexistence-temperature-substitutable.png", width = 8, height = 6)

