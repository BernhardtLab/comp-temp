



R1N <- (m1 * k1N)/ (r1 - m1)
R1P <- (m1 * k1P)/ (r1 - m1)



R2P <- (m2 * k2P)/ (r2 - m2)


(R2N <- (0.5 * 2)/ (1.5 - 0.5)) ## 1, ### this mean r2 = 1.5, r1 = 2, k1P = 1.5, k2P = 2
(R1P <- (0.5 * 1.5)/ (2 - 0.5)) ## 0.5

# a11 <- c1P / (D * (SP - R1P));  c1P = 7, D = 0.5, SP = 1
# a12 <- c2P / (D * (SP - R1P)) c2P = 4.5, D = 0.5, SP = 1
# a21 <- c1N / (D * (SN - R2N)); C1N = 11, SN = 2
# a22 <- c2N / (D * (SN - R2N)); C2N = 13, SN = 2


### alphas should be between 7 and 30
# a11 <- c1P / (D * (1 - R1P))
(a11 <- 7 / (0.5 * (1 - R1P)))
(a12 <-  4.5 / (0.5 * (1 - R1P))) 
(a21 <-  11 / (0.5 * (2 - R2N)))
(a22 <- 13 / (0.5 * (2 - R2N)))

# sqrt((a11*a12)/(a22*a21)) this needs to get smaller

# a11 <- 31
# a12 <- 18
# a21 <- 23
# a22 <- 26


rho <- sqrt((a12*a21)/(a11*a22)) #niche overlap
stabil_potential <- 1 - rho #stabilizing potential
fit_ratio <- sqrt((a11*a12)/(a22*a21)) #fitness ratio 
coexist <- rho < fit_ratio &  fit_ratio < 1/rho
coexist
