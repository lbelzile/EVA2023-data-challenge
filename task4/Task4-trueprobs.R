library(mev)
setwd(this.path::here())
source("Task4-functions.R")
# Below, we assume that clusters are independent
# and compute the joint probability of exceedance 
# for each in turn, so as to compare with our
# estimates and discuss

# Clusters of inverted logistic max-stable
B <- 1e4L
theta <- seq(from = 0.4, to = 0.9, length.out = B)
m2 <- c(7L, 5L)
m5 <- c(7L, 2L)
# Average probability for cluster C2 and C5
logp2C2 <- copula:::lsum(-(sum(m2)^theta)*log(300)) - log(B)
logp2C5 <- copula:::lsum(-(sum(m5)^theta)*log(300)) - log(B)


logp1C2_i <- sapply(theta, function(theta_i){
  -expme(z = qgev(p = c(rep(1/300, m2[1]), rep(12/300, m2[2])), 1,1,1),
         model = "log", par = theta_i)
  })
logp1C5_i <- sapply(theta, function(theta_i){
  -expme(z = qgev(p = c(rep(1/300, m5[1]), rep(12/300, m5[2])), 1,1,1),
         model = "log", par = theta_i)
})
logp1C2 <- copula:::lsum(logp1C2_i) - log(B)
logp1C5 <- copula:::lsum(logp1C5_i) - log(B)


# Clusters of asymptotically dependent data
m1 <- c(3L, 5L)
m4 <- c(5L, 3L)

# The following approximation follows from de Haan and Ferreira, 2006, Theorem 9.3.1
logp1C1_i <- sapply(theta, function(theta_i){
  log(expme_log_min(u = qgev(p = 1 - c(rep(1/300, m1[1]), rep(12/300, m1[2])), 1,1,1), par = theta_i))
})
logp1C4_i <- sapply(theta, function(theta_i){
  log(expme_log_min(u = qgev(p = 1 - c(rep(1/300, m4[1]), rep(12/300, m4[2])), 1,1,1), par = theta_i))
})
logp1C1 <- copula:::lsum(logp1C1_i) - log(B)
logp1C4 <- copula:::lsum(logp1C4_i) - log(B)

# Same probability level
logp2C1_i <- sapply(theta, function(theta_i){
  c(0, sapply(1:sum(m1), function(i){lchoose(sum(m1), i)+ (i^theta_i)*log(1-1/300)}))
})
# The below is a (quite accurate) approximation
# logp2C1_i <- sapply(theta, function(theta_i){
#   log(expme_log_min(u = qgev(p = 1 - rep(1/300, sum(m1)), 1,1,1), par = theta_i))
# })
# copula:::lsum(logp2C1_i) - log(B)

logp2C1 <- logp2C4 <- 
  copula:::lsum(copula:::lssum(logp2C1_i, 
                 signs = rep(c(1,-1), length.out = sum(m1) + 1))) - log(B)

# Cluster of Gaussian
B <- 1000L
m3 <- c(3L, 10L)
rho <- seq(0.1, 0.7, length.out = B)
set.seed(1234)
logp2C3_i <- sapply(rho, function(rho_i){
  Cov <- diag(rho_i, sum(m3)) + matrix(1-rho_i, sum(m3), sum(m3))
  TruncatedNormal::pmvnorm(mu = rep(0, sum(m3)),
                           sigma = Cov, 
                           lb = rep(-Inf, sum(m3)),
                           ub = rep(qnorm(1/300), sum(m3)))
})
logp2C3 <- copula:::lsum(log(logp2C3_i)) - log(B)
                           
logp1C3_i <- sapply(rho, function(rho_i){
  Cov <- diag(rho_i, sum(m3)) + matrix(1-rho_i, sum(m3), sum(m3))
  TruncatedNormal::pmvnorm(mu = rep(0, sum(m3)),
                          sigma = Cov, 
                          lb = rep(-Inf, sum(m3)),
                          ub = c(rep(qnorm(1/300), m3[1]), rep(qnorm(12/300, m3[2]))))
})
logp1C3 <- copula:::lsum(log(logp1C3_i)) - log(B)

# Probability by cluster
probs <- rbind(
  p1 = c(logp1C1, logp1C2, logp1C3, logp1C4, logp1C5),
  p2 = c(logp2C1, logp2C2, logp2C3, logp2C4, logp2C5))
# Probability under independence
save(probs, file = "../outputs/C4_truth_postmortem.RData")
