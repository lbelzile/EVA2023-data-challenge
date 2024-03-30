setwd(this.path::here())
library(mev)
source("Task4-functions.R")


# Load data and change column names to avoid duplicate labels

UtopulaU1_data <- read.csv(file = "../data/UtopulaU1.csv", header = TRUE)
colnames(UtopulaU1_data) <- paste0("I1_", c(paste0(0, 1:9), 10:25))
UtopulaU2_data <- read.csv(file = "../data/UtopulaU2.csv", header = TRUE)
colnames(UtopulaU2_data) <- paste0("I2_", c(paste0(0, 1:9), 10:25))
# Merge databases
Utopula <- data.frame(UtopulaU1_data, UtopulaU2_data)


#####################################################
####  Exploratory data analysis                ######
#####################################################
# Kendall tau correlation
cormat <- pcaPP::cor.fk(Utopula)
# relabel to remove "Y" from the matrix
colnames(cormat) <- rownames(cormat) <- 1:50

# Get order and pick observations from last cluster
clust_order <- corrplot::corrMatOrder(corr = cormat,
                                      order = 'hclust',
                                      hclust.method = "ward.D2")
# Compute hclust separately to sort them by island
hclust_od1 <- cutree(tree = hclust(as.dist(1 - cormat),
                                   method = "ward.D2"),
                     k = 5)
# Label-switching: use same cluster nomenclature as editorial
hclust_od <- dplyr::case_match(.x = hclust_od1,
                               3 ~ 1,
                               4 ~ 2,
                               5 ~ 3,
                               2 ~ 4,
                               1 ~ 5)
clust_order <- order(as.integer(paste0(hclust_od, c(paste0(0, 1:9), 11:50))))
# Transform to uniform scale
unif <- exp(-exp(-Utopula)) # Gumbel to uniform


#####################################################
####     APPROACH 1 - EMPIRICAL EXTRAPOLATION  ######
######         A LA LEDFORD-TAWN               ######
#####################################################
# Compute eta coefficient for all values in the cluster
# at a lower threshold and use Ledford and Tawn (1997) to extrapolate
# We do the extrapolation in exponential scale
# Unequal approach with an extension of Wadsworth and Tawn (2013)
prob1_cl <- eta_p1 <-  numeric(5)
prob2_cl <- eta_p2 <- numeric(5)
qlev <- 0.985 # leaves 200 observations
s1_exp <- qexp(1 - 1 / 300)
s2_exp <- qexp(1 - 12 / 300)
edata <- -log(1 - unif)

for (cl in 1:5) {
  edata_cl <- edata[, which(hclust_od == cl)]
  # Structural variable
  Texp <- apply(edata_cl, 1, min)
  # Pick threshold
  s0_exp <- quantile(Texp, qlev)
  eta_p2[cl] <- mean(Texp[Texp > s0_exp]) - s0_exp
  # Compute empirical probability that the structure variable min
  # is above upper bound
  # p_ub <- mean(apply(unif[ , which(hclust_od == cl)], 1, min) > ub)
  prob2_cl[cl] <-
    exp(-(s1_exp - s0_exp) / pmin(eta_p2[cl], 1)) * (1 - qlev)

  clsize <- ncol(edata_cl)
  # Compute size of cluster and number of variables in U1 and U2
  u1size <- sum((1:50)[hclust_od == cl] <= 25)
  u2size <- clsize - u1size
  # Under exchangeability, data are exchangeable
  combos <- combn(m = u2size, x = clsize)
  # Estimate eta for weighted data for all permutations
  eta_cl <- tail_prob <- numeric(ncol(combos))
  for (j in 1:ncol(combos)) {
    weights <- rep(1, clsize)
    weights[combos[, j]] <- s2_exp / s1_exp
    wedat <- sweep(edata_cl, 2, weights, FUN = "/")
    un <- quantile(apply(wedat, 1, min), qlev)
    Texp <- apply(wedat, 1, min)
    eta_cl[j] <- mean(Texp[Texp > un] - un)
    tail_prob[j] <- exp(-(s1_exp - un) / pmin(1, eta_cl[j]))
  }
  eta_p1[cl] <- mean(eta_cl)
  # Compute weights relative to u1
  prob1_cl[cl] <- mean(tail_prob) * (1 - qlev)
}
log(prob2_cl)
log(prob1_cl)
# Check we respect ordering constraints
isTRUE(all(log(prob1_cl) > log(prob2_cl)))

# Print table with results
semipar_results <- rbind(prob1_cl, prob2_cl)
semipar_tbl <- apply(rbind(log(prob1_cl), log(prob2_cl)),
                     1:2, function(x) {
                       paste("$", sprintf("%.2f", x), "$")
                     })
rownames(semipar_tbl) <- c("$\\log \\widehat{p}_1$",
                           "$\\log \\widehat{p}_2$")

# TABLE 7
cat(
   knitr::kable(
     semipar_tbl,
     format = "latex",
     col.names = paste0("$C_", 1:5, "$"),
     booktabs = TRUE,
     escape = FALSE,
     align = c("rrrrr")),
   file = "../tables/Table7.tex"
)

