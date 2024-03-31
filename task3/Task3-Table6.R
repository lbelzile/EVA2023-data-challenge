setwd(this.path::here())
library(texmex)
library(evd)
library(geometricMVE)
source("Task3-functions.R")
Coputopia <- as.matrix(read.csv("../data/Coputopia.csv", header = TRUE)[,3:5])

# Create a plot of coefficients of tail dependence and probability estimates
HRV <- probs_hrv(quantiles = seq(0.9, 0.995, by = 0.001),
                 data = Coputopia,
                 a1 = 0.96, # 25 points above
                 a2 = 0.97, # 26 points above
                 qlev = "theoretical")
# There are 25 exceedances if we pick 96%, 10 if we pick 98%
pdf("../figures/C3_HRV.pdf", width = 8.5, height = 5)
par(mfrow = c(1,2), mar = c(4,5,1,1), bty = "l")
matplot(HRV$quantiles,
        HRV$eta1,
        type = 'l',
        lty = c(2,1,2),
        col = "black",
        ylab = expression(eta),
        xlab = 'quantile level',
        ylim = c(0,1),
        yaxs = "i")
matplot(x = HRV$quantiles,
        y = HRV$eta2,
        type = 'l',
        lty = c(2,1,2),
        col = "grey",
        add = TRUE)
plot(x = HRV$quantiles, 
     y = 1e6*HRV$p1, 
     pch = 20, 
     ylab = expression("probability ("%*%10^6~")"),
     xlab = "quantile level",
     panel.first = {abline(h = 1e6*HRV$probs[1])},
     ylim = 1e6*c(1e-6, 2e-5))
points(x = HRV$quantiles, 
       y = 1e6*HRV$p2, col = 'grey',
       panel.first = {abline(h = 1e6*HRV$probs[2], col = 'grey')})
dev.off()


# Dependence quantile
thresh <- c(0.9, 0.95, 0.96, 0.97, 0.98)
a1 <- 0.96
a2 <- 0.97

set.seed(2023)
resultsC3 <- list()
for(i in seq_along(thresh)){
   resultsC3[[i]] <- estimate_probs_C3(data = Coputopia, thresh = thresh[i])
   print(i)
   save(resultsC3, file = "../outputs/Task3-results.RData")
}
res_tab <- matrix(
   data = unlist(lapply(
   resultsC3, function(x){
      sprintf(fmt = "%.2f", 1e6*x)})), 
   nrow = 3)[,c(seq(1L, 9L, by = 2L), seq(2L, 10L, by = 2L))]
rownames(res_tab) <- c("conditional", "HRV", "geometric")
cat(
knitr::kable(res_tab,
             booktabs = TRUE, 
             align = "r",
             format = "latex"),
 file = "../tables/Table6.tex")
