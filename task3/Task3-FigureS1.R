setwd(this.path::here())
library(mgcv)
source("Task3-functions.R")
Coputopia <- as.matrix(read.csv("../data/Coputopia.csv", header = TRUE)[,3:5])

# Create a plot of coefficients of tail dependence and probability estimates
qu <- seq(0.9, 0.995, by = 0.001)
HRV <- probs_hrv(quantiles = qu,
                 data = Coputopia,
                 a1 = 0.96, # 25 points above
                 a2 = 0.97, # 26 points above
                 qlev = "theoretical")
set.seed(202404)
bootstrap <- boot::boot(data = Coputopia,
                        R = 999L,
                        statistic = function(data, indices, ...){
 est <- probs_hrv(quantiles = qu,
            data = data[indices,],
            a1 = 0.96, # 25 points above
            a2 = 0.97, # 26 points above
            qlev = "theoretical")
 log(c(est$p1, est$p2))
})
logprobsCI <- t(sapply(
  X = 1:(2*length(qu)),
  FUN = function(i){
    boot::boot.ci(bootstrap, index = i, type = "perc")$perc[4:5]
  }))

# There are 25 exceedances if we pick 96%, 10 if we pick 98%
pdf("../figures/FigureS1.pdf", width = 8.5, height = 5)
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
     y = predict(mgcv::gam(log(HRV$p1) ~ s(HRV$quantiles))),
     type = "l",
     ylab = expression("log prob."),
     xlab = "quantile level",
     ylim = c(-14, -10))
# points(HRV$quantiles, y = log(HRV$p1), pch = 20, cex = 0.5)
matplot(qu,
        y = logprobsCI[1:length(qu),],
        lty = 2,
        add = TRUE, col = 1, type = "l")
rug(side = 2, log(HRV$probs[1]),ticksize = -0.01)
rug(side = 2, log(HRV$probs[2]), col = 'grey',ticksize = -0.01)
lines(x = HRV$quantiles,
      y = predict(loess(log(HRV$p2)~ HRV$quantiles)),
       col = 'grey')
# points(HRV$quantiles, y = log(HRV$p2), pch = 20, cex = 0.5,col = 'grey')
matplot(qu,
        y = logprobsCI[length(qu) + 1:length(qu),],
        lty = 2,
        col = 'grey',
        add = TRUE, type = "l")
dev.off()
