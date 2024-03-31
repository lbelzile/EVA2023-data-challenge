## Load libraries
setwd(this.path::here())
library(mev)
library(Ternary)
library(patchwork)
data <- read.csv("../data/Coputopia.csv", header = TRUE)
d12 <- mev::taildep(data = with(data, cbind(Y1, Y2)[Season == "S1",]),
             u = seq(0.8, 0.995, by = 0.01),
             method = list(eta = "betacop",
                           chi = "betacop"),
             confint = "lrt",
             plot = FALSE)
d13 <- mev::taildep(data = cbind(data$Y1, data$Y3),
                    u = seq(0.8, 0.995, by = 0.01),
                    method = list(eta = "betacop",
                                  chi = "betacop"),
                    confint = "lrt",
                    plot = FALSE)
d23 <- mev::taildep(data = cbind(data$Y2, data$Y3),
                    u = seq(0.8, 0.995, by = 0.01),
                    method = list(eta = "betacop",
                                  chi = "betacop"),
                    confint = "lrt",
                    plot = FALSE)
# Transform observations to unit Frechet scale and plot angles
frechet <- mev::qgev(mev::pgev(as.matrix(data[,3:5])), loc = 1, shape = 1)
radius <- rowSums(frechet)
angs <- frechet[,1:3]/radius


pdf("../figures/Figure6a.pdf", width = 4, height = 4)
par(mar = rep(0.2, 4))
TernaryPlot(
   atip = expression(W[1]), 
   btip = expression(W[2]), 
   ctip = expression(W[3]),grid.minor.lines = 0, 
   axis.labels = seq(0, 1, length.out = 11))
# Colour plot background
# ColourTernary(TernaryDensity(angs[radius > quantile(radius, 0.95),], resolution = 20L))
TernaryPoints(angs[radius > quantile(radius, 0.95),],  pch = 20, cex = 0.5)
dev.off()
pdf("../figures/Figure6b.pdf", width = 5, height = 6)
par(mfrow = c(2,1), bty = "l", mar = c(4,4,0.5,0.1))
matplot(x = d12$u,
        y = d12$eta,
        lty = c(1,2,2),
        col = 1,
        ylim = c(0.5,1),
        yaxs = "i",
        type = "l",
        ylab = expression(eta),
        xlab = "Quantile level",
        lwd = 1.25,
        cex.axis = 1.25,
        cex.lab = 1.25)
matplot(x = d13$u,
        y = d13$eta,
        lty = c(1,2,2),
        col = "gray90",
        add = TRUE,
        type = "l")
matplot(x = d23$u,
        y = d23$eta,
        lty = c(1,2,2),
        add = TRUE,
        col = "gray50",
        type = "l")
matplot(x = d12$u,
        y = d12$chi,
        lty = c(1,2,2),
        lwd = 1.25,
        col = 1,
        ylim = c(0,0.5),
        yaxs = "i",
        type = "l",
        ylab = expression(chi),
        xlab = "Quantile level",
        cex.lab = 1.25,
        cex.axis = 1.25)
matplot(x = d13$u,
        y = d13$chi,
        lty = c(1,2,2),
        col = "gray90",
        add = TRUE,
        type = "l")
matplot(x = d23$u,
        y = d23$chi,
        lty = c(1,2,2),
        add = TRUE,
        col = "gray50",
        type = "l")
dev.off()


