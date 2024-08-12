setwd(this.path::here())
load("../submission/AnswerC1.Rdata")
library(evgam)
library(mev)
load("../outputs/models_and_data_imputed.Rdata")
load("../outputs/AmaurotTestSet_CP.Rdata")
load("../outputs/threshold_levels_lists.Rdata")

source("Task1-functions.R")
test_data <- read.csv("../data/AmaurotTestSet.csv")

scale <- with(test_data,
     exp(-1) + 18.7*sqrt(V2) + 9*(1+ log(V3))^2 + 5.71*WindSpeed^1.5)/10
shape <- -1/pi^2
loc <- with(test_data, ifelse(
  Season == "S1", 112-6*abs(WindDirection), 110-5*abs(WindDirection)^0.9)
)
if(!"outputs/C1_point_estimates.RData" %in% list.files("..",recursive = TRUE)){
n <- nrow(test_data)
qtrue <- ptrue <- numeric(n)
B <- 3e7
for(i in seq_len(n)){
  set.seed(2023)
 samp <- mev::rgp(n = B, loc = 0, scale = scale[i], shape = shape)
 if(test_data$Season[i] == "S1"){
   Y <- na.omit(ifelse(samp > loc[i], samp, ifelse(
     rbeta(B, samp, loc[i]) > runif(B), samp, NA)))
 } else{
   Y <- na.omit(ifelse(samp > loc[i], samp, ifelse(
     rbeta(B, exp(2 + (samp-loc[i])/30), 1) > runif(B), samp, NA)))
 }
 qtrue[i] <- quantile(Y, 0.9999)
 ptrue[i] <- mean(Y >= 112)
}
qtrue2 <- ifelse(ptrue < 1e-4,
                  qtrue,
                  qgp(pmax(0,(ptrue-1e-4)/ptrue),
                      loc = 112,
                      scale = scale-112/pi^2,
                      shape = shape))
 save(qtrue, ptrue, qtrue2, file = "../outputs/C1_point_estimates.RData")
} else{
  load("../outputs/C1_point_estimates.RData")
}
r1 <- range(qtrue)
r2 <- range(AnswerC1)
r3 <- c(min(r1[1], r2[1]), max(r1[2], r2[2]))


pdf("../figures/Figure3.pdf", width = 8.5, height = 4)
par(mfrow = c(1,2), mar = c(4,4,1,1), bty = "l")

cols <- ifelse(qtrue2 > AnswerC1[,2] & qtrue2 < AnswerC1[,3], "black", "gray")
plot(x = qtrue2,
     y = AnswerC1[,1],
     xlab = "True quantile",
     ylab = "Predicted quantile",
     ylim = r3,
     xlim = r3,
     cex.lab = 1.25,
     col = cols,
     pch = 20,
     panel.first = {abline(a=0,b=1)})

for(i in 1:100){
  segments(x0 = qtrue2[i], y1 = AnswerC1[i,2], y0 = AnswerC1[i,3], col = cols[i])
}



zeta <- 0.98
probs <- 0.9999 ## probability of interest
m <- 300 # number of days in a year
theta <- 1
nsim <- 1e4 # number of simulation
N.year <- 1 / (m * (1 - probs))  # number of years for the return level

#---------------------------

results_bestmodel <- fun_evgam_gpd(
   fit_thr_ald = thr.levels.lists[[4]],
   model.gpd = model_gpd.list[[44]],
   data_train = Amaurot.imp.final,
   data_test = AmaurotTestSet,
   m = m,
   zeta = 1 - zeta,
   theta = theta,
   N.year = N.year,
   nsim = nsim
)

AnswerC1 <- results_bestmodel$C50_95_results[, 1:3]

r1 <- range(qtrue)
r2 <- range(AnswerC1)
r3 <- c(min(r1[1], r2[1]), max(r1[2], r2[2]))


cols <- ifelse(qtrue2 > AnswerC1[, 2] & qtrue2 < AnswerC1[, 3],
               "black",
               "gray")
plot(
   x = qtrue2,
   y = AnswerC1[, 1],
   xlab = "True quantile",
   ylab = "Predicted quantile",
   ylim = r3,
   xlim = r3,
   cex.lab = 1.25,
   col = cols,
   pch = 20,
   panel.first = {
      abline(a = 0, b = 1)
   }
)

for (i in 1:100) {
   segments(
      x0 = qtrue2[i],
      y1 = AnswerC1[i, 2],
      y0 = AnswerC1[i, 3],
      col = cols[i]
   )
}
dev.off()
