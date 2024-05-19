setwd(this.path::here())
library(geometricMVE)
library(texmex)
source("Task3-functions.R")
data <- pgumbel(read.csv("../data/Coputopia.csv")[, -(1:2)])
n <- nrow(data)
train_size <- n / 3

### Combination
levs1 <- 0.95
levs2 <- 0.96
thresh <- 0.94
## note that levs1 and levs2 should be larger than thresh in case of HT models

prob.levels <- c(levs1, levs2, 0.5)
# probability levels at which to perform the cross-validation

# Break down observations to perform stratified sampling
# Extract the same number of observations for which data[,3] < 0.5
g1 <- which(data[,3] < 0.5)
g2 <- which(data[,3] >= 0.5)
ng1 <- length(g1)
ng2 <- length(g2)

### check the empirical probabilities and if we have enough number of data points
exceed.p1 <- exceed.p2 <- c()
for (i in 1:10) {
  ran.ind <- sample(1:n, size = train_size, replace = FALSE)
  train_data <- data[ran.ind,]
  test_data <- data[-ran.ind,]
  exceed.p1[i] <- sum(apply(test_data, 1, min) > prob.levels[1])
  test_data <- test_data[test_data[, 3] < prob.levels[3], -3]
  exceed.p2[i] <- sum(test_data[, 1] > prob.levels[2] & test_data[, 2] > prob.levels[2])
}

summary(exceed.p1)
summary(exceed.p2)
# fit model once to obtain suitable initial values

##asymmetric logistic 2d
fit.geom_asyLog.2d <- fit_geometric.2d(
  data = data[data[, 3] < prob.levels[3], 1:2],
  gauge.fun = gauge_rvad,
  initial.val =  c(3,  0.9),
  lower.limit = c(0, 0),
  upper.limit = c(5, 1),
  nsim = 1e5L,
  thresh = thresh
)

##asymmetric logistic 3d
pdf("../figures/Task3_QQplot_geometric.pdf", width = 5, height = 5)
par(mar = c(4.5, 4.5, 0.5, 0.5), bty = "l")
fit.geom_asyLog.3d <- fit_geometric.3d(
  data = data,
  gauge.fun = gauge_rvad_all3d,
  initial.val =  c(3, rep(0.33, length.out = 3), 0.9),
  lower.limit = c(0, rep(0, 4)),
  upper.limit = c(5, rep(1, 4)),
  nsim = 1e5L,
  thresh = thresh,
  diagnostic = TRUE
)
dev.off()

## HT model
fit.HT1 <- fit_HT_fixed_margin(
  data = data,
  cond.var = 1,
  thresh = thresh)
fit.HT2 <- fit_HT_fixed_margin(
  data = data,
  cond.var = 2,
  thresh = thresh)

### calculations of probabilities by doing the three-fold cross-validation repeatedly
B <- 100
emp.probs <-  array(dim = c(B, 3, 2))
## first dim is the number of bootstraps samples,
## second is for the 3-folds,
## third is for the three probabilities
est.probs <- array(dim = c(B, 3, 3, 3))
## first dim is the number of bootstraps samples,
## second is for the 3-folds,
## third is for the three method (HT, HRV, geometics) and
## fourth is for the three probabilities
bin.holds <- array(dim = c(B, 3, 2, 2))
## first dim is the number of bootstraps samples,
## second is for the 3-folds,
## third is for the all the 14000 binary points
## fourth is for the two types of probabilities (p1 and p2)

for (i in 1:B) {
  ## loop over three folds
  set.seed(i + 12345)
  inds1 <- sample(g1, size = ng1)
  inds2 <- sample(g2, size = ng2)
  for (j in 1:3) {
    ## loop over three folds
    indtrain <- c(inds1[1:(ng1/3)+(j-1)*ng1/3], inds2[1:(ng2/3)+(j-1)*ng2/3])
    train_data <- data[indtrain,]
    test_data <-  data[-indtrain,]

    est.probs[i, j, , ] <- estimate_probs_C3(
      data = as.matrix(train_data),
      nsim = 1e5,
      a1 = 0.93,
      a2 = 0.93,
      thresh = thresh,
      init.val.geom.3d = fit.geom_asyLog.3d$mle,
      init.val.geom.2d = fit.geom_asyLog.2d$mle,
      # init.val.HT = fit.HT1$est.param,
      levs = prob.levels,
      repl = 10
    )
    emp.probs[i, j, ] <- emp_C3(
      data = as.matrix(test_data),
      levs = prob.levels)
    bin.holds[i, j, , 1] <-
      as.integer(table(apply(test_data, 1, function(x){ isTRUE(min(x) > prob.levels[1])})))
    bin.holds[i, j, , 2] <-
      as.integer(table(ifelse(test_data[, 1] > prob.levels[2] & test_data[, 2] > prob.levels[2] & test_data[, 3] < prob.levels[3], 1L, 0L)))
  }
  print(i)
}

save(
  est.probs,
  emp.probs,
  bin.holds,
  file = "../outputs/Task3_crossvalidation.RData"
)


load("../outputs/Task3_crossvalidation.RData")

brierscore <- logscore <- array(dim = c(B,3,3))
# Compute Brier and logarithmic scores
B <- dim(emp.probs)[1]
for(b in 1:B){
  for(i in 1:3){ # method (HT, HRV, geometric)
    for(j in 1:3){
        brierscore[b,i,j] <- -sum(bin.holds[b,,1,min(j,2)]*est.probs[b, ,i,j]^2 +
                                    bin.holds[b,,2,min(j,2)]*(1-est.probs[b,,i,j])^2)
        logscore[b,i,j] <- -sum(bin.holds[b,,1,min(j,2)]*log(est.probs[b,,i,j]) +
                                    bin.holds[b,,2,min(j,2)]*log(1-est.probs[b,,i,j]))
    }
  }
}
# Need to maximize the score (positively oriented scores)
apply(brierscore, 2:3, mean)
apply(brierscore, 2:3, sd)
# All within the margin of error

# According to logscore, geometric best for p1 and p2, HT worst
mlogscore <- apply(logscore, 2:3, mean)/nrow(data)
# Standard error (not independent, but...)
apply(logscore, 2:3, sd)/nrow(data)/sqrt(B)
options(knitr.kable.NA = '')
row.names(mlogscore) <- c("conditional","HRV","geometric")
cat(
knitr::kable(mlogscore, 
              col.names = c("$\\widetilde{p}_1$", "$\\widetilde{p}_2$", "$\\widetilde{p}_2^{\\star}$"),
              booktabs = TRUE,
              format = "latex",
              digits = 2,
              escape = FALSElinesep = ""),
file = "../tables/Table6.tex")
