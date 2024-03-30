setwd(this.path::here())
library(evgam)
load("../outputs/threshold_levels_lists.Rdata")
load("../outputs/models_and_data_imputed.Rdata")
load("../outputs/Amaurot_imputed.Rdata")
Amaurot <- read.csv("../data/Amaurot.csv")
load("../outputs/train_test_data_imputed.Rdata")
load("../outputs/fitted-thr-models_and_GP-models.Rdata")
source("Task2-functions.R")
i = 4 # 0.98 threshold
j = 2 # Model 2
fit_thr_ald = thr.levels.lists[[i]]
### output of the fitted thresholds models i
model.gpd = model_gpd.list[[j]]
### the formula for the GPD models j
data_train = Amaurot.train
newdata = Amaurot.train
m = 1 # modified March 29th, 2024
Nyear = 300*200 # modified March 29th, 2024
method = "npboot"
zeta = 1 - fit_thr_ald$tau
ncov = 1e4L
set.seed(2023)
### the estimate of the thresholds based on the covariates
data_train$threshold <- predict(fit_thr_ald, newdata = data_train)$location
data_train$excess <- data_train$Y - data_train$threshold
## put NAs to the non-threshold excesses
data_train$excess <- with(
  data_train,
  ifelse(excess <= 0, NA, excess))
## fit GP distribution to threshold exceedances
fit_gpd <- evgam::evgam(formula = model.gpd,
                        data = data_train,
                        family = "gpd")
# Simulate nsim values from the posterior and compute associated
# thresholds, scale and shape for each
  ndat <- dplyr::slice_sample(newdata, n = ncov)
  # Generate thresholds from fitted ALD model and
  # obtain approximate posterior pred of zeta_u
  nsimu <- 100
  # Sample from approximate posterior of the ALD parameters
    # via a Gaussian approximation on the transformed scale
    u <- predict(fit_thr_ald, newdata = ndat)[,'location']
    post_ald <- simulate(
      fit_thr_ald,
      newdata = ndat,
      nsim = nsimu,
      type = "response")
    # Compute the probability of exceeding u
    # treating the latter as a fixed quantity
    pthresh <- sapply(seq_len(nsimu), function(i){
      pald(q = u,
           location = post_ald$location[,i],
           scale = post_ald$scale[,i],
           tau = 1 - zeta)})
  # Simulate generalized Pareto parameters and compute
  # resulting scale and shape for each combination of covariates
  gp_sim <- simulate(
    fit_gpd,
    newdata = ndat,
    nsim = nsimu,
    type = "response")
  xs <- seq(155, 250, by = 0.1)
 logtailprob <- matrix(0, ncol = nsimu, nrow = length(xs))
  for(j in seq_len(nsimu)){
    for(i in seq_along(xs)){
      x <- xs[i]
  logtailprob[i,j] <- log(mean(exp(log(1-pthresh[,j]) + revdbayes::pgp(
    q = x, loc = u, scale = gp_sim[[1]][,j],
    shape = gp_sim[[2]][,j], log.p = TRUE, lower.tail = FALSE))))
    }
  }
retlevs <- sapply(seq_len(nsimu),
       function(j){
        uniroot(f = function(x){ log(mean(exp(log(1-pthresh[,j]) + revdbayes::pgp(
           q = x, loc = u, scale = gp_sim[[1]][,j],
           shape = gp_sim[[2]][,j], log.p = TRUE, lower.tail = FALSE)))) + log(Nyear)
       }, interval = c(100, 300))$root})
loss <- optimize(f = function(q){
  lossfun(qhat = q, q = retlevs)
}, interval = c(0.9, 1.1)*range(retlevs), tol = 1e-3)$minimum


pdf("../figures/FigureS2.pdf", width = 6, height = 6)
par(bty = "l")
matplot(xs, logtailprob, type = "l", col = scales::alpha(1, 0.1), lty = 1,
        xlab = "quantile",
        ylab = "log tail probability",
        xlim = c(170,220),
        xaxs = "i")
abline(h = -log(60000))
rug(retlevs)
abline(v = loss, lty = 2)
dev.off()
#
# quant_covariates <- apply(
#   pthresh, 1, quantile,
#   probs = c(0.05,0.1, 0.025,0.5,0.75,0.9,0.95))
