## Load libraries
setwd(this.path::here())
library(mev)
library(revdbayes)
library(evgam)
load("../outputs/threshold_levels_lists.Rdata")
load("../outputs/models_and_data_imputed.Rdata")
load("../outputs/Amaurot_imputed.Rdata")
Amaurot <- read.csv("../data/Amaurot.csv")
load("../outputs/train_test_data_imputed.Rdata")
load("../outputs/fitted-thr-models_and_GP-models.Rdata")
source("Task2-functions.R")
## Load data
data <- read.csv(  file = "../data/Amaurot.csv")
data <- data |>
   dplyr::mutate(Season = as.integer(as.factor(Season)))
Y <- data$Y

# Naive explanatory analysis, without covariates
u <- quantile(Y, 0.95)

# Fit generalized Pareto model and compute profile likelihood
mle <- fit.gpd(xdat = Y, threshold = u, show = TRUE)
profile <- mev::gpd.pll(
   psi = seq(170, 240, length.out = 101),
   param = "quant",
   dat = Y,
   threshold = u,
   p = 1/(0.1*60000))
# Obtain return levels
retlev <- function(threshold, scale, shape,
                   zetau, m, nty){
   threshold + scale/shape*((m*nty*zetau)^shape-1)
}
# Maximum likelihood estimator
ret <- retlev(threshold = u,
              scale = coef(mle)['scale'],
              shape = coef(mle)['shape'],
              zetau = 0.1,
              m = 200,
              nty = 300)

# Bayesian analysis using ratio-of-uniform method
bay_gp <- revdbayes::rpost_rcpp(
   n = 1e5,
   model = "bingp",
   data = Y,
   # nrep = 100,
   thresh = u,
   prior = revdbayes::set_prior(prior = "mdi",
                     model = "gp"))
# Extract parameters from the posterior
sigmau <- bay_gp$sim_vals[,"sigma[u]"]
xi <- bay_gp$sim_vals[,"xi"]
zeta <- bay_gp$bin_sim_vals
# Compute posterior distribution of return levels
retlev_gp <- bay_gp$thresh + sigmau/xi*((60000*zeta)^xi-1)
plot(density(retlev_gp),
     main = "",
     xlab = "",
     xlim = c(170,240))

## Alternative formulation with Poisson point process
## Not reported in the paper, as result is VERY similar
# bay_ipp <- revdbayes::rpost_rcpp(n = 1e5,
#                                  model = "pp",
#                                  noy = 70,
#                                  data = Y,
#                                  thresh = u,
#                                  prior = set_prior(prior = "mdi",
#                                                    model = "pp"))
# post_pp <- bay_ipp$sim_vals
# retlev_pp <- revdbayes::qgev(p = 0.364,
#                              loc = post_pp[,'mu'] - post_pp[,'sigma']*
#                                 (1-200^post_pp[,'xi'])/post_pp[,'xi'],
#                              scale = post_pp[,'sigma']*200^post_pp[,'xi'],
#                              shape = post_pp[,'xi'])
# lines(density(retlev_pp), col = "grey")

# Compute loss function pointwise for posterior
loss_uncond <- lossfun(qhat = (qu <- seq(180, 250, by = 0.01)),
                q = retlev_gp)
# Optimal value to report
qu[which.min(loss)]
# 0-1 loss gives back negative log density of posterior
kdens <- density(retlev_gp, adjust = 2, bw = "SJ")

###############################################################################

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
loss_cond <- optimize(f = function(q){
  lossfun(qhat = q, q = retlevs)
}, interval = c(0.9, 1.1)*range(retlevs), tol = 1e-3)$minimum



pdf("../figures/Figure5.pdf", width = 7.5, height = 3.5)
par(bty = "n", mar = c(4,4,1,1), mfrow = c(1,2))
# Profile log-likelihood with confints
# plot(profile,
#      xlim = c(170, 220),
#      ylim = c(-5,0),
#      ylab = "Profile log-likelihood",
#      xlab = "Return level")
# Plot loss function as function of qhat
# UNCONDITIONAL
plot(qu, loss_uncond-min(loss_uncond),
     type = 'l',
     xlab = "Return level",
     ylab = "Loss",
     xlim = c(160,240),
     lwd = 2)
lines(x = kdens$x[kdens$x < 220],
      y = (max(log(kdens$y))-log(kdens$y))[kdens$x < 220],
      xlim = c(160, 220), lty = 2)
abline(v = qu[which.min(loss)], lwd = 0.5, lty = 2)


# CONDITIONAL on covariates, model 2, threshold 0.98
# Since curve is nearly linear, subsample value of return levels to reduce
# plot size
sub <- seq(1, 951, by = 15L)
matplot(xs[sub], logtailprob[sub,], type = "l", col = scales::alpha(1, 0.1), lty = 1,
        xlab = "Return level",
        ylab = "Log tail probability",
        xlim = c(170,220),
        xaxs = "i")
abline(h = -log(60000))
rug(retlevs)
abline(v = loss_cond, lty = 2, lwd = 0.5)
dev.off()
