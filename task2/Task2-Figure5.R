## Load libraries
setwd(this.path::here())
library(mev)
library(revdbayes)
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
   thresh = u,
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
   prior = set_prior(prior = "mdi",
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

# Loss function
lossfun <- function(qhat, q = retlev){
   loss <- function(qhat, q){
      mean(ifelse(0.99*q > qhat,
                  0.99*(0.99*q-qhat),
                  ifelse(1.01*q < qhat,
                         0.1*(qhat-1.01*q),
                         0)))
   }
   retval <- numeric(length = length(qhat))
   for(i in seq_along(retval)){
      retval[i] <- loss(qhat[i], q = q)
   }
   retval
}
# Compute loss function pointwise for posterior
loss <- lossfun(qhat = (qu <- seq(180, 250, by = 0.01)),
                q = retlev_gp)
# Optimal value to report
qu[which.min(loss)]
# 0-1 loss gives back negative log density of posterior
kdens <- density(retlev_gp, adjust = 2, bw = "SJ")


pdf("../figures/C2_unconditional.pdf", width = 8, height = 4)
par(bty = "n", mar = c(4,4,1,1), mfrow = c(1,2))
# Profile log-likelihood with confints
plot(profile,
     xlim = c(170, 220),
     ylim = c(-5,0),
     xlab = "return level")
# Plot loss as function of qhat
plot(qu, loss-min(loss),
     type = 'l',
     xlab = "return level",
     ylab = "loss",
     xlim = c(160,240),
     lwd = 2)
lines(x = kdens$x[kdens$x < 220],
      y = (max(log(kdens$y))-log(kdens$y))[kdens$x < 220],
      xlim = c(160, 220), lty = 2)
abline(v = qu[which.min(loss)], lwd = 0.5, col = "grey", lty = 2)
dev.off()
   
   
