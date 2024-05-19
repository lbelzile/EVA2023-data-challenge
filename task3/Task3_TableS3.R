setwd(this.path::here())
source("Task3-functions.R") ### to load the fit_geometric, quantR, and  pgumbel
data_coputopia <- read.csv("../data/Coputopia.csv")
data <- data_coputopia[, -c(1:2)]
library(geometricMVE)
thresh <- 0.95
y <- 6
qrR <- quantR(data = data, thresh = thresh)
set.seed(1234)
########## Fit the model with different choices of gauge functions
####### Gaussian gauge function
fit.Gauss <- fit_geometric(
  data = data,
  gauge.fun = gauge_gaussian3d,
  initial.val =  c(3, rep(0.33, length.out = 3)),
  lower.limit = c(0, rep(-1, length.out = 3)),
  upper.limit = c(5, rep(1, length.out = 3)),
  nsim = 1e5L,
  thresh = thresh,
  qr = qrR
)
#
# simulated.data <- sim.3d(
#   w = fit.Gauss$wexc,
#   r0w = fit.Gauss$r0w,
#   nsim = 1e5,
#   par = fit.Gauss$mle,
#   gfun = gauge_gaussian3d,
#   k = 3 * y / max(fit.Gauss$r0w)
# )

###### asymmetric logistic
fit.alog <- fit_geometric(
  data = data,
  gauge.fun = gauge_rvad_all3d,
  initial.val =  c(3, rep(0.33, length.out = 3), 0.9),
  lower.limit = c(0, rep(0, 4)),
  upper.limit = c(5, rep(1, 4)),
  nsim = 1e5L,
  thresh = thresh,
  qr = qrR
)


####### Clayton
fit.Clayton <- fit_geometric(
  data = data,
  gauge.fun = gauge_clayton3d,
  initial.val =  c(3, 0.33),
  lower.limit = c(0, 0),
  upper.limit = c(5, Inf),
  nsim = 1e5L,
  thresh = thresh,
  qr = qrR
)


### asymmetric logistic with trivariate dependence parameter set to 0
gauge_rvad_all3d_competition <- function(xyz, par) {
  gauge_rvad_all3d(xyz = xyz,
                   par = par,
                   theta123 = 0)
}

gauge_alog <- function(xyz, par) {
  gauge_rvad_all3d(xyz = xyz,
                   par = c(3, par))
}

fit.alog2 <- fit_geometric(
  data = data,
  gauge.fun = gauge_rvad_all3d_competition,
  initial.val = c(3, rep(0.33, length.out = 3)),
  lower.limit = c(0, rep(0, 4)),
  upper.limit = c(5, rep(1, 4)),
  nsim = 1e5L,
  thresh = thresh,
  qr = qrR
)
fit.alog3 <- fit_geometric(
  data = data,
  gauge.fun = gauge_alog,
  initial.val = rep(0.8, 4),
  lower.limit = rep(0, 4),
  upper.limit = rep(1, 4),
  nsim = 1e5L,
  thresh = thresh,
  qr = qrR
)


##### mixtures of  Gaussian and asymmetric logistics gauge functions
fit.mixt_gauss_alog <- fit_geometric(
  data = data,
  add.gauge = TRUE,
  gauge1 = gauge_gaussian3d,
  gauge2 = gauge_rvad_all3d,
  initial.val = c(3, rep(0.33, length.out = 3), rep(0.33, length.out = 3), 0.9),
  par1ind = 2:4,
  par2ind = 5:8,
  lower.limit = c(0, rep(-1, 3), rep(0, 4)),
  upper.limit = c(5, rep(1, 3), rep(1, 4)),
  nsim = 1e5L,
  thresh = thresh,
  qr = qrR
)


##### mixtures of Clayton and asymmetric logistics gauge functions
fit.mixt_clayton_alog <- fit_geometric(
  data = data,
  add.gauge = TRUE,
  gauge1 = gauge_clayton3d,
  gauge2 = gauge_rvad_all3d,
  initial.val = c(3, 0.33, rep(0.33, length.out = 3), 0.9),
  par1ind = 2,
  par2ind = 3:6,
  lower.limit = c(0, 0, rep(0, 4)),
  upper.limit = c(5, Inf, rep(1, 4)),
  nsim = 1e5L,
  thresh = thresh,
  qr = qrR
)

################## model comparison in terms of AIC values
AIC <-
  c(
    fit.Gauss$AIC,
    fit.Clayton$AIC,
    fit.alog$AIC,
    fit.alog2$AIC,
    fit.mixt_gauss_alog$AIC,
    fit.mixt_clayton_alog$AIC
  )
names(AIC) <-
  c("Gaussian",
    "Clayton",
    "alog",
    "alogm3",
    "mixt_gauss_alog",
    "mixt_clayton_alog")
cat(
  knitr::kable(
    t(AIC),
    align = "cccccc",
    col.names = paste0("$M_", 1:6, "$"),
    escape = FALSE,
    digits = 2,
    booktabs = TRUE,
    format = "latex"
  ),
  file = "../tables/TableS3.tex"
)



comp.gauge_fun <- list(
  Gaussian = fit.Gauss,
  Clayton = fit.Clayton,
  "asymmetric logistic" = fit.alog,
  "restricted asymmetric logistic" = fit.alog2,
  "mixture of Gaussian and asymmetric logistic" = fit.mixt_gauss_alog,
  "mixture of Clayton and asymmetric logistic" = fit.mixt_clayton_alog
)

# Print coefficients of best model
round(comp.gauge_fun[[which.min(AIC)]]$mle, 2)

save(comp.gauge_fun, file = "../outputs/comp.gauge_fun.RData")
