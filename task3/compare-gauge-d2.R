setwd(this.path::here())
source("Task3-functions.R") ### to load the fit_geometric, quantR, and  pgumbel
data_coputopia <- read.csv("../data/Coputopia.csv")
data <- data_coputopia[, -c(1:2)]
data<- data[data[,3]<evd::qgumbel(0.5), 1:2]

library(geometricMVE)
thresh <- 0.95
########## Fit the model with different choices of gauge functions
####### Gaussian gauge function
fit.Gauss <- fit_geometric.2d(
  data = data,
  gauge.fun = gauge_gaussian,
  initial.val =  c(3, 0.5),
  lower.limit = c(0, -1),
  upper.limit = c(5, 1),
  nsim = 1e5L,
  thresh = thresh,
)

######  logistic type
fit.logistic  <- fit_geometric.2d(
  data = data,
  gauge.fun = gauge_rvad,
  initial.val =  c(3, 0.9),
  lower.limit = c(0, 0),
  upper.limit = c(5, 1),
  nsim = 1e5L,
  thresh = thresh,
)


######  inverse logistic type
fit.invlogistic <- fit_geometric.2d(
  data = data,
  gauge.fun = gauge_invlogistic,
  initial.val =  c(3, 0.9),
  lower.limit = c(0, 0),
  upper.limit = c(5, 1),
  nsim = 1e5L,
  thresh = thresh,
)


######  square-type gauge
fit.square <- fit_geometric.2d(
  data = data,
  gauge.fun = gauge_square,
  initial.val =  c(3, 0.9),
  lower.limit = c(0, 0),
  upper.limit = c(5, 1),
  nsim = 1e5L,
  thresh = thresh,
)


##### mixtures of  Gaussian and  logistics type gauge functions
fit.mixt_gauss_logistic <- fit_geometric.2d(
  data = data,
  add.gauge = TRUE,
  gauge1 = gauge_gaussian,
  gauge2 = gauge_rvad,
  initial.val = c(3, 0.5, 0.9),
  par1ind = 2,
  par2ind = 3,
  lower.limit = c(0, -1, 0),
  upper.limit = c(5, 1, 1),
  nsim = 1e5L,
  thresh = thresh,
)


##### mixtures of Clayton and asymmetric logistics gauge functions
fit.mixt_gauss_invlogistic <- fit_geometric.2d(
  data = data,
  add.gauge = TRUE,
  gauge1 = gauge_gaussian,
  gauge2 = gauge_invlogistic,
  initial.val = c(3, 0.5, 0.9),
  par1ind = 2,
  par2ind = 3,
  lower.limit = c(0, -1, 0),
  upper.limit = c(5, 1, 1),
  nsim = 1e5L,
  thresh = thresh,
)


comp.gauge_fun <- list(
  Gaussian = fit.Gauss,
  logistic = fit.logistic,
  "inverse logistic" = fit.invlogistic,
  "square logistic" = fit.square,
  "mixture of Gaussian and logistic" = fit.mixt_gauss_logistic,
  "mixture of Gaussian and inverse logistic" = fit.mixt_gauss_invlogistic
)

save(comp.gauge_fun, file = "../outputs/comp.gauge_fun2d.RData")

################## model comparison in terms of AIC values
AIC <-
  c(
    fit.Gauss$AIC,
    fit.logistic$AIC,
    fit.invlogistic$AIC,
    fit.square$AIC,
    fit.mixt_gauss_logistic$AIC,
    fit.mixt_gauss_invlogistic$AIC
  )
AIC

# 
# cat(
#   knitr::kable(
#     t(AIC),
#     align = "cccccc",
#     col.names = paste0("$M_", 1:6, "$"),
#     escape = FALSE,
#     digits = 2,
#     booktabs = TRUE,
#     format = "latex"
#   ),
#   file = "../tables/TableS1-dim2.tex"
# )
