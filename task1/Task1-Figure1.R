setwd(this.path::here())
library(mgcv)
library(ismev)
library(evgam)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpointdensity)

Amaurot <- read.csv("../data/Amaurot.csv")
data <- Amaurot |>
   dplyr::mutate(Season = as.factor(Season)) |>
   tidyr::drop_na()

model_thr_ald <- list(Y ~ Season + s(V1, bs = "cs") + s(V2, bs = "cs") +
                         s(V3, bs = "cs") + s(V4, bs = "cs") +
                         s(WindDirection , bs = "cc") + s(WindSpeed , bs = "cs") +
                         s(Atmosphere, bs = "cs"),
                      ~ Season + s(V1, bs = "cs") + s(V2, bs = "cs") +
                         s(V3, bs = "cs") + s(V4, bs = "cs") +
                         s(WindDirection , bs = "cc") + s(WindSpeed , bs = "cs") +
                         s(Atmosphere, bs = "cs"))

## fitting asymmetric Laplace for the threshold

# fit_thr_ald_0.95 <- evgam(model_thr_ald,
#                           data = data,
#                           family = "ald",
#                           ald.args = list(tau = 0.95))
fit_thr_ald_0.98 <- evgam(model_thr_ald,
                          data = data,
                          family = "ald",
                          ald.args = list(tau = 0.98))

# threshold_0.95 <- predict(fit_thr_ald_0.95)$location
threshold_0.98 <- predict(fit_thr_ald_0.98)$location

scatter.plotter <- function(X, Y, varname1, varname2, minv = NA){
   
   # exceed_0.95 <- which(Y > threshold_0.95)
   exceed_0.98 <- which(Y > threshold_0.98)
   
   # X_0.95 <- X[setdiff(exceed_0.95, exceed_0.98)]
   # Y_0.95 <- Y[setdiff(exceed_0.95, exceed_0.98)]
   X_0.98 <- X[exceed_0.98]
   Y_0.98 <- Y[exceed_0.98]
   p <- ggplot() +
      geom_hex(aes(x = X, y = Y), bins = 50, alpha = 0.5) +
      # geom_point(aes(x = X_0.95, y = Y_0.95), shape = "cross") +
      geom_point(aes(x = X_0.98, y = Y_0.98),
                 colour = "black", fill = "white", shape = 23) +
     colorspace::scale_fill_continuous_sequential(
       palette = "Purple-Yellow") + 
      labs(x = varname1,
           y = "") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      scale_x_continuous(limits = c(minv, NA),
                         expand = expansion(mult = c(0, 0.1))) +
      theme_classic() +
      theme(legend.position = "none")
   p
}

g1 <- with(data, scatter.plotter(X = V1, Y = Y, varname1 = expression(V[1])))
g2 <- with(data, scatter.plotter(X = V2, Y = Y, varname1 = expression(V[2])))
g3 <- with(data, scatter.plotter(X = V3, Y = Y, varname1 = expression(V[3])), minv = 0)
g4 <- with(data, scatter.plotter(X = V4, Y = Y, varname1 = expression(V[4])), minv = 0)
g5 <- with(data, scatter.plotter(X = WindSpeed,
                                 Y = Y, varname1 = "Wind speed"))
g6 <- with(data, scatter.plotter(X = Atmosphere,
                                 Y = Y, varname1 = "Atmosphere"))
g7 <-  with(data, scatter.plotter(X = WindDirection,
                                  Y = Y, varname1 = "Wind direction")) +
   scale_x_continuous(limits = c(-pi, pi),
                      breaks = seq(-pi, pi/2, by = pi/2),
                      labels =  c(expression(pi),
                                  expression(3*pi/2),
                                  expression(0),
                                  expression(pi/2)
                      )) +
   scale_y_continuous(breaks = NULL) +
   coord_polar(start = +pi/2, direction = -1) + 
  theme(legend.text = element_text(size = 15))
ggsave("../figures/Figure1.pdf",
       width = 10,
       height = 6,
       plot = (g1 | g2 | g3 | g4) / (g6 | g5 | g7),
       dpi = 300)
