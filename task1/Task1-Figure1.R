## Load libraries
setwd(this.path::here())
library(mev)
library(circular)
library(changepoint)
library(energy)
library(ggplot2)
library(patchwork)
## Load data
data <- read.csv(file = "../data/Amaurot.csv")
data <- data |>
   dplyr::mutate(Season = as.integer(as.factor(Season)))
Y <- data$Y

# Testing for independence between groups
set.seed(2023)
groups <- list(
   c("V1", "V2"),
   c("V3", "Season"),
   "V4",
   c("WindSpeed", "xWindDirection", "yWindDirection"),
   "Atmosphere"
)
dataTests <- data
dataTests$Season <- ifelse(data$Season == "S1", 0, 1)
dataTests$xWindDirection <- cos(dataTests$WindDirection)
dataTests$yWindDirection <- sin(dataTests$WindDirection)
indep_tests <- replicate(
   n = 10,
   expr = {
       apply(combn(length(groups), m = 2), 2, function(inds) {
          data_mat <-
             as.matrix(na.omit(dataTests[, unlist(c(groups[inds[1]], groups[inds[2]]))]))
          # Subsample because otherwise there is a memory allocation error
          obsnum <-
             sample.int(n = nrow(data_mat), size = 2500)
          # Extract the length of the first block
          l <- length(unlist(groups[inds[1]]))
          # Energy test of independence with distance covariance
          test <-
             energy::indep.test(x = data_mat[obsnum, 1:l],
                                y = data_mat[obsnum, -(1:l)],
                                R = 399) # number of bootstrap replications
          test$p.value
       })
    })
# reject independence test between Wind and Atmosphere


## Fit changepoint model
changepoint <-
   changepoint::cpt.meanvar(data = na.omit(data$WindDirection),
                            method = "AMOC")
# Manually identify threshold
changepoint@cpts[1] <- 8200
cpos <-
   changepoint@cpts[1] + sum(which(is.na(data$WindDirection)) < changepoint@cpts[1])

deg(circular::mean.circular(data$WindDirection[1:cpos], na.rm = TRUE))

deg(circular::mean.circular(data$WindDirection[-(1:cpos)], na.rm = TRUE))

# Proportion of 1 vs 2
prop_wind <- cpos / length(na.omit(data$WindDirection))
rwinddummy <- function(n) {
   2 - rbinom(n = n, size = 1, prob = prop_wind)
}


# Fit circular density estimator
k1 <- circular::density.circular(
   x =
      circular::as.circular(
         na.omit(data$WindDirection[1:cpos]),
         modulo = "2pi",
         units = "radian"
      ),
   from = -pi,
   to = pi,
   bw = 75
)
k2 <- circular::density.circular(
   x =
      circular::as.circular(
         na.omit(data$WindDirection[-(1:cpos)]),
         modulo = "2pi",
         units = "radian"
      ),
   from = -pi,
   to = pi,
   bw = 75
)
dens <- cbind(k1$x, k1$y, k2$y)

predict_ang <- function(angle, k1, k2) {
   k1dens <- sapply(angle, function(ang) {
      k1$y[which.min((k1$x - ang) %% (2 * pi))]
   })
   k2dens <- sapply(angle, function(ang) {
      k2$y[which.min((k2$x - ang) %% (2 * pi))]
   })
   apply(cbind(k1dens, k2dens), 1, which.min)
}
# Example of prediction
angle <- seq(0, 2 * pi, length.out = 101)
predict_ang(angle, k1, k2)
pi_scales <-
   scales::math_format(
      .x * pi,
      format = function(x)
         x / pi
   )
g1 <- ggplot(
   data = data.frame(
      x = 1:21000 / 300,
      y = data$WindDirection %% (2 * pi),
      cluster = factor(c(rep(1, cpos), rep(0, 21000 - cpos)))
   ),
   mapping = aes(x = x, y = y,
                 group = cluster)
) +
   # ggpointdensity::geom_pointdensity() +
   geom_hex(bins = 40) +
   geom_vline(xintercept = cpos / 300) +
   labs(x = "Year", y = "", subtitle = "Wind direction (radians) ") +
   viridis::scale_fill_viridis() +
   scale_x_continuous(limits = c(0, 70), expand = c(0, 0)) +
   scale_y_continuous(
      limits = c(0, 2 * pi),
      expand = c(0, 0),
      breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
      labels = pi_scales
   ) +
   theme_classic() +
   theme(legend.position = "none")

g2 <- ggplot(
   data = data.frame(
      x = as.numeric(c(k1$x, k2$x)),
      y = as.numeric(c(k1$y, k2$y)),
      group = factor(c(rep(1, length(
         k1$x
      )),
      rep(2, length(
         k2$x
      ))))
   ),
   mapping = aes(
      x = x,
      y = y,
      group = group,
      col = group
   )
) +
   geom_line() +
   scale_color_grey() +
   scale_x_continuous(
      limits = c(0, 2 * pi),
      expand = c(0, 0),
      breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
      labels = pi_scales
   )  +
   scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
   labs(x = "Angles (in radians)", y = "") +
   theme_classic() +
   theme(legend.position = "none")
pdf("../figures/Figure1.pdf",
    width = 7,
    height = 3.5)
g1 + g2
dev.off()
