setwd(this.path::here())
library(ggplot2)
library(patchwork)
library(mev)
source("Task4-functions.R")


# Load data and change column names to avoid duplicate labels

UtopulaU1_data <- read.csv(file = "../data/UtopulaU1.csv", header = TRUE)
colnames(UtopulaU1_data) <- paste0("I1_", c(paste0(0, 1:9), 10:25))
UtopulaU2_data <- read.csv(file = "../data/UtopulaU2.csv", header = TRUE)
colnames(UtopulaU2_data) <- paste0("I2_", c(paste0(0, 1:9), 10:25))
# Merge databases
Utopula <- data.frame(UtopulaU1_data, UtopulaU2_data)


#####################################################
####  Exploratory data analysis                ######
#####################################################
# Kendall tau correlation
cormat <- pcaPP::cor.fk(Utopula)
# relabel to remove "Y" from the matrix
colnames(cormat) <- rownames(cormat) <- 1:50

# Get order and pick observations from last cluster
clust_order <- corrplot::corrMatOrder(corr = cormat,
                                      order = 'hclust',
                                      hclust.method = "ward.D2")
# Compute hclust separately to sort them by island
hclust_od1 <- cutree(tree = hclust(as.dist(1 - cormat),
                                   method = "ward.D2"),
                     k = 5)
# Label-switching: use same cluster nomenclature as editorial
hclust_od <- dplyr::case_match(.x = hclust_od1,
                               3 ~ 1,
                               4 ~ 2,
                               5 ~ 3,
                               2 ~ 4,
                               1 ~ 5)
clust_order <- order(as.integer(paste0(hclust_od, c(paste0(0, 1:9), 11:50))))
# Transform to uniform scale
unif <- exp(-exp(-Utopula)) # Gumbel to uniform

qlev <- 0.05
data <- apply(unif, 2, qlaplace)
v1 <- qlaplace(1 - 1 / 300)
v2 <- qlaplace(1 - 12 / 300)
res <- list()
for(cl in 1:5){
  clsize <- sum(hclust_od == cl)
  u1size <- sum((1:50)[hclust_od == cl] <= 25)
  # Create container for parameter estimates (alpha, beta)
  # and confidence intervals at each threshold
  # Obtain point estimates by maximizing pseudo-likelihood
  opt <- Rsolnp::solnp(
    pars = c(0.2, 0.1, rep(0, 2)),
    fun = function(x) {
      -sum(eht_pll(
        par = x,
        thresh = 1 - qlev,
        data = data[, hclust_od == cl],
        type = "norm"
      ))
    },
    LB = c(-1, rep(-1e3, 3)),
    UB = c(rep(1, 2), rep(1e3, 2)),
    control = list(trace = FALSE)
  )
  alpha <- opt$pars[1]
  beta <- opt$pars[2]
  residuals <- # use function to get residuals
    residuals_eht(
      alpha = alpha,
      beta = beta,
      data = data[, hclust_od == cl],
      thresh = 1-qlev
    )

  res[[cl]] <- c(residuals$res)
  plot(density(res[[cl]], bw = "SJ"))

  fit_sn <- sn::selm(res[[cl]] ~ 1, family = "SN")
  (coef(fit_sn))
  fit_st <- sn::selm(res[[cl]] ~ 1, family = "ST")
  (coef(fit_st))
}

# Illustrate that the function is monotone for given (alpha, beta)
zmin <- min(residuals$res[1,])
v1 <- qlaplace(1 - 1 / 300)
funy0 <- function(y0){
  zmin - (v1 -alpha*y0)/y0^beta
}

nres <- sapply(res, length)
sres <- lapply(res, scale)
sapply(res, moments::kurtosis)
sapply(res, moments::skewness)

g1 <- ggplot() +
  stat_function(fun = dnorm, xlim = c(-10, 10),
                col = "black", n = 1e3L) +
  geom_density(data = data.frame(
    z = unlist(sres),
    cluster = factor(unlist(sapply(1:5, function(i){
      rep(i, length.out = nres[i])})))),
    mapping = aes(x = z, color = cluster, fill = cluster),
    bw = "SJ", alpha = 0.1) +
  scale_x_continuous(limits = c(-5, 8)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(y = "", subtitle = "density", x = "scaled residuals") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.5))

g2 <- ggplot() +
  stat_function(fun = funy0, xlim = c(3, 150)) +
  geom_vline(xintercept = uniroot(funy0, interval = c(3, 200))$root,
             linetype = "dashed") +
  labs(y = "", x = expression(y[0]),
       subtitle = expression(z[min] - (v[1]-alpha*y[0])/y[0]^beta)) +
  theme_classic()
g1 + g2
ggsave("../figures/FigureS1.pdf",
       width = 8,
       height = 4)
