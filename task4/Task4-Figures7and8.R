setwd(this.path::here())
library(mev)
library(ggplot2)
library(patchwork)
source("Task4-functions.R")
generateFigures <- FALSE



# Load data and change column names to avoid duplicate labels

UtopulaU1_data <-
  read.csv(file = "../data/UtopulaU1.csv", header = TRUE)
colnames(UtopulaU1_data) <- paste0("I1_", c(paste0(0, 1:9), 10:25))
UtopulaU2_data <-
  read.csv(file = "../data/UtopulaU2.csv", header = TRUE)
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
table(hclust_od)
clust_order <-
  order(as.integer(paste0(hclust_od, c(paste0(0, 1:9), 10:50))))
# Transform to uniform scale
unif <- exp(-exp(-Utopula)) # Gumbel to uniform

#-------------------------------------------------------------------------------
# Figure 7

if (generateFigures) {
  tailcor <- array(dim = c(50, 50, 3))
  # Compute tail correlation for each pair
  for (i in 2:50) {
    for (j in 1:(i - 1)) {
      chi_ij <- try(mev::taildep(
        data = unif[, c(j, i)],
        depmeas = "chi",
        method = list(chi = "betacop"),
        confint = "lrt",
        u = 0.975,
        plot = FALSE
      )$chi,
      silent = TRUE)
      if (!inherits(chi_ij, "try-error") & length(chi_ij) > 0) {
        tailcor[j, i,] <- tailcor[i, j,] <- chi_ij
      }
    }
  }
  # Extract point estimates
  tailcormat <- tailcor[, , 1]
  diag(tailcormat) <- 1
  # Keep same ordering/names as previous plot
  colnames(tailcormat) <- rownames(tailcormat) <- rownames(cormat)
  grayp <- grey(c(rep(1, 50), 
                  1-c(seq(0,0.1, length.out = 15), 
                 seq(0.1,1, length.out = 36))))
  # Correlation plots
  pdf("../figures/Figure7.pdf",
      width = 11,
      height = 5)
  par(mfrow = c(1, 2), mar = c(0, 7, 1, 5))
  corr <- apply(cormat, 1:2, function(x) {
    pmax(0, x)
  })[clust_order, clust_order]
  corrplot::corrplot(corr,
    order = 'original',
    tl.col = "black",
    tl.cex = 0.7,
    cl.cex = 1.3,
    col = grayp,
    col.lim = c(0, 1),
    addgrid.col = "white",
    tl.pos = "d",
    method = 'shade'
  )

  corrplot::corrplot(
    tailcormat[clust_order, clust_order],
    tl.col = "black",
    order = "original",
    tl.cex = 0.7,
    cl.cex = 1.3,
    col = grayp,
    col.lim = c(0, 1),
    addgrid.col = "white",
    tl.pos = "d",
    method = 'shade'
  )
  dev.off()
}

#-------------------------------------------------------------------------------
# Figure 8
useq <- seq(0.85, 0.995, by = 0.005)
nu <- length(useq)
# Compute confidence intervals via nonparametric bootstrap
B <- 1000L
set.seed(2024)
chi <- array(0, dim = c(5, nu, B))
eta <- array(0, dim = c(5, nu, B))
for(b in seq_len(B)){
  for (cl in seq_len(5)) {
  # Compute tail correlation for each pair
  ncombo <- t(combn(which(hclust_od == cl), 2))
  num_na_eta <- num_na_chi <- rep(0, nu)
  if(b == 1L){
    index <- seq_len(nrow(unif))
  } else{
    # resample indices with replacement
    index <- sample.int(nrow(unif), replace = TRUE)
  }
  for (i in seq_len(nrow(ncombo))) {
    tails <- try(mev::taildep(
      data = unif[index, ncombo[i,]],
      depmeas = c("chi", "eta"),
      method = list(chi = "emp",
                    eta = "emp"),
      confint = "wald",
      # Don't estimate the margins
      empirical.transformation = FALSE,
      u = useq,
      plot = FALSE
    ),
    silent = TRUE)

    if (!inherits(tails, "try-error")) {
      cand_chi <- !is.na(tails$chi[, 1])
      num_na_chi <- num_na_chi + cand_chi
      chi[cl, cand_chi,b] <-
        chi[cl, cand_chi,b] + tails$chi[cand_chi, 1]
      cand_eta <- !is.na(tails$eta[, 1])
      num_na_eta <- num_na_eta + cand_eta
      eta[cl, cand_eta,b] <-
        eta[cl, cand_eta,b] + tails$eta[cand_eta, 1]
    }
  }
  chi[cl,,b] <- chi[cl,,b] / num_na_chi
  eta[cl,,b] <- eta[cl,,b] / num_na_eta
  }
  print(b)
}



# Plot chi and eta
if (generateFigures) {

  depregime <- data.frame(
    u = rep(useq, each = 5),
    chi = c(chi[,,1]),
    eta = c(eta[,,1]),
    chi_lower = c(apply(chi, 1:2, quantile, 0.025)),
    chi_upper = c(apply(chi, 1:2, quantile, 0.975)),
    eta_lower = c(apply(eta, 1:2, quantile, 0.025)),
    eta_upper = c(apply(eta, 1:2, quantile, 0.975)),
    cluster = factor(rep(1:5, length.out = 5 * nu))
  )
  g1 <- ggplot(data = depregime,
               mapping = aes(
                 x = u,
                 y = chi,
                 ymin = chi_lower,
                 ymax = chi_upper,
                 group = cluster,
                 color = cluster
               )) +
    geom_pointrange(size = 0.2) +
    geom_line(linewidth = 1.5) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_x_continuous(
      expand = c(0.01, 0.015),
      limits = c(0.85, 1),
      breaks = c(0.85, 0.9, 0.95, 1),
      labels = c("0.85", "0.9", "0.95", "1")
    ) +
    labs(x = "Quantile level",
         y = "",
         subtitle = expression(chi)) +
    scale_color_viridis_d(name = "Cluster") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20))
  g2 <- ggplot(data = depregime,
               mapping = aes(
                 x = u,
                 y = eta,
                 ymin = eta_lower,
                 ymax = eta_upper,
                 group = cluster,
                 color = cluster
               )) +
    geom_pointrange(size = 0.2) +
    geom_line(linewidth = 1.5) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_x_continuous(
      expand = c(0.01, 0.01),
      limits = c(0.85, 1),
      breaks = c(0.85, 0.9, 0.95, 1),
      labels = c("0.85", "0.9", "0.95", "1")
    ) +
    scale_color_viridis_d(name = "Cluster") +
    labs(x = "Quantile level",
         y = "",
         subtitle = expression(eta)) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20))

  pdf("../figures/Figure8.pdf",
      width = 12,
      height = 6)
  g1 + g2 + plot_layout(guides = "collect") &
    theme(text = element_text(size = 20),
    legend.position = "bottom"
    )
  dev.off()
}
