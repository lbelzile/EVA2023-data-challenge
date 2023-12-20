setwd(this.path::here())
library(knitr)
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



#####################################################
####     APPROACH 3 - ONLY FOR AD CLUSTERS     ######
####       MULTIVARIATE GENERALIZED PARETO     ######
#####################################################
######################################################
# Asymptotically dependent models
# with exchangeable dependence structure
######################################################
for (k in 1:2) {
  cl <- c(1, 4)[k]
  mthresh <- 0.5 # marginal threshold for censoring
  thresh_seq <- seq(0.90, 0.99, by = 0.005)
  # risk threshold for max
  clsize <- sum(hclust_od == cl)
  u1size <- sum((1:50)[hclust_od == cl] <= 25)
  Cl <- unif[, hclust_od == cl]

  par_seq <- seq(0.2, 0.95, by = 0.02)
  alpha_seq <- seq(0.2, 0.95, by = 0.01)

  lev_par <- mev::qgp(
    p = 1 / 300,
    loc = 1,
    scale = 1,
    shape = 1,
    lower.tail = FALSE
  )
  lev_par2 <- mev::qgp(
    p = 12 / 300,
    loc = 1,
    scale = 1,
    shape = 1,
    lower.tail = FALSE
  )

  s1 <- c(rep(lev_par, u1size),
          rep(lev_par2, clsize - u1size))
  s2 <- rep(lev_par, clsize)
  # Compute profile log likelihood
  # Huesler-Reiss model and logistic models

  fit_br <- matrix(nrow = length(thresh_seq),
                   ncol = 11L)
  colnames(fit_br) <- c(
    "threshold",
    "coef",
    "lower",
    "upper",
    "coefu",
    "loweru",
    "upperu",
    "maxll",
    "chi",
    "logp1",
    "logp2"
  )
  fit_log <- fit_br

  for (th in seq_along(thresh_seq)) {
    thresh <- thresh_seq[th]
    data <- Cl[apply(Cl, 1, max) > thresh, ]
    log_cens <- fit_mgp(
      data = data,
      thresh = thresh,
      mthresh = mthresh,
      model = "log",
      confint = "lr",
      type = "censored"
    )
    log_ucens <- fit_mgp(
      data = data,
      thresh = thresh,
      mthresh = mthresh,
      model = "log",
      confint = "lr",
      type = "uncensored"
    )
    br_cens <- fit_mgp(
      data = data,
      thresh = thresh,
      mthresh = mthresh,
      model = "br",
      confint = "wald",
      type = "censored"
    )
    br_ucens <- fit_mgp(
      data = data,
      thresh = thresh,
      mthresh = mthresh,
      model = "br",
      confint = "wald",
      type = "uncensored"
    )
    chi_log <- 2 - mev::expme(
      z = rep(1, 2),
      par = list(alpha = log_cens$mle),
      model = "log"
    )
    chi_br <- 2 - mev::expme(
      z = rep(1, 2),
      par = list(Lambda = Lmat_exch(par = br_cens$mle,
                                    D = 2)),
      model = "br"
    )

    thresh_par <- mev::qgp(
      p = thresh,
      loc = 1,
      scale = 1,
      shape = 1
    )
    Vu_br <- expme(
      z = rep(thresh_par, clsize),
      par = list(Lambda = Lmat_exch(par = br_cens$mle,
                                    D = clsize)),
      model = "br",
      method = "TruncatedNormal"
    )
    vXiu_p1_br <- expme_BR_min(z = s1,
                               Lambda = Lmat_exch(par = br_cens$mle,
                                                  D = clsize))
    vXiu_p2_br <- expme_BR_min(z = s2,
                               Lambda = Lmat_exch(par = br_cens$mle,
                                                  D = clsize))
    prob_under <- mean(rowSums(Cl > thresh) >= 1)

    # The same, but with the logistic model
    Vu_log <- expme(
      z = rep(thresh_par, clsize),
      par = list(alpha = log_cens$mle),
      model = "log"
    )
    vXiu_p1_log <- expme_log_min(s1, par = log_cens$mle)
    vXiu_p2_log <- expme_log_min(s2, par = log_cens$mle)

    fit_br[th, ] <- c(
      thresh,
      br_cens$mle,
      br_cens$confint,
      br_ucens$mle,
      br_ucens$confint,
      br_cens$maxll,
      chi_br,
      log(vXiu_p1_br / Vu_br * prob_under),
      log(vXiu_p2_br / Vu_br * prob_under)
    )
    fit_log[th, ] <- c(
      thresh,
      log_cens$mle,
      log_cens$confint,
      log_ucens$mle,
      log_ucens$confint,
      log_cens$maxll,
      chi_log,
      log(vXiu_p1_log / Vu_log * prob_under),
      log(vXiu_p2_log / Vu_log * prob_under)
    )
  }

  ## Check the formula via Monte Carlo for a toy example
  # d <- 10
  # test <- mev::rparp(n = 1e5,
  #            riskf = "max",
  #            d = d,
  #            param = opt_par_log$mle,
  #            model = "log")
  # mean(apply(test, 1, function(x){min(x) > 1}))
  # expme_log_min(rep(1, d), par = opt_par_log$mle) / d
  # expme(z = rep(1, d), par = list(alpha = 1/opt_par_log$mle),
  #       model = "log")

  ## Check the formula via Monte Carlo for a toy example
  # d <- 5
  # test <- mev::rparp(n = 1e5,
  #            riskf = "max",
  #            d = d,
  #            sigma = Lmat_exch(par = opt_par_br$mle, D = d),
  #            model = "hr")
  # mean(apply(test, 1, function(x){min(x) > 1}))
  # expme_BR_min(rep(2, d),
  #  Lambda = Lmat_exch(par = opt_par_br$mle, D = d)) /
  # expme(z = rep(2, d),par = list(
  # Lambda = Lmat_exch(par = opt_par_br$mle, D = d)),
  #       model = "br", method = "TruncatedNormal")

  results <- list(br = fit_br,
                  log = fit_log)

  save(results,
       file = paste0("../outputs/Task4-results_AD_cl", cl, ".RData"),
       version = 2)

  pdf(
    file = paste0("../figures/Task4-AD_cluster_tstab", cl, ".pdf"),
    width = 10,
    height = 5
  )
  par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
  matplot(
    results$log[, 'threshold'],
    results$log[, c('coefu', 'loweru', 'upperu')],
    type = "l",
    lty = c(1, 2, 2),
    bty = "l",
    ylab = expression(alpha),
    xlab = 'quantile level',
    col = c("grey")
  )
  matplot(
    results$log[, 'threshold'],
    results$log[, c('coef', 'lower', 'upper')],
    add = TRUE,
    type = "l",
    lwd = 1,
    lty = c(1, 2, 2),
    col = c("black")
  )
  matplot(
    results$br[, 'threshold'],
    results$br[, c('coef', 'lower', 'upper')],
    type = "l",
    lty = c(1, 2, 2),
    bty = "l",
    ylab = expression(lambda),
    xlab = 'quantile level',
    col = c("grey")
  )
  matplot(
    results$br[, 'threshold'],
    results$br[, c('coefu', 'loweru', 'upperu')],
    type = "l",
    add = TRUE,
    lty = c(1, 2, 2),
    col = c("black")
  )
  dev.off()
}




# Load results and print a table for LaTeX
load("../outputs/results_AD_cl1.RData")
results_C1 <- results
load("../outputs/results_AD_cl4.RData")
results_C4 <- results

thind <- 11
res_AD <- cbind(
  with(results_C1, c(br[thind, c("coef", "logp1", "logp2", "chi")])),
  with(results_C4, c(br[thind, c("coef", "logp1", "logp2", "chi")])),
  with(results_C1, c(log[thind, c("coef", "logp1", "logp2", "chi")])),
  with(results_C4, c(log[thind, c("coef", "logp1", "logp2", "chi")]))
)
fres_AD <- apply(res_AD, 1:2, function(x) {
  paste0("$", sprintf("%1.3f", x), "$")
})
rownames(fres_AD) <- c(
  "coefficients",
  "$\\log \\widehat{p}_1$",
  "$\\log \\widehat{p}_2$",
  "$\\widehat{\\chi}$"
)
cat(
  print(
    xtable::xtable(fres_AD),
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    include.colnames = FALSE,
    booktabs = TRUE
  ),
  file = "../tables/Table9.tex"
)
