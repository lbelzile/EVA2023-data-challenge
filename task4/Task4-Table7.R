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
hclust_od1 <- cutree(tree = hclust(as.dist(1 - cormat), method = "ward.D2"), k = 5)
# Label-switching: use same cluster nomenclature as editorial
hclust_od <- dplyr::case_match(.x = hclust_od1, 3 ~ 1, 4 ~ 2, 5 ~ 3, 2 ~ 4, 1 ~ 5)
clust_order <- order(as.integer(paste0(hclust_od, c(
  paste0(0, 1:9), 10:50
))))
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

  s1 <- c(rep(lev_par, u1size), rep(lev_par2, clsize - u1size))
  s2 <- rep(lev_par, clsize)
  # Compute profile log likelihood
  # Huesler-Reiss model and logistic models

  fit_br <- matrix(nrow = length(thresh_seq), ncol = 11L)
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
  fit_log_min <- matrix(ncol = 8L, nrow = length(thresh_seq))
  colnames(fit_log_min) <- c("mle",
                             "loglik",
                             "lower",
                             "upper",
                             "nexc",
                             "chi",
                             "logp1",
                             "logp2")
  fit_neglog_min <- fit_hr_min <- fit_log_min
  for (th in seq_along(thresh_seq)) {
    maxCl <- apply(Cl, 1, max)
    thresh <- quantile(maxCl, thresh_seq[th])
    data <- Cl[maxCl > thresh, ]
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
    # neglog_ucens <- fit_mgp(
    #   data = data,
    #   thresh = thresh,
    #   mthresh = mthresh,
    #   model = "neglog",
    #   confint = "lr",
    #   type = "uncensored"
    # )
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
    chi_log <- expme_log_min(u = rep(1, clsize), par = log_cens$mle)
    chi_br <- expme_BR_min(z = rep(1, clsize),
                           Lambda = Lmat_exch(par = br_cens$mle, D = clsize))

    thresh_par <- mev::qgp(
      p = thresh,
      loc = 1,
      scale = 1,
      shape = 1
    )
    Vu_br <- expme(
      z = rep(thresh_par, clsize),
      par = list(Lambda = Lmat_exch(par = br_cens$mle, D = clsize)),
      model = "br",
      method = "TruncatedNormal"
    )
    vXiu_p1_br <- expme_BR_min(z = s1,
                               Lambda = Lmat_exch(par = br_cens$mle, D = clsize))
    vXiu_p2_br <- expme_BR_min(z = s2,
                               Lambda = Lmat_exch(par = br_cens$mle, D = clsize))
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
      log(vXiu_p1_log) - log(Vu_log) + log(prob_under),
      log(vXiu_p2_log) - log(Vu_log) + log(prob_under)
    )

    # Try out with a different risk functional - minimum is large
    minCl <- apply(Cl, 1, min)
    thresh <- quantile(minCl, probs = thresh_seq[th])
    data <- Cl[minCl > thresh, ]
    # Fit three parametric models
    log_min <- fit_mingp(data = data,
                         thresh = thresh,
                         family = "log")
    # We get the parameter on the boundary of the parameter space
    neglog_min <- fit_mingp(data = data,
                            thresh = thresh,
                            family = "neglog")
    hr_min <- fit_mingp(
      data = data,
      thresh = thresh,
      family = "hr",
      lb = 0.01
    )
    # Compute probability of exceedance and chi

    chi_logmin <- expme_log_min(u = rep(1, clsize), par = log_min['mle'])
    vXiu_p1_logmin <- log(expme_log_min(s1, par = log_min['mle'])) -
      log(expme_log_min(rep(1 / (1 - thresh), ncol(data)),
                        par = log_min['mle'])) +
      log(log_min['nexc'] / nrow(Cl))
    vXiu_p2_logmin <- log(expme_log_min(s2, par = log_min['mle'])) -
      log(expme_log_min(rep(1 / (1 - thresh), ncol(data)),
                        par = log_min['mle'])) +
      log(log_min['nexc'] / nrow(Cl))
    fit_log_min[th, ] <- c(log_min, chi_logmin,
                           vXiu_p1_logmin, vXiu_p2_logmin)
    chi_neglogmin <- expme_neglog_min(u = rep(1, clsize),
                                      par = neglog_min['mle'])
    vXiu_p1_neglogmin <- log(expme_neglog_min(s1, par = neglog_min['mle'])) -
      log(expme_neglog_min(rep(1 / (1 - thresh), ncol(data)),
                           par = neglog_min['mle'])) +
      log(neglog_min['nexc'] / nrow(Cl))
    vXiu_p2_neglogmin <- log(expme_neglog_min(s2, par = neglog_min['mle'])) -
      log(expme_neglog_min(rep(1 / (1 - thresh), ncol(data)),
                           par = neglog_min['mle'])) +
      log(neglog_min['nexc'] / nrow(Cl))
    fit_neglog_min[th, ] <- c(neglog_min,
                              chi_neglogmin,
                              vXiu_p1_neglogmin,
                              vXiu_p2_neglogmin)

    Lambda <- Lmat_exch(par = hr_min['mle'], D = clsize)
    chi_hrmin <- expme_BR_min(z = rep(1, clsize), Lambda = Lambda)
    vXiu_p1_hrmin <- log(expme_BR_min(s1, Lambda)) -
      log(expme_BR_min(rep(1 / (1 - thresh), ncol(data)), Lambda)) +
      log(hr_min['nexc'] / nrow(Cl))
    vXiu_p2_hrmin <- log(expme_BR_min(s2, Lambda)) -
      log(expme_BR_min(rep(1 / (1 - thresh), ncol(data)), Lambda)) +
      log(hr_min['nexc'] / nrow(Cl))
    fit_hr_min[th, ] <- c(hr_min, chi_hrmin, vXiu_p1_hrmin, vXiu_p2_hrmin)
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

  results <- list(br = fit_br, log = fit_log, logmin = fit_log_min,
                  neglogmin = fit_neglog_min, hrmin = fit_hr_min)

  save(
    results,
    file = paste0("../outputs/Task4-results_AD_cl", cl, ".RData"),
    version = 2
  )



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
load("../outputs/Task4-results_AD_cl1.RData")
results_C1 <- results
load("../outputs/Task4-results_AD_cl4.RData")
results_C4 <- results

thind <- 11
res_AD <- cbind(
  with(results_C1, c(br[thind, c("coef", "logp1", "logp2", "chi")])),
  with(results_C1, c(hrmin[thind, c("mle", "logp1", "logp2", "chi")])),
  with(results_C1, c(log[thind, c("coef", "logp1", "logp2", "chi")])),
  with(results_C1, c(logmin[thind, c("mle", "logp1", "logp2", "chi")])),
  with(results_C4, c(br[thind, c("coef", "logp1", "logp2", "chi")])),
  with(results_C4, c(hrmin[thind, c("mle", "logp1", "logp2", "chi")])),
  with(results_C4, c(log[thind, c("coef", "logp1", "logp2", "chi")])),
  with(results_C4, c(logmin[thind, c("mle", "logp1", "logp2", "chi")]))
)
fres_AD <- apply(res_AD, 1:2, function(x) {
  paste0("$", sprintf("%1.3f", x), "$")
})
rownames(fres_AD) <- c(
  "coefficients",
  "$\\log \\widehat{p}_{1,c}$",
  "$\\log \\widehat{p}_{2,c}$",
  "$\\widehat{\\chi}_{|C_k|}$"
)
cat(
  print(
    xtable::xtable(fres_AD),
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    include.colnames = FALSE,
    booktabs = TRUE
  ),
  file = "../tables/Table7.tex"
)
