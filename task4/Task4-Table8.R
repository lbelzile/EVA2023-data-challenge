setwd(this.path::here())
library(chandwich)
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
####     APPROACH 2 - CONDITIONAL EXTREMES     ######
####     EXCHANGEABLE HEFFERNAN-TAWN MODEL     ######
#####################################################
# Quantile levels for threshold stability plots
qlevs <- seq(0.005, 0.05, by = 0.005)
data <- apply(unif, 2, qlaplace)
v1 <- qlaplace(1 - 1 / 300)
v2 <- qlaplace(1 - 12 / 300)
# Check transformation is correct
#stopifnot(isTRUE(all(rank(data[1,])== rank(Utopula[1,]))))

# Repeat the following steps for each cluster
for (k in 1:5) {
  cl <- k  #pick cluster
  clsize <- sum(hclust_od == cl)
  u1size <- sum((1:50)[hclust_od == cl] <= 25)
  # Create container for parameter estimates (alpha, beta)
  # and confidence intervals at each threshold
  param_mat <- array(data = NA,
                     dim = c(2, 3, length(qlevs)))
  for (qlev_ind in seq_along(qlevs)) {
    # For each quantile level
    qlev <- qlevs[qlev_ind]
    #####################################################
    ####  Heffernan-Tawn model estimation          ######
    #####################################################
    # Obtain point estimates by maximizing pseudo-likelihood
    opt <- Rsolnp::solnp(
      pars = c(0.2, 0.1, rep(0, 3)),
      fun = function(x) {
        -sum(eht_pll(
          par = x,
          thresh = 1 - qlev,
          data = data[, hclust_od == cl],
          type = "skewnorm"
        ))
      },
      LB = c(-1, rep(-1e3, 3), -1e3),
      UB = c(rep(1, 2), rep(1e3, 2), 1e3),
      control = list(trace = FALSE)
    )
    # Store point estimates
    param_mat[, 1, qlev_ind] <- opt$pars[1:2]

    if (opt$pars[1] < 0.98) {
      which.pars <- c("alpha", "beta")
      # Adjust the likelihood for repeated data
      adjprof <- adjust_loglik(
        loglik = eht_pll,
        p = 5,
        # number of parameters
        thresh = 1 - qlev,
        # marginal threshold
        data = data[, hclust_od == cl],
        mle = opt$pars,
        type = "skewnorm",
        par_names = c("alpha", "beta", "mu", "sigma", "kappa")
      )
      # Compute adjusted confidence intervals
      confints <- chandwich:::conf_intervals(adjprof,
                                             which_pars = which.pars)
      param_mat[(1:2)[c("alpha", "beta") %in% which.pars], -1, qlev_ind] <-
        confints$prof_CI
    }
  }
  # Parameter stability plots
  pdf(
    paste0("../figures/Task4-threshold_stability", cl, ".pdf"),
    width = 10,
    height = 4.5
  )
  par(mfrow = c(1, 2),
      mar = c(4, 4.5, 1, 1),
      bty = "l")
  matplot(
    x = 1 - qlevs,
    y = t(param_mat[1, , ]),
    type = "l",
    lty = c(1, 2, 2),
    col = 1,
    ylim = c(-0.2, 1),
    ylab = expression(alpha),
    xlab = "quantile level"
  )
  matplot(
    x = 1 - qlevs,
    y = t(param_mat[2, , ]),
    type = "l",
    lty = c(1, 2, 2),
    ylim = c(0, 1),
    col = 1,
    ylab = expression(beta),
    xlab = "quantile level"
  )
  dev.off()



  #####################################################
  ####  Extract residuals and estimate prob.     ######
  #####################################################
  # Extract parameters for the selected threshold

  probs <- matrix(NA, nrow = length(qlevs), ncol = 3)
  for (ql in seq_along(qlevs)) {
    alpha <- param_mat[1, 1, ql]
    beta <- param_mat[2, 1, ql]
    th <- 1 - qlevs[ql]
    residuals <- # use function to get residuals
      residuals_eht(
        alpha = alpha,
        beta = beta,
        data = data[, hclust_od == cl],
        thresh = th
      )

    # Test equality of distribution for residuals
    # print(energy::eqdist.etest(
    #   residuals$res,
    #   sizes = residuals$nexc,
    #   distance = FALSE,
    #   method = "original",
    #   R = 999))

    Zmin <- apply(residuals$res, 1, min)
    y0s <- qexp_level(Zmin, v1)
    # Compute Monte Carlo estimator
    log(mean(exp(-y0s)))
    # A more numerically stable version is here
    p2 <- logp(y0s)
    probs[ql, 3] <- p2
    # Using a permutation approach
    ngclust <- table(ifelse(which(hclust_od == cl) <= 25, 1, 2))
    #U2 has lower defense
    combos <- combn(m = ngclust[1], sum(ngclust))
    p1p_s <- vector("numeric", length = ncol(combos))
    for (i in seq_len(ncol(combos))) {
      clustid <- rep(2, sum(ngclust))
      clustid[combos[, i]] <- 1
      residuals2 <- # use function to get residuals
        residuals_eht(
          alpha = alpha,
          beta = beta,
          data = data[, hclust_od == cl],
          thresh = th,
          testIndep = FALSE,
          group = clustid
        )
      Zmin1 <- with(residuals2,
                    apply(res[, seq_len(ng - 1)], 1, min))
      Zmin2 <- with(residuals2,
                    apply(res[, -seq_len(ng - 1)], 1, min))
      # Plot residuals
      # plot(Zmin1, Zmin2)
      y0s1 <- qexp_level(Zmin1, v1)
      y0s2 <- qexp_level(Zmin2, v2)
      p1p_s[i] <- logp(pmax(y0s1, y0s2))
      # print(i)
    }
    p1p <- mean(p1p_s)
    probs[ql, 2] <- p1p
    # Using the identity of the columns
    residuals2 <- # use function to get residuals
      residuals_eht(
        alpha = alpha,
        beta = beta,
        data = data[, hclust_od == cl],
        thresh = th,
        testIndep = FALSE,
        group = ifelse(which(hclust_od == cl) <= 25, 1, 2)
      )
    Zmin1 <- with(residuals2,
                  apply(res[, seq_len(ng - 1)], 1, min))
    Zmin2 <- with(residuals2,
                  apply(res[, -seq_len(ng - 1)], 1, min))
    # Plot residuals
    # plot(Zmin1, Zmin2)
    v2 <- qlaplace(1 - 12 / 300)
    y0s1 <- qexp_level(Zmin1, v1)
    y0s2 <- qexp_level(Zmin2, v2)
    p1 <- logp(pmax(y0s1, y0s2))
    probs[ql, 1] <- p1
  }
  # Inconsistent estimates: the estimated
  #  probability is smaller than before...
  # but we know p1 < p2!
  # Culprit is difference in residuals
  # logp(qexp_level(with(residuals2,
  #                      apply(res, 1, min)), v1))
  # Probability of exceedance under independence
  p2_indep <- -clsize * log(300)
  p1_indep <- -u1size * log(300) - (clsize - u1size) * log(12 / 300)

  # Assign names for ease
  dimnames(param_mat) <- list(
    param = c("alpha", "beta"),
    est = c("estimate", "lower", "upper"),
    qlev = paste0(qlevs)
  )
  colnames(probs) <- c("p1", "p1p", "p2")
  rownames(probs) <- 1 - qlevs
  results <- list(
    prob = probs,
    indep = c(p1 = p1_indep, p2 = p2_indep),
    parameters = param_mat
  )
  save(results,
       file = paste0("../outputs/C4_results_cl", cl, ".RData"),
       version = 2)
}


# Load all results
interm_res <- list.files(pattern = "../outputs/C4_results_cl")
C4 <- matrix(nrow = 5, ncol = 2)
params <- matrix(nrow = 5, ncol = 2)
for (i in seq_along(interm_res)) {
  load(file = interm_res[i])
  C4[i, ] <- results$prob[4, 2:3]
  params[i, ] <- results$parameters[, 1, 4]
}

tab <- apply(rbind(
  table(hclust_od),
  sprintf("%.3f", t(C4[, 1])),
  sprintf("%.3f", t(C4[, 2])),
  paste0("(", apply(params, 1, function(x) {
    paste(sprintf("%.2g", x),  collapse = ", ")
  }), ")")
),
1:2, function(x) {
  paste0("$", x, "$")
})
rownames(tab) <- c(
  "cluster size",
  "$\\log(\\widehat{p}_1)$",
  "$\\log(\\widehat{p}_2)$",
  "$(\\widehat{\\alpha}, \\widehat{\\beta})$"
)
cat(
  knitr::kable(
    tab,
    format = "latex",
    col.names = paste0("$C_", 1:5, "$"),
    booktabs = TRUE,
    escape = FALSE,
    align = c("rrrrr")
  ),
  file = "../tables/Table8.tex"
)



# Return results for qlev = 0.02 (98%)
C4 <- exp(colSums(C4))
# save(C4, file = "../submission/AnswerC4.Rdata", version = 2)

# Plot profile likelihood (adjusted)
# for the Heffernan-Tawn model for AI clusters
for (cl in c(2, 3, 5)) {
  opt <- Rsolnp::solnp(
    pars = c(0.2, 0.1, rep(0, 3)),
    fun = function(x) {
      -sum(eht_pll(
        par = x,
        thresh = .98,
        data = data[, hclust_od == cl],
        type = "skewnorm"
      ))
    },
    LB = c(-1, rep(-1e3, 3), -1e3),
    UB = c(rep(1, 2), rep(1e3, 2), 1e3),
    control = list(trace = FALSE)
  )
  adjprof <- adjust_loglik(
    loglik = eht_pll,
    p = 5,
    # number of parameters
    thresh = 0.98,
    # marginal threshold
    data = data[, hclust_od == cl],
    mle = opt$pars,
    type = "skewnorm",
    par_names = c("alpha", "beta", "mu", "sigma", "kappa")
  )
  assign(
    x = paste0("prof_cl", cl),
    chandwich::conf_region(
      adjprof,
      which_pars = c("alpha", "beta"),
      type = "vertical",
      conf = 99
    )
  )
}

pdf("../figures/Task4-profile_cl.pdf",
    width = 12,
    height = 5)
par(mfrow = c(1, 3),
    mar = c(4, 4, 1, 1),
    bty = "l")
plot(
  prof_cl2,
  conf = c(50, 75, 90, 95, 99),
  xlab = expression(alpha),
  ylab = expression(beta),
  xlim = c(0, 0.2),
  ylim = c(0.2, 0.45),
  xaxs = "i",
  yaxs = "i"
)
points(params[2, 1], params[2, 2], pch = 4)
plot(
  prof_cl3,
  conf = c(50, 75, 90, 95, 99),
  xlab = expression(alpha),
  ylab = expression(beta),
  xlim = c(0.15, 0.45),
  ylim = c(0.35, 0.5),
  xaxs = "i",
  yaxs = "i"
)
points(params[3, 1], params[3, 2], pch = 4)
plot(
  prof_cl5,
  conf = c(50, 75, 90, 95, 99),
  xlab = expression(alpha),
  ylab = expression(beta),
  xlim = c(0.04, 0.15),
  ylim = c(0.1, 0.35),
  xaxs = "i",
  yaxs = "i"
)
points(params[5, 1], params[5, 2], pch = 4)
dev.off()
