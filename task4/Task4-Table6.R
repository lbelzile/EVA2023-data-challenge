setwd(this.path::here())
library(chandwich)
library(mev)
library(ggplot2)
library(patchwork)
library(reshape)
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

# Transform to uniform scale
unif <- exp(-exp(-Utopula)) # Gumbel to uniform

#####################################################
####     APPROACH 1 - EMPIRICAL EXTRAPOLATION  ######
######         A LA LEDFORD-TAWN               ######
#####################################################
# Compute eta coefficient for all values in the cluster
# at a lower threshold and use Ledford and Tawn (1997) to extrapolate
# We do the extrapolation in exponential scale
# Unequal approach with an extension of Wadsworth and Tawn (2013)
prob1_cl <- eta_p1 <-  numeric(5)
prob2_cl <- eta_p2 <- numeric(5)
qlev <- 0.985 # leaves 200 observations
s1_exp <- qexp(1 - 1 / 300)
s2_exp <- qexp(1 - 12 / 300)
edata <- -log(1 - unif)

for (cl in 1:5){
  edata_cl <- edata[, which(hclust_od == cl)]
  clsize <- ncol(edata_cl)
  u1size <- sum((1:50)[hclust_od == cl] <= 25)

  r2 <- tail_emp_exchange(qlev = qlev, edata = edata_cl)
  r1 <- tail_emp_exchange(qlev = qlev, edata = edata_cl,
                          group = c(rep(1, u1size),
                                    rep(2, clsize - u1size)))
  prob1_cl[cl] <- r1$prob
  prob2_cl[cl] <- r2$prob
  eta_p1[cl] <- r1$eta
  eta_p2[cl] <- r2$eta
}
log(prob2_cl)
log(prob1_cl)
# Check we respect ordering constraints
isTRUE(all(log(prob1_cl) > log(prob2_cl)))

# Print table with results
semipar_results <- rbind(prob1_cl, prob2_cl)
semipar_tbl <- apply(rbind(log(prob1_cl), log(prob2_cl)),
                     1:2, function(x) {
                       paste("$", sprintf("%.3f", x), "$")
                     })
rownames(semipar_tbl) <- c("$\\log \\widehat{p}_{1,a}$",
                           "$\\log \\widehat{p}_{2,a}$")


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
for (cl in 1:5) {
  clsize <- sum(hclust_od == cl)
  u1size <- sum((1:50)[hclust_od == cl] <= 25)
  # Create container for parameter estimates (alpha, beta)
  # and confidence intervals at each threshold
  param_mat <- param_mat_u <- array(data = NA,
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
      fun = function(par, ...) {
        -sum(eht_pll(
          par = par,
          thresh = 1 - qlev,
          data = data[, hclust_od == cl],
          type = "skewnorm",
          constrain = FALSE
        ))
      },
      ineqfun = function(par, ...){
        optim_constr_ht(
          par = par,
          constants = data_constr(
            par = par,
            thresh = 1 - qlev,
            data = data[, hclust_od == cl]))
      },
      ineqLB = 0,
      ineqUB = Inf,
      LB = c(-1, rep(-1e3, 3), -1e3),
      UB = c(rep(1, 2) + 1e-10, rep(1e3, 2), 1e3),
      control = list(trace = FALSE)
    )
    opt_unconst <- Rsolnp::solnp(
      pars = c(0.2, 0.1, rep(0, 3)),
      fun = function(x) {
        -sum(eht_pll(
          par = x,
          thresh = 1 - qlev,
          data = data[, hclust_od == cl],
          type = "skewnorm",
          constrain = FALSE
        ))
      },
      LB = c(-1, rep(-1e3, 3), -1e3),
      UB = c(rep(1, 2) + 1e-10, rep(1e3, 2), 1e3),
      control = list(trace = FALSE)
    )
    # Fit the model with asymptotic dependence, because
    # the ridge region won't be visited by the optimization algorithm
    opt2 <- Rsolnp::solnp(
      pars = c(rep(0, 3)),
      fun = function(x) {
        -sum(eht_pll(
          par = c(0.9999,0, x),
          thresh = 1 - qlev,
          data = data[, hclust_od == cl],
          type = "skewnorm",
          constrain = TRUE
        ))
      },
      LB = rep(-1e3, 3),
      UB = rep(1e3, 3),
      control = list(trace = FALSE)
    )
    if(tail(opt$values,1) < tail(opt2$values,1)){
      alphabeta <- opt$pars[1:2]
    } else{
      alphabeta <- c(1,0)
    }
    if(qlev == 0.02){
      # Profile likelihood (not accounting for clustering)
      if(cl == 3L){
      constraint_field <- eht_constraints(
        alpha_grid = seq(-1, 1, length.out = 101),
        beta_grid = seq(-1, 1, length.out = 101),
        thresh = 1 - qlev,
        data = data[, hclust_od == cl]
      )
      }
      profile_skewnorm <- eht_pll_profile(
        alpha_grid = seq(-1, 1, length.out = 101),
        beta_grid = seq(-1, 1, length.out = 101),
        thresh = 1 - qlev,
        data = data[, hclust_od == cl],
        type = "skewnorm")
      # profile_norm <- eht_pll_profile(
      #   alpha_grid = seq(-1, 1, length.out = 101),
      #   beta_grid = seq(-1, 1, length.out = 101),
      #   thresh = 1 - qlev,
      #   data = data[, hclust_od == cl],
      #   type = "norm")
       pdf(paste0("../figures/Task4-profile_", cl,".pdf"), width = 5, height = 5)
       par(mfrow = c(1, 1), mar = c(4.5,4.5,0.5,0.5), bty = "l")
       with(profile_skewnorm,
            image(x = alpha,
                  y = beta,
                  z = pll,
                  col = viridis::viridis(256),
            xlab = expression(alpha),
            ylab = expression(beta)))
       points(x = alphabeta[1], alphabeta[2], pch = 20)
       points(x = opt_unconst$pars[1], opt_unconst$pars[2])
       dev.off()
      if(cl == 3L){
      # Create two plots, and fix dimensions
      rownames(profile_skewnorm$pll) <- profile_skewnorm$beta
      colnames(profile_skewnorm$pll) <- profile_skewnorm$alpha
      rownames(constraint_field$constrain) <- constraint_field$beta
      colnames(constraint_field$constrain) <- constraint_field$alpha
      # Profile pseudo log likelihood, with ggplot2
      g1 <- ggplot() +
        geom_tile(data = reshape::melt(profile_skewnorm$pll, as.is = TRUE),
                  aes(x = X1, y = X2, fill = value)) +
        geom_point(data = data.frame(X1 = alphabeta[1],
                                     X2 = alphabeta[2]),
                   col = "black",
                   aes(x = X1, y = X2),
                   size = 2) +
        geom_point(data = data.frame(X1 = opt_unconst$pars[1],
                                     X2 = opt_unconst$pars[2]),
                   shape = 2,
                   aes(x = X1, y = X2),
                   size = 1.5) +
        # coord_fixed() +
        scale_fill_viridis_c(na.value = "white")  +
        scale_y_continuous(limits = c(-1,1), expand = expansion()) +
        scale_x_continuous(limits = c(-1,1), expand = expansion()) +
        labs(x = expression(alpha), y = expression(beta)) +
        theme_classic() +
        theme(legend.position = "none")
        g2 <- ggplot() +
          geom_tile(data = reshape::melt(
            apply(constraint_field$constrain, 1:2, function(x){min(0,x)})),
                    aes(x = X1, y = X2, fill = value)) +
          # coord_fixed() +
          scale_fill_viridis_c(na.value = "white")  +
          scale_y_continuous(limits = c(-1,1), expand = expansion()) +
          scale_x_continuous(limits = c(-1,1), expand = expansion()) +
          labs(x = expression(alpha), y = expression(beta),fill = "penalty") +
          theme_classic() +
          theme(legend.position = "right")
      g1 + g2
      ggsave(filename = "../figures/FigureS3.pdf", width = 10, height = 5)
      dev.off()
      # with(constraint_field,
      #      image(
      #       x = alpha,
      #       y = beta,
      #       z = apply(constrain, 1:2, function(x){min(0, x)}),
      #       col = viridis::viridis(256),
      #       xlab = expression(alpha),
      #       ylab = expression(beta)))
      # }
      # with(profile_norm,
      #      fields::image.plot(
      #        x = alpha,
      #        y = beta,
      #        z = pll,
      #        col = viridis::viridis(256)))
    }
    # Store point estimates
    param_mat[, 1, qlev_ind] <- alphabeta
    param_mat_u[,1,qlev_ind] <- opt_unconst$pars[1:2]
    if (alphabeta[1] < 0.98) {
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
    parameters = param_mat,
    uparameters = param_mat_u
  )
  save(results,
       file = paste0("../outputs/C4_results_cl", cl, ".RData"),
       version = 2)
}


# Load all results
interm_res <- list.files(path = "../outputs/",
                         full.names = TRUE,
                         pattern = "C4_results_cl")
C4 <- matrix(nrow = 5, ncol = 2)
params <- matrix(nrow = 5, ncol = 2)
for (i in seq_along(interm_res)) {
  load(file = interm_res[i])
  C4[i, ] <- results$prob[4, 2:3]
  params[i, ] <- results$parameters[, 1, 4]
}
load("../outputs/C4_truth_postmortem.RData")
tab <- rbind(
  paste0("$",  as.integer(table(hclust_od)),
         "$ $(", as.integer(table(hclust_od[1:25])),
         " \\mid ", as.integer(table(hclust_od[26:50])), ")$"),
  paste0("$", sprintf("%.3f", probs[1,]), "$"),
  semipar_tbl[1,],
  paste0("$", sprintf("%.3f", t(C4[, 1])), "$"),
  paste0("$", sprintf("%.3f", probs[2,]), "$"),
  semipar_tbl[2,],
  paste0("$", sprintf("%.3f", t(C4[, 2])), "$"),
  paste0("$(", apply(params, 1, function(x) {
    paste(sprintf("%.2g", x),  collapse = ", ")
  }), ")$")
)
rownames(tab) <- c(
  "cluster size",
  "$\\log p_{1}$",
  "$\\log \\widehat{p}_{1,a}$",
  "$\\log \\widehat{p}_{1,b}$",
  "$\\log p_{2}$",
  "$\\log \\widehat{p}_{2,a}$",
  "$\\log \\widehat{p}_{2,b}$",
  "$(\\widehat{\\alpha}, \\widehat{\\beta})$"
)
cat(
  knitr::kable(
    tab,
    format = "latex",
    col.names = paste0("$C_", 1:5, "$"),
    booktabs = TRUE,
    escape = FALSE,
    align = c("rrrrr"),
    linesep = ""
  ),
  file = "../tables/Table6.tex"
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
        type = "skewnorm",
        constrain = FALSE
      ))
    },
    LB = c(-1, rep(-1e3, 3), -1e3),
    UB = c(rep(1, 2), rep(1e3, 2), 1e3),
    control = list(trace = FALSE)
  )
  params[cl,] <- opt$pars[1:2]
  adjprof <- adjust_loglik(
    loglik = eht_pll,
    constrain = FALSE,
    p = 5,  # number of parameters
    thresh = 0.98,  # marginal threshold
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
    mar = c(5, 4, 1, 1),
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
  yaxs = "i",
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
  yaxs = "i",
)
points(params[5, 1], params[5, 2], pch = 4)
dev.off()
