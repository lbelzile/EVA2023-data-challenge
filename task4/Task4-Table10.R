setwd(this.path::here())
library(mev)
library(mvPot)
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
####     CROSS-VALIDATION FOR AD CLUSTERS      ######
#####################################################
for (k in 1:2){
  cl <- c(1, 4)[k]
  mthresh <- 0.5 # marginal threshold for left-censoring
  # risk threshold for max
  clsize <- sum(hclust_od == cl)
  u1size <- sum((1:50)[hclust_od == cl] <= 25)
  Cl <- unif[, hclust_od == cl]
  pairs <- combn(clsize, 2)
  th <- 0.97 # threshold level
  vlev <- 0.99 # validation level
  vlev2 <- 0.96
  # Containers for l2 distance
  dchi_log_p <- rep(0, ncol(pairs))
  dchi_br_p <- rep(0, ncol(pairs))
  dchi_ht_p <- rep(0, ncol(pairs))
  # Containers for w2 distance
  dw2_log_p <- rep(0, ncol(pairs))
  dw2_br_p <- rep(0, ncol(pairs))
  dw2_ht_p <- rep(0, ncol(pairs))
  for (i in seq_len(ncol(pairs))) {
    set.seed(i)
    # extract data from the cluster, ignoring validation data
    sdat <- Cl[, -pairs[, i]]
    # subsample to retrieve exceedances
    sdat <- sdat[apply(sdat, 1, max) > th, ]
    # compute empirical measure (eq. 7)
    emp_chi <- mean(rowSums(Cl[, pairs[, i]] > vlev) == 2) / (1 - vlev)
    # Consider different probabilities for w2 (eq. 8)
    emp_w2 <- mean((Cl[, pairs[1, i]] > vlev &
                       Cl[, pairs[2, i]] > vlev2) |
                   (Cl[, pairs[1, i]] > vlev2 &
                         Cl[, pairs[2, i]] > vlev)) /
      (2 * (1 - th))
    # Fit the multivariate generalized Pareto models
    opt_par_log_cv <- fit_mgp(
      data = sdat,
      model = "log",
      mthresh = mthresh,
      thresh = th
    )
    opt_par_br_cv <- fit_mgp(
      data = sdat,
      model = "br",
      mthresh = mthresh,
      thresh = th
    )
    # Compute l2 distance (pointwise)
    dchi_log_p[i] <-
      expme_log_min(u = rep(1, 2),
                    par = opt_par_log_cv$mle) - emp_chi
    dchi_br_p[i] <-
      expme_BR_min(z = rep(1, 2),
                   Lambda = Lmat_exch(opt_par_br_cv$mle, D = 2)) -
      emp_chi
    # Compute w2 distance (pointwise)
    dw2_log_p[i] <- expme_log_min(
      u = c(vlev, vlev2),
      par = opt_par_log_cv$mle) /
      expme(z = rep(th, 2), par = opt_par_log_cv$mle, model = "log") - emp_w2
    dw2_log_p[i] <-
      expme_log_min(u = rep(1, 2),
                    par = opt_par_log_cv$mle) - emp_chi
    dw2_br_p[i] <- expme_BR_min(
      z = c(vlev, vlev2),
      Lmat_exch(opt_par_br_cv$mle, D = 2)) /
      expme(z = rep(th, 2),
            par = list(Lambda = Lmat_exch(opt_par_br_cv$mle, D = 2)),
            model = "br") - emp_w2
   # fit the exchangeable conditional extremes model
    opt <- Rsolnp::solnp(
      pars = c(0.2, 0.1, rep(0, 3)),
      fun = function(x) {
        -sum(eht_pll(
          par = x,
          thresh = th,
          data = apply(Cl, 2, qlaplace)[, -pairs[, i]],
          type = "skewnorm"
        ))
      },
      LB = c(-1, rep(-1e3, 3), -1e3),
      UB = c(rep(1, 2), rep(1e3, 2), 1e3),
      control = list(trace = FALSE)
    )
    alpha <- opt$pars[1]
    beta <- opt$pars[2]
    residuals <- # use function to get residuals
      c(residuals_eht(
        alpha = alpha,
        beta = beta,
        data = apply(Cl, 2, qlaplace)[, -pairs[, i]],
        thresh = th
      )$res)
    y0 <- rexp(n = length(residuals)) + qlaplace(vlev)
    dchi_ht_p[i] <-
      mean(alpha * y0 + y0 ^ beta * residuals > qlaplace(vlev)) -
      emp_chi
   y0 <- rexp(n = length(residuals)) + qlaplace(vlev2)
   dw2_ht_p[i] <-
    mean(alpha * y0 + y0 ^ beta * residuals > qlaplace(vlev))/
    mean(y0 > qlaplace(th) & (alpha * y0 + y0 ^ beta * residuals > qlaplace(th))) -
    emp_w2
  }

  dchi_l2_log_p <- sqrt(sum(dchi_log_p^2)) / ncol(pairs)
  dchi_l2_br_p <- sqrt(sum(dchi_br_p^2)) / ncol(pairs)
  dchi_l2_ht_p <- sqrt(sum(dchi_ht_p^2)) / ncol(pairs)
  dw2_log_p <- sqrt(sum(dw2_log_p^2)) / ncol(pairs)
  dw2_br_p <- sqrt(sum(dw2_br_p^2)) / ncol(pairs)
  dw2_ht_p <- sqrt(sum(dw2_ht_p^2)) / ncol(pairs)
  # Ratio of l2 measures (smaller is better)
  dchi_l2_log_p / dchi_l2_br_p
  dchi_l2_log_p / dchi_l2_ht_p
  # Ratio of w2 measures (smaller is better)
  dw2_log_p / dw2_br_p
  dw2_log_p / dw2_ht_p
  # Test for equality using paired samples
  t.test(x = dchi_log_p^2,
         y = dchi_br_p^2,
         paired = TRUE)
  t.test(x = dchi_log_p^2,
         y = dchi_ht_p^2,
         paired = TRUE)

  triples <- combn(clsize, 3)
  th <- 0.97
  vlev <- 0.99
  dchi_log_t <- rep(0, ncol(triples))
  dchi_br_t <- rep(0, ncol(triples))
  dchi_ht_t <- rep(0, ncol(triples))
  for (i in seq_along(dchi_log_t)) {
    set.seed(i)
    sdat <- Cl[, -triples[, i]]
    sdat <- sdat[apply(sdat, 1, max) > th, ]
    emp_chi <- mean(rowSums(Cl[, triples[, i]] > vlev) == 3) / (1 - vlev)
    opt_par_log_cv <- fit_mgp(
      data = sdat,
      model = "log",
      mthresh = mthresh,
      thresh = th
    )
    dchi_log_t[i] <-
      expme_log_min(u = rep(1, 3),
                    par = opt_par_log_cv$mle) -
      emp_chi
    opt_par_br_cv <- fit_mgp(
      data = sdat,
      model = "br",
      mthresh = mthresh,
      thresh = th
    )
    dchi_br_t[i] <-
      expme_BR_min(z = rep(1, 3),
                   Lambda = Lmat_exch(opt_par_br_cv$mle, D = 3)) -
      emp_chi

    opt <- Rsolnp::solnp(
      pars = c(0.2, 0.1, rep(0, 3)),
      fun = function(x) {
        -sum(eht_pll(
          par = x,
          thresh = th,
          data = apply(Cl, 2, qlaplace)[, -triples[, i]],
          type = "skewnorm"
        ))
      },
      LB = c(-1, rep(-1e3, 3), -1e3),
      UB = c(rep(1, 2), rep(1e3, 2), 1e3),
      control = list(trace = FALSE)
    )
    alpha <- opt$pars[1]
    beta <- opt$pars[2]
    residuals <- # use function to get residuals
      residuals_eht(
        alpha = alpha,
        beta = beta,
        data = apply(Cl, 2, qlaplace)[, -triples[, i]],
        thresh = th
      )$res

    dchi_ht_t[i] <- mean(replicate(expr = {
      y0 <- rexp(n = length(residuals)) + qlaplace(vlev)

      mean(sapply(1:nrow(residuals),
                  function(j) {
                    sum(alpha * y0[j] + y0[j] ^ beta *
                          residuals[j, sample(ncol(residuals), 2)] >
                          qlaplace(vlev)) == 2L
                  }))
    }, n = 100))  -  emp_chi
  }

  dchi_l2_log_t <- sqrt(sum(dchi_log_t^2)) / ncol(triples)
  dchi_l2_br_t <- sqrt(sum(dchi_br_t^2)) / ncol(triples)
  dchi_l2_ht_t <- sqrt(sum(dchi_ht_t^2)) / ncol(triples)
  # Ratio of l2-distance
  dchi_l2_log_t / dchi_l2_br_t
  dchi_l2_log_t / dchi_l2_ht_t
  # Test for equality using paired samples
  t.test(x = dchi_log_t^2,
         y = dchi_br_t^2,
         paired = TRUE)
  t.test(x = dchi_log_t^2,
         y = dchi_ht_t^2,
         paired = TRUE)

results_gof <- rbind(
    c(dchi_l2_log_p, dchi_l2_br_p, dchi_l2_ht_p),
    c(dw2_log_p, dw2_br_p, dw2_ht_p),
    c(dchi_l2_log_t, dchi_l2_br_t, dchi_l2_ht_t)
  )



 save(
    results_gof,
    file = paste0("../outputs/Task4-results_AD_gof_cl", cl, ".RData"),
    version = 2
  )
}

load("../outputs/Task4-results_AD_gof_cl1.RData")
gof1 <- apply(results_gof, 1:2, function(x){sprintf(1000*x, fmt = "%.1f")})

load("../outputs/Task4-results_AD_gof_cl4.RData")
gof2 <- apply(results_gof, 1:2, function(x){sprintf(1000*x, fmt = "%.1f")})
  # Print table with distances
rownames(gof1) <- rownames(gof2) <- c(
  "pairs",
  "triples",
  "unequal")
colnames(gof1) <- colnames(gof2) <- c("logistic","Hüsler--Rei{\\ss}","Heffernan--Tawn")

# Create table
tab <- xtable::xtable(x = cbind(t(gof1[c(1,3,2),]),
                               t(gof2[c(1,3,2),])),
                   align = c("lrrrrrr"))


addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- paste0(
  paste0('& \\multicolumn{3}{c}{cluster ', c(1,4), '}',
         collapse=''),
  '\\\\',
  "\\cmidrule(l){2-4} \\cmidrule(l){5-7}",
  "Model &",
  paste0(rep(c("pairs",
               "triples",
               "unequal"), length.out = 6L),
         collapse = " & ",
         sep = " "),"\\\\")
print(tab,
      booktabs = TRUE,
      floating = FALSE,
      include.colnames = FALSE,
     #file = "../tables/Table10.tex",
      add.to.row = addtorow,
      sanitize.rownames.function = identity)
