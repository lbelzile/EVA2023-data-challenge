setwd(this.path::here())
load("../outputs/results_AD_cl4.RData")
library(ggplot2)
library(patchwork)
library(mev)

qlev <- qgev(seq(0.8, 0.99, by = 0.01), 1,1,1)
d <- 8
nmixt <- 100
pars <- seq(0.4, 0.9, length.out = nmixt)
nsim <- 1e6
nobs <- nsim / nmixt
if(! "outputs/C4_tstab_proptbl.RData" %in% list.files("..", recursive = TRUE)){
  set.seed(2023)
prop_tbl <- matrix(0, nrow = length(qlev), ncol = nmixt)
for (nrep in 1:1000) {
  samp <- matrix(NA, nrow = nsim, ncol = d)
  for (i in seq_along(pars)) {
    samp[(i - 1) * nobs + 1:nobs, ] <- rmev(
      n = nobs,
      d = d,
      model = "log",
      param = pars[i]
    )
  }

  for (j in seq_along(qlev)) {
    prop_tbl[j, ] <- prop_tbl[j, ]  + table(factor(ceiling(which(
      apply(samp, 1, max) > qlev[j]
    ) / nobs),
    levels = 1:nmixt),
    exclude = FALSE)
  }
}

prop_tbl <- prop_tbl / rowSums(prop_tbl, na.rm = TRUE)
save(prop_tbl, file = "../outputs/C4_tstab_proptbl.RData")
} else{
  load(file = "../outputs/C4_tstab_proptbl.RData")
}

df1 <- data.frame(
  pars = c(pars, pars, pars),
  prop = c(prop_tbl[1,], prop_tbl[11,], prop_tbl[16,]),
  thresh = factor(rep(c("0.8", "0.9", "0.95"), each = length(pars))))
g1 <- ggplot(data = data.frame(
  thresh = results$log[, 'threshold'],
  results$log[, c('coef', 'lower', 'upper')]),
  mapping = aes(x = thresh, y = coef, ymin = lower, ymax = upper)) +
  geom_pointrange(fatten = 1) +
  scale_x_continuous(breaks = seq(0.9, 1, by = 0.025),
                     labels = c("0.9", "0.925", "0.95", "0.975", "1")) +
  labs(x = "Quantile level",
       y = "",
       subtitle = expression("Logistic dependence parameter " * alpha)) +
  theme_classic()
g2 <- ggplot(data = df1,
             mapping = aes(x = pars,
                           y = prop*100,
                           col = thresh)) +
  geom_hline(yintercept = 1) +
  geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x), show.legend = TRUE) +
  # geom_point() +
  scale_color_grey() +
  labs(x = expression(alpha),
       y = "",
       color = "Quantile level",
       subtitle = "Proportion of exceedances (%)") +
  theme_classic() +
  theme(legend.position = c(0.3, 0.8))


pdf("../figures/FigureS4.pdf", width = 6, height = 3)
g1 + g2
dev.off()
