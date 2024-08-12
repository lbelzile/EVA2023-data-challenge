setwd(this.path::here())
library(evgam)
# Load new simulated data and true quantiles
load("../outputs/C1_check_coverage.RData") # load set of 10K new data points for external validation
load("../outputs/threshold_levels_lists.Rdata") # get fitted ALD models
source("../task2/Task2-functions.R") # quantile function of ALD
data <- read.csv("../data/Amaurot.csv")

nsimu <- 1000L
nobs <- nrow(qthresh)
summary_th <- matrix(nrow = 4, ncol = 5L)
for (mod_id in 1:4) {
  mod <- thr.levels.lists[[mod_id]]
  # Compute the MLE for the quantile at level tau
  predquant <- predict(mod, newdata = newTest)[seq_len(nobs), 1]
  # Draw from the ALD posterior to get the probability of true quantiles
  qqplot(predquant, qthresh[, mod_id], panel.first = {
    abline(a = 0, b = 1)
  })
  post_ald <- simulate(mod,
                       newdata = newTest,
                       nsim = nsimu,
                       type = "response")
  quants <- matrix(nrow = nobs, ncol = 5L)
  for (i in seq_len(nobs)) {
    quants[i, ] <- quantile(
      pald(
        q = qthresh[i, mod_id],
        location = post_ald$location[i, ],
        scale = post_ald$scale[i, ],
        tau = mod$tau
      ),
      probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
    )
  }
  summary_th[mod_id,] <-
    c(mod$tau,
      mean(quants[, 3]) - mod$tau, # bias
      sd(quants[, 3])/sqrt(10000), # std. error
      1 - mean(quants[, 2] > mod$tau | quants[, 4] < mod$tau), # coverage 50%
      1 - mean(quants[, 1] > mod$tau | quants[, 5] < mod$tau)  # coverage 80%
    )
}
cat(
knitr::kable(
  summary_th*100,
  col.names = c("level", "bias", "std. error", "coverage (50\\%)", "coverage (80\\%)"),
  booktabs = TRUE,
  format = "latex",
  digits = c(0,2,3,1,1),
  escape = FALSE)
    file = "../tables/TableS3.tex")
