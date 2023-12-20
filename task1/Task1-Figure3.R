library(evgam)
library(mev)
setwd(this.path::here())
load("../outputs/Amaurot_imputed.Rdata") ## imputed Amaurot data
load("../outputs/C1_models.Rdata")  ### GP models and threshold model list


data.train <- Amaurot.imp.final
data.train$Season <- ifelse(data.train$Season == "S1", 1, 0)

prob.choice <- 0.95
ny <- 70
npy <- 300
n <- ny * npy


### fit a threshold model
fit_thr_ald <- evgam(
  model_thr_ald,
  data.train,
  family = "ald",
  ald.args = list(tau = prob.choice)
)
data.train$threshold <- predict(fit_thr_ald)$location
data.train$excess <- data.train$Y - data.train$threshold
is.na(data.train$excess[data.train$excess <= 0]) <- TRUE

### extract exceedances
data_excess <- data.train[!is.na(data.train$excess), ]
####  split the data in three folds
k.fold <- 3L
nmod <- length(model_gpd.list)
ne <- nrow(data_excess) ## number of exceedances

### permutations over all the possibility of train and test datasets
arrangement <- arrangements::permutations(
    x = c("fold1", "fold2", "fold3"),
    replace = FALSE)

### Score function
score_interv <- function(l, u, y, alpha = 0.5){
  (u - l) + 2/alpha*(l - y)*I(y < l) + 2/alpha*(y - u)*I(y > u)
}

final_scores <- matrix(
  nrow = nrow(arrangement),
  ncol = nmod)
final_cover <- matrix(
   nrow = nrow(arrangement),
   ncol = nmod)
# Number of replications
R <- 100L
mod_id <- 2
results <- matrix(nrow = R, ncol = nmod)
coverage <- matrix(nrow = R, ncol = nmod)
preds <- list()
for(r in seq_len(R)){
  set.seed(r + 2023)
  ind_sample <- rep(1:k.fold, length.out = ne)[sample(
    x = 1:ne,
    size = ne,
    replace = FALSE)] ### randomly divide the data into three fold
  fold_data <- list(
    "fold1" = data_excess[ind_sample == 1L,],
    "fold2" = data_excess[ind_sample == 2L,],
    "fold3" = data_excess[ind_sample == 3L,])
## fit all the 7 GP models to folds
list.fit_gpd <- list(
  "fold1" = lapply(model_gpd.list, function(model_gpd) {
    evgam(model_gpd, fold_data$fold1, family = "gpd")
  }),
  "fold2" = lapply(model_gpd.list, function(model_gpd) {
    evgam(model_gpd, fold_data$fold2, family = "gpd")
  }),
  "fold3" = lapply(model_gpd.list, function(model_gpd) {
    evgam(model_gpd, fold_data$fold3, family = "gpd")
  })
)

for (ll in 1:nrow(arrangement)) {
  ### loop over the 6 choices of  exchanging the role over folds
  fold1 <- fold_data[[arrangement[ll, 1]]]
  fold2 <- fold_data[[arrangement[ll, 2]]]
  fold3 <- fold_data[[arrangement[ll, 3]]]


  ### predicting the probability level in train 1
  p.star <- matrix(nrow = nrow(fold3), ncol = nmod)
  for (i in seq_len(nmod)) {
    pars2 <- predict(
      object = list.fit_gpd[[arrangement[ll,1]]][[i]],
      newdata = fold3)
    p.star[, i] <- revdbayes::pgp(
        q = fold3$excess,
        loc = 0,
        scale = exp(pars2[,'logscale']),
        shape = pars2[,'shape']
      )
  }
  ## Some probabilities can be predicted to 1
  p.star <- apply(p.star, 1:2, function(x){ pmin(x, 1-1e-6)})


  ### predicting the quantile level in train 2 along with the 50% CIs
  ### simulate scale and shape parameter of the GPD for train 2
  nsim <- 1000L
  simulated_pars <-
    lapply(list.fit_gpd[[arrangement[ll,2]]],
           function(gpd_model)
      simulate(
        gpd_model,
        newdata = fold3,
        nsim = nsim,
        type = "response"
      ))
  
  y.star <- array(dim = c(nrow(fold3), nmod, nsim))
  ### loop over the models
  for (i in seq_len(nmod)) {
    ### loop over the set of fold3 data
      for (j in seq_len(nsim)) {
        ## loop over the simulated samples
        y.star[, i, j] <- revdbayes::qgp(
          p.star[, i],
          loc = 0,
          scale = simulated_pars[[i]]$scale[, j],
          shape = simulated_pars[[i]]$shape[, j]
        )
      }
  }

  #### 50% CI  for all the 7 models and all the observations in fold3 sets
  l_CI <- apply(y.star,
                MARGIN = c(1, 2),
                FUN = quantile,
                probs = 0.25)
  u_CI <- apply(y.star,
                MARGIN = c(1, 2),
                FUN = quantile,
                probs = 0.75)
  if(r == 1L & ll %in% c(1L, 4L, 5L)){
     preds[[match(ll, c(1L, 4L, 5L))]] <- cbind(
        truth = fold3$excess, 
        pred = apply(y.star[,mod_id,],1, median),
        lower = l_CI[,mod_id],
        upper = u_CI[,mod_id])
  }

  ### calculating the scores for all the 7 models
  ### and for all the excesses in fold3 data
  scores <- matrix(nrow = nrow(fold3), ncol = nmod)
  for (i in 1:length(model_gpd.list)) {
      scores[, i] <- score_interv(
        l = l_CI[,i],
        u = u_CI[,i],
        y = fold3$excess)

    }
  final_cover[ll, ] <- colSums(l_CI < fold3$excess & fold3$excess < u_CI)
  final_scores[ll, ] <- apply(scores, MARGIN = 2, FUN = sum)

}
results[r,] <- colSums(final_scores) / (2*ne)
coverage[r,] <- colSums(final_cover)/(2*ne)
print(r)
}
save(results, coverage,
     file = "../outputs/cross-validation-interval-score-95.RData",
     version = 2)

# Create quantile-quantile plot for model 2
library(ggplot2)
library(patchwork)


g1 <- ggplot(
  data = data.frame(score = as.vector(results),
                    model = factor(rep(1:7, each = R))),
  mapping = aes(x = model, y = score)) +
  scale_y_continuous(limits = c(0, 16.2), expand = c(0,0)) +
  ggdist::stat_slabinterval(point_interval = "median_qi",.width = 0.5) +
  labs(x = "Model", y = "", subtitle = "Interval score") +
  theme_classic()

g2 <- ggplot(
   data = data.frame(score = 100*as.vector(coverage),
                     model = factor(rep(1:7, each = R))),
   mapping = aes(x = model, y = score)) +
   scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
   ggdist::stat_slabinterval(point_interval = "median_qi",.width = 0.5) +
   labs(x = "Model", y = "", subtitle = "Coverage") +
   theme_classic()
preds_mat <- as.data.frame(rbind(preds[[1]], preds[[2]], preds[[3]]))
cols <- with(preds_mat, ifelse(truth > lower & truth < upper, 1, 2))
mean(cols) - 1
g3 <- ggplot(data = preds_mat, 
             aes(x = truth, 
                 y = pred, 
                 ymin = lower,
                 ymax = upper)) +
   geom_abline(slope = 1, intercept = 0) + 
   geom_point(size = 0.5, color = c("black","grey")[cols]) +
   geom_errorbar(color = c("black","grey")[cols]) +
   scale_x_continuous(limits = c(0, max(preds_mat))) +
   scale_y_continuous(limits = c(0, max(preds_mat))) +
   labs(x = "True exceedance", y = "", subtitle = "Prediction with 50% interval") +
   theme_classic()


g1 + g2 + g3
ggsave("../figures/Figure3.pdf", width = 7.5, height = 2.5)
