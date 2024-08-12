setwd(this.path::here())
devtools::install_github("nicolagnecco/erf")
library(erf)
library(quantregForest)
library(evgam)

load("../outputs/Amaurot_imputed.Rdata")

Y <- Amaurot.imp.final$Y
V1 <- Amaurot.imp.final$V1
V2 <- Amaurot.imp.final$V2
V3 <- Amaurot.imp.final$V3
V4 <- Amaurot.imp.final$V4
WindDirection <- Amaurot.imp.final$WindDirection
WindSpeed <- Amaurot.imp.final$WindSpeed
Atmosphere <- Amaurot.imp.final$Atmosphere
Season <- ifelse(Amaurot.imp.final$Season == "S2", 1, 0)
CP <- ifelse(Amaurot.imp.final$CP == "1", 1, 0)

data <- data.frame(Y = Y,
                   V1 = V1,
                   V2 = V2,
                   V3 = V3,
                   V4 = V4,
                   Season = Season,
                   WindDirection = WindDirection,
                   WindSpeed = WindSpeed,
                   Atmosphere = Atmosphere,
                   CP = CP)
rm(Amaurot.imp.final)
#-------------------------------------------------------------------------------

prob.choices <- seq(0.95, 0.99, by = 0.01)
X <- cbind(V1, V2, V3, V4, WindDirection, WindSpeed, Atmosphere, Season, CP)

# Cross-validation
set.seed(2023)
K <- 10L
groups <- split(
   x = sample.int(nrow(X), nrow(X)),
   f = factor(rep(seq_len(K), length.out = nrow(X))))
predictions <- array(dim = c(nrow(X), 3L, length(prob.choices)))
proportions <- array(dim = c(K, 3L, length(prob.choices)))
for(k in seq_len(K)){
   print(Sys.time())
# quantregForest
quantregforest.object <- quantregForest(
   x = X[-groups[[k]],], 
   y = Y[-groups[[k]]], 
   nthreads = 4, 
   keep.inbag = TRUE)

predictions[groups[[k]], 1, ] <- predict(quantregforest.object, 
                                         newdata = X[groups[[k]],],
                                         what = prob.choices)

rm(quantregforest.object)
# ERF
erf.object <- erf(X = X[-groups[[k]],],
                  Y = Y[-groups[[k]]], 
                  intermediate_quantile = 0.8)

predictions[groups[[k]], 2, ] <- predict(erf.object, 
                          newdata = X[groups[[k]],], 
                          quantiles = prob.choices) 
rm(erf.object)

# evgam
model_thr_ald <- list(Y ~ Season + CP + s(V1, bs = "cs") + 
                         s(V2, bs = "cs") + s(V3, bs = "cs") + s(V4, bs = "cs") +
                         s(WindDirection, bs = "cs") + s(WindSpeed, bs = "cs") + 
                         s(Atmosphere, bs = "cs"),
                      ~ Season + CP + s(V1, bs = "cs") + 
                         s(V2, bs = "cs") + s(V3, bs = "cs") + s(V4, bs = "cs") +
                         s(WindDirection, bs = "cs") + s(WindSpeed, bs = "cs") + 
                         s(Atmosphere, bs = "cs"))
for(i in seq_along(prob.choices)){
   fit_thr_ald <- evgam(model_thr_ald,
                        data = data[-groups[[k]],],
                        family = "ald", 
                        ald.args = list(tau = prob.choices[i]))
   predictions[groups[[k]], 3, i] <- predict(
      object = fit_thr_ald,
      newdata = data[groups[[k]],],
      type = "response")$location
}
proportions[k,,] <- apply(predictions[groups[[k]], , ] , 2:3, function(yhat){
   mean(yhat > Y[groups[[k]]])
})
}
# End of cross-validation loop
props <- apply(predictions, 2:3, function(yhat){
   mean(yhat > Y)
})
sd_props <- apply(proportions, 2:3, sd)

save(props, sd_props,
     file = "../outputs/Task1_threshold_model_comparison.Rdata")

#-------------------------------------------------------------------------------

props.table <- apply(props, 1:2, function(x){sprintf(100*x, fmt = "%.2f")})

rownames(props.table) <- c("\\texttt{quantregForest}", "\\texttt{erf}", "\\texttt{evgam}")
colnames(props.table) <- c("$q_{0.95}$","$q_{0.96}$","$q_{0.97}$","$q_{0.98}$","$q_{0.99}$")
cat(knitr::kable(props.table, 
                 format = "latex", 
                 escape = FALSE, 
                 booktabs = TRUE),
    file = "../tables/Table1.tex")
