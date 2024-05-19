setwd(this.path::here())
library(texmex)
library(geometricMVE)
source("Task3-functions.R")
Coputopia <- pgumbel(as.matrix(read.csv("../data/Coputopia.csv", header = TRUE)[,3:5]))

# Dependence quantile
thresh <- c(0.9, 0.95, 0.96, 0.97, 0.98)
a1 <- 0.96
a2 <- 0.97

set.seed(2023)
resultsC3 <- list()
for(i in seq_along(thresh)){
   resultsC3[[i]] <- estimate_probs_C3(data = Coputopia, thresh = thresh[i])
   print(i)
   save(resultsC3, file = "../outputs/Task3-results.RData")
}

res_tab_log <- log(rbind(
  sapply(resultsC3, function(x){x[,1]}), 
  sapply(resultsC3, function(x){x[,2]}), 
  sapply(resultsC3, function(x){x[,3]}))[-5,]) # remove columns of NAs
res_tab <- apply(res_tab_log, 1:2, function(x){
      paste0("$", sprintf(fmt = "%.2f", x), "$")})
res_tab <- cbind(c("", "$\\log p_1$", "", "","", "$\\log p_2$", "",""),
                 rep(c("conditional", "HRV", "geometric", 
                 "conditional", "geometric", "conditional (2)", "HRV (2)", "geometric (2)")),
                 res_tab)
cat(
knitr::kable(res_tab,
             booktabs = TRUE,
             align = "llrrrrrr",
             row.names = FALSE,
             escape = FALSE,
             col.names = c("", "method", "$q_{0.90}$", "$q_{0.95}$", "$q_{0.96}$",
                           "$q_{0.97}$","$q_{0.98}$"),
             format = "latex",
             linesep = ""),
 file = "../tables/Table4.tex")
