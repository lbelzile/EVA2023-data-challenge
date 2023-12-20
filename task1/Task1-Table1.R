setwd(this.path::here())
library(missForest)
library(mice)
remotes::install_github("cran/imputeMissings")
library(imputeMissings)

# Compute the distance between angles using sine-cosine
rmse_circular <- function(y, yhat){
  mean(sqrt((cos(y)-cos(yhat))^2 + (sin(y)-sin(yhat))^2))
}
# Load data and keep only complete cases (since MCAR)
Amaurot <- read.csv("../data/Amaurot.csv")
Amaurot$Season <- factor(Amaurot$Season)
AmaurotComplete <- na.omit(Amaurot)
# Compute the fraction of missing values in each variable
nmiss <- apply(Amaurot, 2, function(x){sum(is.na(x))})
perc_miss <- nmiss/sum(nmiss)

# Create a dataset with MCAR values
set.seed(2023)
holes <- mice::ampute(data = AmaurotComplete,
                      prop = 0.1,
                      freq = perc_miss,
                      mech = "MCAR")$amp
forms <- mice::make.formulas(holes)
forms$V1 <- update(forms$V1, ~ Y + V2)
forms$V2 <- update(forms$V2, ~ Y + V1)
forms$V3 <- update(forms$V3, ~ Y + Season)
forms$V4 <- update(forms$V4, ~ Y)
forms$WindDirection <- update(forms$WindDirection, ~ Y + Season + WindSpeed)
forms$WindSpeed <- update(forms$WindSpeed, ~ Y + Season + cos(WindDirection) + sin(WindDirection))

imp_mice <- mice::mice(holes, m = 10, formulas = forms)

imp_forest <- missForest(holes, maxiter = 10)

imp_median <- imputeMissings::impute(holes)

V1mod <- mgcv::gam(V1 ~ s(Y) + s(V2), data = holes)
V2mod <- mgcv::gam(V2 ~ s(Y) + s(V1), data = holes)
V3mod <- mgcv::gam(V3 ~ s(Y) + Season, data = holes)
V4mod <- mgcv::gam(V4 ~ s(Y), data = holes)
WindSpeedmod <- mgcv::gam(WindSpeed ~ s(WindDirection, bs = "cc", knots = c(0, 2*pi)) + s(Atmosphere), data = holes)
WindDirectionmod <- mgcv::gam(WindDirection ~ s(WindSpeed) + s(Atmosphere), data = holes)


missCols <- which(nmiss > 0)
missIndices <- apply(holes[,missCols], 2, function(x){which(is.na(x))})
colSD <- apply(AmaurotComplete[,missCols], MARGIN = 2, FUN = sd)
rmse <- matrix(0, ncol = length(missCols), nrow = 4)
for(i in seq_along(missCols)){
  if(names(missCols)[i] != "WindDirection"){
  for(j in 1:10){
   rmse[1,i] <- rmse[1,i]  + mean((AmaurotComplete[missIndices[[i]],missCols[i]] - complete(imp_mice, j)[missIndices[[i]],missCols[i]])^2)
  }
  rmse[1,i] <- sqrt(rmse[1,i]/10)
  rmse[2,i] <- sqrt(mean((AmaurotComplete[missIndices[[i]],missCols[i]] - imp_forest$ximp[missIndices[[i]],missCols[i]])^2))
  rmse[3,i] <- sqrt(mean((AmaurotComplete[missIndices[[i]],missCols[i]] - imp_median[missIndices[[i]],missCols[i]])^2))
  } else{
    for(j in 1:10){
      rmse[1,i] <- rmse[1,i]  + rmse_circular(
        AmaurotComplete[missIndices[[i]],missCols[i]],
        complete(imp_mice, j)[missIndices[[i]],missCols[i]])
    }
    rmse[1,i] <- rmse[1,i]/10
    rmse[2,i] <- rmse_circular(AmaurotComplete[missIndices[[i]],missCols[i]], imp_forest$ximp[missIndices[[i]],missCols[i]])
    rmse[3,i] <- rmse_circular(AmaurotComplete[missIndices[[i]],missCols[i]], imp_median[missIndices[[i]],missCols[i]])
  }
}


rmse[4,1] <- sqrt(mean((predict(V1mod, holes[is.na(holes$V1),]) - AmaurotComplete[missIndices$V1, missCols['V1']])^2))
rmse[4,2] <- sqrt(mean((predict(V2mod, holes[is.na(holes$V2),]) - AmaurotComplete[missIndices$V2, missCols['V2']])^2))
rmse[4,3] <- sqrt(mean((predict(V3mod, holes[is.na(holes$V3),]) - AmaurotComplete[missIndices$V3, missCols['V3']])^2))
rmse[4,4] <- sqrt(mean((predict(V4mod, holes[is.na(holes$V4),]) - AmaurotComplete[missIndices$V4, missCols['V4']])^2))
rmse[4,5] <- sqrt(mean((predict(WindSpeedmod, holes[is.na(holes$WindSpeed),]) - AmaurotComplete[missIndices$WindSpeed, missCols['WindSpeed']])^2))
rmse[4,6] <- rmse_circular(predict(WindDirectionmod, holes[is.na(holes$WindDirection),]),
                           AmaurotComplete[missIndices$WindDirection, missCols['WindDirection']])

rmse_tab <- apply(rmse[c(4,3,1,2),], 1:2, function(x){sprintf(x, fmt = "%.2f")})
rownames(rmse_tab) <- c("\\texttt{gam}","\\texttt{impute}","\\texttt{mice}","\\texttt{missForest}")
colnames(rmse_tab) <- c("$V_1$","$V_2$","$V_3$","$V_4$","Wind speed","Wind direction")
cat(knitr::kable(rmse_tab, 
                 format = "latex", 
                 escape = FALSE, 
                 booktabs = TRUE),
     file = "../tables/Table1.tex")
