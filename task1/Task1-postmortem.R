
# Generate new data using instructions from section 4.1 of the Data challenge editorial
data <- read.csv("../data/Amaurot.csv")
nB <- 1e4
set.seed(2023)
W <- ifelse(runif(nB) < 0.4,
            rnorm(nB, mean = 30, sd = 3),
            rnorm(nB, mean = 36, sd = sqrt(6.25)))
V2 <- W + rgamma(nB, rate = 0.3, shape = 1.2)
V1 <- W + rnorm(nB, sd = 2)
Season <- rep(c("S1", "S2"), each = nB/2)


V3 <- ifelse(Season == "S1",
             sn::rsn(nB, xi = 4, omega = 4, alpha = 5),
             mev::qgev(mev::pgev(0, 4, 4, 0.2) + runif(n = nB)*mev::pgev(0, 4, 4, 0.2, lower.tail = FALSE), loc = 4, scale = 4, shape = 0.2))
V4 <- pmax(0, mev::rgev(nB, 0,1, 0.1))
WindDirection <- as.numeric(circular::rmixedvonmises(
   n = nB,
   mu1 = circular::circular(atan2(x = 0.3,y = 0.5)),
   mu2 = circular::circular(atan2(x = -0.4,y = -0.4)),
   kappa1 = 2, kappa2 = 4, prop = 0.4))
WindDirection <- round(ifelse(WindDirection > pi, -2*pi + WindDirection, WindDirection),3)
WindSpeed <- TruncatedNormal::rtnorm(
   n = 1,
   mu = 3 + 2*sin(WindDirection),
   lb = rep(0, nB),
   ub = rep(Inf, nB),
   sd = 1)

t.test(V1, data$V1)
t.test(abs(V1 - mean(V1)), abs(data$V1 - mean(V1)))
t.test(V2, data$V2)
t.test(abs(V2 - mean(V2)), abs(data$V2 - mean(V2)))
t.test(V3, data$V3)
t.test(V4, data$V4)
t.test(WindDirection, data$WindDirection)
t.test(WindSpeed, data$WindSpeed)

newTest <- data.frame(V1, V2, V3, V4, Season = factor(Season), WindDirection, WindSpeed, Atmosphere = sample(data$Atmosphere, nB))
if(!"../outputs/C1_check_coverage.RData" %in% list.files(recursive = TRUE)){
   qtrue <- ptrue <- numeric(nB)
   B <- 3e7
   scale <- with(newTest,
                 exp(-1) + 18.7*sqrt(V2) + 9*(1+ log(V3))^2 + 5.71*WindSpeed^1.5)/10
   shape <- -1/pi^2
   loc <- with(newTest, ifelse(
      Season == "S1", 112-6*abs(WindDirection), 110-5*abs(WindDirection)^0.9)
   )
   for(i in seq_len(nB)){
      print(i)
      set.seed(2023)
      samp <- mev::rgp(n = B, loc = 0, scale = scale[i], shape = shape)
      if(newTest$Season[i] == "S1"){
         Y <- na.omit(ifelse(samp > loc[i], samp, ifelse(
            rbeta(B, samp, loc[i]) > runif(B), samp, NA)))
      } else{
         Y <- na.omit(ifelse(samp > loc[i], samp, ifelse(
            rbeta(B, exp(2 + (samp-loc[i])/30), 1) > runif(B), samp, NA)))
      }
      qtrue[i] <- quantile(Y, 0.9999)
      ptrue[i] <- mean(Y >= 112)
   }
   qtrue2 <- ifelse(ptrue < 1e-4,
                    qtrue,
                    qgp(pmax(0,(ptrue-1e-4)/ptrue),
                        loc = 112,
                        scale = pmax(1e-6, scale-112/pi^2),
                        shape = shape))
   newTest$Qtruth <- qtrue2
   save(newTest, file = "../outputs/C1_check_coverage.RData")
}