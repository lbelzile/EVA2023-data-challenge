setwd(this.path::here())
remotes::install_github("samperochkin/tautest")
library(tautests)

UtopulaU1_data <- read.csv(file = "../data/UtopulaU1.csv", header = TRUE)
colnames(UtopulaU1_data) <- paste0("I1_", c(paste0(0, 1:9), 10:25))
UtopulaU2_data <- read.csv(file = "../data/UtopulaU2.csv", header = TRUE)
colnames(UtopulaU2_data) <- paste0("I2_", c(paste0(0, 1:9), 10:25))
# Merge databases
Utopula <- data.frame(UtopulaU1_data, UtopulaU2_data)


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

clust_order <- order(as.integer(paste0(hclust_od1, c(paste0(0, 1:9), 11:50))))
# Use block structured data
X <- as.matrix(Utopula[,clust_order])
n <- nrow(X)
d <- ncol(X)
p <- choose(d, 2)

# Function from Samuel Perreault
constrainSigma <- function(Sh, clus, debug = FALSE) {
   #### fixed parameters ####
   d <- length(clus)
   p <- d * (d - 1) / 2
   K <- length(unique(clus))
   
   #### vectorizing ####
   
   ij <- tautests::R2IJ(1:p)
   kl <- matrix(clus[c(ij)], p, 2) |> apply(1, sort) |> t()
   kl0 <- unique(kl)
   L <- nrow(kl0)
   
   s_kl <- kl[, 1] + (K + 1) * kl[, 2]
   s_kl0 <- unique(s_kl) # should be in same order as u_kl
   s_kl <- match(s_kl, s_kl0)
   if (length(s_kl0) != nrow(kl0))
      cat("problem in generating unique block ids\n")
   
   # loop to make sure we don't run into memory problems
   for (u in 1:L) {
      for (v in u:L) {
         kl_u <- kl0[u, ]
         kl_v <- kl0[v, ]
         
         if (debug)
            cat("----------\n case", kl_u, "and", kl_v, "\n")
         
         r_u <- which(s_kl == u)
         s_v <- which(s_kl == v)
         
         if (length(intersect(kl_u, kl_v)) == 0) {
            # no overlap possible
            Sh[r_u, s_v] <- Sh[s_v, r_u] <- mean(Sh[r_u, s_v])
            next
         }
         
         # collect entries from both blocks
         ij_u <- ij[r_u, ]
         ij_v <- ij[s_v, ]
         p_u <- nrow(ij_u)
         p_v <- nrow(ij_v)
         
         psi <- mapply(
            FUN = function(x, y) {
               psi <- clus[intersect(ij_u[x, ], ij_v[y, ])]
               c(rep(0, 2 - length(psi)), sort(psi))
            },
            x = rep(1:p_u, times = p_v),
            y = rep(1:p_v, each = p_u)
         ) |>
            t() |> c() |> array(dim = c(p_u, p_v, 2))
         
         s_psi <- psi[, , 1] + (L + 1) * psi[, , 2]
         groups <- split(1:length(s_psi), s_psi)
         
         for (g in groups) {
            s <- ceiling(g / p_u)
            r <- g - (s - 1) * p_u
            Sh[cbind(r_u[r], s_v[s])] <-
               Sh[cbind(s_v[s], r_u[r])] <- mean(Sh[cbind(r_u[r], s_v[s])])
            # image(t(Sh)[rev(r_u),s_v])
         }
         
         if (debug) {
            l1 <- min(length(unique(kl_u)), length(unique(kl_v)))
            l2 <- max(length(unique(kl_u)), length(unique(kl_v)))
            l12 <- length(intersect(kl_u, kl_v))
            
            category <- which(
               c(
                  l1 == 1 & l2 == 1 & l12 == 1,
                  l1 == 1 & l2 == 2 & l12 == 1,
                  l1 == 1 & l2 == 1 & l12 == 0,
                  l1 == 2 & l2 == 2 & l12 == 2,
                  l1 == 1 & l2 == 2 & l12 == 0,
                  l1 == 2 & l2 == 2 & l12 == 1,
                  l1 == 2 & l2 == 2 & l12 == 0
               )
            )
            cat("produces", length(groups), "groups \n")
            cat("DIFF", length(groups) - c(3, 2, 1, 4, 1, 2, 1)[category], "\n")
            cat("ACTUAL DIFF ****",
                length(unique(c(Sh[r_u, s_v]))) - c(3, 2, 1, 4, 1, 2, 1)[category],
                "****\n")
         }
         
      }
   }
   
   return(Sh)
}



# compute tau and its jackknife variance in O(n*log(n))
TSh <- tautests::tau_and_jack(X)
Sh <- TSh$jack_var
Th <- TSh$tau
image(t(Th)[, d:1])
th <- Th[tautests::R2IJ(1:p)]

hc <- hclust(as.dist(1 - abs(Th)), "ward.D2")
plot(hc)
K <- 5
clus <- cutree(hc, K)

# structure Sigma and compute inverse and principal square root
St <- constrainSigma(Sh, clus)
image(t(St)[, p:1])

eig <- eigen(St)
plot(eig$values)
Sti <- eig$vectors %*% diag(1 / eig$values) %*% t(eig$vectors)
Sti2 <- eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)

# checkup (when all clusters are greater than 3)
length(unique(c(St))) == (K ^ 4 + 6 * K ^ 3 + 11 * K ^ 2 + 6 * K) / 8


# construct membership matrix B
ij <- tautests::R2IJ(1:p)
A <- matrix(0, d, K)
A[cbind(1:d, clus)] <- 1
A <- t(apply(ij, 1, function(x)
   colSums(A[x, ])))
A_u <- unique(A)
B <- matrix(0, nrow(A), nrow(A_u))
B[cbind(1:p, apply(A, 1, function(a)
   which.max(colSums(a == t(
      A_u
   )))))] <- 1
B <- B[, colSums(B) != 0]
L <- ncol(B)

# Alternative B construction, pooling all off-diagonal entries to zero
B2 <- matrix(0, nrow = p, ncol = 6)
for (i in 1:5) {
   B2[, i] <- tcrossprod(clus == i)[upper.tri(Th)]
}
B2[, 6] <- rowSums(B2) == 0
L <- 6L

# construct Gamma matrix and covariance matrix for Monte Carlo
G <- B2 %*% MASS::ginv(B2)
IG <- diag(p) - G
tt <- G %*% th
Tt <- diag(d)
Tt[ij] <- Tt[ij[, 2:1]] <- tt
par(mfrow = c(1, 2), mar = c(2, 2, 1, 1))
image(t(Th)[, d:1], zlim = range(c(Th, Tt)))
image(t(Tt)[, d:1], zlim = range(c(Th, Tt)))
St_dag <- Sti2 %*% IG %*% St %*% IG %*% Sti2

# p-values
mc_z <- mvtnorm::rmvnorm(1e4, sigma = St_dag)

loss_E <- n * mahalanobis(th, tt, Sti, TRUE)
pval_E <- pchisq(loss_E, p - L, lower.tail = FALSE)
mc_loss_E <- apply(mc_z, 1, function(x)
   c(crossprod(x)))
hist(mc_loss_E)
abline(v = loss_E,
       col = 2,
       lty = 3,
       lwd = 2)
# c(pval_E, mean(loss_E <= mc_loss_E)) # just for checking

loss_M <- sqrt(n) * max(abs(Sti2 %*% (th - tt)))
mc_loss_M <- apply(mc_z, 1, function(x)
   max(abs(x)))
pval_M <- mean(loss_M <= mc_loss_M)
hist(mc_loss_M)
abline(v = loss_M,
       col = 2,
       lty = 3,
       lwd = 2)
pval_M

c(pval_E, pval_M)

# Case for the I/n as scaling factor
St_dag <- (IG %*% St %*% IG)

# p-values
mc_z <- mvtnorm::rmvnorm(1e4, sigma = St_dag)

loss_E <- n * crossprod(th - tt)[1, 1]
mc_loss_E <- apply(mc_z, 1, function(x)
   c(crossprod(x)))
hist(mc_loss_E)
abline(v = loss_E,
       col = 2,
       lty = 3,
       lwd = 2)
pval_E <- mean(loss_E <= mc_loss_E) # just for checking

loss_M <- sqrt(n) * max(abs(th - tt))
mc_loss_M <- apply(mc_z, 1, function(x)
   max(abs(x)))
pval_M <- mean(loss_M <= mc_loss_M)
hist(mc_loss_M)
abline(v = loss_M,
       col = 2,
       lty = 3,
       lwd = 2)
pval_M

c(pval_E, pval_M)
