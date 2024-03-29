
library(GPfit)
library(fields)
library(lhs)

# load("rda/hydro-SA-FF3-nPC4_all64.rda")
# 
# dim(bs[[1]]$bss.model$l2)
# head(bs[[1]]$bss.model$l2)
# 
# dim(bms[[1]]$bss.model$l2)
# head(bms[[1]]$bss.model$l2)
# 
# dim(bs[[1]]$bss.model$r)
# head(bs[[1]]$bss.model$r)
# 
# my_order <- c()
# for (j in 1:p) {
#   my_order[j] <- 
#   sum(-colMeans(bms[[1]]$bss.model$l2)
#       [which(bms[[1]]$bss.model$term1 == j | bms[[1]]$bss.model$term2 == j)])
# }
# 
# bms[[1]]$order
# order(my_order)
# 
# for (j in 1:p) {
#   my_order[j] <- 
#     sum(-colMeans(bs[[1]]$bss.model$l2)
#         [which(bs[[1]]$bss.model$term1 == j | bs[[1]]$bss.model$term2 == j)])
# }
# order(my_order)
# bs[[1]]$order


tic <- proc.time()[3]

ntrain <- 100
N <- 1000 # npred

# toy example
# prediction grid/true grid
l1 <- l2 <- 11
x1 <- seq(0,10,length=l1)/10
x2 <- seq(0,pi/2, length=l2)/(pi/2)

y1true <- y2true <- matrix(NA, l1, l2)
for (i in 1:l1) for (j in 1:l2) y1true[i,j] <-  10*x1[i] * cos(pi/2*x2[j])
for (i in 1:l1) for (j in 1:l2) y2true[i,j] <-  10*x1[i] * sin(pi/2*x2[j])

image.plot(y1true, main = "ytrue", asp=1); contour(y1true, add=T, col=0)
image.plot(y2true, main = "ytrue", asp=1); contour(y2true, add=T, col=0)

B <- 10
d <- 2
indices <- matrix(NA, B, d + choose(d,2))
summies <- c()
for (b in 1:B) {
  print(b)
  
  xtrain <- randomLHS(n = ntrain, k = 2)
  y1train <- y2train <- c()
  for (i in 1:nrow(xtrain)) y1train[i] <- 10*xtrain[i,1] * cos(pi/2*xtrain[i,2])
  for (i in 1:nrow(xtrain)) y2train[i] <- 10*xtrain[i,1] * sin(pi/2*xtrain[i,2])
  
  # compare sensitivity estimates with values in Gamboa (pg 4)
  xpred <- randomLHS(n = N, k = 2)
  a1  <- GP_fit(X = xtrain, Y = y1train)
  ap1 <- predict(a1, xpred)
  # y1predmat <- matrix(ap1$Y_hat, l1, l2)
  # image.plot(y1predmat, main = "ypred")
  # contour(y1predmat, add=T,col="gray")
  # contour(y1true, add=T,col=0)
  
  a2  <- GP_fit(X = xtrain, Y = y2train)
  ap2 <- predict(a2, xpred)
  # y2predmat <- matrix(ap2$Y_hat, l1, l2)
  # image.plot(y2predmat)
  # contour(y2predmat, add=T,col="gray")
  # contour(y2true, add=T,col=0)
  
  Ey <- mean(ap1$Y_hat)
  Vary <- (t(ap1$Y_hat) %*% ap1$Y_hat)/nrow(xpred) - Ey^2
  
  # first order sensitivity indices
  d <- 2
  Mprime <- randomLHS(N, d)
  S <- EE2j <- rep(NA, d)
  for (j in 1:d) {
    Mjprime <- Mprime
    Mjprime[,j] <- xpred[,j] #change later perhaps...
    pMprime <- predict(a1, Mjprime)
    EE2j[j] <- (t(ap1$Y_hat) %*% pMprime$Y_hat)/(N-1)
    S[j] <- (EE2j[j] - Ey^2)/Vary
  }
  
  S
  
  # second order sensitivity indices
  S_2wis <- EE2j <- t1 <- t2 <- rep(NA, choose(d,2))
  it <- 1
  for (j in 1:(d-1)) {
    for (k in (j+1):d) {
      Mjprime <- Mprime
      Mjprime[,c(j,k)] <- xpred[,c(j,k)] #change later perhaps...
      pMprime <- predict(a1, Mjprime)
      EE2j[it] <- (t(ap1$Y_hat) %*% pMprime$Y_hat)/(N-1)
      S_2wis[it] <- (EE2j[it] - Ey^2)/Vary
      t1 <- j; t2 <- k
      it <- it + 1
    }
  }
  
  (S_2wi_adjs <- S_2wis - sum(S))
  
  indices[b,] <- c(S, S_2wi_adjs)
  summies[b] <- sum(c(S, S_2wi_adjs))
}

par(mfrow=c(2,2))
hist(indices[,1]); abline(v=0.52, col="blue")
hist(indices[,2]); abline(v=0.36, col="blue")
hist(indices[,3]); abline(v=0.12, col="blue")
hist(summies); abline(v=1, col="blue")

toc <- proc.time()[3]

toc - tic
(toc - tic)/3600
