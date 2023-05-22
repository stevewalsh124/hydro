# Obtain main effects estimates in another way
# Right now, *only* trying to get a plot for x1's main effect
# Will extend to x2, x3, and x4 afterwards.

PDF <- F
if(PDF) pdf("pdf/UQ_MEs.pdf")

library(fields)
library(LaplacesDemon)

# This loads results from a GP_fit for each FPC
load("rda/hydro-SA-FF3-nv9-nPC4_all64.rda")

# this subsets the part of the full factorial design that varies x1
# by nv = 9 (e.g.) values, while x2,x3,x4 are varied by n_ff (e.g., 3)
facDes <- facDes[1:((p-1)^3*nv),]

# p,      the number of input parameters (4)
# n_pc,   the number of FPCs used to model the functional output (4)
# facDes, the prediction grid (FF padded with nv elements for each x_i)
#         facDes depends on p=4, nv (e.g., 9) and n_ff (e.g., 3)
# as,     a list of GP_fits (with length = n_pc)
# bases,  the FPCs (each column of bases is an FPC)
# smvals, stellar mass values; the x axis of the functional output

# This function outputs the estimated main effect for x1 (the first ME)
# for each of the n_pc = 4 functional principal components (FPCs)
# The GP fits estimate the weights for each corresponding FPC

get_first_ME <- function(as, plot = F){
  MEs_by_FPC <- list()
  if(p != length(as)) stop("dimension mismatch :(")
  
  # for each FPC (ie, for each column of "bases")...
  for (i_pc in 1:n_pc) {
    
    # Get the GP fit for the i'th FPC (i_pc from 1 to n_pc)
    a1 <- as[[i_pc]] 
    sig2 <- a1$sig2 #marginal variance
    Delta <- diag(rep(a1$delta, nrow(plgp::distance(a1$X)))) #nugget mtx
    
    # create covariance matrix Sigma for the training data
    # based on separable fit (here, p goes from 1 to 4)
    R1 <- function(D) exp(-10^(a1$beta[1]) * D)
    R2 <- function(D) exp(-10^(a1$beta[2]) * D)
    R3 <- function(D) exp(-10^(a1$beta[3]) * D)
    R4 <- function(D) exp(-10^(a1$beta[4]) * D)
    
    # squared distances of the design points for each dim
    Dx1 <- plgp::distance(a1$X[,1])
    Dx2 <- plgp::distance(a1$X[,2])
    Dx3 <- plgp::distance(a1$X[,3])
    Dx4 <- plgp::distance(a1$X[,4])
    
    # the ntrain x ntrain covariance matrices for each dim (64x64)
    R1x1 <- R1(Dx1)
    R2x2 <- R2(Dx2)
    R3x3 <- R3(Dx3)
    R4x4 <- R4(Dx4)
    
    # combine to get the separable/anisotropic covariance matrix
    # for the training data
    Sigma <- (sig2 * R1x1*R2x2*R3x3*R4x4) + Delta
    if(plot) image.plot(Sigma, main = paste("train, PC", i_pc))
    
    # Now for the prediction grid...
    # facDes is the factorial design for prediction locations
    npred <- nrow(facDes) # npred = n1pred = ... = n4pred
    
    # Get squared distances for each column of x in the factorial design
    Dx1p <- plgp:::distance(facDes[,1])
    Dx2p <- plgp:::distance(facDes[,2])
    Dx3p <- plgp:::distance(facDes[,3])
    Dx4p <- plgp:::distance(facDes[,4])
    
    # and construct the corresponding covariance matrix for the pred locns
    R1x1p <- R1(Dx1p)
    R2x2p <- R2(Dx2p)
    R3x3p <- R3(Dx3p)
    R4x4p <- R4(Dx4p)
    
    # npred x npred covariance matrix for the prediction (facDes) locations
    # **problem** the diagonal elements != sig2 when I have the means
    # of the other dimensions included...
    Sigma_X1M <- sig2 * R1x1p #* mean(R2x2p) * mean(R3x3p) * mean(R4x4p)
    range(Sigma_X1M)
    if(plot) image.plot(Sigma_X1M, main = "pred", asp=1)
    # dim(Sigma_X1M)
    
    # Final component: the cross covariance matrix between
    # training and prediction locations
    # Same **problem** as above; if I divide each auxiliary dimension
    # by n_ff or npred, etc, the results look very wrong,
    # and the max value will not equal sig2
    CX1MX1 <- sig2 * R1(plgp:::distance(facDes[,1], a1$X[,1])) *
      (1/1 * (R2(plgp:::distance(facDes[,2], a1$X[,2])))) *
      (1/1 * (R3(plgp:::distance(facDes[,3], a1$X[,3])))) *
      (1/1 * (R4(plgp:::distance(facDes[,4], a1$X[,4]))))
    range(CX1MX1)
    if(plot) image.plot(CX1MX1, asp=1, main = "cross")
    # dim(CX1MX1)
    
    # make the complete covariance matrix for the joint [ypred, ytrain]
    # recall, y here is the weights for the corresponding FPC
    BigSig <- rbind( cbind(Sigma_X1M, CX1MX1),
                     cbind(t(CX1MX1), Sigma) )
    
    if(plot) image.plot(t(BigSig), main = "BigSig")
    
    # **problem** I subtracted 0.5, otherwise results didn't seem sensible...
    condl_mean <- CX1MX1 %*% solve(Sigma) %*% a1$X[,1] - 0.5
    
    # this subsets only the nv unique values for the weights
    # using the part of the facDes where x1 is varied by nv
    condl_avg <- c()
    for (i in 1:nv) condl_avg[i] <- mean(condl_mean[((i-1)*n_ff^(p-1))+1:(n_ff^(p-1))])
    
    # plot the weights (black), vs the average weight for each setting of x1 (red)
    if(plot) {plot(condl_mean, type="l"); plot(condl_mean[1:(nv*(p-1)^3)], type="l")}
    if(plot) for(i in 1:nv) lines(((i-1)*n_ff^(p-1))+1:(n_ff^(p-1)), 
                                  rep(condl_avg[i],length(((i-1)*n_ff^(p-1))+1:(n_ff^(p-1)))),
                                  col="red")
    if(plot) plot(condl_avg, col="red")
    
    # covariance matrix from kriging equations
    condl_Sig <-  Sigma_X1M - CX1MX1 %*% solve(Sigma) %*% t(CX1MX1)
    
    if(plot) matplot(t(outer(c(condl_mean), bases[,i_pc])), type="l", 
                     ylim=range(t(outer(c(condl_mean), bases[,i_pc]))))
    if(plot) matplot(t(outer(c(condl_avg), bases[,i_pc])), type="l")
    
    MEs_by_FPC[[i_pc]] <- t(outer(c(condl_avg), bases[,i_pc]))
  }
  return(MEs_by_FPC)
}

# Function returns a list, where the i'th list element
# is the outer product of the nv weights 
# with the corresponding (i'th) basis function
firstME <- get_first_ME(as, plot = T)
# for (i_pc in 1:n_pc) matplot(firstME[[i_pc]], type="l", main = i_pc)

gr_colors <- paste0("gray",round(seq(20,80, length=nv)))
matplot(smvals, Reduce("+",firstME), type="l", lty=1, 
        col = gr_colors)

if(PDF) dev.off()

