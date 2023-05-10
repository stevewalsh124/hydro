# learn about BSS-ANOVA output

## This is from the bss-anova function:
## yhat - a nrow(X.new) vector of the posterior mean of predicted y's
## yreal - a (nreal) x (nrow(X.new)) matrix:
##          each row gives X predictions for a given posterior realization 
## curves - a (# functional components) x (nreal) x (nrow(X.new)) array:
##           e.g., curves[j,r,] provides predictions for j-th functional component
##	     for the r-th posterior realization


# run new32_hydro.R first to get a bss-anova fit (b1_pred, ..., b4_pred is four PCs)
load("/projects/precipit/hydro/new32_hydro.rda")

# by default, 800 MCMC draws
# there are n_sm = 11 values of stellar mass, 
# and m = 32 runs in this case
all_hats <- array(NA, c(800,n_sm,m))


# this is for a main effects only (meo) fit, so just the four main effects
# cu is which batch of posterior curves we are looking at
for(cu in 1:800){
  # combine the four main effects for PC1
  my_hats1 <- outer(bases[,1], b1_pred_meo$curves[1,cu,]) + 
    outer(bases[,1], b1_pred_meo$curves[2,cu,]) + 
    outer(bases[,1], b1_pred_meo$curves[3,cu,]) + 
    outer(bases[,1], b1_pred_meo$curves[4,cu,])
  # combine the four main effects for PC2
  my_hats2 <- outer(bases[,2], b2_pred_meo$curves[1,cu,]) + 
    outer(bases[,2], b2_pred_meo$curves[2,cu,]) + 
    outer(bases[,2], b2_pred_meo$curves[3,cu,]) + 
    outer(bases[,2], b2_pred_meo$curves[4,cu,])  
  # combine the four main effects for PC3
  my_hats3 <- outer(bases[,3], b3_pred_meo$curves[1,cu,]) + 
    outer(bases[,3], b3_pred_meo$curves[2,cu,]) + 
    outer(bases[,3], b3_pred_meo$curves[3,cu,]) + 
    outer(bases[,3], b3_pred_meo$curves[4,cu,])
  # combine the four main effects for PC4
  my_hats4 <- outer(bases[,4], b4_pred_meo$curves[1,cu,]) + 
    outer(bases[,4], b4_pred_meo$curves[2,cu,]) + 
    outer(bases[,4], b4_pred_meo$curves[3,cu,]) + 
    outer(bases[,4], b4_pred_meo$curves[4,cu,])
  # combine all of these into the 32 predictions (there are 32 training, but also 32 testing)
  all_hats[cu,,] <- my_hats1 + my_hats2 + my_hats3 + my_hats4
}

# average over the MCMC draws
my_created_preds = apply(all_hats,c(2,3),mean)

# compare this prediction with the output we have
matplot(my_created_preds, type="l", col="gray", main = "compare yhat vs handmade-curves-preds", lwd=2, lty=1,
        ylim = range(c(my_created_preds, new32_preds_bm - mean0mat))) # compare this with 
matplot(new32_preds_bm - mean0mat, type="l", col="green", lty=3, lwd=2, add=T)
legend("bottomleft", c("yhat","curves-pred"), lty=1, col=c("gray","green"))
all.equal(my_created_preds, new32_preds_bm - mean0mat)

# this illustrates 
matplot(b1_pred_meo$yreal, type="n")
for (i in 1:m) {
  abline(h=b1_pred_meo$yhat[i], col=i)
  lines(b1_pred_meo$yreal[,i], col=i)
}

# this is one of the GPs we are fitting
x1_pred <- des_pred[,1]
plot(x1_pred, b1_pred$curves[1,cu,])
points(x1_pred, b1_pred_meo$curves[1,cu,], col="blue")

