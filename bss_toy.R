# toy function to see if BSS-ANOVA code uses a 
# higher order interaction component f_0

a <- 1
nburn <- 1000
ntot <- 2000

PDF <- T
if(PDF) pdf(paste0("pdf/bss_toy_a",a,"_ntot",ntot,".pdf"))

# generate training grid
nx <- 25
x <- seq(0,1,length.out=nx)
xmat <- expand.grid(x,x,x)

# generate prediction grid
np <- nx - 1
x_pred <- seq(0.01,0.99,length.out=np)
xmat_pred <- expand.grid(x_pred, x_pred, x_pred)

# deterministic output
y <- 1.5*cos(pi*xmat[,2]) + a*(xmat[,1]-.5)*cos(pi*xmat[,2])

par(mfrow=c(1,1))
y_arr <- array(y, c(nx,nx))
image.plot(y_arr, main = "true surface", xlab = expression(x[1]), ylab = expression(x[2]))
# persp(y_arr, theta = 45, phi = 50, col = c(1:5))

# plot each main effect
par(mfrow=c(1,3))
plot(xmat[,1],y, col="gray", xlab = expression(x[1])); lines(x, rowMeans(y_arr), col="blue")
plot(xmat[,2],y, col="gray", xlab = expression(x[2])); lines(x, colMeans(y_arr), col="blue")
plot(xmat[,3],y, col="gray", xlab = expression(x[3])); lines(x, rep(0, length(x)), col="blue")

source("../bssanova.R")

#######################################
# MAIN EFFECTS AND TWO-WAY INITERACNS #
#######################################
b1 = bssanova(xmat,y, BTE = c(nburn,ntot,10))
b1_pred = predict.bssanova(xmat_pred, b1)
ypred_arr <- array(b1_pred$yhat,c(np,np))
image.plot(ypred_arr, main = "prediction (MEs and 2WIs)")

# create prediction results based on "curves"; (trying to find higher order f_0)
nmcmc <- nrow(b1_pred$yreal)
my_hats <- my_hats_meo <- matrix(NA, nrow = nrow(xmat_pred), ncol = nmcmc)
for(cu in 1:nmcmc) my_hats[,cu] <- 
  b1_pred$curves[1,cu,]+b1_pred$curves[2,cu,]+
  b1_pred$curves[3,cu,]+b1_pred$curves[4,cu,]+
  b1_pred$curves[5,cu,]+b1_pred$curves[6,cu,]

# average over the MCMC draws
my_preds = rowMeans(my_hats)

# plot the results
mypred_arr <- array(my_preds,c(np,np))
image.plot(mypred_arr, main = "my prediction")
image.plot(mypred_arr - ypred_arr, main = "prediction difference")
hist(mypred_arr - ypred_arr)

#####################
# MAIN EFFECTS ONLY #
#####################
b1_meo = bssanova(xmat,y,int.order = 1, BTE = c(nburn,ntot,10))
b1_pred_meo = predict.bssanova(xmat_pred, b1_meo)
ypred_arr_meo <- array(b1_pred_meo$yhat,c(np,np))
image.plot(ypred_arr_meo, main = "prediction (main effects only)")

# this is for a main effects only (meo) fit, so just the three main effects
# cu is which batch of posterior curves we are looking at
for(cu in 1:nmcmc) my_hats_meo[,cu] <- 
  b1_pred_meo$curves[1,cu,]+
  b1_pred_meo$curves[2,cu,]+
  b1_pred_meo$curves[3,cu,]

# average over the MCMC draws
my_preds_meo = rowMeans(my_hats_meo)

# plot the results
mypred_arr_meo <- array(my_preds_meo,c(np,np))
image.plot(mypred_arr_meo, main = "my prediction (meo)")
image.plot(mypred_arr_meo - ypred_arr_meo, main = "prediction difference (meo)")
hist(mypred_arr_meo - ypred_arr_meo)

if(PDF) dev.off()

save.image(paste0("/projects/precipit/hydro/bss_toy_a",a,"_ntot",ntot,".rda"))

