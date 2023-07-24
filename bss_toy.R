# toy function to see if BSS-ANOVA code uses a 
# higher order interaction component f_0

library(fields)

a <- 1
nburn <- 1000
ntot <- 2000

WRITE <- F
PDF <- F
if(PDF) pdf(paste0("pdf/bss_toy_a",a,"_ntot",ntot,".pdf"))

# generate training grid
nx <- 11
x <- seq(0,1,length.out=nx)
xmat <- expand.grid(x,x,x)

# generate prediction grid
np <- nx - 1
x_pred <- seq(0.01,0.99,length.out=np)
xmat_pred <- expand.grid(x_pred, x_pred, x_pred)

# deterministic output
y <- 1.5*cos(pi*xmat[,2]) + a*(xmat[,1]-.5)*cos(pi*xmat[,2]) + 5

# make y mean zero first
ybar <- mean(y)
y <- y - mean(y)

yr <- range(y)

par(mfrow=c(1,1))
y_arr <- array(y + ybar, c(nx,nx))
image.plot(y_arr, main = "true surface", xlab = expression(x[1]), ylab = expression(x[2]))
persp(y_arr, theta = 230, phi = 35)

# obtain main effect estimates
# average out all other params
ME_x1 <- ME_x2 <- ME_x3 <- c()
for (i in 1:nx) {
  ME_x1[i] <- mean(y[which(xmat[,1]==x[i])]) #- ybar
  ME_x2[i] <- mean(y[which(xmat[,2]==x[i])]) #- ybar
  ME_x3[i] <- mean(y[which(xmat[,3]==x[i])]) #- ybar
}

# convert vector to a matrix
# each of these assume that variable is on the x axis
# transpose to have as the y axis
ME_x1_mat <- matrix(rep(ME_x1, nx), nx, nx)
ME_x2_mat <- matrix(rep(ME_x2, nx), nx, nx)
ME_x3_mat <- matrix(rep(ME_x3, nx), nx, nx)

# obtain two-way interactions
TWI_12 <- TWI_13 <- TWI_23 <- matrix(NA, nx, nx)
for (i in 1:nx) {
  for (j in 1:nx) {
    TWI_12[i,j] <- mean(y[which(xmat[,1]==x[i] & xmat[,2]==x[j])]) -
      ME_x1[i] - ME_x2[j] #+ ybar
    TWI_13[i,j] <- mean(y[which(xmat[,1]==x[i] & xmat[,3]==x[j])]) -
      ME_x1[i] - ME_x3[j] #+ ybar
    TWI_23[i,j] <- mean(y[which(xmat[,2]==x[i] & xmat[,3]==x[j])]) -
      ME_x2[i] - ME_x3[j] #+ ybar
  }
}

image.plot(TWI_12, zlim = yr)
image.plot(TWI_13, zlim = yr)
image.plot(TWI_23, zlim = yr)

persp(TWI_12, zlim = yr, theta = 230, phi = 35,
      xlab = "x1", ylab = "x2", zlab = "y")
persp(TWI_13, zlim = yr, theta = 230, phi = 35,
           xlab = "x1", ylab = "x3", zlab = "y")
persp(TWI_23, zlim = yr, theta = 230, phi = 35,
           xlab = "x2", ylab = "x3", zlab = "y")

# reconstruct y on a 2D grid with x1,x2
y_rec <- ybar + ME_x1_mat + t(ME_x2_mat) + ME_x3_mat + TWI_12 + TWI_13 + TWI_23

par(mfrow=c(1,2))
image.plot((y_rec), main = "recreated")
image.plot(y_arr, main = "true surface")
persp((y_rec), main = "recreated", theta = 230, phi = 35,
      xlab = "x1", ylab = "x2", zlab = "y")
persp(y_arr, main = "true surface", theta = 230, phi = 35,
      xlab = "x1", ylab = "x2", zlab = "y")

all.equal(y_rec, y_arr)

# plot each main effect
par(mfrow=c(1,3))
plot(xmat[,1],y, col="gray", xlab = expression(x[1])); lines(x, rowMeans(y_arr - ybar), col="blue")
plot(xmat[,2],y, col="gray", xlab = expression(x[2])); lines(x, colMeans(y_arr - ybar), col="blue")
plot(xmat[,3],y, col="gray", xlab = expression(x[3])); lines(x, rep(0, length(x)), col="blue")

source("bssanova.R")

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

if(WRITE) save.image(paste0("/projects/precipit/hydro/bss_toy_a",a,"_ntot",ntot,".rda"))

