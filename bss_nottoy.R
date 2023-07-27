# only look at unique combos
# there is duplication when nx == nv
# Using des=="FFS" fixes this problem

# p <- 4
# n_ff <- 11

library(fields)

load("/projects/precipit/hydro/hydro-SA-lo-trainac-predFFS21-nPC4.rda")
nx <- n_ff

PDF <- T
if(PDF) pdf(paste0("pdf/bss_nottoy_",des,"_",n_ff,".pdf"))

x <- seq(0,1,length.out=nx)
xmat <- expand.grid(x,x,x,x)

# generate prediction grid
# np <- nx - 1
# x_pred <- seq(0.01,0.99,length.out=np)
# xmat_pred <- expand.grid(x_pred, x_pred, x_pred)

# for starters, only look at first stellar mass value
# deterministic output
y <-  etaEmu[1,1:nx^p]
# deterministic output
# y <- 1.5*cos(pi*xmat[,2]) + (xmat[,1]-.5)*cos(pi*xmat[,2]) +
#   xmat[,2]*xmat[,1]^2*cos(pi*(xmat[,3])) + 5 +
#   xmat[,2]^2*xmat[,3]*cos(pi*(xmat[,4]))

# make y mean zero first
ybar <- mean(y)
y <- y - mean(y)

par(mfrow=c(1,1))
y_arr <- array(y + ybar, c(nx,nx))
image.plot(y_arr, main = "estimated surface from FF des", xlab = expression(x[1]), ylab = expression(x[2]))
persp(y_arr, theta = 230, phi = 35)

# obtain main effect estimates
# average out all other params
ME_x1 <- ME_x2 <- ME_x3 <- ME_x4 <- c()
for (i in 1:nx) {
  ME_x1[i] <- mean(y[which(xmat[,1]==x[i])]) #- ybar
  ME_x2[i] <- mean(y[which(xmat[,2]==x[i])]) #- ybar
  ME_x3[i] <- mean(y[which(xmat[,3]==x[i])]) #- ybar
  ME_x4[i] <- mean(y[which(xmat[,4]==x[i])]) #- ybar
}

# convert vector to a matrix
# each of these assume that variable is on the x axis
# transpose to have as the y axis
ME_x1_mat <- matrix(rep(ME_x1, nx), nx, nx)
ME_x2_mat <- matrix(rep(ME_x2, nx), nx, nx)
ME_x3_mat <- matrix(rep(ME_x3, nx), nx, nx)
ME_x4_mat <- matrix(rep(ME_x4, nx), nx, nx)

yr <- range(c(ME_x1, ME_x2, ME_x3, ME_x4))

par(mfrow=c(2,2))
image.plot(ME_x1_mat, zlim = yr)
image.plot(ME_x2_mat, zlim = yr)
image.plot(ME_x3_mat, zlim = yr)
image.plot(ME_x4_mat, zlim = yr)

ME_x1_full <- matrix(rep(ME_x1, nx^(p-1))[1:nx^2],nx,nx)
ME_x2_full <- matrix(rep(rep(ME_x2, each = nx), nx^(p-2))[1:nx^2],nx,nx)
ME_x3_full <- matrix(rep(rep(ME_x3, each = nx^(p-2)), nx)[1:nx^2],nx,nx)
ME_x4_full <- matrix(rep(ME_x4, each = nx^(p-1))[1:nx^2],nx,nx)

# obtain two-way interactions
TWI_12 <- TWI_13 <- TWI_14 <- TWI_23 <- 
  TWI_24 <- TWI_34 <- c()#matrix(NA, nx, nx)
kk <- 1
TWI_inds <- NULL
for (i in 1:nx) {
  for (j in 1:nx) {
    TWI_12[kk] <- mean(y[which(xmat[,1]==x[i] & xmat[,2]==x[j])]) -
      ME_x1[i] - ME_x2[j] #+ ybar
    TWI_13[kk] <- mean(y[which(xmat[,1]==x[i] & xmat[,3]==x[j])]) -
      ME_x1[i] - ME_x3[j] #+ ybar
    TWI_14[kk] <- mean(y[which(xmat[,1]==x[i] & xmat[,4]==x[j])]) -
      ME_x1[i] - ME_x4[j] #+ ybar
    TWI_23[kk] <- mean(y[which(xmat[,2]==x[i] & xmat[,3]==x[j])]) -
      ME_x2[i] - ME_x3[j] #+ ybar
    TWI_24[kk] <- mean(y[which(xmat[,2]==x[i] & xmat[,4]==x[j])]) -
      ME_x2[i] - ME_x4[j] #+ ybar
    TWI_34[kk] <- mean(y[which(xmat[,3]==x[i] & xmat[,4]==x[j])]) -
      ME_x3[i] - ME_x4[j] #+ ybar
    TWI_inds <- rbind(TWI_inds, c(i,j))
    kk <- kk + 1
  }
}

par(mfrow=c(2,3))
image.plot(array(TWI_12,c(nx,nx)), zlim = yr)
image.plot(array(TWI_13,c(nx,nx)), zlim = yr)
image.plot(array(TWI_14,c(nx,nx)), zlim = yr)
image.plot(array(TWI_23,c(nx,nx)), zlim = yr)
image.plot(array(TWI_24,c(nx,nx)), zlim = yr)
image.plot(array(TWI_34,c(nx,nx)), zlim = yr)

# get subsets of 2WIs based on the 2D plot
# for starters, fix x3 and x4 at 0
# since these are all zeros in facDes[1:nx^2,]

# column selection depends on "xyz": TWI_xyz_sub
# any locations of "xyz" that are 3 or 4 (for x3 and x4)
# we want them ==1, since these are x3=0 and x4=0
TWI_12_sub <- TWI_12
TWI_13_sub <- TWI_13[which(TWI_inds[,2]==1)]
TWI_14_sub <- TWI_14[which(TWI_inds[,2]==1)]
TWI_23_sub <- TWI_23[which(TWI_inds[,2]==1)]
TWI_24_sub <- TWI_24[which(TWI_inds[,2]==1)]
TWI_34_sub <- TWI_34[which(TWI_inds[,1]==1 & TWI_inds[,2]==1)]

# repeat the TWIs as appropriate based on the design of facDes[1:nx^2,]
TWI_12_full <- t(array(TWI_12_sub, c(nx,nx)))
TWI_13_full <- t(array(rep(TWI_13_sub, nx), c(nx,nx)))
TWI_14_full <- t(array(rep(TWI_14_sub, nx), c(nx,nx)))
TWI_23_full <- t(array(rep(TWI_23_sub, each = nx), c(nx,nx)))
TWI_24_full <- t(array(rep(TWI_24_sub, each = nx), c(nx,nx)))
TWI_34_full <- t(array(rep(TWI_34_sub, nx^2), c(nx,nx)))

# persp(TWI_12, zlim = yr, theta = 230, phi = 35,
#       xlab = "x1", ylab = "x2", zlab = "y")
# persp(TWI_13, zlim = yr, theta = 230, phi = 35,
#       xlab = "x1", ylab = "x3", zlab = "y")
# persp(TWI_23, zlim = yr, theta = 230, phi = 35,
#       xlab = "x2", ylab = "x3", zlab = "y")

# obtain three-way interactions
TWI_123 <- TWI_124 <- TWI_134 <- TWI_234 <- c()
kk <- 1
HWI_inds <- NULL
for (i in 1:nx) {
  for (j in 1:nx) {
    for (k in 1:nx) {
      TWI_123[kk] <- mean(y[which(xmat[,1]==x[i] & xmat[,2]==x[j] & xmat[,3]==x[k])]) -
        TWI_12[which(TWI_inds[,1]==i & TWI_inds[,2]==j)] - 
        TWI_13[which(TWI_inds[,1]==i & TWI_inds[,2]==k)] - 
        TWI_23[which(TWI_inds[,1]==j & TWI_inds[,2]==k)] - 
        ME_x1[i] - ME_x2[j] - ME_x3[k]
      
      TWI_124[kk] <- mean(y[which(xmat[,1]==x[i] & xmat[,2]==x[j] & xmat[,4]==x[k])]) -
        TWI_12[which(TWI_inds[,1]==i & TWI_inds[,2]==j)] - 
        TWI_14[which(TWI_inds[,1]==i & TWI_inds[,2]==k)] - 
        TWI_24[which(TWI_inds[,1]==j & TWI_inds[,2]==k)] - 
        ME_x1[i] - ME_x2[j] - ME_x4[k]
      
      TWI_134[kk] <- mean(y[which(xmat[,1]==x[i] & xmat[,3]==x[j] & xmat[,4]==x[k])]) -
        TWI_13[which(TWI_inds[,1]==i & TWI_inds[,2]==j)] - 
        TWI_14[which(TWI_inds[,1]==i & TWI_inds[,2]==k)] - 
        TWI_34[which(TWI_inds[,1]==j & TWI_inds[,2]==k)] - 
        ME_x1[i] - ME_x3[j] - ME_x4[k]
      
      TWI_234[kk] <- mean(y[which(xmat[,2]==x[i] & xmat[,3]==x[j] & xmat[,4]==x[k])]) -
        TWI_23[which(TWI_inds[,1]==i & TWI_inds[,2]==j)] - 
        TWI_24[which(TWI_inds[,1]==i & TWI_inds[,2]==k)] - 
        TWI_34[which(TWI_inds[,1]==j & TWI_inds[,2]==k)] - 
        ME_x2[i] - ME_x3[j] - ME_x4[k]
      HWI_inds <- rbind(HWI_inds, c(i,j,k))
      kk <- kk + 1
    }
  }
}

par(mfrow=c(2,2))
hist(TWI_123);abline(v=mean(TWI_123), col="blue")
hist(TWI_124);abline(v=mean(TWI_124), col="blue")
hist(TWI_134);abline(v=mean(TWI_134), col="blue")
hist(TWI_234);abline(v=mean(TWI_234), col="blue")

# get subsets of 3WIs based on the 2D plot
# for starters, fix x1 and x4 at 0
# since these are all zeros in facDes[1:nx^2,]

# column selection depends on "xyz": TWI_xyz_sub
# any locations of "xyz" that are 1 or 4 (for x1 and x4)
# we want them ==1, since these are x1=0 and x4=0
TWI_123_sub <- TWI_123[which(HWI_inds[,3]==1)]
TWI_124_sub <- TWI_124[which(HWI_inds[,3]==1)]
TWI_134_sub <- TWI_134[which(HWI_inds[,2]==1 & HWI_inds[,3]==1)]
TWI_234_sub <- TWI_234[which(HWI_inds[,2]==1 & HWI_inds[,3]==1)]

# repeat the values as appropriately, based on facDes[1:nx^2,]
TWI_123_full <- t(array(TWI_123_sub,c(nx,nx)))
TWI_124_full <- t(array(TWI_124_sub,c(nx,nx)))
TWI_134_full <- t(array(rep(TWI_134_sub, nx),c(nx,nx)))
TWI_234_full <- t(array(rep(TWI_234_sub, each = nx),c(nx,nx)))

# obtain four-way interaction
FWI_1234 <- c()
kk <- 1
FWI_inds <- NULL
for (i in 1:nx) {
  for (j in 1:nx) {
    for (k in 1:nx) {
      for (l in 1:nx) {
        FWI_1234[kk] <- mean(y[which(xmat[,1]==x[i] & xmat[,2]==x[j] & 
                                       xmat[,3]==x[k] & xmat[,4]==x[l])]) -
          TWI_123[which(HWI_inds[,1]==i & HWI_inds[,2]==j & HWI_inds[,3]==k)] - 
          TWI_124[which(HWI_inds[,1]==i & HWI_inds[,2]==j & HWI_inds[,3]==l)] - 
          TWI_134[which(HWI_inds[,1]==i & HWI_inds[,2]==k & HWI_inds[,3]==l)] - 
          TWI_234[which(HWI_inds[,1]==j & HWI_inds[,2]==k & HWI_inds[,3]==l)] - 
          TWI_12[which(TWI_inds[,1]==i & TWI_inds[,2]==j)] - 
          TWI_13[which(TWI_inds[,1]==i & TWI_inds[,2]==k)] - 
          TWI_14[which(TWI_inds[,1]==i & TWI_inds[,2]==l)] - 
          TWI_23[which(TWI_inds[,1]==j & TWI_inds[,2]==k)] - 
          TWI_24[which(TWI_inds[,1]==j & TWI_inds[,2]==l)] - 
          TWI_34[which(TWI_inds[,1]==k & TWI_inds[,2]==l)] - 
          ME_x1[i] - 
          ME_x2[j] - 
          ME_x3[k] - 
          ME_x4[l]
        
        FWI_inds <- rbind(FWI_inds, c(i,j,k,l))
        kk <- kk + 1
      }
    }
  }
}

FWI_1234_full <- t(array(FWI_1234[which(FWI_inds[,3]==1 & FWI_inds[,4]==1)],c(nx,nx)))

# look at the 2D projections of each component with respect to x1 and x2
par(mfrow=c(2,2))
image.plot(ME_x1_full, main ="x1", zlim=yr)
image.plot(ME_x2_full, main ="x2", zlim=yr)
image.plot(ME_x3_full, main ="x3", zlim=yr)
image.plot(ME_x4_full, main ="x4", zlim=yr)

par(mfrow=c(2,3))
image.plot(TWI_12_full, main ="x1,x2", zlim=yr)
image.plot(TWI_13_full, main ="x1,x3", zlim=yr)
image.plot(TWI_14_full, main ="x1,x4", zlim=yr)
image.plot(TWI_23_full, main ="x2,x3", zlim=yr)
image.plot(TWI_24_full, main ="x2,x4", zlim=yr)
image.plot(TWI_34_full, main ="x3,x4", zlim=yr)


par(mfrow=c(2,2))
image.plot(TWI_123_full, main = "x1,x2,x3", zlim=yr)
image.plot(TWI_124_full, main = "x1,x2,x4", zlim=yr)
image.plot(TWI_134_full, main = "x1,x3,x4", zlim=yr)
image.plot(TWI_234_full, main = "x2,x3,x4", zlim=yr)

image.plot(FWI_1234_full, main = "x1,x2,x3,x4", zlim=yr)

# reconstruct y on a 2D grid with x1,x2
# y_rec <- 
#   ybar + 
#   matrix(ME_x1_full,nx,nx) + matrix(ME_x2_full,nx,nx) + 
#   matrix(ME_x3_full,nx,nx) + matrix(ME_x4_full,nx,nx)  +
#   (matrix(TWI_12_full,nx,nx)) + matrix(TWI_13_full,nx,nx) +
#   matrix(TWI_14_full,nx,nx) + (matrix(TWI_23_full,nx,nx)) +
#   t(matrix(TWI_24_full,nx,nx)) + t(matrix(TWI_34_full,nx,nx)) +
#   (matrix(TWI_123_full,nx,nx)) + (matrix(TWI_124_full,nx,nx)) +
#   (matrix(TWI_134_full,nx,nx)) + (matrix(TWI_234_full,nx,nx)) +
#   matrix(FWI_1234_full,nx,nx)
good_boys <- NULL
for (i1 in 0:1) {
  for (i2 in 0:1) {
    for (i3 in 0:1) {
      for (i4 in 0:1) {
        for(i5 in 0:1){
          for (i6 in 0:1) {
            for (i7 in 0:1) {
              for (i8 in 0:1) {
                for (i9 in 0:1) {
                  for (i10 in 0:1) {
                    for (i11 in 0:1) {

                        if(i1) i1p <- TWI_12_full else i1p <- t(TWI_12_full)
                        if(i2) i2p <- TWI_13_full else i2p <- t(TWI_13_full)
                        if(i3) i3p <- TWI_14_full else i3p <- t(TWI_14_full)
                        if(i4) i4p <- TWI_23_full else i4p <- t(TWI_23_full)
                        if(i5) i5p <- TWI_24_full else i5p <- t(TWI_24_full)
                        if(i6) i6p <- TWI_34_full else i6p <- t(TWI_34_full)
                        if(i7) i7p <- TWI_123_full else i7p <- t(TWI_123_full)
                        if(i8) i8p <- TWI_124_full else i8p <- t(TWI_124_full)
                        if(i9) i9p <- TWI_134_full else i9p <- t(TWI_134_full)
                        if(i10) i10p <- TWI_234_full else i10p <- t(TWI_234_full)
                        if(i11) i11p <- FWI_1234_full else i11p <- t(FWI_1234_full)
                        
                        y_rec <- 
                          ybar + 
                          ME_x1_full + ME_x2_full + ME_x3_full + ME_x4_full +
                          i1p + i2p + i3p + i4p + i5p + i6p + 
                          i7p + i8p + i9p + i10p + i11p
                if(sum(abs(y_rec - y_arr)) < 1e-10){
                  skrt <- (c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11))
                  good_boys <- rbind(good_boys, skrt)
                }
                    
              }}}}
            }
          }
        }
      }
    }
  }
}

colnames(good_boys) <- c("x1x2","x1x3","x1x4","x2x3","x2x4","x3x4","x123","x124","x134","x234","x1234")
my_is <- good_boys[1,]
i1 <- my_is[1]
i2 <- my_is[2]
i3 <- my_is[3]
i4 <- my_is[4]
i5 <- my_is[5]
i6 <- my_is[6]
i7 <- my_is[7]
i8 <- my_is[8]
i9 <- my_is[9]
i10 <- my_is[10]
i11 <- my_is[11]
if(i1) i1p <- TWI_12_full else i1p <- t(TWI_12_full)
if(i2) i2p <- TWI_13_full else i2p <- t(TWI_13_full)
if(i3) i3p <- TWI_14_full else i3p <- t(TWI_14_full)
if(i4) i4p <- TWI_23_full else i4p <- t(TWI_23_full)
if(i5) i5p <- TWI_24_full else i5p <- t(TWI_24_full)
if(i6) i6p <- TWI_34_full else i6p <- t(TWI_34_full)
if(i7) i7p <- TWI_123_full else i7p <- t(TWI_123_full)
if(i8) i8p <- TWI_124_full else i8p <- t(TWI_124_full)
if(i9) i9p <- TWI_134_full else i9p <- t(TWI_134_full)
if(i10) i10p <- TWI_234_full else i10p <- t(TWI_234_full)
if(i11) i11p <- FWI_1234_full else i11p <- t(FWI_1234_full)

y_rec <- 
  ybar + 
  ME_x1_full + ME_x2_full + ME_x3_full + ME_x4_full +
  i1p + i2p + i3p + i4p + i5p + i6p + 
  i7p + i8p + i9p + i10p + i11p
  
  sum(abs(y_rec - y_arr))
hist(y_rec - y_arr)
all.equal(y_rec, y_arr)
image.plot(y_rec - y_arr, main = "difference")

par(mfrow=c(1,2))
image.plot(y_rec, main = "recreated", zlim = range(y_rec, y_arr))
image.plot(y_arr, main = "original FF surface", zlim = range(y_rec, y_arr))
persp((y_rec), main = "recreated", theta = 230, phi = 35,
      xlab = "x1", ylab = "x2", zlab = "y")
persp(y_arr, main = "original FF surface", theta = 230, phi = 35,
      xlab = "x1", ylab = "x2", zlab = "y")

all.equal(y_rec, y_arr)
round(y_rec - y_arr, 10)

# par(mfrow=c(1,3))
# plot(xmat[,1],y, col="gray", xlab = expression(x[1])); lines(x, rowMeans(y_arr - ybar), col="blue")
# plot(xmat[,2],y, col="gray", xlab = expression(x[2])); lines(x, colMeans(y_arr - ybar), col="blue")
# plot(xmat[,3],y, col="gray", xlab = expression(x[3])); lines(x, rep(0, length(x)), col="blue")
# plot(xmat[,4],y, col="gray", xlab = expression(x[4])); lines(x, rep(0, length(x)), col="blue")

if(PDF) dev.off()