# make a 2D heat map for each smval to show two-way interactions

library(fields)

PNG <- T

p = 4
n_pc = 4  # number of principal components to model
des <- "FF"
n_ff <- 11

load(paste0("/projects/precipit/hydro/hydro-SA-",des,if(des=="FF"){n_ff},
            if(des=="unif"){paste0("-",m)},"-nPC",n_pc,"_all64.rda"))

if(nv != n_ff) stop("this dim mismatch is not good for image plots")

# get ranges so they are consistent across all plots
zr <- c(0,0)
for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    for(sm in 1:n_sm){
      two_way_ints <-two_way_ints1 <-two_way_ints2 <- matrix(NA, nv, nv)
      for (ii in 1:nv) {
        for (jj in 1:nv) {
          two_way_ints[ii,jj] <- mean((bssEmu_meo[,which(facDes[,i] == prange[ii] &
                                                           facDes[,j] == prange[jj])])[sm,]) - avg_mean[sm]
          two_way_ints1[ii,jj] <- mean((bssEmu[,which(facDes[,i] == prange[ii] &
                                                        facDes[,j] == prange[jj])])[sm,]) - avg_mean[sm]
          two_way_ints2[ii,jj] <- mean((etaEmu[,which(facDes[,i] == prange[ii] &
                                                        facDes[,j] == prange[jj])])[sm,]) - avg_mean[sm]
          
        }
      }
      zr <- range(c(zr,two_way_ints, two_way_ints1, two_way_ints2))
    }
  }
}

## BSS-ANOVA (MEs only)
if(PNG) png("png/two_way_ints_BSS1.png", width=2000,height=2000,res=300)
par(mfrow=c(6,11), pty="s",oma=c(3,3,3,5),mar=c(0,0,0,0))
for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    # par(mfrow=c(2,3),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
    for(sm in 1:n_sm){
      two_way_ints <- matrix(NA, nv, nv)
      for (ii in 1:nv) {
        for (jj in 1:nv) {
          two_way_ints[ii,jj] <- mean((bssEmu_meo[,which(facDes[,i] == prange[ii] &
                                                           facDes[,j] == prange[jj])])[sm,]) - avg_mean[sm]
        }
      }
      image(two_way_ints, col=heat.colors(100),
            zlim = zr,ylab="", yaxt="n",xlab="", xaxt="n", asp=1)
      contour(two_way_ints, levels = c(-0.5,-0.25,0,.25,.5), add=T)
      if(sm == 1) mtext(paste(i,j,sep=","), side=2)
    }
  }
}
par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(3,15,3,0.5))
image.plot(matrix(rep(zr,2),2,2), legend.only = T, col=heat.colors(100))
par(mfrow=c(1,1),oma=c(0,0,0,0))
mtext("Stellar mass value", side=1, line=-2)
dev.off()

### BSS-ANOVA: MEs and 2WIs
# get ranges so they are consistent across all plots
if(PNG) png("png/two_way_ints_BSS2.png", width=2000,height=2000,res=300)
par(mfrow=c(6,11), pty="s",oma=c(3,3,3,5),mar=c(0,0,0,0))
for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    # par(mfrow=c(2,3),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
    for(sm in 1:n_sm){
      two_way_ints <- matrix(NA, nv, nv)
      for (ii in 1:nv) {
        for (jj in 1:nv) {
          two_way_ints[ii,jj] <- mean((bssEmu[,which(facDes[,i] == prange[ii] &
                                                       facDes[,j] == prange[jj])])[sm,]) - avg_mean[sm]
        }
      }
      image(two_way_ints, col=heat.colors(100),
            zlim = zr,ylab="", yaxt="n",xlab="", xaxt="n", asp=1)
      contour(two_way_ints, levels = c(-0.5,-0.25,0,.25,.5), add=T)
      if(sm == 1) mtext(paste(i,j,sep=","), side=2)
    }
  }
}
par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(3,15,3,0.5))
image.plot(matrix(rep(zr,2),2,2), legend.only = T, col=heat.colors(100))
par(mfrow=c(1,1),oma=c(0,0,0,0))
mtext("Stellar mass value", side=1, line=-2)
dev.off()

# GP-PC
# get ranges so they are consistent across all plots
if(PNG) png("png/two_way_ints_GPPC.png", width=2000,height=2000,res=300)
par(mfrow=c(6,11), pty="s",oma=c(3,3,3,5),mar=c(0,0,0,0))
for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    # par(mfrow=c(2,3),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
    for(sm in 1:n_sm){
      two_way_ints <- matrix(NA, nv, nv)
      for (ii in 1:nv) {
        for (jj in 1:nv) {
          two_way_ints[ii,jj] <- mean((etaEmu[,which(facDes[,i] == prange[ii] &
                                                       facDes[,j] == prange[jj])])[sm,]) - avg_mean[sm]
        }
      }
      image(two_way_ints, col=heat.colors(100),
            zlim = zr,ylab="", yaxt="n",xlab="", xaxt="n", asp=1)
      contour(two_way_ints, levels = c(-0.5,-0.25,0,.25,.5), add=T)
      if(sm == 1) mtext(paste(i,j,sep=","), side=2)
    }
  }
}
par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(3,15,3,0.5))
image.plot(matrix(rep(zr,2),2,2), legend.only = T, col=heat.colors(100))
par(mfrow=c(1,1),oma=c(0,0,0,0))
mtext("Stellar mass value", side=1, line=-2)
dev.off()

