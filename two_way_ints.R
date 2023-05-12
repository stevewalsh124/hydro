# make a 2D heat map for each smval to show two-way interactions

library(fields)

p = 4
n_pc = 4  # number of principal components to model
des <- "FF"
n_ff <- 11

load(paste0("/projects/precipit/hydro/hydro-SA-",des,if(des=="FF"){n_ff},
            if(des=="unif"){paste0("-",m)},"-nPC",n_pc,"_all64.rda"))

if(nv != n_ff) stop("this dim mismatch is not good for image plots")

pdf(paste0("pdf/two_way_ints_",nv,".pdf"))

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
      image.plot(two_way_ints, xlab = paste("param",i), ylab = paste("param",j),
                 main = paste("bss-anova (MEs only), smval", sm),
                 zlim = zr)
    }
  }
}


### BSS-ANOVA: MEs and 2WIs
# get ranges so they are consistent across all plots
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
      image.plot(two_way_ints, xlab = paste("param",i), ylab = paste("param",j),
                 main = paste("bss-anova (ME and 2WIs), smval", sm),
                 zlim = zr)
    }
  }
}

# GP-PC
# get ranges so they are consistent across all plots
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
      image.plot(two_way_ints, xlab = paste("param",i), ylab = paste("param",j),
                 main = paste("GP-PC, smval", sm),
                 zlim = zr)
    }
  }
}

dev.off()

