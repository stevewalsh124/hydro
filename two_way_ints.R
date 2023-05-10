# Set up the basic function to produce simulated (fake)
# spectra on which we'll fit an emulator.

p=4
n_pc = 4  # number of principal components to model
des <- "FF"

if(des =="unif") m <- 100 else m <- 1
if(des == "FF"){
  n_ff <- 7
  m <- n_ff^(p-1)
} 

load(paste0("rda/hydro-SA-",des,if(des=="FF"){n_ff},
            if(des=="unif"){paste0("-",m)},"-nPC",n_pc,"_all64.rda"))

if(nv != n_ff) stop("this dim mismatch is not good for image plots")

pdf(paste0("pdf/two_way_ints_",nv,".pdf"))

# # make images in 2x3 window, since 11 smvals
# for (i in 1:(p-1)) {
#   for (j in (i+1):p) {
#     print(c(i,j))
#     
#     unique(c(which(c(T,F,F,F,T)), which(c(T,F,F,T,F))))
#     
#     # make a 2D heat map for each smval to show interactions
#     par(mfrow=c(2,3),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
#     yr = range(mainEffs_bss - avg_mean)
#     for(sm in 1:n_sm){
#       two_way_ints <- matrix(NA, nv, nv)
#       for (ii in 1:nv) {
#         for (jj in 1:nv) {
#           print(c(ii,jj))
#           two_way_ints[ii,jj] <- mean((bssEmu[,which(facDes[,i] == prange[ii] &
#                                                       facDes[,j] == prange[jj])])[sm,]) - avg_mean[sm]
#         }
#       }
#       image(two_way_ints,axes=F); box()
#       # text(-3,yr[1]+0.01,pnames[sm],adj=c(0,0))
#       if(sm %in% c(1,4,7,10)) axis(2)
#       if(sm %in% c(4:6, 9:11)) axis(1)
#       if(sm %in% c(1,7)) {
#         mtext(paste("param",i),side=1,line=2.5,outer=T)
#         mtext(paste("param",j),side=2,line=2.5,outer=T)
#         mtext("bss-anova (ME and 2WIs) main effect estimates", side=3,outer = T)
#       }
#     }
#   }
# }

# get ranges so they are consistent across all plots
zr <- c(0,0)
for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    for(sm in 1:n_sm){
      two_way_ints <- matrix(NA, nv, nv)
      for (ii in 1:nv) {
        for (jj in 1:nv) {
          two_way_ints[ii,jj] <- mean((bssEmu[,which(facDes[,i] == prange[ii] &
                                                       facDes[,j] == prange[jj])])[sm,]) - avg_mean[sm]
          
        }
      }
      zr <- range(c(zr,two_way_ints))
    }
  }
}

for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    # make a 2D heat map for each smval to show interactions
    # par(mfrow=c(2,3),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
    yr = range(mainEffs_bss - avg_mean)
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

dev.off()
