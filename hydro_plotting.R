# Plots for hydro sensitivity analysis

PNG <- F
# need 8GB to load the 500 unif RData
# get unif results from bssanovaHydro.R first
load("rda/smvals_lo.rda")
load("rda/modRuns_lo_ac.rda")

# plot the 32 hydro runs
if(PNG) png("png/32runs.png",height = 1000, width=1500,res = 300)
par(mar=c(4.5,4.5,1,1)+0.1)
matplot(smvals, modRuns[,order(modRuns[1,1:32])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32),
        ylim = range(modRuns), xlab = "log10(stellar mass)", ylab="log10(GSMF)")
if(PNG) dev.off()

if(PNG) png("png/64runs.png",height = 1000, width=1500,res = 300)
par(mar=c(4.5,4.5,1,1)+0.1)
matplot(smvals, modRuns[,order(modRuns[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(64),
        ylim = range(modRuns), xlab = "log10(stellar mass)", ylab="log10(GSMF)")
if(PNG) dev.off()

load("/projects/precipit/hydro/hydro-SA-lo-trainac-predFF7-nPC4.rda")
# look at estimated main effects; w/ average removed only (w/ and w/o below)
if(PNG) png("png/main_effects_bss_meo.png", width = 2400, height = 800, res = 300)
# Plot bss-anova (MEs only) main effects
yr = range(mainEffs_bss_meo - avg_mean)
par(mfrow=c(1,4),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
for(k in 1:p){
  matplot(smvals,mainEffs_bss_meo[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)', if(rm_avg){" [avg removed]"}),side=2,line=2,outer=T)
mtext("bss-anova (MEs only) main effect estimates", side=3,outer = T)
if(PNG) dev.off()

# look at estimated main effects w/ and w/o average removed
if(PNG) png("png/main_effects_bss_meo_yesandnoavg.png", width = 2400, height = 1200, res = 300)
pnames = c(paste0("param",1:4))#c('omega_m','omega_b','sigma_8','h','n_s','w_0','w_a','omega_nu','z')
par(mfrow=c(2,4),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
grcolors = paste0("grey",round(seq(90,25,length=11)),sep='')

yr = range(mainEffs_bss_meo - avg_mean)
for(k in 1:p){
  matplot(smvals,mainEffs_bss_meo[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  # if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)', if(rm_avg){" [average removed top row]"}),side=2,line=2.5,outer=T)
mtext("bss-anova (MEs only) main effect estimates", side=3,outer = T)

yr = range(mainEffs_bss_meo)
for(k in 1:p){
  matplot(smvals,mainEffs_bss_meo[,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
if(PNG) dev.off()

# sum of squares on sqrt scale for the FF design (for all 3 methods)
if(PNG) png("png/sqrt_ss_FF.png", width = 2400, height = 800, res = 300)
{
  par(mfrow=c(1,3),oma=c(0,0,0,0),mar=c(4,4,1.8,1))
  if(ERROR) plot(1:ntest,sqrt(ss),xlab='effect',ylab='sqrt SS',log='')
  if(ERROR) mtext('true model',side=3,line=.2,cex=.9)
  plot(1:2^p,sqrt(ssEmu),xlab='effect',ylab='emu sqrt SS',col='blue',log='',ylim=ss_lim)
  mtext('GP-PC emulator',side=3,line=.2,cex=.9)
  abline(v=c(1.5, 5.5, 11.5, 15.5), col="gray", lty=2)
  if(ERROR) plot(1:2^p,sqrt(sseEmu),xlab='effect',ylab='emu sqrt SS',col='blue',log='',ylim=ss_lim)
  if(ERROR) mtext('emu trained on y+e',side=3,line=.2,cex=.9)
  if(ERROR) plot(1:256,sqrt(ssEmuError),xlab='effect',ylab='emu+error sqrt SS',col='red',log='',ylim=ss_lim)
  if(ERROR) mtext('emu + e',side=3,line=.2,cex=.9)
  plot(1:2^p,sqrt(ssBssEmu),xlab='effect',ylab='bss emu sqrt SS',col='red',log='',ylim=ss_lim)
  mtext('bss emu',side=3,line=.2,cex=.9)
  abline(v=c(1.5, 5.5, 11.5, 15.5), col="gray", lty=2)
  plot(1:2^p,sqrt(ssBssEmu_meo),xlab='effect',ylab='bss emu (MEs only) sqrt SS',col='red',log='',ylim=ss_lim)
  mtext('bss emu (MEs only)',side=3,line=.2,cex=.9)
  abline(v=c(1.5, 5.5, 11.5, 15.5), col="gray", lty=2)
}
if(PNG) dev.off()


load("rda/modRuns_hydro.rda")
load("rda/LOOCV_preds.rda")
load("rda/LOOCV_preds_b.rda")
load("rda/LOOCV_preds_bm.rda")
load("rda/modRuns_pred.rda")
load("rda/new32_preds.rda")
load("rda/new32_preds_b.rda")
load("rda/new32_preds_bm.rda")

# Older version of new32_compare without axis labels
# par(mfrow=c(1,3),oma=c(1,1,1,1),mar=c(0,0,0,0))
# ylim_all <- range(c(modRuns, new32_preds, new32_preds_b, new32_preds_bm))
# 
# matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
# lwd=2, col=viridisLite::viridis(32), axes=F,
#         ylim = ylim_all, xlab = "log10(smass)", ylab="log10(GSMF)")
# box();text(10.6,max(ylim_all)-.06,"GP",adj=c(0,0), cex=1.25)
# matplot(smvals, new32_preds[,order(modRuns_pred[1,])], type="l", lty=1, 
# lwd=2, col=viridisLite::viridis(32), axes=F, add = T)
# 
# matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
# lwd=2, col=viridisLite::viridis(32), axes=F, 
#         ylim = ylim_all, xlab = "log10(smass)", ylab="log10(GSMF)")
# box();text(10.6,max(ylim_all)-.06,"BSS (2nd order)",adj=c(0,0), cex=1.25)
# matplot(smvals, new32_preds_b[,order(modRuns_pred[1,])], type="l", lty=1, 
# lwd=2, col=viridisLite::viridis(32), axes=F, add = T)
# 
# matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
# lwd=2, col=viridisLite::viridis(32), axes=F, 
#         ylim = ylim_all,xlab = "log10(smass)", ylab="log10(GSMF)")
# box();text(10.6,max(ylim_all)-.06,"BSS (1st order)",adj=c(0,0), cex=1.25)
# matplot(smvals, new32_preds_bm[,order(modRuns_pred[1,])], type="l", lty=1, 
# lwd=2, col=viridisLite::viridis(32), axes=F, add = T)



# Compare the estimates for the new 32 runs for each of the 3 methods
if(PNG) png("png/new32_compare.png",width = 1200, height = 400, res = 150)
yr = range(mainEffs_bss_meo - avg_mean)
par(mfrow=c(1,3),oma=c(4,4,1,1),mar=c(0,0,0,0))

ylim_all <- range(c(modRuns, new32_preds, new32_preds_b, new32_preds_bm))

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
        lwd=2, col=viridisLite::viridis(32), axes=F,
        ylim = ylim_all, xlab = "log10(smass)", ylab="log10(GSMF)")
box();text(10.6,max(ylim_all)-.06,"GP",adj=c(0,0), cex=1.25)
matplot(smvals, new32_preds[,order(modRuns_pred[1,])], type="l", lty=1, 
        lwd=2, col=viridisLite::viridis(32), axes=F, add = T)
axis(1); axis(2)

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
        lwd=2, col=viridisLite::viridis(32), axes=F, 
        ylim = ylim_all, xlab = "log10(smass)", ylab="log10(GSMF)")
box();text(10.6,max(ylim_all)-.06,"BSS (2nd order)",adj=c(0,0), cex=1.25)
matplot(smvals, new32_preds_b[,order(modRuns_pred[1,])], type="l", lty=1, 
        lwd=2, col=viridisLite::viridis(32), axes=F, add = T)
axis(1)
matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
        lwd=2, col=viridisLite::viridis(32), axes=F, 
        ylim = ylim_all,xlab = "log10(smass)", ylab="log10(GSMF)")
box();text(10.6,max(ylim_all)-.06,"BSS (1st order)",adj=c(0,0), cex=1.25)
matplot(smvals, new32_preds_bm[,order(modRuns_pred[1,])], type="l", lty=1, 
        lwd=2, col=viridisLite::viridis(32), axes=F, add = T)
axis(1)

mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)'),side=2,line=2.5,outer=T)
# mtext("bss-anova (MEs only) main effect estimates", side=3,outer = T)
if(PNG) dev.off()



if(PNG) png("png/LOOCV_predictions.png",width = 1200, height = 400, res = 150)
yr = range(mainEffs_bss_meo - avg_mean)
par(mfrow=c(1,3),oma=c(4,4,1,1),mar=c(0,0,0,0))

ylim_all <- range(c(modRuns, LOOCV_preds, LOOCV_preds_b, LOOCV_preds_bm))

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
        lwd=2, col=viridisLite::viridis(32), axes=F,
        ylim = ylim_all, xlab = "log10(smass)", ylab="log10(GSMF)")
box();text(10.6,max(ylim_all)-.06,"GP",adj=c(0,0), cex=1.25)
matplot(smvals, LOOCV_preds[,order(modRuns_pred[1,])], type="l", lty=1, 
        lwd=2, col=viridisLite::viridis(32), axes=F, add = T)
axis(1); axis(2)

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
        lwd=2, col=viridisLite::viridis(32), axes=F, 
        ylim = ylim_all, xlab = "log10(smass)", ylab="log10(GSMF)")
box();text(10.6,max(ylim_all)-.06,"BSS (2nd order)",adj=c(0,0), cex=1.25)
matplot(smvals, LOOCV_preds_b[,order(modRuns_pred[1,])], type="l", lty=1, 
        lwd=2, col=viridisLite::viridis(32), axes=F, add = T)
axis(1)
matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, 
        lwd=2, col=viridisLite::viridis(32), axes=F, 
        ylim = ylim_all,xlab = "log10(smass)", ylab="log10(GSMF)")
box();text(10.6,max(ylim_all)-.06,"BSS (1st order)",adj=c(0,0), cex=1.25)
matplot(smvals, LOOCV_preds_bm[,order(modRuns_pred[1,])], type="l", lty=1, 
        lwd=2, col=viridisLite::viridis(32), axes=F, add = T)
axis(1)

mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)'),side=2,line=2.5,outer=T)
# mtext("bss-anova (MEs only) main effect estimates", side=3,outer = T)
if(PNG) dev.off()



# RMSE as a function of stellar mass (sm)
if(PNG) png("png/RMSE_fn_of_sm_batch.png", width = 2400, height = 1200, res = 300)
par(mfrow=c(1,1), mar=c(4, 4, 2, 2) + 0.1)
rmseXsm <- rmseXsm_b <- rmseXsm_bm <- c()
for (i in 1:nrow(modRuns_pred)) {
  rmseXsm[i] <- sqrt(mean((modRuns_pred[i,] - new32_preds[i,])^2))
  rmseXsm_b[i] <- sqrt(mean((modRuns_pred[i,] - new32_preds_b[i,])^2))
  rmseXsm_bm[i] <- sqrt(mean((modRuns_pred[i,] - new32_preds_bm[i,])^2))
}

plot(smvals,rmseXsm, type = "l", ylim=range(rmseXsm, rmseXsm_b, rmseXsm_bm), 
     xlab="log10(stellar mass)", ylab="RMSE")
lines(smvals,rmseXsm_b, col="red", lty=2)
lines(smvals,rmseXsm_bm, col="blue", lty=3)
legend(legend = c("GP-PC", "BSS (2nd order)", "BSS (1st order)"), "topleft", 
       lty = 1:3, col=c("black","red","blue"))
if(PNG) dev.off()

if(PNG) png("png/RMSE_fn_of_sm_loocv.png", width = 2400, height = 1200, res = 300)
par(mfrow=c(1,1), mar=c(4, 4, 2, 2) + 0.1)
rmseXsm <- rmseXsm_b <- rmseXsm_bm <- c()
for (i in 1:nrow(modRuns)) {
  rmseXsm[i] <- sqrt(mean((modRuns[i,] - LOOCV_preds[i,])^2))
  rmseXsm_b[i] <- sqrt(mean((modRuns[i,] - LOOCV_preds_b[i,])^2))
  rmseXsm_bm[i] <- sqrt(mean((modRuns[i,] - LOOCV_preds_bm[i,])^2))
}

plot(smvals,rmseXsm, type = "l", ylim=range(rmseXsm, rmseXsm_b, rmseXsm_bm), 
     xlab="log10(stellar mass)", ylab="RMSE")
lines(smvals,rmseXsm_b, col="red", lty=2)
lines(smvals,rmseXsm_bm, col="blue", lty=3)
legend(legend = c("GP-PC", "BSS (2nd order)", "BSS (1st order)"), "topleft", 
       lty = 1:3, col=c("black","red","blue"))
if(PNG) dev.off()



# illustrate decomposition: mean, PCs, and how GP or similar can model weights
# run new32_hydro.R, or at least get a b1_pred from there
if(PNG) png("png/mean_PCs_oneW.png", width = 2400, height = 800, res=300)
par(mfrow=c(1,3), mar=c(4,4,1.8,1))
plot(smvals, mean0mat[,1], type="l", xlab= "log10(stellar mass)", ylab = "log10(GSMF)")
matplot(bases, type="l", xlab="log10(stellar mass)", ylab = "log10(GSMF)")

x1_pred <- facDes[,1]
plot(x1_pred[order(x1_pred)], bps[[1]]$curves[1,800,][order(x1_pred)], type="l", xlab = expression(kappa), ylab = expression(W[1]))
# points(x1_pred, b1_pred_meo$curves[1,cu,], col="blue")
if(PNG) dev.off()
