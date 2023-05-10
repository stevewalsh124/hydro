library(GPfit)
source("../bssanova.R")

# main effects plots from hydro runs
nruns <- 32
p <- 4
nv <- 9
m <- 1
PDF <- T

if(PDF) pdf("pdf/LOOCV_hydro_complete.pdf")

rm_avg <- T

des1 = read.table("hydro.design")
colnames(des1) <- c("kappa", "EGW", "NPERH_AGN", "SeedMass")
# m = nrow(des1)
maxvals = c(8,4,4,4e6)
minvals = c(2,.25,1.5,4e5)

kappas <- des1[,1]
EGWs <- des1[,2]
NPERH_AGNs <- des1[,3]
SeedMasss <- des1[,4]

hydro_runs <- list()
for (i in 1:nruns) {
  n_c <- nchar(format(SeedMasss[i], scientific = T))
  seed_sc <- paste0(substr(format(SeedMasss[i], scientific = T),1, n_c - 3), 
                    substr(format(SeedMasss[i], scientific = T), n_c, n_c))
  
  # match values to the 128MPC folder names
  
  if(nchar(kappas[i])==4) kappas[i] <- paste0(kappas[i], "0")
  if(nchar(EGWs[i])==4) EGWs[i] <- paste0(EGWs[i], "0")
  if(nchar(EGWs[i])==6) EGWs[i] <- round(as.numeric(EGWs[i]), 3)
  if(nchar(NPERH_AGNs[i])==3) NPERH_AGNs[i] <- paste0(NPERH_AGNs[i], "00")
  if(nchar(NPERH_AGNs[i])==4) NPERH_AGNs[i] <- paste0(NPERH_AGNs[i], "0")
  if(nchar(seed_sc)==6) seed_sc <- paste0(substr(seed_sc, 1, 4),"0",substr(seed_sc, 5, 6))
  
  run_dir <- paste0("KAPPA_", kappas[i],
                    "_EGW_", EGWs[i],
                    "_NPERH_AGN_", NPERH_AGNs[i],
                    "_SEED_",seed_sc)
  
  hydro_runs[[i]] <- read.table(paste0("ProfileData/SCIDAC_RUNS/128MPC_RUNS/", 
                                       run_dir, 
                                       "/analysis_pipeline/profiles/GalStellarMassFunction_624.txt"))
}

sub_rg <-  which(log10(hydro_runs[[1]][,1]) >= 10 &  log10(hydro_runs[[1]][,1]) <= 11.5)
if(!exists("smvals")) smvals <- log10(hydro_runs[[1]][,1])[sub_rg]

n_sm <- length(smvals)
modRuns <- matrix(NA, n_sm, nruns)
for (i in 1:nruns) modRuns[,i] <- log10(hydro_runs[[i]][sub_rg,2])

plot(smvals, modRuns[,1], type="l", 
     xlab = "log10(Stellar Mass)", ylab = "log10(GSMF_Apperture)", main = "Subset x in [10,11.5]", ylim = range(modRuns))
for (i in 2:nruns) { lines(smvals, modRuns[,i], col=i)}

# plot(smvals, 10^(modRuns[,1]), type="l", 
#      xlab = "log10(Stellar Mass)", ylab = "GSMF_Apperture", main = "Subset x in [10,11.5]; y not log10", ylim = 10^(range(modRuns)))
# for (i in 2:nruns) { lines(smvals, 10^(modRuns[,i]), col=i)}

des001 = matrix(NA,nrow=nrow(des1),ncol=p)
for(k in 1:p) des001[,k] = (des1[,k]-minvals[k])/(maxvals[k]-minvals[k])

# LOOCV: skip the first, train on the rest
LOOCV_preds <- LOOCV_preds_b <- LOOCV_preds_bm <- matrix(NA, nrow = n_sm, ncol = nruns)

for (wts in 1:nruns) {
  print(wts) #which run to skip
  eta <- modRuns[,-wts]
  des_train <- des001[-wts,]
  eta_pred <- modRuns[,wts]
  des_pred <- t(des001[wts,])
  
  # # write it to a file; each column a simulation output
  # write(t(eta),file='csv/eta-100x111.txt',ncol=111)
  # # matplot(smvals,matrix(scan('eta-100x111.txt'),nrow=100,byrow=T),type='l')
  
  # for fun, let's use SVD and see the dimensionality in eta
  mean0mat = matrix(apply(eta,1,mean),nrow=n_sm,ncol=ncol(eta))
  eta0 = eta - mean0mat
  a = svd(eta0)
  # plot(a$d);abline(h=0)   # looks like 3 PCs should work
  # plot the basis elements
  # matplot(smvals,a$u%*%sqrt(diag(a$d)),type='l')
  
  # look at coefficients for each basis function
  # coef1 = a$v[,1]
  # hist(coef1)
  # scale the coefficients so they have variance = 1
  coef = a$v*sqrt(m)
  # accordingly, scale the bases so when multiplied by
  # the coef's, we'll get the spectra back
  bases = a$u%*%diag(a$d)/sqrt(m)
  
  spectraFull = bases%*%t(coef)
  # par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1.5,1))
  # matplot(smvals,bases%*%t(coef),type='l',ylab='P(k)')
  # mtext('all bases',side=3,lin=.2)
  # # one basis
  # matplot(smvals,bases[,1:1]%*%t(coef[,1:1]),type='l',ylab='P(k)')
  # mtext('1 basis',side=3,lin=.2)
  # # two basis
  # matplot(smvals,bases[,1:2]%*%t(coef[,1:2]),type='l',ylab='P(k)')
  # mtext('2 bases',side=3,lin=.2)
  # # 3 basis
  # matplot(smvals,bases[,1:3]%*%t(coef[,1:3]),type='l',ylab='P(k)')
  # mtext('3 bases',side=3,lin=.2)
  # matplot(smvals,bases%*%t(coef),type='l',ylab='P(k)')
  # mtext('all bases',side=3,lin=.2)
  # # 4 basis
  # matplot(smvals,bases[,1:4]%*%t(coef[,1:4]),type='l',ylab='P(k)')
  # mtext('4 bases',side=3,lin=.2)
  # # 5 basis
  # matplot(smvals,bases[,1:5]%*%t(coef[,1:5]),type='l',ylab='P(k)')
  # mtext('5 bases',side=3,lin=.2)
  # # 6 basis
  # matplot(smvals,bases[,1:6]%*%t(coef[,1:6]),type='l',ylab='P(k)')
  # mtext('6 bases',side=3,lin=.2)
  # # # 7 basis
  # # matplot(smvals,bases[,1:7]%*%t(coef[,1:7]),type='l',ylab='P(k)')
  # # mtext('7 bases',side=3,lin=.2)
  # 
  # # note, the spectra can be recovered by adding back in the mean
  # # one basis
  # matplot(smvals,bases[,1:1]%*%t(coef[,1:1])+mean0mat,type='l',ylab='P(k)')
  # mtext('1 basis',side=3,lin=.2)
  # matplot(smvals,bases[,1:2]%*%t(coef[,1:2])+mean0mat,type='l',ylab='P(k)')
  # mtext('2 bases',side=3,lin=.2)
  # matplot(smvals,bases[,1:3]%*%t(coef[,1:3])+mean0mat,type='l',ylab='P(k)')
  # mtext('3 bases',side=3,lin=.2)
  # matplot(smvals,bases[,1:4]%*%t(coef[,1:4])+mean0mat,type='l',ylab='P(k)')
  # mtext('4 bases',side=3,lin=.2)
  
  # try fitting GPs to the coefficients

  a1 = GP_fit(des_train,coef[,1])
  a2 = GP_fit(des_train,coef[,2])
  a3 = GP_fit(des_train,coef[,3])
  a4 = GP_fit(des_train,coef[,4])
  a5 = GP_fit(des_train,coef[,5])
  a6 = GP_fit(des_train,coef[,6])
  a7 = GP_fit(des_train,coef[,7])
  a8 = GP_fit(des_train,coef[,8])
  a9 = GP_fit(des_train,coef[,9])
  a10 = GP_fit(des_train,coef[,10])
  a11 = GP_fit(des_train,coef[,11])
  
  # make a 2^8 design
  # defined outside of the if/else statement
  # estimate main effect for parameter 1
  # repeat the udesign matrix nv times
  
  # make the predictions for the 2^8 design padded with nv variations for each of the p params
  a1_pred = predict(a1,des_pred)
  a2_pred = predict(a2,des_pred)
  a3_pred = predict(a3,des_pred)
  a4_pred = predict(a4,des_pred)
  a5_pred = predict(a5,des_pred)
  a6_pred = predict(a6,des_pred)
  a7_pred = predict(a7,des_pred)
  a8_pred = predict(a8,des_pred)
  a9_pred = predict(a9,des_pred)
  a10_pred = predict(a10,des_pred)
  a11_pred = predict(a11,des_pred)
  
  mean0pred = matrix(mean0mat[,1],nrow=n_sm,ncol=1)
  eta1pred = outer(bases[,1],a1_pred$Y_hat)
  eta2pred = outer(bases[,2],a2_pred$Y_hat)
  eta3pred = outer(bases[,3],a3_pred$Y_hat)
  eta4pred = outer(bases[,4],a4_pred$Y_hat)
  eta5pred = outer(bases[,5],a5_pred$Y_hat)
  eta6pred = outer(bases[,6],a6_pred$Y_hat)
  eta7pred = outer(bases[,7],a7_pred$Y_hat)
  eta8pred = outer(bases[,8],a8_pred$Y_hat)
  eta9pred = outer(bases[,9],a9_pred$Y_hat)
  eta10pred = outer(bases[,10],a10_pred$Y_hat)
  eta11pred = outer(bases[,11],a11_pred$Y_hat)
  LOOCV_preds[,wts] = eta1pred+eta2pred+eta3pred+eta4pred+
    eta5pred+eta6pred+eta7pred+eta8pred+eta9pred+eta10pred+mean0pred
  
  
  if(nrow(des_pred)==1) des_pred <- rbind(des_pred, des_pred)
  # try fitting bssanova to the coefficients
  b1 = bssanova(des_train,coef[,1])
  b2 = bssanova(des_train,coef[,2])
  b3 = bssanova(des_train,coef[,3])
  b4 = bssanova(des_train,coef[,4])
  b5 = bssanova(des_train,coef[,5])
  b6 = bssanova(des_train,coef[,6])
  b7 = bssanova(des_train,coef[,7])
  b8 = bssanova(des_train,coef[,8])
  b9 = bssanova(des_train,coef[,9])
  b10 = bssanova(des_train,coef[,10])
  b11 = bssanova(des_train,coef[,11])
  
  
  b1_meo = bssanova(des_train,coef[,1], int.order = 1)
  b2_meo = bssanova(des_train,coef[,2], int.order = 1)
  b3_meo = bssanova(des_train,coef[,3], int.order = 1)
  b4_meo = bssanova(des_train,coef[,4], int.order = 1)
  b5_meo = bssanova(des_train,coef[,5], int.order = 1)
  b6_meo = bssanova(des_train,coef[,6], int.order = 1)
  b7_meo = bssanova(des_train,coef[,7], int.order = 1)
  b8_meo = bssanova(des_train,coef[,8], int.order = 1)
  b9_meo = bssanova(des_train,coef[,9], int.order = 1)
  b10_meo = bssanova(des_train,coef[,10], int.order = 1)
  b11_meo = bssanova(des_train,coef[,11], int.order = 1)
  
  # now for the bssanova fits
  # des_pred <- rbind(des_pred, des_pred) #quick hack when there was only 1 pred location
  b1_pred = predict.bssanova(des_pred, b1)
  b2_pred = predict.bssanova(des_pred, b2)
  b3_pred = predict.bssanova(des_pred, b3)
  b4_pred = predict.bssanova(des_pred, b4)
  b5_pred = predict.bssanova(des_pred, b5)
  b6_pred = predict.bssanova(des_pred, b6)
  b7_pred = predict.bssanova(des_pred, b7)
  b8_pred = predict.bssanova(des_pred, b8)
  b9_pred = predict.bssanova(des_pred, b9)
  b10_pred = predict.bssanova(des_pred, b10)
  b11_pred = predict.bssanova(des_pred, b11)
  
  # now for the bssanova fits (main effects only)
  b1_pred_meo = predict.bssanova(des_pred, b1_meo)
  b2_pred_meo = predict.bssanova(des_pred, b2_meo)
  b3_pred_meo = predict.bssanova(des_pred, b3_meo)
  b4_pred_meo = predict.bssanova(des_pred, b4_meo)
  b5_pred_meo = predict.bssanova(des_pred, b5_meo)
  b6_pred_meo = predict.bssanova(des_pred, b6_meo)
  b7_pred_meo = predict.bssanova(des_pred, b7_meo)
  b8_pred_meo = predict.bssanova(des_pred, b8_meo)
  b9_pred_meo = predict.bssanova(des_pred, b9_meo)
  b10_pred_meo = predict.bssanova(des_pred, b10_meo)
  b11_pred_meo = predict.bssanova(des_pred, b11_meo)
  
  # now bss
  bss1pred = outer(bases[,1],b1_pred$yhat)
  bss2pred = outer(bases[,2],b2_pred$yhat)
  bss3pred = outer(bases[,3],b3_pred$yhat)
  bss4pred = outer(bases[,4],b4_pred$yhat)
  bss5pred = outer(bases[,5],b5_pred$yhat)
  bss6pred = outer(bases[,6],b6_pred$yhat)
  bss7pred = outer(bases[,7],b7_pred$yhat)
  bss8pred = outer(bases[,8],b8_pred$yhat)
  bss9pred = outer(bases[,9],b9_pred$yhat)
  bss10pred = outer(bases[,10],b10_pred$yhat)
  bss11pred = outer(bases[,11],b11_pred$yhat)
  
  bss1pred_meo = outer(bases[,1],b1_pred_meo$yhat)
  bss2pred_meo = outer(bases[,2],b2_pred_meo$yhat)
  bss3pred_meo = outer(bases[,3],b3_pred_meo$yhat)
  bss4pred_meo = outer(bases[,4],b4_pred_meo$yhat)
  bss5pred_meo = outer(bases[,5],b5_pred_meo$yhat)
  bss6pred_meo = outer(bases[,6],b6_pred_meo$yhat)
  bss7pred_meo = outer(bases[,7],b7_pred_meo$yhat)
  bss8pred_meo = outer(bases[,8],b8_pred_meo$yhat)
  bss9pred_meo = outer(bases[,9],b9_pred_meo$yhat)
  bss10pred_meo = outer(bases[,10],b10_pred_meo$yhat)
  bss11pred_meo = outer(bases[,11],b11_pred_meo$yhat)
  

  LOOCV_preds_b[,wts] <- (bss1pred+bss2pred+bss3pred+bss4pred+bss5pred+
                            bss6pred+bss7pred+bss8pred+bss9pred+bss10pred+bss11pred)[,1] + mean0pred
  LOOCV_preds_bm[,wts] <- (bss1pred_meo+bss2pred_meo+bss3pred_meo+bss4pred_meo+bss5pred_meo+
                             bss6pred_meo+bss7pred_meo+bss8pred_meo+bss9pred_meo+bss10pred_meo+bss11pred_meo)[,1] + mean0pred
  
}

# pdf("LOOCV_compare.pdf")
par(mfrow=c(1,1))
ylim_all <- range(c(modRuns, LOOCV_preds, LOOCV_preds_b, LOOCV_preds_bm))

matplot(smvals, modRuns[,order(modRuns[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "GP-PC (solid) vs simulated (dashed) log10(GMSF)",
        xlab = "log10(smass)", ylab="log10(GSMF)")

matplot(smvals, modRuns[,order(modRuns[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "GP-PC (solid) vs simulated (dashed) log10(GMSF)",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, LOOCV_preds[,order(modRuns[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32), add = T)

matplot(smvals, modRuns[,order(modRuns[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "BSS (MEs & 2WIs) (solid) vs simulated (dashed) log10(GMSF)",
        xlab = "log10(smass)", ylab="log10(GSMF)")

matplot(smvals, modRuns[,order(modRuns[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "BSS (MEs & 2WIs) (solid) vs simulated (dashed) log10(GMSF)",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, LOOCV_preds_b[,order(modRuns[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32), add = T)

matplot(smvals, modRuns[,order(modRuns[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "BSS (MEs only) (solid) vs simulated (dashed) log10(GMSF)",
        xlab = "log10(smass)", ylab="log10(GSMF)")

matplot(smvals, modRuns[,order(modRuns[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "BSS (MEs only) (solid) vs simulated (dashed) log10(GMSF)",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, LOOCV_preds_bm[,order(modRuns[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32), add = T)

matplot(smvals, LOOCV_preds[,order(modRuns[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "GP-PC",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, LOOCV_preds_b[,order(modRuns[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32),
        ylim = ylim_all, main = "BSS (MEs and 2WIs)",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, LOOCV_preds_bm[,order(modRuns[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32),
        ylim = ylim_all, main = "BSS (MEs only)",
        xlab = "log10(smass)", ylab="log10(GSMF)")


if(PDF) dev.off()

save(modRuns, file = "rda/modRuns_hydro.rda")
save(LOOCV_preds, file = "rda/LOOCV_preds.rda")
save(LOOCV_preds_b, file = "rda/LOOCV_preds_b.rda")
save(LOOCV_preds_bm, file = "rda/LOOCV_preds_bm.rda")
save.image("/projects/precipit/hydro/LOOCV_hydro.rda")

