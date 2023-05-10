library(GPfit)
source("bssanova.R")

# main effects plots from hydro runs
m <- 32 # number of runs
p <- 4
nv <- 9

PDF <- T
if(PDF) pdf("pdf/new32_hydro_rev.pdf")

rm_avg <- T

des1 = read.table("hydro.design_c")
des1_pred = read.table("hydro.design")
colnames(des1) <- c("kappa", "EGW", "NPERH_AGN", "SeedMass")
# m = nrow(des1)
maxvals = c(8,4,4,4e6)
minvals = c(2,.25,1.5,4e5)

kappas <- des1[,1]
EGWs <- des1[,2]
NPERH_AGNs <- des1[,3]
SeedMasss <- des1[,4]

hydro_runs <- list()
for (i in 1:m) {
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
  
  hydro_runs[[i]] <- read.table(paste0("128MPC_RUNS_NEW/", 
                                       run_dir, 
                                       "/analysis_pipeline/extract/GalStellarMassFunction_624.txt"))
}

sub_rg <-  which(log10(hydro_runs[[1]][,1]) >= 10 &  log10(hydro_runs[[1]][,1]) <= 11.5)
if(!exists("smvals")) smvals <- log10(hydro_runs[[1]][,1])[sub_rg]

n_sm <- length(smvals)
modRuns <- matrix(NA, n_sm, m)
for (i in 1:m) modRuns[,i] <- log10(hydro_runs[[i]][sub_rg,2])

plot(smvals, modRuns[,1], type="l", 
     xlab = "log10(Stellar Mass)", ylab = "log10(GSMF_Apperture)", main = "Subset x in [10,11.5]", ylim = range(modRuns))
for (i in 2:m) { lines(smvals, modRuns[,i], col=i)}





colnames(des1_pred) <- c("kappa", "EGW", "NPERH_AGN", "SeedMass")
maxvals = c(8,4,4,4e6)
minvals = c(2,.25,1.5,4e5)

kappas <- des1_pred[,1]
EGWs <- des1_pred[,2]
NPERH_AGNs <- des1_pred[,3]
SeedMasss <- des1_pred[,4]

hydro_runs <- list()
for (i in 1:m) {
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
modRuns_pred <- matrix(NA, n_sm, m)
for (i in 1:m) modRuns_pred[,i] <- log10(hydro_runs[[i]][sub_rg,2])

plot(smvals, modRuns_pred[,1], type="l", 
     xlab = "log10(Stellar Mass)", ylab = "log10(GSMF_Apperture)", main = "Subset x in [10,11.5]", ylim = range(modRuns_pred))
for (i in 2:m) { lines(smvals, modRuns_pred[,i], col=i)}



# plot(smvals, 10^(modRuns[,1]), type="l", 
#      xlab = "log10(Stellar Mass)", ylab = "GSMF_Apperture", main = "Subset x in [10,11.5]; y not log10", ylim = 10^(range(modRuns)))
# for (i in 2:m) { lines(smvals, 10^(modRuns[,i]), col=i)}

des001 = des001_pred = matrix(NA,nrow=nrow(des1),ncol=p)
for(k in 1:p) des001[,k] = (des1[,k]-minvals[k])/(maxvals[k]-minvals[k])
for(k in 1:p) des001_pred[,k] = (des1_pred[,k]-minvals[k])/(maxvals[k]-minvals[k])

# use the first 32 to train, the second 32 to test
new32_preds <- new32_preds_b <- new32_preds_bm <- matrix(NA, nrow = n_sm, ncol = m)

eta <- modRuns
des_train <- des001

des_pred <- des001_pred

# # write it to a file; each column a simulation output
# write(t(eta),file='csv/eta-100x111.txt',ncol=111)
# # matplot(smvals,matrix(scan('eta-100x111.txt'),nrow=100,byrow=T),type='l')

# for fun, let's use SVD and see the dimensionality in eta
mean0mat = matrix(apply(eta,1,mean),nrow=n_sm,ncol=ncol(eta))
eta0 = eta - mean0mat
a = svd(eta0)
plot(a$d);abline(h=0)   # looks like 3 PCs should work
# plot the basis elements
matplot(smvals,a$u%*%sqrt(diag(a$d)),type='l')

# look at coefficients for each basis function
# coef1 = a$v[,1]
# hist(coef1)
# scale the coefficients so they have variance = 1
coef = a$v*sqrt(m)
# accordingly, scale the bases so when multiplied by
# the coef's, we'll get the spectra back
bases = a$u%*%diag(a$d)/sqrt(m)

spectraFull = bases%*%t(coef)

# try fitting GPs to the coefficients

a1 = GP_fit(des_train,coef[,1])
a2 = GP_fit(des_train,coef[,2])
a3 = GP_fit(des_train,coef[,3])
a4 = GP_fit(des_train,coef[,4])
# a5 = GP_fit(des_train,coef[,5])
# a6 = GP_fit(des_train,coef[,6])
# a7 = GP_fit(des_train,coef[,7])
# a8 = GP_fit(des_train,coef[,8])
# a9 = GP_fit(des_train,coef[,9])
# a10 = GP_fit(des_train,coef[,10])


# make a 2^8 design
# defined outside of the if/else statement
# estimate main effect for parameter 1
# repeat the udesign matrix nv times

# make the predictions for the 2^8 design padded with nv variations for each of the p params
a1_pred = predict(a1,des_pred)
a2_pred = predict(a2,des_pred)
a3_pred = predict(a3,des_pred)
a4_pred = predict(a4,des_pred)
# a5_pred = predict(a5,des_pred)
# a6_pred = predict(a6,des_pred)
# a7_pred = predict(a7,des_pred)
# a8_pred = predict(a8,des_pred)
# a9_pred = predict(a9,des_pred)
# a10_pred = predict(a10,des_pred)

mean0pred = matrix(mean0mat[,1],nrow=n_sm,ncol=1)
eta1pred = outer(bases[,1],a1_pred$Y_hat)
eta2pred = outer(bases[,2],a2_pred$Y_hat)
eta3pred = outer(bases[,3],a3_pred$Y_hat)
eta4pred = outer(bases[,4],a4_pred$Y_hat)
# eta5pred = outer(bases[,5],a5_pred$Y_hat)
# eta6pred = outer(bases[,6],a6_pred$Y_hat)
# eta7pred = outer(bases[,7],a7_pred$Y_hat)
# eta8pred = outer(bases[,8],a8_pred$Y_hat)
# eta9pred = outer(bases[,9],a9_pred$Y_hat)
# eta10pred = outer(bases[,10],a10_pred$Y_hat)
new32_preds = eta1pred+eta2pred+eta3pred+eta4pred+mean0mat #+eta6pred+eta7pred+eta8pred+eta9pred+eta10pred


# try fitting bssanova to the coefficients
b1 = bssanova(des_train,coef[,1])
b2 = bssanova(des_train,coef[,2])
b3 = bssanova(des_train,coef[,3])
b4 = bssanova(des_train,coef[,4])

b1_meo = bssanova(des_train,coef[,1], int.order = 1)
b2_meo = bssanova(des_train,coef[,2], int.order = 1)
b3_meo = bssanova(des_train,coef[,3], int.order = 1)
b4_meo = bssanova(des_train,coef[,4], int.order = 1)

# now for the bssanova fits
# des_pred <- rbind(des_pred, des_pred) #quick hack when there was only 1 pred location
b1_pred = predict.bssanova(des_pred, b1)
b2_pred = predict.bssanova(des_pred, b2)
b3_pred = predict.bssanova(des_pred, b3)
b4_pred = predict.bssanova(des_pred, b4)
# now for the bssanova fits (main effects only)
b1_pred_meo = predict.bssanova(des_pred, b1_meo)
b2_pred_meo = predict.bssanova(des_pred, b2_meo)
b3_pred_meo = predict.bssanova(des_pred, b3_meo)
b4_pred_meo = predict.bssanova(des_pred, b4_meo)

# now bss
bss1pred = outer(bases[,1],b1_pred$yhat)
bss2pred = outer(bases[,2],b2_pred$yhat)
bss3pred = outer(bases[,3],b3_pred$yhat)
bss4pred = outer(bases[,4],b4_pred$yhat)
bssEmu = bss1pred+bss2pred+bss3pred+bss4pred+mean0mat

bss1pred_meo = outer(bases[,1],b1_pred_meo$yhat)
bss2pred_meo = outer(bases[,2],b2_pred_meo$yhat)
bss3pred_meo = outer(bases[,3],b3_pred_meo$yhat)
bss4pred_meo = outer(bases[,4],b4_pred_meo$yhat)

new32_preds_b <- bssEmu
new32_preds_bm <- bss1pred_meo+bss2pred_meo+bss3pred_meo+bss4pred_meo+mean0mat

# pdf("pdf/new32_compare.pdf")
par(mfrow=c(1,1))
ylim_all <- range(c(modRuns, new32_preds, new32_preds_b, new32_preds_bm))

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "GP-PC (solid) vs simulated (dashed)",
        xlab = "log10(smass)", ylab="log10(GSMF)")

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "GP-PC (solid) vs simulated (dashed)",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, new32_preds[,order(modRuns_pred[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32), add = T)

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "BSS (MEs & 2WIs) (solid) vs simulated (dashed)",
        xlab = "log10(smass)", ylab="log10(GSMF)")

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "BSS (MEs & 2WIs) (solid) vs simulated (dashed)",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, new32_preds_b[,order(modRuns_pred[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32), add = T)

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "BSS (MEs only) (solid) vs simulated (dashed)",
        xlab = "log10(smass)", ylab="log10(GSMF)")

matplot(smvals, modRuns_pred[,order(modRuns_pred[1,])], type="l", lty=2, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "BSS (MEs only) (solid) vs simulated (dashed)",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, new32_preds_bm[,order(modRuns_pred[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32), add = T)

matplot(smvals, new32_preds[,order(modRuns_pred[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32), 
        ylim = ylim_all, main = "GP-PC",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, new32_preds_b[,order(modRuns_pred[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32),
        ylim = ylim_all, main = "BSS (MEs and 2WIs)",
        xlab = "log10(smass)", ylab="log10(GSMF)")
matplot(smvals, new32_preds_bm[,order(modRuns_pred[1,])], type="l", lty=1, lwd=2, col=viridisLite::viridis(32),
        ylim = ylim_all, main = "BSS (MEs only)",
        xlab = "log10(smass)", ylab="log10(GSMF)")


if(PDF) dev.off()

save(modRuns_pred, file = "rda/modRuns_pred_rev.rda")
save(new32_preds, file = "rda/new32_preds_rev.rda")
save(new32_preds_b, file = "rda/new32_preds_b_rev.rda")
save(new32_preds_bm, file = "rda/new32_preds_bm_rev.rda")
