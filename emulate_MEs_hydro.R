# main effects plots from hydro runs
nruns <- 32
p <- 4
nv <- 9
m <- 1
PDF <- T

if(PDF) pdf("emulate_MEs_hydro.pdf")

rm_avg <- T

des0 = read.table("hydro.design")
colnames(des0) <- c("kappa", "EGW", "NPERH_AGN", "SeedMass")
# m = nrow(des0)
des1 = des0
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
  
  hydro_runs[[i]] <- read.table(paste0("128MPC_RUNS/", 
                                       run_dir, 
                                       "/analysis_pipeline/profiles/GalStellarMassFunction_624.txt"))
}

sub_rg <-  which(log10(hydro_runs[[1]][,1]) >= 10 &  log10(hydro_runs[[1]][,1]) <= 11.5)
if(!exists("smvals")) smvals <- log10(hydro_runs[[1]][,1])[sub_rg]

n_sm <- length(smvals)
modRuns <- matrix(NA, n_sm, nruns)
for (i in 1:nruns) modRuns[,i] <- log10(hydro_runs[[i]][sub_rg,2])
modRuns_hydro <- modRuns
save(modRuns_hydro, file = "modRuns_hydro.rda")

plot(smvals, modRuns[,1], type="l", 
     xlab = "log10(Stellar Mass)", ylab = "log10(GSMF_Apperture)", main = "Subset x in [10,11.5]", ylim = range(modRuns))
for (i in 2:nruns) { lines(smvals, modRuns[,i], col=i)}

# plot(smvals, 10^(modRuns[,1]), type="l", 
#      xlab = "log10(Stellar Mass)", ylab = "GSMF_Apperture", main = "Subset x in [10,11.5]; y not log10", ylim = 10^(range(modRuns)))
# for (i in 2:nruns) { lines(smvals, 10^(modRuns[,i]), col=i)}

des01 = matrix(NA,nrow=nrow(des1),ncol=p)
for(k in 1:p) des01[,k] = (des1[,k]-minvals[k])/(maxvals[k]-minvals[k])

# LOOCV: skip the first, train on the rest
# wts <- 1 #which run to skip
eta <- modRuns#[,-wts]
des01 <- des01#[-wts,]

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
coef1 = a$v[,1]
hist(coef1)
# scale the coefficients so they have variance = 1
coef = a$v*sqrt(m)
# accordingly, scale the bases so when multiplied by
# the coef's, we'll get the spectra back
bases = a$u%*%diag(a$d)/sqrt(m)

spectraFull = bases%*%t(coef)
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1.5,1))
matplot(smvals,bases%*%t(coef),type='l',ylab='P(k)')
mtext('all bases',side=3,lin=.2)
# one basis
matplot(smvals,bases[,1:1]%*%t(coef[,1:1]),type='l',ylab='P(k)')
mtext('1 basis',side=3,lin=.2)
# two basis
matplot(smvals,bases[,1:2]%*%t(coef[,1:2]),type='l',ylab='P(k)')
mtext('2 bases',side=3,lin=.2)
# 3 basis
matplot(smvals,bases[,1:3]%*%t(coef[,1:3]),type='l',ylab='P(k)')
mtext('3 bases',side=3,lin=.2)
matplot(smvals,bases%*%t(coef),type='l',ylab='P(k)')
mtext('all bases',side=3,lin=.2)
# 4 basis
matplot(smvals,bases[,1:4]%*%t(coef[,1:4]),type='l',ylab='P(k)')
mtext('4 bases',side=3,lin=.2)
# 5 basis
matplot(smvals,bases[,1:5]%*%t(coef[,1:5]),type='l',ylab='P(k)')
mtext('5 bases',side=3,lin=.2)
# 6 basis
matplot(smvals,bases[,1:6]%*%t(coef[,1:6]),type='l',ylab='P(k)')
mtext('6 bases',side=3,lin=.2)
# # 7 basis
# matplot(smvals,bases[,1:7]%*%t(coef[,1:7]),type='l',ylab='P(k)')
# mtext('7 bases',side=3,lin=.2)

# note, the spectra can be recovered by adding back in the mean
# one basis
matplot(smvals,bases[,1:1]%*%t(coef[,1:1])+mean0mat,type='l',ylab='P(k)')
mtext('1 basis',side=3,lin=.2)
matplot(smvals,bases[,1:2]%*%t(coef[,1:2])+mean0mat,type='l',ylab='P(k)')
mtext('2 bases',side=3,lin=.2)
matplot(smvals,bases[,1:3]%*%t(coef[,1:3])+mean0mat,type='l',ylab='P(k)')
mtext('3 bases',side=3,lin=.2)
matplot(smvals,bases[,1:4]%*%t(coef[,1:4])+mean0mat,type='l',ylab='P(k)')
mtext('4 bases',side=3,lin=.2)

# try fitting GPs to the coefficients
library(GPfit)

print("fit GP1")
a1 = GP_fit(des01,coef[,1])
print("fit GP2")
a2 = GP_fit(des01,coef[,2])
print("fit GP3")
a3 = GP_fit(des01,coef[,3])
print("fit GP4")
a4 = GP_fit(des01,coef[,4])
print("fit GP5")
a5 = GP_fit(des01,coef[,5])
# print("fit GP6")
# a6 = GP_fit(des01,coef[,6])
# print("fit GP7")
# a7 = GP_fit(des01,coef[,7])
# print("fit GP8")
# a8 = GP_fit(des01,coef[,8])
# print("fit GP9")
# a9 = GP_fit(des01,coef[,9])
# print("fit GP10")
# a10 = GP_fit(des01,coef[,10])


# make a 2^8 design
# defined outside of the if/else statement
# estimate main effect for parameter 1
# repeat the udesign matrix nv times
facDes = expand.grid(0:1,0:1,0:1,0:1)
repfacDes <- facDes[rep(seq_len(nrow(facDes)), each = nv), ]

# set up a vector that varies parameter value with each copy of udesign
prange = seq(0,1,length=nv)
repPrange = rep(prange, 2^(p-1))

# replace the value for each parameter with repPrange
list_facDes <- list()
for (i in 1:p) {
  des2 = repfacDes
  des2[,i] <- repPrange
  list_facDes[[i]] <- des2[!duplicated(des2),] # since the facdes had 0 and 1 for each param, there are duplicate columns
}

facDes_all <- as.matrix(dplyr::bind_rows(list_facDes))
facDes_all <- unname(facDes_all)

# make the predictions for the 2^8 design padded with nv variations for each of the p params
a1_pred = predict(a1,facDes_all)
a2_pred = predict(a2,facDes_all)
a3_pred = predict(a3,facDes_all)
a4_pred = predict(a4,facDes_all)
a5_pred = predict(a5,facDes_all)
# a6_pred = predict(a6,facDes_all)
# a7_pred = predict(a7,facDes_all)
# a8_pred = predict(a8,facDes_all)
# a9_pred = predict(a9,facDes_all)
# a10_pred = predict(a10,facDes_all)

mean0pred = matrix(mean0mat[,1],nrow=n_sm,ncol=2^(p-1)*p*nv)
eta1pred = outer(bases[,1],a1_pred$Y_hat)
eta2pred = outer(bases[,2],a2_pred$Y_hat)
eta3pred = outer(bases[,3],a3_pred$Y_hat)
eta4pred = outer(bases[,4],a4_pred$Y_hat)
eta5pred = outer(bases[,5],a5_pred$Y_hat)
# eta6pred = outer(bases[,6],a6_pred$Y_hat)
# eta7pred = outer(bases[,7],a7_pred$Y_hat)
# eta8pred = outer(bases[,8],a8_pred$Y_hat)
# eta9pred = outer(bases[,9],a9_pred$Y_hat)
# eta10pred = outer(bases[,10],a10_pred$Y_hat)
etafacDes = eta1pred+eta2pred+eta3pred+eta4pred+eta5pred+mean0pred #+eta6pred+eta7pred+eta8pred+eta9pred+eta10pred

# do this for all p parameters
mainEffs = array(NA,c(length(smvals),nv,p))
for(k in 1:p){
  # collect sumulations
  sims1 = etafacDes[,(k-1)*(m*nv*2^(p-1))+1:(m*nv*2^(p-1))]
  # reshape into a 3d array by kval,rep, paramval
  sims1a = array(sims1,c(n_sm,m,nv))
  # compute main effect
  mainEffs[,,k] = apply(sims1a,c(1,3),mean)
}

if(rm_avg) avg_mean <- rowMeans(modRuns) else avg_mean <- rep(0, nrow(modRuns))

# make a plot - you could do this with ggplot if you'd rather
pnames = c(paste0("param",1:4))#c('omega_m','omega_b','sigma_8','h','n_s','w_0','w_a','omega_nu','z')

par(mfrow=c(2,4),oma=c(4,4,1,1),mar=c(0,0,0,0))
yr = range(mainEffs - avg_mean)
grcolors = paste0("grey",round(seq(90,25,length=11)),sep='')
for(k in 1:p){
  matplot(smvals,mainEffs[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)', if(rm_avg){" [average removed top row]"}),side=2,line=2.5,outer=T)

mainEffs_real <- mainEffs
mainEffs_true <- array(dim = dim(mainEffs_real))

yr = range(mainEffs_real)
for(k in 1:p){
  matplot(smvals,mainEffs_real[,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
# mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
# mtext('log10(GSMF_A)',side=2,line=2.5,outer=T)


nvm1 <- nv - 1
par(mfrow=c(1,1))
for (it in 1:p) {
  matplot(smvals, mainEffs_real[,,it], type = "l", lty=1, col = grcolors)
  for (pp in (0:nvm1)/nvm1) {
    mainEffs_true[,nvm1*pp+1,it] <- pp*mainEffs_real[,1,it]+(1-pp)*mainEffs_real[,nv,it]
    lines(smvals, mainEffs_true[,nvm1*pp+1,it], col ="red", lty=2, lwd=2)
  }
}

if(PDF) dev.off()

avg_mean_hydro <- avg_mean
mainEffs_true_hydro <- mainEffs_true
save(avg_mean_hydro, file = "avg_mean_hydro.rda")
save(mainEffs_real, file = "mainEffs_real.rda")
save(mainEffs_true_hydro, file = "mainEffs_true_hydro.rda")
