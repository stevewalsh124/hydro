# Set up the basic function to produce simulated (fake)
# spectra on which we'll fit an emulator.

tic <- proc.time()[3]

library(GPfit)
source("bssanova.R")

p = 4     # number of parameters to vary for sensitivity analysis
n_pc = 4  # number of principal components to model
if(n_pc < 1 | n_pc > 11) stop("number of PCs has to be between 1 and 11")

ERROR=F  # try putting error into the model output to see what happens
PDF=T    # do we make pdfs?
WRITE=T  # write csvs and rda files?

# which testing design to use?
# FF: full factorial design padded with nv levels for each parameter; ntest = 2^(p-1) * (p*nv)
# unif: uniform sampling design; ntest = m * (p*nv)
# condl: conditional, each of nv settings for each param is 0.5 elsewhere (ntest = 1 * (p*nv))
des <- "FF"
if(!(des %in% c("FF","unif","condl"))) stop('des must be in c("FF","unif","condl")')

if(des =="unif") m <- 100 else m <- 1
if(des == "FF"){
  n_ff <- 3
  m <- n_ff^(p-1)
} 

wts <- 1

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

if(PDF) pdf(paste0("pdf/LOOCV_MEs.pdf"))

nv = 9
prange = seq(0,1,length=nv)

load("rda/smvals_full.rda")
smvals = smvals_full[which(smvals_full >= 10 & smvals_full <= 11.5)]

# make a function that produces spectrum-like output, depending
# on the 4 model parameters
# simModel <- function(pmat){
#   # makes a simulated output (over [-3,1]) for each row of pmat
#   # pmat holds the 8-d parameter inputs.  The output is 100 x nrow(pmat)
#   # smvals = seq(8.5, 12.5, length.out = 29)
#   # xlocs = seq(8, 13, length.out = 8)         # knot locs
#   # v0 = c(3.5,5,8.5,8.6,10.0,7.2,10.9,15.3) # knot values
#   # dv1 = c(.35, .3, .21, .15, .15, .45, .2, .5)
#   # dv2 = c(0,0,0,.1,.2,.2,.15,.1)
#   # K = matrix(0,nrow=29,ncol=8)
#   # for(i in 1:8) K[,i] = dnorm(smvals,mean=xlocs[i],sd=.7)
#   eta0 = outer(avg_mean_hydro, rep(1, m))
#   eta1 = outer(mainEffs_real[,1,1],pmat[,1])+outer(mainEffs_real[,nv,1],1-pmat[,1]) 
#   eta2 = outer(mainEffs_real[,1,2],pmat[,2])+outer(mainEffs_real[,nv,2],1-pmat[,2]) 
#   eta3 = outer(mainEffs_real[,1,3],pmat[,3])+outer(mainEffs_real[,nv,3],1-pmat[,3]) 
#   eta4 = outer(mainEffs_real[,1,4],pmat[,4])+outer(mainEffs_real[,nv,4],1-pmat[,4]) 
#   # return a 100xnrow(pmat) matrix of output
#   return(eta0+eta1+eta2+eta3+eta4)
# }

# make some emulator runs - each column is a model run
# modRuns = simModel(udesign)
load("rda/modRuns_hydro.rda")
modRuns_first32 <- modRuns
load("rda/modRuns_hydro_new32.rda")
modRuns <- cbind(modRuns_first32, modRuns_hydro)[,-wts]
nruns <- ncol(modRuns)
# plot them
matplot(smvals, modRuns, type="l", 
        xlab = "Stellar Mass", ylab = "GSMF_Apperture", main = "both x and y log10")

matplot(smvals, 10^modRuns, type="l", 
        xlab = "Stellar Mass", ylab = "GSMF_Apperture", main = "only x log10")
# read in the space filling LHC design from the Mira-Titan paper
des0 = read.table("hydro.design")
des00 = read.table("hydro.design_c")
des1 <- rbind(des0, des00)[-wts,]
colnames(des1) <- c("kappa", "EGW", "NPERH_AGN", "SeedMass")
# m = nrow(des0)
maxvals = c(8,4,4,4e6)
minvals = c(2,.25,1.5,4e5)

des01 = matrix(NA,nrow=nrow(des1),ncol=p)
for(k in 1:p) des01[,k] = (des1[,k]-minvals[k])/(maxvals[k]-minvals[k])

# plot to check
# if(PDF) pdf('pdf/designnrunsx8.pdf',width=7,height=7)
pairs(des1,pch=c(rep('.',32),rep("*",32)))
# if(PDF) dev.off()
pairs(des01,pch=c(rep('.',32),rep("*",32)))

# write out the 01 design easily readable by matlab
# if(WRITE) write(t(des01),file='hydro_design_out_64.txt',ncol=p,sep=' ')

# create the simulation output for this design
eta = modRuns
matplot(smvals,eta,type='l')

# write it to a file; each column a simulation output
# if(WRITE) write(t(eta),file='eta-hydro_64.txt',ncol=nruns)
# matplot(smvals,matrix(scan('eta-100xnruns.txt'),nrow=100,byrow=T),type='l')

n_sm <- length(smvals)

# for fun, let's use SVD and see the dimensionality in eta
mean0mat = matrix(apply(eta,1,mean),nrow=n_sm,ncol=nruns)
eta0 = eta - mean0mat
a = svd(eta0)
plot(a$d)   # looks like 3 pc's should work
# plot the basis elements
matplot(smvals,a$u%*%sqrt(diag(a$d)),type='l')

# look at coefficients for each basis function
coef1 = a$v[,1]
hist(coef1)
# scale the coefficients so they have variance = 1
coef = a$v*sqrt(nruns)
# accordingly, scale the bases so when multiplied by
# the coef's, we'll get the spectra back
bases = a$u%*%diag(a$d)/sqrt(nruns)

spectraFull = bases%*%t(coef)
par(mfrow=c(2,3),oma=c(0,0,0,0),mar=c(4,4,1.5,1))
matplot(smvals,bases%*%t(coef),type='l',ylab='P(k)')
mtext(paste0('all (',ncol(bases),') bases'),side=3,lin=.2)
# one basis
matplot(smvals,bases[,1:1]%*%t(coef[,1:1]),type='l',ylab='P(k)')
mtext('1 basis',side=3,lin=.2)
# 3 basis
matplot(smvals,bases[,1:3]%*%t(coef[,1:3]),type='l',ylab='P(k)')
mtext('3 bases',side=3,lin=.2)
# 5 basis
matplot(smvals,bases[,1:5]%*%t(coef[,1:5]),type='l',ylab='P(k)')
mtext('5 bases',side=3,lin=.2)
# 6 basis
matplot(smvals,bases[,1:6]%*%t(coef[,1:6]),type='l',ylab='P(k)')
mtext('6 bases',side=3,lin=.2)
# 9 basis
matplot(smvals,bases[,1:9]%*%t(coef[,1:9]),type='l',ylab='P(k)')
mtext('9 bases',side=3,lin=.2)

# note, the spectra can be recovered by adding back in the mean
# one basis
matplot(smvals,bases%*%t(coef)+mean0mat,type='l',ylab='P(k)')
mtext(paste0('all (',ncol(bases),') bases'),side=3,lin=.2)
# one basis
matplot(smvals,bases[,1:1]%*%t(coef[,1:1])+mean0mat,type='l',ylab='P(k)')
mtext('1 basis',side=3,lin=.2)
# 3 basis
matplot(smvals,bases[,1:3]%*%t(coef[,1:3])+mean0mat,type='l',ylab='P(k)')
mtext('3 bases',side=3,lin=.2)
# 5 basis
matplot(smvals,bases[,1:5]%*%t(coef[,1:5])+mean0mat,type='l',ylab='P(k)')
mtext('5 bases',side=3,lin=.2)
# 6 basis
matplot(smvals,bases[,1:6]%*%t(coef[,1:6])+mean0mat,type='l',ylab='P(k)')
mtext('6 bases',side=3,lin=.2)
# 9 basis
matplot(smvals,bases[,1:9]%*%t(coef[,1:9])+mean0mat,type='l',ylab='P(k)')
mtext('9 bases',side=3,lin=.2)

# try fitting gp's to the coefficients
as <- list()
for (i in 1:n_pc) {
  print(paste("GP",i))
  as[[i]] = GP_fit(des01,coef[,i])
}

# try fitting bssanova to the coefficients
bs <- list()
for (i in 1:n_pc) {
  print(paste("BSS-ANOVA",i))
  bs[[i]] = bssanova(des01,coef[,i])
}

# MEO: Main effects only
bms <- list()
for (i in 1:n_pc) {
  print(paste("BSS-ANOVA-MEO",i))
  bms[[i]] = bssanova(des01,coef[,i], int.order = 1)
}

if(ERROR){
  aes <- list()
  for (i in 1:n_pc) {
    print(paste("GP w/ error",i))
    aes[[i]] = GP_fit(des01,coef[,i]+rnorm(m,mean=0,sd=.001))
  }
}

# Create testing design (facDes)
if(des == "FF"){
  # make a n_ff^(p-1) design
  x_ff <- seq(0,1,len=n_ff)
  facDes=do.call(expand.grid,rep(list(x_ff), p-1))
  # repeat it nv times
  repfacDes = matrix(rep(t(facDes),nv),ncol=p-1,byrow=TRUE)
  
  # set up a vector that varies parameter value with each copy of FF design
  prange = seq(0,1,length=nv)
  repPrange = rep(prange, each =n_ff^(p-1))
  
  # replace the value for each parameter with repPrange
  list_facDes <- list()
  list_facDes[[1]] <- unname(cbind(repPrange, repfacDes))
  for (kk in 2:(p-1)) {
    list_facDes[[kk]] <- unname(cbind(repfacDes[,1:(kk-1)], repPrange, repfacDes[,kk:(p-1)]))
  }
  list_facDes[[p]] <- unname(cbind(repfacDes, repPrange))
  for(kk in 1:p) colnames(list_facDes[[kk]]) <- paste0("p",1:p)
  facDes_all <- do.call("rbind", list_facDes)
  facDes <- facDes_all
  
} else if(des == "condl") {
  
  xx <- matrix(0.5, nv, p)
  
  grid <- seq(0, 1, length.out=nv)
  
  list_facDes <- list()
  for (i in 1:p) {
    xx_temp <- xx
    xx_temp[,i] <- grid
    names(xx_temp) <- paste0("Var", 1:p)
    list_facDes[[i]] <- xx_temp
  }
  
  facDes_all <- as.matrix(dplyr::bind_rows(list_facDes))
  facDes_all <- unname(facDes_all)
  facDes <- facDes_all
  
} else if(des =="unif") {
  
  # generate parameter settings that are uniform over [0,1]^8
  udesign = matrix(runif(p*m),ncol=p) #lhs::randomLHS(m,p) #
  # estimate main effect for parameter 1
  # repeat the udesign matrix nv times
  repUdes = matrix(rep(t(udesign),nv),ncol=p,byrow=TRUE)
  # set up a vector that varies parameter value with each copy of udesign
  repPrange = matrix(prange,ncol=m,nrow=nv)
  repPrange = as.vector(t(repPrange))
  
  # do this for all 8 parameters
  list_facDes <- list()
  for(k in 1:p){
    des1 = repUdes
    des1[,k] = repPrange
    names(des1) <- paste0("Var", 1:p)
    list_facDes[[k]] <- des1
  }
  
  facDes_all <- as.matrix(dplyr::bind_rows(list_facDes))
  facDes_all <- unname(facDes_all)
  facDes <- facDes_all
} else {stop("Value for 'des' argument not recognized")}

ntest <- nrow(facDes)

# make the predictions for the given facDes design
aps <- list()
for (i in 1:n_pc) {
  print(paste("GP pred",i))
  aps[[i]] = predict(as[[i]],facDes)
}

# now for the bssanova fits
bps <- list()
for (i in 1:n_pc) {
  print(paste("BSS-ANOVA pred",i))
  bps[[i]] = predict.bssanova(facDes, bs[[i]])
}

# now for the bssanova fits (main effects only)
bmps <- list()
for (i in 1:n_pc) {
  print(paste("BSS-ANOVA-MEO pred",i))
  bmps[[i]] = predict.bssanova(facDes, bms[[i]])
}

if(ERROR){
  aeps <- list()
  for (i in 1:n_pc) {
    print(paste("GP w/ error pred",i))
    aeps[[i]] = predict(aes[[i]],facDes)
  }
}

mean0pred = matrix(mean0mat[,1],nrow=n_sm,ncol=nrow(facDes))

# create GP-PC predictions
eta_preds <- list()
for (i in 1:n_pc) eta_preds[[i]] = outer(bases[,i],aps[[i]]$Y_hat)
etaEmu <- mean0pred + Reduce('+', eta_preds)

# now bss
bss_preds <- list()
for (i in 1:n_pc) bss_preds[[i]] = outer(bases[,i],bps[[i]]$yhat)
bssEmu = mean0pred + Reduce('+', bss_preds)

# now bss (main effects only)
bssm_preds <- list()
for (i in 1:n_pc) bssm_preds[[i]] = outer(bases[,i],bmps[[i]]$yhat)
bssEmu_meo = mean0pred + Reduce('+', bssm_preds)

if(ERROR){
  etae_preds <- list()
  for (i in 1:n_pc) etae_preds[[i]] = outer(bases[,i],aeps[[i]]$Y_hat)
  etaeEmu <- mean0pred + Reduce('+', etae_preds)
}

# show emulator output compared to training data
if(ntest < 1000){
  par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1.5,1))
  matplot(smvals,eta,type='l',ylab='GSMF_A')
  mtext('32-run training set',side=3,line=.1,cex=.9)
  
  matplot(smvals,etaEmu,type='l',ylab='emulator')
  mtext(paste0(ntest,'-run design, GP-PC'),side=3,line=.1,cex=.9)
  
  matplot(smvals,bssEmu,type='l',ylab='bss emulator')
  mtext(paste0(ntest,'-run design, bss-anova (MEs and 2WIs)'),side=3,line=.1,cex=.9)
  
  matplot(smvals,bssEmu_meo,type='l',ylab='bss (main effects only) emulator')
  mtext(paste0(ntest,'-run design, bss-anova (MEs only)'),side=3,line=.1,cex=.9)
}

# etatrue = simModel(facDes)
# matplot(smvals,etatrue,type='l',ylab='GSMF_A'); mtext('2^4-run design',side=3,line=.1,cex=.9)

# compare 2^k model fits
# grab 10 levels of smass and make the corresponding 
igrab = round(seq(1,n_sm,length=n_sm))
# etakFac = etatrue[igrab,]
# # plot check ... again
# matplot(1:10,etakFac,type='l')
mainEffs_bss = mainEffs_gppc = mainEffs_bss_meo = array(NA,c(n_sm,nv,p))
for(k in 1:p){
  # collect bss-anova (MEs and 2WIs) predictions
  sims1 = bssEmu[,(k-1)*(ntest/p)+1:(ntest/p)]
  # reshape into a 3d array by kval,rep, paramval
  sims1a = array(sims1,c(n_sm,m,nv))
  # compute main effect
  mainEffs_bss[,,k] = apply(sims1a,c(1,3),mean)
  
  # collect bss-anova (main effects only) predictions
  sims1 = bssEmu_meo[,(k-1)*(ntest/p)+1:(ntest/p)]
  # reshape into a 3d array by kval,rep, paramval
  sims1a = array(sims1,c(n_sm,m,nv))
  # compute main effect
  mainEffs_bss_meo[,,k] = apply(sims1a,c(1,3),mean)
  
  # collect GP-PC predictions
  sims1 = etaEmu[,(k-1)*(ntest/p)+1:(ntest/p)]
  # reshape into a 3d array by kval,rep, paramval
  sims1a = array(sims1,c(n_sm,m,nv))
  # compute main effect
  mainEffs_gppc[,,k] = apply(sims1a,c(1,3),mean)
}

rm_avg <- T

if(rm_avg) avg_mean <- rowMeans(modRuns) else avg_mean <- rep(0, nrow(modRuns))

# make a plot - you could do this with ggplot if you'd rather
pnames = c(paste0("param",1:4))#c('omega_m','omega_b','sigma_8','h','n_s','w_0','w_a','omega_nu','z')

par(mfrow=c(2,4),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
grcolors = paste0("grey",round(seq(90,25,length=11)),sep='')

# Plot bss-anova (MEs and 2WIs) main effects estimates
yr = range(mainEffs_bss - avg_mean)
for(k in 1:p){
  matplot(smvals,mainEffs_bss[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)', if(rm_avg){" [average removed top row]"}),side=2,line=2.5,outer=T)
mtext("bss-anova (ME and 2WIs) main effect estimates", side=3,outer = T)

yr = range(mainEffs_bss)
for(k in 1:p){
  matplot(smvals,mainEffs_bss[,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}

# Plot bss-anova (MEs only) main effects
yr = range(mainEffs_bss_meo - avg_mean)
for(k in 1:p){
  matplot(smvals,mainEffs_bss_meo[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
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

# Plot GP-PC main effects
yr = range(mainEffs_gppc - avg_mean)
for(k in 1:p){
  matplot(smvals,mainEffs_gppc[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)', if(rm_avg){" [average removed top row]"}),side=2,line=2.5,outer=T)
mtext("GP-PC main effect estimates", side=3,outer = T)

yr = range(mainEffs_gppc)
for(k in 1:p){
  matplot(smvals,mainEffs_gppc[,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}


# nvm1 <- nv - 1
# par(mfrow=c(1,1))
# for (it in 1:p) {
#   matplot(smvals, mainEffs_real[,,it], type = "l", lty=1, col = grcolors)
#   for (pp in (0:nvm1)/nvm1) {
#     tempy <- pp*mainEffs_real[,1,it]+(1-pp)*mainEffs_real[,nv,it]
#     lines(smvals, tempy, col ="red", lty=2, lwd=2)
#   }
# }



# make the corresponding design
facDesK <- matrix(NA, nrow = ntest*n_sm, ncol = p+1)
for (i in 1:ntest) {
  for (j in 1:n_sm) {
    facDesK[(i-1)*(n_sm) + j,] <- c(smvals[j], facDes[i,])
  }
}
# facDesK = expand.grid(smvals,0:1,0:1,0:1,0:1)
# appropriately reshape the model output
# yfacDesK = as.vector(etakFac)
# ytrue = as.vector(etatrue[igrab,])
yEmu = as.vector(etaEmu[igrab,])
if(ERROR) yeEmu = as.vector(etaeEmu[igrab,])
if(ERROR) yEmuError = yEmu + rnorm(length(yEmu),mean=0,sd=.001)
ybssEmu = as.vector(bssEmu[igrab,])
ybssEmu_meo = as.vector(bssEmu_meo[igrab,])
# convert columns of facDesK to factors
for(k in 1:ncol(facDesK)) facDesK[,k] = as.factor(facDesK[,k])
colnames(facDesK) = c("sm","p1","p2","p3","p4")

# make a dataframe for lm()
if(ERROR){
  facSAdf = data.frame(y1=yEmu,y2=yeEmu,y3=yEmuError,y4=ybssEmu,y5=ybssEmu_meo,facDesK) #y0=ytrue,
} else {
  facSAdf = data.frame(y1=yEmu,y4=ybssEmu,y5=ybssEmu_meo,facDesK) #y0=ytrue,
}

# run the linear model to estimate effects; everything interacts
# with k.
if(ERROR) a = lm(y0 ~ sm+sm:(p1+p2+p3+p4)^3,data=facSAdf)
aEmu = lm(y1 ~ sm+sm:(p1+p2+p3+p4)^3,data=facSAdf)
if(ERROR) aeEmu = lm(y2 ~ sm+sm:(p1+p2+p3+p4)^3,data=facSAdf)
if(ERROR) aEmuError = lm(y3 ~ sm+sm:(p1+p2+p3+p4)^3,data=facSAdf)
aBssEmu = lm(y4 ~ sm+sm:(p1+p2+p3+p4)^3,data=facSAdf)
aBssEmu_meo = lm(y5 ~ sm+sm:(p1+p2+p3+p4)^3,data=facSAdf)

# get the Sum of Squares for each term using the anova() function
if(ERROR) b = anova(a)
bEmu = anova(aEmu)
if(ERROR) beEmu = anova(aeEmu)
if(ERROR) bEmuError = anova(aEmuError)
bBssEmu = anova(aBssEmu)
bBssEmu_meo = anova(aBssEmu_meo)

# make a plot of the sums of squares
ssnames = row.names(bEmu)
if(ERROR) ss = b$`Sum Sq`
ss = ssEmu = bEmu$`Sum Sq`
if(ERROR) sseEmu = beEmu$`Sum Sq`
if(ERROR) ssEmuError = bEmuError$`Sum Sq`
ssBssEmu = bBssEmu$`Sum Sq`
ssBssEmu_meo = bBssEmu_meo$`Sum Sq`

ss_lim <- range(sqrt(c(ss,ssBssEmu, ssBssEmu_meo)))
if(des=="condl"){
  par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1.8,1))
  if(ERROR) plot(1:ntest,sqrt(ss),xlab='effect',ylab='sqrt SS',log='y')
  if(ERROR) mtext('true model',side=3,line=.2,cex=.9)
  plot(1:(p+2),sqrt(ssEmu),xlab='effect',ylab='emu sqrt SS',col='blue',log='y',ylim=ss_lim)
  mtext('GP-PC emulator',side=3,line=.2,cex=.9)
  if(ERROR) plot(1:(p+2),sqrt(sseEmu),xlab='effect',ylab='emu sqrt SS',col='blue',log='y',ylim=ss_lim)
  if(ERROR) mtext('emu trained on y+e',side=3,line=.2,cex=.9)
  if(ERROR) plot(1:(p+2),sqrt(ssEmuError),xlab='effect',ylab='emu+error sqrt SS',col='red',log='y',ylim=ss_lim)
  if(ERROR) mtext('emu + e',side=3,line=.2,cex=.9)
  plot(1:(p+2),sqrt(ssBssEmu),xlab='effect',ylab='bss emu sqrt SS',col='red',log='y',ylim=ss_lim)
  mtext('bss emu',side=3,line=.2,cex=.9)
  plot(1:(p+2),sqrt(ssBssEmu_meo),xlab='effect',ylab='bss emu (MEs only) sqrt SS',col='red',log='y',ylim=ss_lim)
  mtext('bss emu (MEs only)',side=3,line=.2,cex=.9)
} else {
  par(mfrow=c(1,3),oma=c(0,0,0,0),mar=c(4,4,1.8,1))
  if(ERROR) plot(1:2^p,sqrt(ss),xlab='effect',ylab='sqrt SS',log='y')
  if(ERROR) mtext('true model',side=3,line=.2,cex=.9)
  plot(1:2^p,sqrt(ssEmu),xlab='effect',ylab='emu sqrt SS',col='blue',log='y',ylim=ss_lim)
  mtext('GP-PC emulator',side=3,line=.2,cex=.9)
  abline(v=c(1.5, 5.5, 11.5, 15.5), col="gray", lty=2)
  if(ERROR) plot(1:2^p,sqrt(sseEmu),xlab='effect',ylab='emu sqrt SS',col='blue',log='y',ylim=ss_lim)
  if(ERROR) mtext('emu trained on y+e',side=3,line=.2,cex=.9)
  if(ERROR) plot(1:2^p,sqrt(ssEmuError),xlab='effect',ylab='emu+error sqrt SS',col='red',log='y',ylim=ss_lim)
  if(ERROR) mtext('emu + e',side=3,line=.2,cex=.9)
  plot(1:2^p,sqrt(ssBssEmu),xlab='effect',ylab='bss emu sqrt SS',col='red',log='y',ylim=ss_lim)
  mtext('bss emu',side=3,line=.2,cex=.9)
  abline(v=c(1.5, 5.5, 11.5, 15.5), col="gray", lty=2)
  plot(1:2^p,sqrt(ssBssEmu_meo),xlab='effect',ylab='bss emu (MEs only) sqrt SS',col='red',log='y',ylim=ss_lim)
  mtext('bss emu (MEs only)',side=3,line=.2,cex=.9)
  abline(v=c(1.5, 5.5, 11.5, 15.5), col="gray", lty=2)
}


# identify large effects on the plot
# plot(1:2^p,sqrt(ssEmu),xlab='effect',ylab='emu+error sqrt SS',col='red',log='y')
# identify(1:2^p,sqrt(ssEmu), labels = ssnames)
# ssnames[c(1,2,11,16)] # residuals is important

toc <- proc.time()[3]

toc - tic

if(PDF) dev.off()
if(WRITE){
  save(ssEmu, file = paste0("/projects/precipit/hydro/rda/ssEmu",wts))
  save(ssBssEmu, file = paste0("/projects/precipit/hydro/rda/ssBssEmu",wts))
  save(ssBssEmu_meo, file = paste0("/projects/precipit/hydro/rda/ssBssEmu_meo",wts))
  save(mainEffs_gppc, file = paste0("/projects/precipit/hydro/rda/mainEffs_gppc",wts))
  save(mainEffs_bss, file = paste0("/projects/precipit/hydro/rda/mainEffs_bss",wts))
  save(mainEffs_bss_meo, file = paste0("/projects/precipit/hydro/rda/mainEffs_bss_meo",wts))
}
# if(WRITE) save.image(paste0("/projects/precipit/hydro/hydro-SA-",des,if(des=="FF"){n_ff},
#                             if(des=="unif"){paste0("-",m)},"-nPC",n_pc,"_all64_wts",wts,".rda"))

