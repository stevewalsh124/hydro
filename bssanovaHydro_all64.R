# For a given training dataset, obtain batch predictions and 
# main effects estimates, and ANOVA decomps of the variance
# For either held out data or a full factorial design, etc

tic <- proc.time()[3]

do_GPMSA <- T
do_BSS   <- F
# do_BSSM  <- T

if(do_GPMSA) library(GPfit)
if(do_BSS | do_BSSM) source("bssanova.R")

####################################
# Which data do will you train on? #
####################################

# Use the high res runs to train?
hi <- T

# Which design do you want to train with? a, c, or both?
atf <- T
ctf <- F

##############################
# Load the preprocessed data #
##############################

# load the training design (des01) 
load(paste0("rda/des_",if(atf){"a"},if(ctf){"c"},".rda")) 
des_train <- des01
p = ncol(des_train)  # number of input parameters

# load the stellar mass vals
load(paste0("rda/smvals_",ifelse(hi,"hi","lo"),".rda"))

# Load the pre-processed model runs (run preprocessing.R first)
load(paste0("rda/modRuns_",ifelse(hi,"hi","lo"),"_",if(atf){"a"},if(ctf){"c"},".rda"))
eta = modRuns
nruns <- ncol(eta)

#########################################
# Choose # of PCs, and if to do pdf/rda #
#########################################

n_pc = 4         # number of principal components to model
if(n_pc < 1 | n_pc > 11) stop("number of PCs has to be between 1 and 11")

ERROR=F  # try putting error into the model output to see what happens
PDF=F    # do we make pdfs?
WRITE=T  # write csvs and rda files?

################################
# Which testing design to use? #
################################

# FF: full factorial design padded with nv levels for each parameter; ntest = 2^(p-1) * (p*nv)
# unif: uniform sampling design; ntest = m * (p*nv)
# condl: conditional, each of nv settings for each param is 0.5 elsewhere (ntest = 1 * (p*nv))
des <- "lhsS"
if(!(des %in% c("FF","unif","lhs","lhsS","unifS","condl","a","c","ac","lhs"))) 
  stop('des must be in c("FF","unif","lhs","lhsS","unifS","condl","a","c","ac")')

if(des %in% c("unif","lhs","lhsS","unifS")) m <- 10^4 else m <- 1
if(des %in% c("lhsS","unifS")) seed <- 1

if(des == "FF"){
  n_ff <- 2
  m <- n_ff^(p-1)
} 

if(PDF) pdf(paste0("pdf/hydro-SA-",ifelse(hi,"hi-","lo-"),
                   "train",if(atf){"a"},if(ctf){"c"},"-",
                   "pred",des,if(des=="FF"){n_ff},
                   if(des%in%c("unif","lhs","lhsS","unifS")){m},
                   "-nPC",n_pc,".pdf"))

nv = 11
prange = seq(0,1,length=nv)

# plot them
matplot(smvals, eta, type="l", 
        xlab = "log10(Stellar Mass)", ylab = "log10(GSMF_Apperture)", 
        main = paste0(nruns, " training runs, ",
                      ifelse(hi,"hi ","lo "),
                      "res, designs ",if(atf){"a "},if(ctf){"c"}))

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
matplot(bases, type="l")

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
  as[[i]] = GP_fit(des_train,coef[,i])
}

if(do_BSS){
  # try fitting bssanova to the coefficients
  bs <- list()
  for (i in 1:n_pc) {
    print(paste("BSS-ANOVA",i))
    bs[[i]] = bssanova(des_train,coef[,i])
  }
  
  # MEO: Main effects only
  bms <- list()
  for (i in 1:n_pc) {
    print(paste("BSS-ANOVA-MEO",i))
    bms[[i]] = bssanova(des_train,coef[,i], int.order = 1)
  }
}

if(ERROR){
  aes <- list()
  for (i in 1:n_pc) {
    print(paste("GP w/ error",i))
    aes[[i]] = GP_fit(des_train,coef[,i]+rnorm(m,mean=0,sd=.001))
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
  
} else if(des %in% c("unif","lhs")) {
  # generate parameter settings that are uniform over [0,1]^8
  if(des=="unif") udesign = matrix(runif(p*m),ncol=p) 
  if(des=="lhs") udesign <- lhs::randomLHS(m,p) 
  
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
  
} else if(des %in% c("lhsS","unifS")) {
  if(des=="unifS") facDes = matrix(runif(p*m),ncol=p) 
  if(des=="lhsS") facDes <- lhs::randomLHS(m,p) 
  
} else if(des =="a") {
  load("rda/des_a.rda")
  facDes <- des01
} else if(des =="c") {
  load("rda/des_a.rda")
  facDes <- des01
} else if(des =="ac") {
  load("rda/des_a.rda")
  facDesa <- des01
  load("rda/des_c.rda")
  facDesc <- des01
  facDes <- rbind(facDesa, facDesc)
}

ntest <- nrow(facDes)

# make the predictions for the given facDes design
aps <- list()
for (i in 1:n_pc) {
  print(paste("GP pred",i))
  aps[[i]] = predict(as[[i]],facDes)
}

if(do_BSS){
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

if(do_BSS){
  # now bss
  bss_preds <- list()
  for (i in 1:n_pc) bss_preds[[i]] = outer(bases[,i],bps[[i]]$yhat)
  bssEmu = mean0pred + Reduce('+', bss_preds)
  
  # now bss (main effects only)
  bssm_preds <- list()
  for (i in 1:n_pc) bssm_preds[[i]] = outer(bases[,i],bmps[[i]]$yhat)
  bssEmu_meo = mean0pred + Reduce('+', bssm_preds)
}


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

if(des %in% c("a","c")){
  load(paste0("rda/modRuns_",ifelse(hi,"hi","lo"),"_",des,".rda"))
  par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1.5,1))
  matplot(smvals,modRuns,type='l',ylab='GSMF_A')
  mtext('true testing set',side=3,line=.1,cex=.9)
  
  matplot(smvals,etaEmu,type='l',ylab='emulator')
  mtext(paste0('GP-PC predictions'),side=3,line=.1,cex=.9)
  
  matplot(smvals,bssEmu,type='l',ylab='bss emulator')
  mtext(paste0('bss-anova (MEs and 2WIs) predictions'),side=3,line=.1,cex=.9)
  
  matplot(smvals,bssEmu_meo,type='l',ylab='bss (main effects only) emulator')
  mtext(paste0('bss-anova (MEs only) predictions'),side=3,line=.1,cex=.9)
  
  matplot(smvals,modRuns,type='l',ylab='GSMF_A')
  mtext('true testing set',side=3,line=.1,cex=.9)
  
  matplot(smvals,etaEmu - modRuns,type='l',ylab='emulator')
  mtext(paste0('GP-PC errors'),side=3,line=.1,cex=.9)
  
  matplot(smvals,bssEmu - modRuns,type='l',ylab='bss emulator')
  mtext(paste0('bss-anova (MEs and 2WIs) errors'),side=3,line=.1,cex=.9)
  
  matplot(smvals,bssEmu_meo - modRuns,type='l',ylab='bss (main effects only) emulator')
  mtext(paste0('bss-anova (MEs only) errors'),side=3,line=.1,cex=.9)
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
if(WRITE) save.image(paste0("/projects/precipit/hydro/hydro-SA-",
                            ifelse(hi,"hi-","lo-"),"train",if(atf){"a"},if(ctf){"c"},"-",
                            "pred",des,if(des=="FF"){n_ff},
                            if(des%in%c("unif","lhs","lhsS","unifS")){m},
                            if(des %in% c("lhsS","unifS")){paste0("s",seed)},
                            "-nPC",n_pc,".rda"))

