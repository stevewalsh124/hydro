# Quantify uncertainty of main effect estimates via a LOOCV approach

mainEffs_bss_meos <- mainEffs_bsss <- mainEffs_gppcs <- list()
for (i in 1:64) {
  load(paste0("/projects/precipit/hydro/rda/mainEffs_bss_meo",i))
  mainEffs_bss_meos[[i]] <- mainEffs_bss_meo
  load(paste0("/projects/precipit/hydro/rda/mainEffs_bss",i))
  mainEffs_bsss[[i]] <- mainEffs_bss
  load(paste0("/projects/precipit/hydro/rda/mainEffs_gppc",i))
  mainEffs_gppcs[[i]] <- mainEffs_gppc
}

p=4
rm_avg <- T

load("rda/smvals_full.rda")
smvals <- smvals_full[which(smvals_full <= 11.5 & smvals_full >= 10)]
load("rda/modRuns_hydro_all64.rda")
if(rm_avg) avg_mean <- rowMeans(modRuns) else avg_mean <- rep(0, nrow(modRuns))

pdf("pdf/LOOCV_MEs_eval.pdf")

pnames = c(paste0("param",1:4))
par(mfrow=c(1,4),oma=c(4,4,1.5,1),mar=c(0,0,0,0))
grcolors = paste0("grey",round(seq(90,25,length=11)),sep='')

# BSS-ANOVA (main effects only) UQ for MEs
yr = Reduce("range", lapply(mainEffs_bss_meos, function(x) range(x - avg_mean)))
for(k in 1:p){
  matplot(smvals,mainEffs_bss_meo[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  for (i in 1:64) {
    matplot(smvals,mainEffs_bss_meos[[i]][,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F,add=T); box() 
  }
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)', if(rm_avg){" [average removed]"}),side=2,line=2.5,outer=T)
mtext("bss-anova (MEs only) main effect estimates", side=3,outer = T)

yr = Reduce("range", lapply(mainEffs_bss_meos, range)) 
for(k in 1:p){
  matplot(smvals,mainEffs_bss_meo[,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  for (i in 1:64) {
    matplot(smvals,mainEffs_bss_meos[[i]][,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F,add=T); box() 
  }
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)'),side=2,line=2.5,outer=T)
mtext("bss-anova (MEs only) main effect estimates", side=3,outer = T)

# BSS-ANOVA (MEs and 2WIs) UQ for MEs
yr = Reduce("range", lapply(mainEffs_bsss, function(x) range(x - avg_mean))) #range(mainEffs_bss - avg_mean)
for(k in 1:p){
  matplot(smvals,mainEffs_bss[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  for (i in 1:64) {
    matplot(smvals,mainEffs_bsss[[i]][,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F,add=T); box() 
  }
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)', if(rm_avg){" [average removed]"}),side=2,line=2.5,outer=T)
mtext("bss-anova (MEs and 2WIs) main effect estimates", side=3,outer = T)

yr = Reduce("range", lapply(mainEffs_bsss, range)) 
for(k in 1:p){
  matplot(smvals,mainEffs_bss[,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  for (i in 1:64) {
    matplot(smvals,mainEffs_bsss[[i]][,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F,add=T); box() 
  }
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)'),side=2,line=2.5,outer=T)
mtext("bss-anova (MEs and 2WIs) main effect estimates", side=3,outer = T)

# GP-PC UQ for MEs
yr = Reduce("range", lapply(mainEffs_gppcs, function(x) range(x - avg_mean))) #range(mainEffs_gppc - avg_mean)
for(k in 1:p){
  matplot(smvals,mainEffs_gppc[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  for (i in 1:64) {
    matplot(smvals,mainEffs_gppcs[[i]][,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F,add=T); box() 
  }
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)', if(rm_avg){" [average removed]"}),side=2,line=2.5,outer=T)
mtext("GP-PC main effect estimates", side=3,outer = T)

yr = Reduce("range", lapply(mainEffs_gppcs, range)) 
for(k in 1:p){
  matplot(smvals,mainEffs_gppc[,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  for (i in 1:64) {
    matplot(smvals,mainEffs_gppcs[[i]][,,k],type='l',col=grcolors,lty=1,ylim=yr,axes=F,add=T); box() 
  }
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1)) axis(2)
  if(k %in% c(1:4)) axis(1)
}
mtext('log10(Stellar Mass)',side=1,line=2.5,outer=T)
mtext(paste0('log10(GSMF_A)'),side=2,line=2.5,outer=T)
mtext("GP-PC main effect estimates", side=3,outer = T)

dev.off()

