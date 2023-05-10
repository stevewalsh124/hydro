library(GPfit)
library(parallel)

ncores <- 1

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

print(ncores)

p = 4
nv = 9
prange = seq(0,1,length=nv)

load("rda/smvals_full.rda")
smvals = smvals_full[which(smvals_full >= 10 & smvals_full <= 11.5)]

# make some emulator runs - each column is a model run
# modRuns = simModel(udesign)
load("rda/modRuns_hydro.rda")
modRuns_first32 <- modRuns
load("rda/modRuns_hydro_new32.rda")
modRuns <- cbind(modRuns_first32, modRuns_hydro)
nruns <- ncol(modRuns)

# read in the space filling LHC design from the Mira-Titan paper
des0 = read.table("hydro.design")
des00 = read.table("hydro.design_c")
des1 <- rbind(des0, des00)
colnames(des1) <- c("kappa", "EGW", "NPERH_AGN", "SeedMass")
# m = nrow(des0)
maxvals = c(8,4,4,4e6)
minvals = c(2,.25,1.5,4e5)

des01 = matrix(NA,nrow=nrow(des1),ncol=p)
for(k in 1:p) des01[,k] = (des1[,k]-minvals[k])/(maxvals[k]-minvals[k])

# create the simulation output for this design
eta = modRuns

n_sm <- length(smvals)

# for fun, let's use SVD and see the dimensionality in eta
mean0mat = matrix(apply(eta,1,mean),nrow=n_sm,ncol=nruns)
eta0 = eta - mean0mat
a = svd(eta0)

# look at coefficients for each basis function
coef1 = a$v[,1]
# scale the coefficients so they have variance = 1
coef = a$v*sqrt(nruns)
# accordingly, scale the bases so when multiplied by
# the coef's, we'll get the spectra back
bases = a$u%*%diag(a$d)/sqrt(nruns)

spectraFull = bases%*%t(coef)

# try fitting gp's to the coefficients
system.time({
  a1 = GP_fit(des01,coef[,1])
  a2 = GP_fit(des01,coef[,2])
  a3 = GP_fit(des01,coef[,3])
  a4 = GP_fit(des01,coef[,4])
  a5 = GP_fit(des01,coef[,5])
  a6 = GP_fit(des01,coef[,6])
  a7 = GP_fit(des01,coef[,7])
  a8 = GP_fit(des01,coef[,8])
  # a9 = GP_fit(des01,coef[,9])
  # a10 = GP_fit(des01,coef[,10])
  # a11 = GP_fit(des01,coef[,11])
})

system.time({
GP_fn <- function(i) GP_fit(des01,coef[,i])
r <- mclapply(1:8, GP_fn, mc.cores = ncores-1)
})

all.equal(r[[1]], a1)
all.equal(r[[2]], a2)
all.equal(r[[3]], a3)
all.equal(r[[4]], a4)
all.equal(r[[5]], a5)
all.equal(r[[6]], a6)
all.equal(r[[7]], a7)
all.equal(r[[8]], a8)
# all.equal(r[[9]], a9)
# all.equal(r[[10]], a10)
# all.equal(r[[11]], a11)

