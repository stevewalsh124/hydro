# Pre-processing hydro data
# 1) isolate the relevant stellar mass (x) and model run output (y)
# 2) rescale the designs to be in the unit hypercube [0,1]^p
# 3) if doing hi-res, interpolate y for -Inf cases

plot.it <- F # plot the runs?

# Hydro runs are done in batches of 32 
# (orthogonal symmetric LHSs; first 16 symmetric 
# w the second 16, 180* rotation)

# Which set of runs do you want to train on?
# TRUE/FALSE for including designs a and/or c
# TRUE/FALSE for using hires runs
atf <- T
ctf <- T
if(!atf & !ctf) stop("you didn't pick any data")

# do you want to use hires (TRUE; 32MPC) or lowres (FALSE; 128MPC) runs?
hi <- T

# Read in the design information for the inputs
if(atf) des = read.table("hydro.design") # design a, first batch of 32
if(ctf) des = read.table("hydro.design_c") # design c, second bath of 32
if(atf & ctf) des = rbind(read.table("hydro.design"),read.table("hydro.design_c"))

# m = number of runs total you are training on
m = nrow(des)
p = ncol(des)

# four parameters are in a [0,1]^4 grid; have to scale params first
colnames(des) <- c("kappa", "EGW", "NPERH_AGN", "SeedMass")
maxvals = c(8,  4,  4,4e6)
minvals = c(2,.25,1.5,4e5)

kappas <- des[,1]
EGWs <- des[,2]
NPERH_AGNs <- des[,3]
SeedMasss <- des[,4]

# scale the design to the unit hypercube [0,1]^p
des01 = matrix(NA,nrow=nrow(des),ncol=p)
colnames(des01) <- colnames(des)
for(k in 1:p) des01[,k] = (des[,k]-minvals[k])/(maxvals[k]-minvals[k])

# plot to check
if(plot.it){
  pairs(des,pch=c(rep('.',32),rep("*",32)))
  pairs(des01,pch=c(rep('.',32),rep("*",32)))
}

# Do the first batch of m=32 runs first
hydro_runs <- list()

for (i in 1:m) {
  des_dir <- res_end <- NULL
  
  # all hires runs have "extract" directory
  # lowres design c runs have "extract" directory
  # lowres design a runs have "profile" directory  
  if(hi) des_dir <- rep("extract", m) else {
    if(atf) des_dir <- c(des_dir, rep("profiles", 32))
    if(ctf) des_dir <- c(des_dir, rep("extract", 32))
  }
  
  # design c has "_NEW" appended, design a does not
  res_dir <- ifelse(hi, "32MPC_RUNS", "128MPC_RUNS")
  if(atf) res_end <- c(res_end, rep("/",32))
  if(ctf) res_end <- c(res_end, rep("_NEW/",32))

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
  
  hydro_runs[[i]] <- read.table(paste0(res_dir, res_end[i],
                                       run_dir, 
                                       "/analysis_pipeline/",des_dir[i],"/GalStellarMassFunction_624.txt"))
}

# plot all smvals
sm_lb <- ifelse(hi, 8.5, 10); sm_ub <- ifelse(hi, 11, 11.5) 
sub_rg <-  which(log10(hydro_runs[[1]][,1]) >= sm_lb &  log10(hydro_runs[[1]][,1]) <= sm_ub)
smvals <- log10(hydro_runs[[1]][,1])[sub_rg]

n_sm <- length(smvals)
modRuns <- matrix(NA, n_sm, m)
for (i in 1:m) modRuns[,i] <- log10(hydro_runs[[i]][sub_rg,2])

if(hi){
  bad_boys <- which(is.infinite(modRuns[nrow(modRuns),]))
  if(plot.it) matplot(modRuns[,bad_boys], type="l",  main = "which runs have -Inf for highest smass")
  
  # only the highest value of smass for hi res runs has a couple -Infs
  for (j in which(is.infinite(modRuns[nrow(modRuns),]))) {
    smvals01 <- (smvals - min(smvals))/(max(smvals) - min(smvals))
    # interpy <- GP_fit(smvals01[-nrow(modRuns)],modRuns[-nrow(modRuns),j])
    # modRuns[nrow(modRuns),j] <- predict(interpy, smvals01[nrow(modRuns)])$Y_hat
    fitty <- smooth.spline(smvals01[-nrow(modRuns)],modRuns[-nrow(modRuns),j])
    modRuns[nrow(modRuns),j] <- predict(fitty, smvals01[nrow(modRuns)])$y
  }
  if(plot.it) matplot(modRuns[,bad_boys], type="l", main = "interpolated the end")
}


if(plot.it){
  plot(smvals, modRuns[,1], type="l", 
       xlab = "log10(Stellar Mass)", ylab = "log10(GSMF_Apperture)", 
       main = paste0("Subset x in [",sm_lb,",",sm_ub,"]"), 
       ylim = range(modRuns[which(is.finite(modRuns))]))
  for (i in 2:m) { lines(smvals, modRuns[,i], col=i)}
}

if(atf | ctf){
  save(des01, file = paste0("rda/des_",if(atf){"a"},if(ctf){"c"},".rda"))
  save(smvals, file = paste0("rda/smvals_",ifelse(hi,"hi","lo"),".rda"))
  save(modRuns, file = paste0("rda/modRuns_",ifelse(hi,"hi","lo"),"_",
                              if(atf){"a"},if(ctf){"c"},".rda"))
}

