# evaluate LOOCV hydro predictions

pdf("pdf/eval_hydro.pdf")

# load sim info and corresponding preds from LOOCV
load("rda/modRuns_hydro.rda")
load("rda/LOOCV_preds.rda")
load("rda/LOOCV_preds_b.rda")
load("rda/LOOCV_preds_bm.rda")
load("rda/modRuns_pred.rda")
load("rda/new32_preds.rda")
load("rda/new32_preds_b.rda")
load("rda/new32_preds_bm.rda")

# overall root mean squared error (RMSEs) for each approach
sqrt(mean((modRuns - LOOCV_preds)^2))
sqrt(mean((modRuns - LOOCV_preds_b)^2))
sqrt(mean((modRuns - LOOCV_preds_bm)^2))

hist((modRuns - LOOCV_preds)^2, main = paste("LOOCV RMSEs for GP-PC:", round(sqrt(mean((modRuns - LOOCV_preds)^2)),4)))
hist((modRuns - LOOCV_preds_b)^2, main = paste("LOOCV RMSEs for BSS-ANOVA:", round(sqrt(mean((modRuns - LOOCV_preds_b)^2)),4)))
hist((modRuns - LOOCV_preds_bm)^2, main = paste("LOOCV RMSEs for BSS-ANOVA (MEs only)", round(sqrt(mean((modRuns - LOOCV_preds_bm)^2)),4)))

# plot((modRuns - LOOCV_preds)^2, (modRuns - LOOCV_preds_b)^2); abline(0,1)
# plot((modRuns - LOOCV_preds)^2, (modRuns - LOOCV_preds_bm)^2); abline(0,1)
# plot((modRuns - LOOCV_preds_b)^2, (modRuns - LOOCV_preds_bm)^2); abline(0,1)

# RMSE as a function of stellar mass (sm)
rmseXsm <- rmseXsm_b <- rmseXsm_bm <- c()
for (i in 1:nrow(modRuns)) {
  rmseXsm[i] <- sqrt(mean((modRuns[i,] - LOOCV_preds[i,])^2))
  rmseXsm_b[i] <- sqrt(mean((modRuns[i,] - LOOCV_preds_b[i,])^2))
  rmseXsm_bm[i] <- sqrt(mean((modRuns[i,] - LOOCV_preds_bm[i,])^2))
}

plot(rmseXsm, type = "l", ylim=range(rmseXsm, rmseXsm_b, rmseXsm_bm), xlab="stellar mass index", ylab="RMSE", main = "LOOCV")
lines(rmseXsm_b, col="red", lty=2)
lines(rmseXsm_bm, col="blue", lty=3)
legend(legend = c("GP-PC", "BSS (both)", "BSS (MEs only)"), "topleft", lty = 1:3, col=c("black","red","blue"))


# overall root mean squared error (RMSEs) for each approach
sqrt(mean((modRuns_pred - new32_preds)^2))
sqrt(mean((modRuns_pred - new32_preds_b)^2))
sqrt(mean((modRuns_pred - new32_preds_bm)^2))

hist((modRuns_pred - new32_preds)^2, main = paste("batch32 RMSEs for GP-PC:", round(sqrt(mean((modRuns_pred - new32_preds)^2)),4)))
hist((modRuns_pred - new32_preds_b)^2, main = paste("batch32 RMSEs for BSS-ANOVA:", round(sqrt(mean((modRuns_pred - new32_preds_b)^2)),4)))
hist((modRuns_pred - new32_preds_bm)^2, main = paste("batch32 RMSEs for BSS-ANOVA (MEs only)", round(sqrt(mean((modRuns_pred - new32_preds_bm)^2)),4)))

# plot((modRuns_pred - new32_preds)^2, (modRuns_pred - new32_preds_b)^2); abline(0,1)
# plot((modRuns_pred - new32_preds)^2, (modRuns_pred - new32_preds_bm)^2); abline(0,1)
# plot((modRuns_pred - new32_preds_b)^2, (modRuns_pred - new32_preds_bm)^2); abline(0,1)

# RMSE as a function of stellar mass (sm)
rmseXsm <- rmseXsm_b <- rmseXsm_bm <- c()
for (i in 1:nrow(modRuns_pred)) {
  rmseXsm[i] <- sqrt(mean((modRuns_pred[i,] - new32_preds[i,])^2))
  rmseXsm_b[i] <- sqrt(mean((modRuns_pred[i,] - new32_preds_b[i,])^2))
  rmseXsm_bm[i] <- sqrt(mean((modRuns_pred[i,] - new32_preds_bm[i,])^2))
}

plot(rmseXsm, type = "l", ylim=range(rmseXsm, rmseXsm_b, rmseXsm_bm), xlab="stellar mass index", ylab="RMSE", main = "batch32")
lines(rmseXsm_b, col="red", lty=2)
lines(rmseXsm_bm, col="blue", lty=3)
legend(legend = c("GP-PC", "BSS (both)", "BSS (MEs only)"), "topleft", lty = 1:3, col=c("black","red","blue"))

dev.off()