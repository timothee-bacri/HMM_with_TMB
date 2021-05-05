# packages
require("xlsx")
require("microbenchmark")

# data and functions
load("tyt_dataset.RData")
load("Functions_HMM_2012_03_28.RData")
dat <- gimme_data

# aux variables
no.sub <- length(dat)
tab.res <- as.data.frame(matrix(NA, nrow = no.sub, ncol = 12))
names(tab.res) <- c("delta1", "delta2", "gam11", "gam22","lambda1",
                    "lambda2", "llk.iid", "llk.hmm", "AIC.iid", "AIC.hmm", 
                    "BIC.iid", "BIC.hmm")
tab.res.mood <- tab.res
tab.res.arousal <- tab.res
t.width <- 25
t.height <- 21


# mood variable
# -------------

# graphical output to pdf
pdf("Figures tyt data/mood.pdf", height = t.height, width = t.width)
layout(matrix(c(1, 3, 2, 3), 2))
par(las = 1, mai=c(.6, .6, .3, .1),  mgp=c(3, 1, 0), mar=c(4.5, 7.5, 3, 2), xaxs = "i", yaxs = "i", 
    cex.lab=2, cex.axis=2, cex.main = 2.5)      

for (i in 1 : no.sub) {
  # fit models to mood variable
  obs <- dat[[i]][, "mood"] * 8
  gam.ini <- matrix(c(0.95, 0.05, 0.05, 0.95), 2)
  lambda.ini <- c(quantile(obs, 0.05), quantile(obs, 0.95))
  res <- estimate.pois.HMM(obs, gam = gam.ini, lambda = lambda.ini)
  res.fwbw <- FWBW.HMM(obs, delta = res$delta, gam = res$gam, lambda = res$lambda, distn = "Poisson")
  state.seq <- apply(res.fwbw, 2, which.max)
  res
  # histogram
  hist(obs, prob = TRUE, breaks = seq(min(obs) - 0.5, max(obs) + 0.5, by = 1))
  t.x <- seq(min(obs), max(obs), by = 1)  
  t.y <- res$delta[1] * dpois(t.x, lambda = res$lambda[1]) + 
         res$delta[2] * dpois(t.x, lambda = res$lambda[2])
  lines(t.x, t.y)
  # acf
  acf(obs)
  # scatterplot
  plot(obs, type = "l")
  t.mean <- rep(NA, length(obs))
  t.mean[state.seq == 1] <- res$lambda[1]
  t.mean[state.seq == 2] <- res$lambda[2]
  lines(t.mean, lwd = 2, col = "blue")
  # llk, aic, bic
  llk.iid <- sum(log(dpois(obs, lambda = mean(obs))))
  llk.hmm <- res$llkmax
  AIC.iid <- -2 * llk.iid + 2 * 1
  AIC.hmm <- -2 * llk.hmm + 2 * 4
  BIC.iid <- -2 * llk.iid + log(length(obs)) * 1
  BIC.hmm <- -2 * llk.hmm + log(length(obs)) * 4
  # output to table
  tab.res.mood[i, "delta1"] <- res$delta[1]   
  tab.res.mood[i, "delta2"] <- res$delta[2]   
  tab.res.mood[i, "gam11"] <- res$gam[1, 1]   
  tab.res.mood[i, "gam22"] <- res$gam[2, 2]  
  tab.res.mood[i, "lambda1"] <- res$lambda[1]   
  tab.res.mood[i, "lambda2"] <- res$lambda[2]  
  tab.res.mood[i, "llk.iid"]  <- llk.iid 
  tab.res.mood[i, "llk.hmm"]  <- llk.hmm
  tab.res.mood[i, "AIC.iid"]  <- AIC.iid
  tab.res.mood[i, "AIC.hmm"]  <- AIC.hmm
  tab.res.mood[i, "BIC.iid"]  <- BIC.iid
  tab.res.mood[i, "BIC.hmm"]  <- BIC.hmm
}

dev.off()

# export of results
write.xlsx(tab.res.mood, "tyt res fitting HMM to mood data.xlsx")




# arousal variable
# ----------------

# graphical output to pdf
pdf("Figures tyt data/arousal.pdf", height = t.height, width = t.width)
layout(matrix(c(1, 3, 2, 3), 2))
par(las = 1, mai=c(.6, .6, .3, .1),  mgp=c(3, 1, 0), mar=c(4.5, 7.5, 3, 2), xaxs = "i", yaxs = "i", 
    cex.lab=2, cex.axis=2, cex.main = 2.5)      

for (i in 1 : no.sub) {
  # fit models to arousal variable
  obs <- dat[[i]][, "arousal"] * 8
  gam.ini <- matrix(c(0.95, 0.05, 0.05, 0.95), 2)
  lambda.ini <- c(quantile(obs, 0.05), quantile(obs, 0.95))
  res <- estimate.pois.HMM(obs, gam = gam.ini, lambda = lambda.ini)
  res.fwbw <- FWBW.HMM(obs, delta = res$delta, gam = res$gam, lambda = res$lambda, distn = "Poisson")
  state.seq <- apply(res.fwbw, 2, which.max)
  res
  # histogram
  hist(obs, prob = TRUE, breaks = seq(min(obs) - 0.5, max(obs) + 0.5, by = 1))
  t.x <- seq(min(obs), max(obs), by = 1)  
  t.y <- res$delta[1] * dpois(t.x, lambda = res$lambda[1]) + 
    res$delta[2] * dpois(t.x, lambda = res$lambda[2])
  lines(t.x, t.y)
  # acf
  acf(obs)
  # scatterplot
  plot(obs, type = "l")
  t.mean <- rep(NA, length(obs))
  t.mean[state.seq == 1] <- res$lambda[1]
  t.mean[state.seq == 2] <- res$lambda[2]
  lines(t.mean, lwd = 2, col = "blue")
  # llk, aic, bic
  llk.iid <- sum(log(dpois(obs, lambda = mean(obs))))
  llk.hmm <- res$llkmax
  AIC.iid <- -2 * llk.iid + 2 * 1
  AIC.hmm <- -2 * llk.hmm + 2 * 4
  BIC.iid <- -2 * llk.iid + log(length(obs)) * 1
  BIC.hmm <- -2 * llk.hmm + log(length(obs)) * 4
  # output to table
  tab.res.arousal[i, "delta1"] <- res$delta[1]   
  tab.res.arousal[i, "delta2"] <- res$delta[2]   
  tab.res.arousal[i, "gam11"] <- res$gam[1, 1]   
  tab.res.arousal[i, "gam22"] <- res$gam[2, 2]  
  tab.res.arousal[i, "lambda1"] <- res$lambda[1]   
  tab.res.arousal[i, "lambda2"] <- res$lambda[2]  
  tab.res.arousal[i, "llk.iid"]  <- llk.iid 
  tab.res.arousal[i, "llk.hmm"]  <- llk.hmm
  tab.res.arousal[i, "AIC.iid"]  <- AIC.iid
  tab.res.arousal[i, "AIC.hmm"]  <- AIC.hmm
  tab.res.arousal[i, "BIC.iid"]  <- BIC.iid
  tab.res.arousal[i, "BIC.hmm"]  <- BIC.hmm
}

dev.off()

# export of results
write.xlsx(tab.res.arousal, "tyt res fitting HMM to arousal data.xlsx")



# some testing with independent mixture models
# --------------------------------------------

# simulated obs
obs2 <- c(rpois(250, lambda = 2), rpois(150, lambda = 10))
hist(obs2)

# flexmix
require("flexmix")
res.mix1 <- stepFlexmix(obs2 ~ 1, k = 2 : 3, nrep = 50, model = FLXMCmvpois())
summary(res.mix1)
bmod <- getModel(res.mix1, which = "BIC")
summary(bmod)
bmod.refit <- refit(bmod)
summary(bmod)
parameters(bmod)
fitted(bmod)
detach("package:flexmix", unload = TRUE)

# gamlss.mx
require("gamlss.mx")
res.mix2 <- gamlssMX(obs2 ~ 1, family = PO, K=2, data = as.data.frame(obs2))
res.mix2
AIC(res.mix2)
formula(res.mix2)
summary(res.mix2)
coef(res.mix2) # gives only first coef.
exp(res.mix2$models[[1]]$mu.coefficients)
exp(res.mix2$models[[2]]$mu.coefficients) 
detach("package:gamlss.mx", unload = TRUE)



# speed check
# -----------

i <- 12

# fit models to mood variable
obs <- dat[[i]][, "mood"] * 8
gam.ini <- matrix(c(0.95, 0.05, 0.05, 0.95), 2)
lambda.ini <- c(quantile(obs, 0.05), quantile(obs, 0.95))
estimate.pois.HMM(obs, gam = gam.ini, lambda = lambda.ini)
res <- microbenchmark(estimate.pois.HMM(obs, gam = gam.ini, lambda = lambda.ini), times = 1000)
res
quantile(res$time, c(0.025, 0.975)) / 10 ^ 7

myf <- function()
{
  #set.seed(56738)
  x <- runif(10^4)
  sin(x)
  return(x)
}
microbenchmark(myf(), times = 1000)








