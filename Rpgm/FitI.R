root <- "C:/Users/Dimitris/Documents/Papers/Paper9/"
rootOut <- paste(root, "Results/Application/", sep = "")
library(R2WinBUGS)

source(paste(root, "Rpgm/Application/PrepareData.R", sep = ""))
source(paste(root, "Rpgm/Application/Models.R", sep = ""))
setwd(rootOut)


Data$nb <- with(Data, ncZ1 + ncZ2 + ncZ3)
Data$mu0 <- rep(0, Data$nb)
Data$priorR.D <- diag(rep(0.25, Data$nb))
DD <- Data[names(Data) %in% namesPwConst.DP.I]

write.model(PwConst.DP.I, file.path(getwd(), "PwConst.DP.I.txt"))
inits <- list(
    betas1 = c(19.5474516, -0.3092732, 2.0592278, 0.4997537, -2.0044231, 5.4755351, 
               -6.2852801),
    betas2 = c(23.70442388, 0.06231521, 1.56502518, 0.02787711, -0.60827309, 18.99244089, 
               4.85040687),
    betas3 = c(-2.26867314, -0.02474365, -0.12157263, 0.02175284, 0.54916985, 
               -0.22119048, 1.19965498),
    gammas = c(-0.01976846, 0.22264399, 0.01374726),
    alphas = c(-5, -1, -5),
    tau1 = 0.05500224,
    tau2 = 0.1647081,
    rho = 1,
    latent = sample(DD$R, DD$N, TRUE)
)

parms <- c("betas1", "betas2", "betas3", "tau1", "tau2", "gammas", "alphas", "xi", "b")

fit <- bugs(DD, list(inits), parms, model.file = "PwConst.DP.I.txt", n.chains = 1, 
            n.iter = 120000, n.burnin = 20000, n.thin = 50, working.directory = getwd(), 
            clearWD = TRUE)

root. <- paste(rootOut, "fitI.RData", sep = "")
save(fit, file = root.)


###############################################################################################

library(R2WinBUGS)
root <- "C:/Users/Dimitris/Documents/Papers/Paper9/"
load(paste(root, "Results/Application/ParamI/fitI.RData", sep = ""))

fitCoda <- as.mcmc.list(fit)
fitCoda[[1]] <- fitCoda[[1]][, !colnames(fitCoda[[1]]) %in% "deviance"]
plot(fitCoda, ask = TRUE)
summary(fitCoda)

nams <- c("alphas[1]", "alphas[2]", "alphas[3]", "betas1[1]", "betas1[2]", "betas1[3]", 
          "betas1[4]", "betas1[5]", "betas1[6]", "betas1[7]", 
          "betas2[1]", "betas2[2]", "betas2[3]", "betas2[4]", "betas2[5]", "betas2[6]", 
          "betas2[7]", "betas3[1]", "betas3[2]", "betas3[3]", "betas3[4]",
          "betas3[5]", "betas3[6]", "betas3[7]", "deviance", "gammas[1]", "gammas[2]", 
          "gammas[3]", "tau1", "tau2", "xi[1]", "xi[2]", "xi[3]", "xi[4]")
fitCoda. <- fitCoda[[1]][, nams]


str(summary(fitCoda.))



sigma1 <- sqrt(1 / fitCoda[[1]][, "tau1"])
summary(sigma1)

sigma2 <- sqrt(1 / fitCoda[[1]][, "tau2"])
summary(sigma2)



gc()
