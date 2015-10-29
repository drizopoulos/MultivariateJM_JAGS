PwConst.DP.I <- function () {
    for (i in 1:N) {
        # Longitudinal Part -- Outcome 1
        for (j in offset1[i]:(offset1[i+1] - 1)) {
            muy1[j] <- inprod(betas1[1:ncX1], X1[j, 1:ncX1]) + 
                inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1])
            y1[j] ~ dnorm(muy1[j], tau1)
        }
        # Longitudinal Part -- Outcome 2
        for (j in offset2[i]:(offset2[i+1] - 1)) {
            muy2[j] <- inprod(betas2[1:ncX1], X2[j, 1:ncX1]) + 
                inprod(b[i, (ncZ1 + 1):(ncZ1 + ncZ2)], Z2[j, 1:ncZ2])
            y2[j] ~ dnorm(muy2[j], tau2)
        }
        # Longitudinal Part -- Outcome 3
        for (j in offset3[i]:(offset3[i+1] - 1)) {
            logit(PP[j]) <- inprod(betas3[1:ncX3], X3[j, 1:ncX3]) + 
                inprod(b[i, (ncZ1 + ncZ2 + 1):(ncZ1 + ncZ2 + ncZ3)], Z3[j, 1:ncZ3])
            pr[j] <- max(0.00000001, min(0.99999999, PP[j]))
            y3[j] ~ dbern(pr[j])
        }
        # Survival Part
        etaBaseline[i] <- inprod(gammas[1:ncW], W[i, 1:ncW])
        f1.T[i] <- inprod(betas1[1:ncX1], X1t[i, 1:ncX1]) + 
            inprod(b[i, 1:ncZ1], Z1t[i, 1:ncZ1])
        f2.T[i] <- inprod(betas2[1:ncX2], X2t[i, 1:ncX1]) + 
            inprod(b[i, (ncZ1 + 1):(ncZ1 + ncZ2)], Z2t[i, 1:ncZ2])
        f3.T[i] <- inprod(betas3[1:ncX3], X3t[i, 1:ncX1]) + 
            inprod(b[i, (ncZ1 + ncZ2 + 1):(ncZ1 + ncZ2 + ncZ3)], Z3t[i, 1:ncZ3])
        for (q in 1:Q) {
            log.hazard[i, q] <- D[i, q] * (log(xi[q]) + etaBaseline[i] + 
                                               alphas[1] * f1.T[i] + alphas[2] * f2.T[i] + 
                                               alphas[3] * f3.T[i])
            # loop over the Gauss-Kronrod quadrature points
            for (k in 1:15) {
                f1.s[i, q, k] <- inprod(betas1[1:ncX1], X1s[Q*K*(i - 1) + K*(q - 1) + k, 1:ncX1]) + 
                    inprod(b[i, 1:ncZ1], Z1s[Q*K*(i - 1) + K*(q - 1) + k, 1:ncZ1])
                f2.s[i, q, k] <- inprod(betas2[1:ncX2], X2s[Q*K*(i - 1) + K*(q - 1) + k, 1:ncX2]) + 
                    inprod(b[i, (ncZ1 + 1):(ncZ1 + ncZ2)], Z2s[Q*K*(i - 1) + K*(q - 1) + k, 1:ncZ2])
                f3.s[i, q, k] <- inprod(betas3[1:ncX3], X3s[Q*K*(i - 1) + K*(q - 1) + k, 1:ncX3]) + 
                    inprod(b[i, (ncZ1 + ncZ2 + 1):(ncZ1 + ncZ2 + ncZ3)], 
                           Z3s[Q*K*(i - 1) + K*(q - 1) + k, 1:ncZ3])
                SurvLong[i, q, k] <- wk[k] * exp(alphas[1] * f1.s[i, q, k] + 
                                                     alphas[2] * f2.s[i, q, k] + 
                                                     alphas[3] * f3.s[i, q, k])
            }
            approxIntegral[i, q] <- P[i, q] * sum(SurvLong[i, q, ])
            log.survival[i, q] <- - xi[q] * exp(etaBaseline[i]) * approxIntegral[i, q]
        }
        phi[i] <- C - sum(log.hazard[i, ]) - sum(log.survival[i, ])
        zeros[i] ~ dpois(phi[i])
        # Random Effects Part
        latent[i] ~ dcat(pis[1:R])
        for (l in 1:nb) {
            b[i, l] <- lambda[latent[i], l]
        }
        for (k in 1:R) {
            memb[i, k] <- equals(latent[i], k)
        }
    }
    # Priors
    # Longitudinal Part -- Outcome 1
    betas1[1:ncX1] ~ dmnorm(priorMean.betas1[], priorTau.betas1[, ])
    tau1 ~ dgamma(priorA.tau1, priorB.tau1)
    # Longitudinal Part -- Outcome 2
    betas2[1:ncX2] ~ dmnorm(priorMean.betas2[], priorTau.betas2[, ])
    tau2 ~ dgamma(priorA.tau2, priorB.tau2)
    # Longitudinal Part -- Outcome 2
    betas3[1:ncX3] ~ dmnorm(priorMean.betas3[], priorTau.betas3[, ])
    # Survival Part
    gammas[1:ncW] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, ])
    alphas[1:3] ~ dmnorm(priorMean.alphas[], priorTau.alphas[, ])
    for (q in 1:Q) {
        xi[q] ~ dgamma(priorA.xi[q], priorB.xi[q])
    }
    # Random Effects Part
    for (k in 1:R) {
        Tmemb[k] <- sum(memb[, k])
        Fmemb[k] <- step(Tmemb[k] - 1)
    }
    Tcluster <- sum(Fmemb[])
    # constructive DP
    for (k in 1:R) {
        r[k] ~ dbeta(1, rho)
        lambda[k, 1:nb] ~ dmnorm(mu0[], inv.D[, ])
    }
    # stick-breaking
    pis[1] <- r[1]
    for (k in 2:(R-1)) {
        log(pis[k]) <- log(r[k]) + sum(RR[k, 1:k-1])
        for (l in 1:k-1) {
            RR[k, l] <- log(1 - r[l])
        }
    }
    pis[R] <- 1 - sum(pis[1:(R-1)])
    inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)   
    rho ~ dgamma(priorA.rho, priorB.rho)
}

namesPwConst.DP.I <- c(
    "N", "C", 
    "offset1", "y1", "X1", "X1t", "X1s", "Z1", "Z1t", "Z1s", "ncX1", "ncZ1",
    "offset2", "y2", "X2", "X2t", "X2s", "Z2", "Z2t", "Z2s", "ncX2", "ncZ2",
    "offset3", "y3", "X3", "X3t", "X3s", "Z3", "Z3t", "Z3s", "ncX3", "ncZ3",
    "D", "Q", "P", "wk", "W", "ncW", "K", "zeros",
    "nb", "mu0", "R",
    "priorMean.betas1", "priorTau.betas1", "priorA.tau1", "priorB.tau1",
    "priorMean.betas2", "priorTau.betas2", "priorA.tau2", "priorB.tau2",
    "priorMean.betas3", "priorTau.betas3",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
    "priorA.xi", "priorB.xi", "priorR.D", "priorK.D", "priorA.rho", "priorB.rho"
)

###############################################################################################################
###############################################################################################################

PwConst.DP.III <- function () {
    for (i in 1:N) {
        # Longitudinal Part -- Outcome 1
        for (j in offset1[i]:(offset1[i+1] - 1)) {
            muy1[j] <- inprod(betas1[1:ncX1], X1[j, 1:ncX1]) + 
                inprod(b[i, 1:ncZ1], Z1[j, 1:ncZ1])
            y1[j] ~ dnorm(muy1[j], tau1)
        }
        # Longitudinal Part -- Outcome 2
        for (j in offset2[i]:(offset2[i+1] - 1)) {
            muy2[j] <- inprod(betas2[1:ncX1], X2[j, 1:ncX1]) + 
                inprod(b[i, (ncZ1 + 1):(ncZ1 + ncZ2)], Z2[j, 1:ncZ2])
            y2[j] ~ dnorm(muy2[j], tau2)
        }
        # Longitudinal Part -- Outcome 3
        for (j in offset3[i]:(offset3[i+1] - 1)) {
            logit(PP[j]) <- inprod(betas3[1:ncX3], X3[j, 1:ncX3]) + 
                inprod(b[i, (ncZ1 + ncZ2 + 1):(ncZ1 + ncZ2 + ncZ3)], Z3[j, 1:ncZ3])
            pr[j] <- max(0.00000001, min(0.99999999, PP[j]))
            y3[j] ~ dbern(pr[j])
        }
        # Survival Part
        etaBaseline[i] <- inprod(gammas[1:ncW], W[i, 1:ncW]) + b[i, nb]
        for (q in 1:Q) {
            mu[i, q] <- xi[q] * exp(etaBaseline[i]) * T[i, q]
            D[i, q] ~ dpois(mu[i, q])
        }
        # Random Effects Part
        latent[i] ~ dcat(pis[1:R])
        for (l in 1:nb) {
            b[i, l] <- lambda[latent[i], l]
        }
        for (k in 1:R) {
            memb[i, k] <- equals(latent[i], k)
        }
    }
    # Priors
    # Longitudinal Part -- Outcome 1
    betas1[1:ncX1] ~ dmnorm(priorMean.betas1[], priorTau.betas1[, ])
    tau1 ~ dgamma(priorA.tau1, priorB.tau1)
    # Longitudinal Part -- Outcome 2
    betas2[1:ncX2] ~ dmnorm(priorMean.betas2[], priorTau.betas2[, ])
    tau2 ~ dgamma(priorA.tau2, priorB.tau2)
    # Longitudinal Part -- Outcome 3
    betas3[1:ncX3] ~ dmnorm(priorMean.betas3[], priorTau.betas3[, ])
    # Survival Part
    gammas[1:ncW] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, ])
    for (q in 1:Q) {
        xi[q] ~ dgamma(priorA.xi[q], priorB.xi[q])
    }
    # Random Effects Part
    for (k in 1:R) {
        Tmemb[k] <- sum(memb[, k])
        Fmemb[k] <- step(Tmemb[k] - 1)
    }
    Tcluster <- sum(Fmemb[])
    # constructive DP
    for (k in 1:R) {
        r[k] ~ dbeta(1, rho)
        lambda[k, 1:nb] ~ dmnorm(mu0[], inv.D[, ])
    }
    # stick-breaking
    pis[1] <- r[1]
    for (k in 2:(R-1)) {
        log(pis[k]) <- log(r[k]) + sum(RR[k, 1:k-1])
        for (l in 1:k-1) {
            RR[k, l] <- log(1 - r[l])
        }
    }
    pis[R] <- 1 - sum(pis[1:(R-1)])
    inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)   
    rho ~ dgamma(priorA.rho, priorB.rho)
}

namesPwConst.DP.III <- c(
    "N", 
    "offset1", "y1", "X1", "Z1", "ncX1", "ncZ1",
    "offset2", "y2", "X2", "Z2", "ncX2", "ncZ2",
    "offset3", "y3", "X3", "Z3", "ncX3", "ncZ3",
    "D", "T", "Q", "W", "ncW",
    "nb", "mu0", "R",
    "priorMean.betas1", "priorTau.betas1", "priorA.tau1", "priorB.tau1",
    "priorMean.betas2", "priorTau.betas2", "priorA.tau2", "priorB.tau2",
    "priorMean.betas3", "priorTau.betas3", "priorMean.gammas", "priorTau.gammas", 
    "priorA.xi", "priorB.xi", "priorR.D", "priorK.D", "priorA.rho", "priorB.rho"
)

###################################################################################################################
###################################################################################################################

modelI <- function () {
    for (i in 1:N) {
        # Longitudinal Part
        # Outcome 1
        for (j in offset1[i]:(offset1[i + 1] - 1)){
            y1[j] ~ dnorm(mu1[j], tau1)            
            mu1[j] <- inprod(beta1.1[1:ncX1], X1[j, 1:ncX1]) + b[i, 1] + 
                inprod(beta1.2[1:(nk1-1)], Z1[j, 1:(nk1-1)]) + 
                inprod(b[i, 2:nk1], Z1[j, 1:(nk1-1)])
            Z1[j, 1] <- Times1[j]
            for (k in 2:(nk1-1)) {
                P11[j, k] <- pow(max(Times1[j] - knots1[k], 0), 3)
                P21[j, k] <- (knots1[nk1 - 1] - knots1[k]) * 
                    pow(max(Times1[j] - knots1[nk1], 0), 3)
                P31[j, k] <- (knots1[nk1] - knots1[k]) * 
                    pow(max(Times1[j] - knots1[nk1 - 1], 0), 3)
                Z1[j, k] <- P11[j, k] + 
                    (P21[j, k] - P31[j, k]) / (knots1[nk1] - knots1[nk1 - 1])
            }
        }
        # Outcome 2
        for (j in offset2[i]:(offset2[i + 1] - 1)){
            y2[j] ~ dnorm(mu2[j], tau2)
            mu2[j] <- inprod(beta2.1[1:ncX2], X2[j, 1:ncX2]) + b[i, nk1 + 1] + 
                inprod(beta2.2[1:(nk2-1)], Z2[j, 1:(nk2-1)]) + 
                inprod(b[i, (nk1 + 2):(nk1 + nk2)], Z2[j, 1:(nk2-1)])
            Z2[j, 1] <- Times2[j]
            for (k in 2:(nk2-1)) {
                P12[j, k] <- pow(max(Times2[j] - knots2[k], 0), 3)
                P22[j, k] <- (knots2[nk2 - 1] - knots2[k]) * pow(max(Times2[j] - knots2[nk2], 0), 3)
                P32[j, k] <- (knots2[nk2] - knots2[k]) * pow(max(Times2[j] - knots2[nk2 - 1], 0), 3)
                Z2[j, k] <- P12[j, k] + (P22[j, k] - P32[j, k]) / (knots2[nk2] - knots2[nk2 - 1])
            }
        }
        # Outcome 3
        for (j in offset3[i]:(offset3[i + 1] - 1)){
            y3[j] ~ dbern(pr[j])
            pr[j] <- max(0.0000000001, min(0.9999999999, P[j]))
            logit(P[j]) <- inprod(beta3.1[1:ncX3], X3[j, 1:ncX3]) + b[i, nk1 + nk2 + 1] + 
                inprod(beta3.2[1:(nk3-1)], Z3[j, 1:(nk3-1)]) + 
                inprod(b[i, (nk1 + nk2 + 2):(nk1 + nk2 + nk3)], Z3[j, 1:(nk3-1)])
            Z3[j, 1] <- Times3[j]
            for (k in 2:(nk3-1)) {
                P13[j, k] <- pow(max(Times3[j] - knots3[k], 0), 3)
                P23[j, k] <- (knots3[nk3 - 1] - knots3[k]) * pow(max(Times3[j] - knots3[nk3], 0), 3)
                P33[j, k] <- (knots3[nk3] - knots3[k]) * pow(max(Times3[j] - knots3[nk3 - 1], 0), 3)
                Z3[j, k] <- P13[j, k] + (P23[j, k] - P33[j, k]) / (knots3[nk3] - knots3[nk3 - 1])
            }
        }
        # Survival Part
        etaBaseline[i] <- inprod(gamma[], W[i, ])
        f1.T[i] <- beta1.1[1] + beta1.1[2]*W[i, 1] + beta1.1[3]*W[i, 2] + beta1.1[4]*W[i, 3] + b[i, 1] + 
            inprod(beta1.2[1:(nk1-1)], ZT1[i, 1:(nk1-1)]) + inprod(b[i, 2:nk1], ZT1[i, 1:(nk1-1)])
        ZT1[i, 1] <- Time[i]
        for (k in 2:(nk1-1)) {
            PT11[i, k] <- pow(max(Time[i] - knots1[k], 0), 3)
            PT21[i, k] <- (knots1[nk1 - 1] - knots1[k]) * pow(max(Time[i] - knots1[nk1], 0), 3)
            PT31[i, k] <- (knots1[nk1] - knots1[k]) * pow(max(Time[i] - knots1[nk1 - 1], 0), 3)
            ZT1[i, k] <- PT11[i, k] + (PT21[i, k] - PT31[i, k]) / (knots1[nk1] - knots1[nk1 - 1])
        }    
        f2.T[i] <- beta2.1[1] + beta2.1[2]*W[i, 1] + beta2.1[3]*W[i, 2] + 
            beta2.1[4]*W[i, 3] + b[i, nk1 + 1] + 
            inprod(beta2.2[1:(nk2-1)], ZT2[i, 1:(nk2-1)]) + 
            inprod(b[i, (nk1 + 2):(nk1 + nk2)], ZT2[i, 1:(nk2-1)])
        ZT2[i, 1] <- Time[i]
        for (k in 2:(nk2-1)) {
            PT12[i, k] <- pow(max(Time[i] - knots2[k], 0), 3)
            PT22[i, k] <- (knots2[nk2 - 1] - knots2[k]) * pow(max(Time[i] - knots2[nk2], 0), 3)
            PT32[i, k] <- (knots2[nk2] - knots2[k]) * pow(max(Time[i] - knots2[nk2 - 1], 0), 3)
            ZT2[i, k] <- PT12[i, k] + (PT22[i, k] - PT32[i, k]) / (knots2[nk2] - knots2[nk2 - 1])
        }
        f3.T[i] <- beta3.1[1] + beta3.1[2]*W[i, 1] + beta3.1[3]*W[i, 2] + 
            beta3.1[4]*W[i, 3] + b[i, nk1 + nk2 + 1] + 
            inprod(beta3.2[1:(nk3-1)], ZT3[i, 1:(nk3-1)]) + 
            inprod(b[i, (nk1 + nk2 + 2):(nk1 + nk2 + nk3)], ZT3[i, 1:(nk3-1)])
        ZT3[i, 1] <- Time[i]
        for (k in 2:(nk3-1)) {
            PT13[i, k] <- pow(max(Time[i] - knots3[k], 0), 3)
            PT23[i, k] <- (knots3[nk3 - 1] - knots3[k]) * pow(max(Time[i] - knots3[nk3], 0), 3)
            PT33[i, k] <- (knots3[nk3] - knots3[k]) * pow(max(Time[i] - knots3[nk3 - 1], 0), 3)
            ZT3[i, k] <- PT13[i, k] + (PT23[i, k] - PT33[i, k]) / (knots3[nk3] - knots3[nk3 - 1])
        }
        for (q in 1:Q) {
            log.hazard[i, q] <- D[i, q] * (log(xi[q]) + etaBaseline[i] + 
                                               alpha[1]*f1.T[i] + alpha[2]*f2.T[i] + 
                                               alpha[3]*f3.T[i])
            P1[i, q] <- (Up[i, q] + Lo[i, q]) / 2
            P2[i, q] <- (Up[i, q] - Lo[i, q]) / 2
            for (k in 1:15) {
                su[i, q, k] <- P2[i, q] * sk[k] + P1[i, q]
                f1.s[i, q, k] <- beta1.1[1] + beta1.1[2]*W[i, 1] + beta1.1[3]*W[i, 2] + 
                    beta1.1[4]*W[i, 3] + b[i, 1] + 
                            inprod(beta1.2[1:(nk1-1)], ZTs1[i, q, k, 1:(nk1-1)]) + 
                    inprod(b[i, 2:nk1], ZTs1[i, q, k, 1:(nk1-1)])
                ZTs1[i, q, k, 1] <- su[i, q, k]
                for (r in 2:(nk1-1)) {
                    PTs11[i, q, k, r] <- pow(max(su[i, q, k] - knots1[r], 0), 3)
                    PTs21[i, q, k, r] <- (knots1[nk1 - 1] - knots1[r]) * 
                        pow(max(su[i, q, k] - knots1[nk1], 0), 3)
                    PTs31[i, q, k, r] <- (knots1[nk1] - knots1[r]) * 
                        pow(max(su[i, q, k] - knots1[nk1 - 1], 0), 3)
                    ZTs1[i, q, k, r] <- PTs11[i, q, k, r] + 
                        (PTs21[i, q, k, r] - PTs31[i, q, k, r]) / (knots1[nk1] - knots1[nk1 - 1])
                }
                f2.s[i, q, k] <- beta2.1[1] + beta2.1[2]*W[i, 1] + beta2.1[3]*W[i, 2] + 
                    beta2.1[4]*W[i, 3] + b[i, nk1 + 1] + 
                            inprod(beta2.2[1:(nk2-1)], ZTs2[i, q, k, 1:(nk2-1)]) + 
                    inprod(b[i, (nk1 + 2):(nk1 + nk2)], ZTs2[i, q, k, 1:(nk2-1)])
                ZTs2[i, q, k, 1] <- su[i, q, k]
                for (r in 2:(nk2-1)) {
                    PTs12[i, q, k, r] <- pow(max(su[i, q, k] - knots2[r], 0), 3)
                    PTs22[i, q, k, r] <- (knots2[nk2 - 1] - knots2[r]) * 
                        pow(max(su[i, q, k] - knots2[nk2], 0), 3)
                    PTs32[i, q, k, r] <- (knots2[nk2] - knots2[r]) * 
                        pow(max(su[i, q, k] - knots2[nk2 - 1], 0), 3)
                    ZTs2[i, q, k, r] <- PTs12[i, q, k, r] + (PTs22[i, q, k, r] - 
                                                                 PTs32[i, q, k, r]) / (knots2[nk2] - knots2[nk2 - 1])
                }
                f3.s[i, q, k] <- beta3.1[1] + beta3.1[2]*W[i, 1] + beta3.1[3]*W[i, 2] + 
                    beta3.1[4]*W[i, 3] + b[i, nk1 + nk2 + 1] + 
                            inprod(beta3.2[1:(nk3-1)], ZTs3[i, q, k, 1:(nk3-1)]) + 
                    inprod(b[i, (nk1 + nk2 + 2):(nk1 + nk2 + nk3)], ZTs3[i, q, k, 1:(nk3-1)])
                ZTs3[i, q, k, 1] <- su[i, q, k]
                for (r in 2:(nk3-1)) {
                    PTs13[i, q, k, r] <- pow(max(su[i, q, k] - knots3[r], 0), 3)
                    PTs23[i, q, k, r] <- (knots3[nk3 - 1] - knots3[r]) * 
                        pow(max(su[i, q, k] - knots3[nk3], 0), 3)
                    PTs33[i, q, k, r] <- (knots3[nk3] - knots3[r]) * 
                        pow(max(su[i, q, k] - knots3[nk3 - 1], 0), 3)
                    ZTs3[i, q, k, r] <- PTs13[i, q, k, r] + 
                        (PTs23[i, q, k, r] - PTs33[i, q, k, r]) / (knots3[nk3] - knots3[nk3 - 1])
                }
                SurvLong[i, q, k] <- wk[k] * exp(alpha[1]*f1.s[i, q, k] + 
                                                     alpha[2]*f2.s[i, q, k] + 
                                                     alpha[3]*f3.s[i, q, k])
            }
            approxIntegral[i, q] <- P2[i, q] * sum(SurvLong[i, q, ])
            log.survival[i, q] <- - xi[q] * exp(etaBaseline[i]) * approxIntegral[i, q]
        }
        zeros[i] <- 0
        phi[i] <- - (sum(log.hazard[i, ]) + sum(log.survival[i, ])) + C
        zeros[i] ~ dpois(phi[i])
        # Random Effects Part
        latent1[i] ~ dcat(pi1[1:L1])
        for (l in 1:nb) {
            b[i, l] <- lambda1[latent1[i], l]
        }
        # membership
        for (k in 1:L1) {
            memb1[i, k] <- equals(latent1[i], k)
        }
    }
    # Priors
    # Longitudinal Part -- Outcome 1
    for (k in 1:ncX1) {
        beta1.1[k] ~ dnorm(0, 0.8)
    }
    for (k in 1:(nk1-1)) {
        beta1.2[k] ~ dnorm(0, 0.9)
    }
    tau1 ~ dgamma(1.01, 1.01)
    # Longitudinal Part -- Outcome 2
    for (k in 1:ncX2) {
        beta2.1[k] ~ dnorm(0, 0.8)
    }
    for (k in 1:(nk2-1)) {
        beta2.2[k] ~ dnorm(0, 0.9)
    }
    tau2 ~ dgamma(1.01, 1.01)
    # Longitudinal Part -- Outcome 3
    for (k in 1:ncX3) {
        beta3.1[k] ~ dnorm(0, 0.8)
    }
    for (k in 1:(nk3-1)) {
        beta3.2[k] ~ dnorm(0, 0.9)
    }
    # Survival Part
    for (k in 1:ncW) {
        gamma[k] ~ dnorm(0, 0.8)
    }
    for (q in 1:Q) {
        xi[q] ~ dgamma(0.2, 0.2)
    }
    for (k in 1:3) {
        alpha[k] ~ dnorm(0, 1.2)
    }
    # Random Effects Part
    for (k in 1:L1) {
        Tmemb1[k] <- sum(memb1[, k])
        Fmemb1[k] <- step(Tmemb1[k] - 1)
    }
    Tcluster1 <- sum(Fmemb1[])
    # constructive DP
    for (k in 1:L1) {
        r1[k] ~ dbeta(1, rho)
        lambda1[k, 1:nb] ~ dmnorm(nu0[], Taunu[1:nb, 1:nb])
    }
    # stick-breaking
    pi1[1] <- r1[1]
    for (k in 2:(L1-1)) {
        log(pi1[k]) <- log(r1[k]) + sum(R1[k, 1:k-1])
        for (l in 1:k-1) {
            R1[k, l] <- log(1 - r1[l])
        }
    }
    pi1[L1] <- 1 - sum(pi1[1:(L1-1)])
    for (k in 1:nb) {
        nu0[k] <- 0
    }
    Taunu[1:nb, 1:nb] ~ dwish(RR[, ], 25)
    ssigmab[1:nb, 1:nb] <- inverse(Taunu[, ])   
    rho ~ dgamma(2, 0.01)
}

namesI <- c(
    "N", "L1", "C", "nb",
    "ncX1", "ncX2", "ncX3",
    "nk1", "nk2", "nk3",
    "knots1", "knots2", "knots3",
    "offset1", "offset2", "offset3",
    "y1", "y2", "y3",    
    "X1", "X2", "X3",
    "Times1", "Times2", "Times3",
    "Time", "D",
    "Lo", "Up", "W", "ncW",
    "Q", "sk", "wk", "RR"
)

#########################################################################

modelIV <- function () {
    for (i in 1:N) {
        # Longitudinal Part
        # Outcome 1
        for (j in offset1[i]:(offset1[i + 1] - 1)){
            y1[j] ~ dnorm(mu1[j], tau1)            
            mu1[j] <- inprod(beta1.1[1:ncX1], X1[j, 1:ncX1]) + b[i, 1] + 
                inprod(beta1.2[1:(nk1-1)], Z1[j, 1:(nk1-1)]) + 
                inprod(b[i, 2:nk1], Z1[j, 1:(nk1-1)])
            Z1[j, 1] <- Times1[j]
            for (k in 2:(nk1-1)) {
                P11[j, k] <- pow(max(Times1[j] - knots1[k], 0), 3)
                P21[j, k] <- (knots1[nk1 - 1] - knots1[k]) * pow(max(Times1[j] - knots1[nk1], 0), 3)
                P31[j, k] <- (knots1[nk1] - knots1[k]) * pow(max(Times1[j] - knots1[nk1 - 1], 0), 3)
                Z1[j, k] <- P11[j, k] + (P21[j, k] - P31[j, k]) / (knots1[nk1] - knots1[nk1 - 1])
            }
        }
        # Outcome 2
        for (j in offset2[i]:(offset2[i + 1] - 1)){
            y2[j] ~ dnorm(mu2[j], tau2)
            mu2[j] <- inprod(beta2.1[1:ncX2], X2[j, 1:ncX2]) + b[i, nk1 + 1] + 
                inprod(beta2.2[1:(nk2-1)], Z2[j, 1:(nk2-1)]) + 
                inprod(b[i, (nk1 + 2):(nk1 + nk2)], Z2[j, 1:(nk2-1)])
            Z2[j, 1] <- Times2[j]
            for (k in 2:(nk2-1)) {
                P12[j, k] <- pow(max(Times2[j] - knots2[k], 0), 3)
                P22[j, k] <- (knots2[nk2 - 1] - knots2[k]) * pow(max(Times2[j] - knots2[nk2], 0), 3)
                P32[j, k] <- (knots2[nk2] - knots2[k]) * pow(max(Times2[j] - knots2[nk2 - 1], 0), 3)
                Z2[j, k] <- P12[j, k] + (P22[j, k] - P32[j, k]) / (knots2[nk2] - knots2[nk2 - 1])
            }
        }
        # Outcome 3
        for (j in offset3[i]:(offset3[i + 1] - 1)){
            y3[j] ~ dbern(pr[j])
            pr[j] <- max(0.0000000001, min(0.9999999999, P[j]))
            logit(P[j]) <- inprod(beta3.1[1:ncX3], X3[j, 1:ncX3]) + b[i, nk1 + nk2 + 1] + 
                inprod(beta3.2[1:(nk3-1)], Z3[j, 1:(nk3-1)]) + 
                inprod(b[i, (nk1 + nk2 + 2):(nk1 + nk2 + nk3)], Z3[j, 1:(nk3-1)])
            Z3[j, 1] <- Times3[j]
            for (k in 2:(nk3-1)) {
                P13[j, k] <- pow(max(Times3[j] - knots3[k], 0), 3)
                P23[j, k] <- (knots3[nk3 - 1] - knots3[k]) * pow(max(Times3[j] - knots3[nk3], 0), 3)
                P33[j, k] <- (knots3[nk3] - knots3[k]) * pow(max(Times3[j] - knots3[nk3 - 1], 0), 3)
                Z3[j, k] <- P13[j, k] + (P23[j, k] - P33[j, k]) / (knots3[nk3] - knots3[nk3 - 1])
            }
        }
        # Survival Part
        etaBaseline[i] <- inprod(gamma[], W[i, ])
        f1.T[i] <- beta1.1[1] + beta1.1[2]*W[i, 1] + beta1.1[3]*W[i, 2] + beta1.1[4]*W[i, 3] + b[i, 1] + 
            inprod(beta1.2[1:(nk1-1)], ZT1[i, 1:(nk1-1)]) + inprod(b[i, 2:nk1], ZT1[i, 1:(nk1-1)])
        ZT1[i, 1] <- Time[i]
        for (k in 2:(nk1-1)) {
            PT11[i, k] <- pow(max(Time[i] - knots1[k], 0), 3)
            PT21[i, k] <- (knots1[nk1 - 1] - knots1[k]) * pow(max(Time[i] - knots1[nk1], 0), 3)
            PT31[i, k] <- (knots1[nk1] - knots1[k]) * pow(max(Time[i] - knots1[nk1 - 1], 0), 3)
            ZT1[i, k] <- PT11[i, k] + (PT21[i, k] - PT31[i, k]) / (knots1[nk1] - knots1[nk1 - 1])
        }    
        f2.T[i] <- beta2.1[1] + beta2.1[2]*W[i, 1] + beta2.1[3]*W[i, 2] + beta2.1[4]*W[i, 3] + b[i, nk1 + 1] + 
            inprod(beta2.2[1:(nk2-1)], ZT2[i, 1:(nk2-1)]) + inprod(b[i, (nk1 + 2):(nk1 + nk2)], ZT2[i, 1:(nk2-1)])
        ZT2[i, 1] <- Time[i]
        for (k in 2:(nk2-1)) {
            PT12[i, k] <- pow(max(Time[i] - knots2[k], 0), 3)
            PT22[i, k] <- (knots2[nk2 - 1] - knots2[k]) * pow(max(Time[i] - knots2[nk2], 0), 3)
            PT32[i, k] <- (knots2[nk2] - knots2[k]) * pow(max(Time[i] - knots2[nk2 - 1], 0), 3)
            ZT2[i, k] <- PT12[i, k] + (PT22[i, k] - PT32[i, k]) / (knots2[nk2] - knots2[nk2 - 1])
        }
        f3.T[i] <- beta3.1[1] + beta3.1[2]*W[i, 1] + beta3.1[3]*W[i, 2] + beta3.1[4]*W[i, 3] + b[i, nk1 + nk2 + 1] + 
            inprod(beta3.2[1:(nk3-1)], ZT3[i, 1:(nk3-1)]) + 
            inprod(b[i, (nk1 + nk2 + 2):(nk1 + nk2 + nk3)], ZT3[i, 1:(nk3-1)])
        ZT3[i, 1] <- Time[i]
        for (k in 2:(nk3-1)) {
            PT13[i, k] <- pow(max(Time[i] - knots3[k], 0), 3)
            PT23[i, k] <- (knots3[nk3 - 1] - knots3[k]) * pow(max(Time[i] - knots3[nk3], 0), 3)
            PT33[i, k] <- (knots3[nk3] - knots3[k]) * pow(max(Time[i] - knots3[nk3 - 1], 0), 3)
            ZT3[i, k] <- PT13[i, k] + (PT23[i, k] - PT33[i, k]) / (knots3[nk3] - knots3[nk3 - 1])
        }
        for (q in 1:Q) {
            log.hazard[i, q] <- D[i, q] * (log(xi[q]) + etaBaseline[i] + 
                                               alpha[1]*f1.T[i] + alpha[2]*f2.T[i] + 
                                               alpha[3]*f3.T[i])
            P1[i, q] <- (Up[i, q] + Lo[i, q]) / 2
            P2[i, q] <- (Up[i, q] - Lo[i, q]) / 2
            for (k in 1:15) {
                su[i, q, k] <- P2[i, q] * sk[k] + P1[i, q]
                f1.s[i, q, k] <- beta1.1[1] + beta1.1[2]*W[i, 1] + beta1.1[3]*W[i, 2] + 
                    beta1.1[4]*W[i, 3] + b[i, 1] + 
                            inprod(beta1.2[1:(nk1-1)], ZTs1[i, q, k, 1:(nk1-1)]) + 
                    inprod(b[i, 2:nk1], ZTs1[i, q, k, 1:(nk1-1)])
                ZTs1[i, q, k, 1] <- su[i, q, k]
                for (r in 2:(nk1-1)) {
                    PTs11[i, q, k, r] <- pow(max(su[i, q, k] - knots1[r], 0), 3)
                    PTs21[i, q, k, r] <- (knots1[nk1 - 1] - knots1[r]) * 
                        pow(max(su[i, q, k] - knots1[nk1], 0), 3)
                    PTs31[i, q, k, r] <- (knots1[nk1] - knots1[r]) * 
                        pow(max(su[i, q, k] - knots1[nk1 - 1], 0), 3)
                    ZTs1[i, q, k, r] <- PTs11[i, q, k, r] + 
                        (PTs21[i, q, k, r] - PTs31[i, q, k, r]) / (knots1[nk1] - knots1[nk1 - 1])
                }
                f2.s[i, q, k] <- beta2.1[1] + beta2.1[2]*W[i, 1] + beta2.1[3]*W[i, 2] + 
                    beta2.1[4]*W[i, 3] + b[i, nk1 + 1] + 
                            inprod(beta2.2[1:(nk2-1)], ZTs2[i, q, k, 1:(nk2-1)]) + 
                    inprod(b[i, (nk1 + 2):(nk1 + nk2)], ZTs2[i, q, k, 1:(nk2-1)])
                ZTs2[i, q, k, 1] <- su[i, q, k]
                for (r in 2:(nk2-1)) {
                    PTs12[i, q, k, r] <- pow(max(su[i, q, k] - knots2[r], 0), 3)
                    PTs22[i, q, k, r] <- (knots2[nk2 - 1] - knots2[r]) * 
                        pow(max(su[i, q, k] - knots2[nk2], 0), 3)
                    PTs32[i, q, k, r] <- (knots2[nk2] - knots2[r]) * 
                        pow(max(su[i, q, k] - knots2[nk2 - 1], 0), 3)
                    ZTs2[i, q, k, r] <- PTs12[i, q, k, r] + 
                        (PTs22[i, q, k, r] - PTs32[i, q, k, r]) / (knots2[nk2] - knots2[nk2 - 1])
                }
                f3.s[i, q, k] <- beta3.1[1] + beta3.1[2]*W[i, 1] + beta3.1[3]*W[i, 2] + 
                    beta3.1[4]*W[i, 3] + b[i, nk1 + nk2 + 1] + 
                            inprod(beta3.2[1:(nk3-1)], ZTs3[i, q, k, 1:(nk3-1)]) + 
                    inprod(b[i, (nk1 + nk2 + 2):(nk1 + nk2 + nk3)], ZTs3[i, q, k, 1:(nk3-1)])
                ZTs3[i, q, k, 1] <- su[i, q, k]
                for (r in 2:(nk3-1)) {
                    PTs13[i, q, k, r] <- pow(max(su[i, q, k] - knots3[r], 0), 3)
                    PTs23[i, q, k, r] <- (knots3[nk3 - 1] - knots3[r]) * pow(max(su[i, q, k] - knots3[nk3], 0), 3)
                    PTs33[i, q, k, r] <- (knots3[nk3] - knots3[r]) * pow(max(su[i, q, k] - knots3[nk3 - 1], 0), 3)
                    ZTs3[i, q, k, r] <- PTs13[i, q, k, r] + 
                        (PTs23[i, q, k, r] - PTs33[i, q, k, r]) / (knots3[nk3] - knots3[nk3 - 1])
                }
                SurvLong[i, q, k] <- wk[k] * exp(alpha[1]*f1.s[i, q, k] + 
                                                     alpha[2]*f2.s[i, q, k] + 
                                                     alpha[3]*f3.s[i, q, k])
            }
            approxIntegral[i, q] <- P2[i, q] * sum(SurvLong[i, q, ])
            log.survival[i, q] <- - xi[q] * exp(etaBaseline[i]) * approxIntegral[i, q]
        }
        zeros[i] <- 0
        phi[i] <- - (sum(log.hazard[i, ]) + sum(log.survival[i, ])) + C
        zeros[i] ~ dpois(phi[i])
        # Random Effects Part
        latent1[i] ~ dcat(pi1[1:L1])
        for (l in 1:nb) {
            b[i, l] <- lambda1[latent1[i], l]
        }
        # membership
        for (k in 1:L1) {
            memb1[i, k] <- equals(latent1[i], k)
        }
    }
    # Priors
    # Longitudinal Part -- Outcome 1
    for (k in 1:ncX1) {
        beta1.1[k] ~ dnorm(0, 0.8)
    }
    for (k in 1:(nk1-1)) {
        beta1.2[k] ~ dnorm(0, 0.9)
    }
    knots1[1] <- p.knots1[1] * (Diff1 - (nk1 + 1) * gap) + (start1 + gap)
    for (r in 2:nk1) {
        knots1[r] <- p.knots1[r] * (Diff1 - (nk1 + 1) * gap) + (knots1[r-1] + gap)
    }
    p.knots1[1] <- q.knots1[1]
    q.knots1[1] ~ dbeta(a.knots1[1], b.knots1[1])
    for (r in 2:nk1) {
        b.knots1[r-1] <- a.knots1[r] + b.knots1[r] - 1
        p.knots1[r] <- q.knots1[r] * (1 - q.knots1[r-1]) * p.knots1[r-1] / q.knots1[r-1]
        q.knots1[r] ~ dbeta(a.knots1[r], b.knots1[r])
    }
    b.knots1[nk1] <- a.knots1[nk1 + 1]
    for (r in 1:(nk1+1)) {
        a.knots1[r] <- aeq.tau1
    }
    tau1 ~ dgamma(1.01, 1.01)
    # Longitudinal Part -- Outcome 2
    for (k in 1:ncX2) {
        beta2.1[k] ~ dnorm(0, 0.8)
    }
    for (k in 1:(nk2-1)) {
        beta2.2[k] ~ dnorm(0, 0.9)
    }
    knots2[1] <- p.knots2[1] * (Diff2 - (nk2 + 1) * gap) + (start2 + gap)
    for (r in 2:nk2) {
        knots2[r] <- p.knots2[r] * (Diff2 - (nk2 + 1) * gap) + (knots2[r-1] + gap)
    }
    p.knots2[1] <- q.knots2[1]
    q.knots2[1] ~ dbeta(a.knots2[1], b.knots2[1])
    for (r in 2:nk2) {
        b.knots2[r-1] <- a.knots2[r] + b.knots2[r] - 1
        p.knots2[r] <- q.knots2[r] * (1 - q.knots2[r-1]) * p.knots2[r-1] / q.knots2[r-1]
        q.knots2[r] ~ dbeta(a.knots2[r], b.knots2[r])
    }
    b.knots2[nk2] <- a.knots2[nk2 + 1]
    for (r in 1:(nk2+1)) {
        a.knots2[r] <- aeq.tau2
    }
    tau2 ~ dgamma(1.01, 1.01)
    # Longitudinal Part -- Outcome 3
    for (k in 1:ncX3) {
        beta3.1[k] ~ dnorm(0, 0.8)
    }
    for (k in 1:(nk3-1)) {
        beta3.2[k] ~ dnorm(0, 0.9)
    }
    knots3[1] <- p.knots3[1] * (Diff3 - (nk3 + 1) * gap) + (start3 + gap)
    for (r in 2:nk3) {
        knots3[r] <- p.knots3[r] * (Diff3 - (nk3 + 1) * gap) + (knots3[r-1] + gap)
    }
    p.knots3[1] <- q.knots3[1]
    q.knots3[1] ~ dbeta(a.knots3[1], b.knots3[1])
    for (r in 2:nk3) {
        b.knots3[r-1] <- a.knots3[r] + b.knots3[r] - 1
        p.knots3[r] <- q.knots3[r] * (1 - q.knots3[r-1]) * p.knots3[r-1] / q.knots3[r-1]
        q.knots3[r] ~ dbeta(a.knots3[r], b.knots3[r])
    }
    b.knots3[nk3] <- a.knots3[nk3 + 1]
    for (r in 1:(nk3+1)) {
        a.knots3[r] <- aeq.tau3
    }
    # Survival Part
    for (k in 1:ncW) {
        gamma[k] ~ dnorm(0, 0.8)
    }
    for (q in 1:Q) {
        xi[q] ~ dgamma(0.2, 0.2)
    }
    for (k in 1:3) {
        alpha[k] ~ dnorm(0, 1.2)
    }
    # Random Effects Part
    for (k in 1:L1) {
        Tmemb1[k] <- sum(memb1[, k])
        Fmemb1[k] <- step(Tmemb1[k] - 1)
    }
    Tcluster1 <- sum(Fmemb1[])
    # constructive DP
    for (k in 1:L1) {
        r1[k] ~ dbeta(1, rho)
        lambda1[k, 1:nb] ~ dmnorm(nu0[], Taunu[1:nb, 1:nb])
    }
    # stick-breaking
    pi1[1] <- r1[1]
    for (k in 2:(L1-1)) {
        log(pi1[k]) <- log(r1[k]) + sum(R1[k, 1:k-1])
        for (l in 1:k-1) {
            R1[k, l] <- log(1 - r1[l])
        }
    }
    pi1[L1] <- 1 - sum(pi1[1:(L1-1)])
    for (k in 1:nb) {
        nu0[k] <- 0
    }
    Taunu[1:nb, 1:nb] ~ dwish(RR[, ], 25)
    ssigmab[1:nb, 1:nb] <- inverse(Taunu[, ])   
    rho ~ dgamma(2, 0.01)
}

namesIV <- c(
    "N", "L1", "gap", "C", "nb",
    "ncX1", "ncX2", "ncX3",
    "nk1", "nk2", "nk3",
    "aeq.tau1", "aeq.tau2", "aeq.tau3",
    "Diff1", "Diff2", "Diff3",
    "start1", "start2", "start3",
    "offset1", "offset2", "offset3",
    "y1", "y2", "y3",    
    "X1", "X2", "X3",
    "Times1", "Times2", "Times3",
    "Time", "D",
    "Lo", "Up", "W", "ncW",
    "Q", "sk", "wk", "RR"
)

#########################################################################

modelVI <- function () {
    for (i in 1:N) {
        # Longitudinal Part
        # Outcome 1
        for (j in offset1[i]:(offset1[i + 1] - 1)){
            y1[j] ~ dnorm(mu1[j], tau1)            
            mu1[j] <- inprod(beta1.1[1:ncX1], X1[j, 1:ncX1]) + b[i, 1] + 
                inprod(beta1.2[1:(nk1-1)], Z1[j, 1:(nk1-1)]) + 
                inprod(b[i, 2:nk1], Z1[j, 1:(nk1-1)])
            Z1[j, 1] <- Times1[j]
            for (k in 2:(nk1-1)) {
                P11[j, k] <- pow(max(Times1[j] - knots1[k], 0), 3)
                P21[j, k] <- (knots1[nk1 - 1] - knots1[k]) * pow(max(Times1[j] - knots1[nk1], 0), 3)
                P31[j, k] <- (knots1[nk1] - knots1[k]) * pow(max(Times1[j] - knots1[nk1 - 1], 0), 3)
                Z1[j, k] <- P11[j, k] + (P21[j, k] - P31[j, k]) / (knots1[nk1] - knots1[nk1 - 1])
            }
        }
        # Outcome 2
        for (j in offset2[i]:(offset2[i + 1] - 1)){
            y2[j] ~ dnorm(mu2[j], tau2)
            mu2[j] <- inprod(beta2.1[1:ncX2], X2[j, 1:ncX2]) + b[i, nk1 + 1] + 
                inprod(beta2.2[1:(nk2-1)], Z2[j, 1:(nk2-1)]) + 
                inprod(b[i, (nk1 + 2):(nk1 + nk2)], Z2[j, 1:(nk2-1)])
            Z2[j, 1] <- Times2[j]
            for (k in 2:(nk2-1)) {
                P12[j, k] <- pow(max(Times2[j] - knots2[k], 0), 3)
                P22[j, k] <- (knots2[nk2 - 1] - knots2[k]) * pow(max(Times2[j] - knots2[nk2], 0), 3)
                P32[j, k] <- (knots2[nk2] - knots2[k]) * pow(max(Times2[j] - knots2[nk2 - 1], 0), 3)
                Z2[j, k] <- P12[j, k] + (P22[j, k] - P32[j, k]) / (knots2[nk2] - knots2[nk2 - 1])
            }
        }
        # Outcome 3
        for (j in offset3[i]:(offset3[i + 1] - 1)){
            y3[j] ~ dbern(pr[j])
            pr[j] <- max(0.0000000001, min(0.9999999999, P[j]))
            logit(P[j]) <- inprod(beta3.1[1:ncX3], X3[j, 1:ncX3]) + b[i, nk1 + nk2 + 1] + 
                inprod(beta3.2[1:(nk3-1)], Z3[j, 1:(nk3-1)]) + 
                inprod(b[i, (nk1 + nk2 + 2):(nk1 + nk2 + nk3)], Z3[j, 1:(nk3-1)])
            Z3[j, 1] <- Times3[j]
            for (k in 2:(nk3-1)) {
                P13[j, k] <- pow(max(Times3[j] - knots3[k], 0), 3)
                P23[j, k] <- (knots3[nk3 - 1] - knots3[k]) * pow(max(Times3[j] - knots3[nk3], 0), 3)
                P33[j, k] <- (knots3[nk3] - knots3[k]) * pow(max(Times3[j] - knots3[nk3 - 1], 0), 3)
                Z3[j, k] <- P13[j, k] + (P23[j, k] - P33[j, k]) / (knots3[nk3] - knots3[nk3 - 1])
            }
        }
        # Survival Part
        etaBaseline[i] <- inprod(gammas[], W[i, ]) + b[i, nk1 + nk2 + nk3 + 1]
        for (q in 1:Q) {
            mu[i, q] <- xi[q] * exp(etaBaseline[i]) * T[i, q]
            D[i, q] ~ dpois(mu[i, q])
        }
        # Random Effects Part
        latent1[i] ~ dcat(pi1[1:L1])
        for (l in 1:nb) {
            b[i, l] <- lambda1[latent1[i], l]
        }
        # membership
        for (k in 1:L1) {
            memb1[i, k] <- equals(latent1[i], k)
        }
    }
    # Priors
    # Longitudinal Part -- Outcome 1
    for (k in 1:ncX1) {
        beta1.1[k] ~ dnorm(0, 0.01)
    }
    for (k in 1:(nk1-1)) {
        beta1.2[k] ~ dnorm(0, 0.01)
    }
    knots1[1] <- p.knots1[1] * (Diff1 - (nk1 + 1) * gap) + (start1 + gap)
    for (r in 2:nk1) {
        knots1[r] <- p.knots1[r] * (Diff1 - (nk1 + 1) * gap) + (knots1[r-1] + gap)
    }
    p.knots1[1] <- q.knots1[1]
    q.knots1[1] ~ dbeta(a.knots1[1], b.knots1[1])
    for (r in 2:nk1) {
        b.knots1[r-1] <- a.knots1[r] + b.knots1[r] - 1
        p.knots1[r] <- q.knots1[r] * (1 - q.knots1[r-1]) * p.knots1[r-1] / q.knots1[r-1]
        q.knots1[r] ~ dbeta(a.knots1[r], b.knots1[r])
    }
    b.knots1[nk1] <- a.knots1[nk1 + 1]
    for (r in 1:(nk1+1)) {
        a.knots1[r] <- aeq.tau1
    }
    tau1 ~ dgamma(1.01, 1.01)
    # Longitudinal Part -- Outcome 2
    for (k in 1:ncX2) {
        beta2.1[k] ~ dnorm(0, 0.01)
    }
    for (k in 1:(nk2-1)) {
        beta2.2[k] ~ dnorm(0, 0.01)
    }
    knots2[1] <- p.knots2[1] * (Diff2 - (nk2 + 1) * gap) + (start2 + gap)
    for (r in 2:nk2) {
        knots2[r] <- p.knots2[r] * (Diff2 - (nk2 + 1) * gap) + (knots2[r-1] + gap)
    }
    p.knots2[1] <- q.knots2[1]
    q.knots2[1] ~ dbeta(a.knots2[1], b.knots2[1])
    for (r in 2:nk2) {
        b.knots2[r-1] <- a.knots2[r] + b.knots2[r] - 1
        p.knots2[r] <- q.knots2[r] * (1 - q.knots2[r-1]) * p.knots2[r-1] / q.knots2[r-1]
        q.knots2[r] ~ dbeta(a.knots2[r], b.knots2[r])
    }
    b.knots2[nk2] <- a.knots2[nk2 + 1]
    for (r in 1:(nk2+1)) {
        a.knots2[r] <- aeq.tau2
    }
    tau2 ~ dgamma(1.01, 1.01)
    # Longitudinal Part -- Outcome 3
    for (k in 1:ncX3) {
        beta3.1[k] ~ dnorm(0, 0.01)
    }
    for (k in 1:(nk3-1)) {
        beta3.2[k] ~ dnorm(0, 0.01)
    }
    knots3[1] <- p.knots3[1] * (Diff3 - (nk3 + 1) * gap) + (start3 + gap)
    for (r in 2:nk3) {
        knots3[r] <- p.knots3[r] * (Diff3 - (nk3 + 1) * gap) + (knots3[r-1] + gap)
    }
    p.knots3[1] <- q.knots3[1]
    q.knots3[1] ~ dbeta(a.knots3[1], b.knots3[1])
    for (r in 2:nk3) {
        b.knots3[r-1] <- a.knots3[r] + b.knots3[r] - 1
        p.knots3[r] <- q.knots3[r] * (1 - q.knots3[r-1]) * p.knots3[r-1] / q.knots3[r-1]
        q.knots3[r] ~ dbeta(a.knots3[r], b.knots3[r])
    }
    b.knots3[nk3] <- a.knots3[nk3 + 1]
    for (r in 1:(nk3+1)) {
        a.knots3[r] <- aeq.tau3
    }
    # Survival Part
    for (k in 1:ncW) {
        gammas[k] ~ dnorm(0, 2.5)
    }
    for (q in 1:Q) {
        xi[q] ~ dgamma(0.2, 0.2)
    }
    # Random Effects Part
    for (k in 1:L1) {
        Tmemb1[k] <- sum(memb1[, k])
        Fmemb1[k] <- step(Tmemb1[k] - 1)
    }
    Tcluster1 <- sum(Fmemb1[])
    # constructive DP
    for (k in 1:L1) {
        r1[k] ~ dbeta(1, rho)
        lambda1[k, 1:nb] ~ dmnorm(nu0[], Taunu[1:nb, 1:nb])
    }
    # stick-breaking
    pi1[1] <- r1[1]
    for (k in 2:(L1-1)) {
        log(pi1[k]) <- log(r1[k]) + sum(R1[k, 1:k-1])
        for (l in 1:k-1) {
            R1[k, l] <- log(1 - r1[l])
        }
    }
    pi1[L1] <- 1 - sum(pi1[1:(L1-1)])
    for (k in 1:nb) {
        nu0[k] <- 0
    }
    Taunu[1:nb, 1:nb] ~ dwish(RR[, ], 25)
    ssigmab[1:nb, 1:nb] <- inverse(Taunu[, ])   
    rho ~ dgamma(2, 0.01)
}

namesVI <- c(
    "N", "L1", "gap", "nb",
    "ncX1", "ncX2", "ncX3",
    "nk1", "nk2", "nk3",
    "aeq.tau1", "aeq.tau2", "aeq.tau3",
    "Diff1", "Diff2", "Diff3",
    "start1", "start2", "start3",
    "offset1", "offset2", "offset3",
    "y1", "y2", "y3",    
    "X1", "X2", "X3",
    "Times1", "Times2", "Times3",
    "T", "D",
    "W", "ncW", "Q", "RR"
)
