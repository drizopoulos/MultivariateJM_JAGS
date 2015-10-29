haem <- read.csv(paste(root, "Data/haem.csv", sep = ""))
haem <- haem[c("id", "haematocrit", "yearse", "gender", "weight", "age", "failure", 
               "Time")]
names(haem) <- c("id", "haematocrit", "yearseHAEM", "gender", "weight", "age", "failure", 
                 "Time")
haem$haematocrit <- 100 * haem$haematocrit
#
prot <- read.csv(paste(root, "Data/PROT-New.csv", sep = ""))
prot$proteinuriaN <- as.numeric(prot$proteinuria >= 1)
prot$proteinuria <- factor(prot$proteinuria >= 1, levels = c(TRUE, FALSE), 
                           labels = c("prot", "no prot"))
prot <- prot[c("PATID", "proteinuria", "proteinuriaN", "yearse", "failure", "followyear")]
names(prot) <- c("id", "proteinuria", "proteinuriaN", "yearsePROT", "failure", "Time")
#
gfr <- read.csv(paste(root, "Data/GFR-New.csv", sep = ""))
names(gfr) <- c("id", "failure", "failure10", "yearseGFR", "gfr", "Time")
gfr <- gfr[c("id", "gfr", "yearseGFR", "failure", "Time")]


ind.gfr <- gfr$id %in% unique(haem$id)
gfr <- gfr[ind.gfr, ]
gfr <- gfr[order(gfr$id, gfr$yearseGFR), ]
ind.prot <- prot$id %in% unique(haem$id)
prot <- prot[ind.prot, ]
prot <- prot[order(prot$id, prot$yearsePROT), ]
haem <- haem[order(haem$id, haem$yearseHAEM), ]

# weight
sp <- sapply(split(haem$weight, haem$id), "[", 1)
prot$weight <- rep(sp, tapply(prot$id, prot$id, length))
gfr$weight <- rep(sp, tapply(gfr$id, gfr$id, length))
# age
sp <- sapply(split(haem$age, haem$id), "[", 1)
prot$age <- rep(sp, tapply(prot$id, prot$id, length))
gfr$age <- rep(sp, tapply(gfr$id, gfr$id, length))
# gender
sp <- sapply(split(haem$gender, haem$id), "[", 1)
prot$gender <- rep(sp, tapply(prot$id, prot$id, length))
gfr$gender <- rep(sp, tapply(gfr$id, gfr$id, length))
# Time
sp <- sapply(split(gfr$Time, gfr$id), "[", 1)
haem$Time <- rep(sp, tapply(haem$id, haem$id, length))

gfr$id <- factor(gfr$id)
haem$id <- factor(haem$id)
prot$id <- factor(prot$id)
gfr <- gfr[!is.na(gfr$gfr), ]
haem <- haem[!is.na(haem$haematocrit), ]
prot <- prot[!is.na(prot$proteinuria), ]

sp <- split(haem, haem$id)
lsp <- lapply(sp, function (x) x[x$yearseHAEM < x$Time, ]) 
haem <- do.call("rbind", lsp)
sp <- split(gfr, gfr$id)
lsp <- lapply(sp, function (x) x[x$yearseGFR < x$Time, ]) 
gfr <- do.call("rbind", lsp)
sp <- split(prot, prot$id)
lsp <- lapply(sp, function (x) x[x$yearsePROT < x$Time, ]) 
prot <- do.call("rbind", lsp)

row.names(prot) <- seq_len(nrow(prot))
row.names(haem) <- seq_len(nrow(haem))
row.names(gfr) <- seq_len(nrow(gfr))

surv <- data.frame(
    id = names(tapply(gfr$id, gfr$id, "[", 1)),
    Time = tapply(gfr$Time, gfr$id, "[", 1),
    failure = tapply(gfr$failure, gfr$id, "[", 1),
    weight = tapply(gfr$weight, gfr$id, "[", 1),
    age = tapply(gfr$age, gfr$id, "[", 1),
    gender = factor(tapply(gfr$gender, gfr$id, "[", 1), levels = 1:2, 
                    labels = c("female", "male"))
)

rm(ind.gfr, ind.prot, lsp, sp)

# Gauss-Kronrod
sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, 
        -0.405845151377397166906606412076961, 0, 0.405845151377397166906606412076961, 
        0.741531185599394439863864773280788, 0.949107912342758524526189684047851, 
        -0.991455371120812639206854697526329, -0.864864423359769072789712788640926, 
        -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 
        0.207784955007898467600689403773245, 0.586087235467691130294144838258730, 
        0.864864423359769072789712788640926, 0.991455371120812639206854697526329)

wk <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 
        0.190350578064785409913256402421014, 0.209482141084727828012999174891714, 
        0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 
        0.063092092629978553290700663189204, 0.022935322010529224963732008058970, 
        0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 
        0.204432940075298892414161999234649, 0.204432940075298892414161999234649, 
        0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 
        0.022935322010529224963732008058970)

###########################################################################################################
###########################################################################################################
###########################################################################################################

library("splines")
N <- length(unique(gfr$id))
K <- 15
Q <- 4
nk1 <- 3
nk2 <- 3
nk3 <- 3
offset1 <- with(gfr, as.vector(c(1, 1 + cumsum(tapply(id, id, length)))))
offset2 <- with(haem, as.vector(c(1, 1 + cumsum(tapply(id, id, length)))))
offset3 <- with(prot, as.vector(c(1, 1 + cumsum(tapply(id, id, length)))))
y1 <- gfr$gfr
y2 <- haem$haematocrit
y3 <- prot$proteinuriaN
#
Time <- as.vector(surv$Time)
event <- as.vector(surv$failure)
W <- cbind(surv$age, 1*(surv$gender == "male"), surv$weight); ncW <- ncol(W)
qs <- quantile(Time[event == 1], seq(0, 1, len = Q + 1), names = FALSE)[-c(1, Q + 1)]
qs <- c(0, qs, max(Time) + 1)
ind <- findInterval(Time, qs, rightmost.closed = TRUE)
D <- matrix(0, length(ind), Q)
D[cbind(seq_along(ind), ind)] <- 1
D <- D * event
Tiq <- outer(Time, qs, pmin)
Lo <- Tiq[, 1:Q]
Up <- Tiq[, 2:(Q+1)]
T <- Up - Lo
P <- T / 2
P1 <- (Up + Lo) / 2
st <- matrix(0, N, K*Q)
skQ <- rep(sk, Q)
for (i in 1:N) {
    st[i, ] <- rep(P[i, ], each = K) * skQ + rep(P1[i, ], each = K)
}
#
Times1 <- gfr$yearseGFR
Times2 <- haem$yearseHAEM
Times3 <- prot$yearsePROT
NS1 <- ns(Times1, nk1)
NS2 <- ns(Times2, nk1)
NS3 <- ns(Times3, nk1)
X1 <- cbind(1, gfr$age, 1*(gfr$gender == "male"), gfr$weight, NS1); ncX1 <- ncol(X1)
X2 <- cbind(1, haem$age, 1*(haem$gender == "male"), haem$weight, NS2); ncX2 <- ncol(X2)
X3 <- cbind(1, prot$age, 1*(prot$gender == "male"), prot$weight, NS3); ncX3 <- ncol(X3)
Z1 <- cbind(1, NS1); ncZ1 <- ncol(Z1)
Z2 <- cbind(1, NS2); ncZ2 <- ncol(Z2)
Z3 <- cbind(1, NS3); ncZ3 <- ncol(Z3)
NSt1 <- ns(Time, knots = attr(NS1, "knots"), Boundary.knots = attr(NS1, "Boundary.knots"))
NSt2 <- ns(Time, knots = attr(NS2, "knots"), Boundary.knots = attr(NS2, "Boundary.knots"))
NSt3 <- ns(Time, knots = attr(NS3, "knots"), Boundary.knots = attr(NS3, "Boundary.knots"))
X1t <- cbind(1, surv$age, 1*(surv$gender == "male"), surv$weight, NSt1)
X2t <- cbind(1, surv$age, 1*(surv$gender == "male"), surv$weight, NSt2)
X3t <- cbind(1, surv$age, 1*(surv$gender == "male"), surv$weight, NSt3)
Z1t <- cbind(1, NSt1)
Z2t <- cbind(1, NSt2)
Z3t <- cbind(1, NSt3)

NSs1 <- ns(c(t(st)), knots = attr(NS1, "knots"), Boundary.knots = attr(NS1, "Boundary.knots"))
NSs2 <- ns(c(t(st)), knots = attr(NS2, "knots"), Boundary.knots = attr(NS2, "Boundary.knots"))
NSs3 <- ns(c(t(st)), knots = attr(NS3, "knots"), Boundary.knots = attr(NS3, "Boundary.knots"))
X1s <- cbind(1, surv$age, 1*(surv$gender == "male"), surv$weight, NSs1)
X2s <- cbind(1, surv$age, 1*(surv$gender == "male"), surv$weight, NSs2)
X3s <- cbind(1, surv$age, 1*(surv$gender == "male"), surv$weight, NSs3)
Z1s <- cbind(1, NSs1)
Z2s <- cbind(1, NSs2)
Z3s <- cbind(1, NSs3)

knots1 <- quantile(Times1, seq(0.05, 0.95, length = nk1))
knots2 <- quantile(Times2, seq(0.05, 0.95, length = nk2))
knots3 <- quantile(Times3, seq(0.05, 0.95, length = nk3))

Data <- list(
    qs = qs,
    N = N, R = 8, C = 5000, Q = Q, gap = 2,
    ncX1 = ncX1, ncX2 = ncX2, ncX3 = ncX3,
    ncZ1 = ncZ1, ncZ2 = ncZ2, ncZ3 = ncZ3,
    nk1 = nk1, nk2 = nk2, nk3 = nk3,
    aeq.tau1 = 2, aeq.tau2 = 2, aeq.tau3 = 2,
    offset1 = offset1, offset2 = offset2, offset3 = offset3,
    y1 = y1, y2 = y2, y3 = y3,
    X1 = X1, X2 = X2, X3 = X3,
    X1t = X1t, X2t = X2t, X3t = X3t,
    X1s = X1s, X2s = X2s, X3s = X3s,
    Z1 = Z1, Z2 = Z2, Z3 = Z3,
    Z1t = Z1t, Z2t = Z2t, Z3t = Z3t,
    Z1s = Z1s, Z2s = Z2s, Z3s = Z3s,
    Times1 = Times1, Times2 = Times2, Times3 = Times3,
    knots1 = knots1, knots2 = knots2, knots3 = knots3,
    Diff1 = diff(range(Times1)), Diff2 = diff(range(Times2)), Diff3 = diff(range(Times3)),
    start1 = min(Times1), start2 = min(Times2), start3 = min(Times3),
    Time = Time, T = T, D = D, Lo = Lo, Up = Up, P = P, zeros = rep(0, N),
    W = W, ncW = ncW, K = 15,
    sk = sk, wk = wk,
    priorMean.betas1 = rep(0, ncX1), 
    priorTau.betas1 = diag(rep(0.001, ncX1)), priorA.tau1 = 0.1, priorB.tau1 = 0.1,
    priorMean.betas2 = rep(0, ncX2), 
    priorTau.betas2 = diag(rep(0.001, ncX2)), priorA.tau2 = 0.1, priorB.tau2 = 0.1,
    priorMean.betas3 = rep(0, ncX2), 
    priorTau.betas3 = diag(rep(0.001, ncX3)),
    priorMean.gammas = rep(0, ncW),
    priorTau.gammas = diag(rep(0.01, ncW)),
    priorMean.alphas = rep(0, 3), priorTau.alphas = diag(rep(0.001, 3)),
    priorA.xi = rep(0.1, Q), priorB.xi = rep(0.1, Q), priorK.D = 15,
    priorA.rho = 2, priorB.rho = 0.1
)


grb <- ls()
grb <- grb[-match(c("Data", "root"), grb)]
rm(list = c(grb, "grb"))














gc()
