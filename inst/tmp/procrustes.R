library(EDMAinR)
library(geomorph)

## -------- functions for rotating/translating -----------

## degree to radian
deg2rad <- function(deg) deg * pi /180

## make a rotation matrix based on degrees
Rmat2d <- function(deg) {
    rad <- deg2rad(deg)
    matrix(c(cos(rad), sin(rad), -sin(rad), cos(rad)), 2, 2)
}

## rotate M around its centroid (i.e. not origin)
rotate_2d <- function(M, deg=0, center=TRUE) {
    tr <- colMeans(M)
    Mpr <- t(t(M) - tr)
    out <- Mpr %*% Rmat2d(deg)
    t(t(out) + tr)
}

## translate M
translate_2d <- function(M, tx=0, ty=0) {
    t(t(M) + c(tx, ty))
}

## rotate and translate M
rot_trans_2d <- function(M, deg=0, tx=0, ty=0) {
    translate_2d(rotate_2d(M, deg), tx, ty)
}

## randomly rotate and translate M within bounds of +/- range
rnd_rot_trans_2d <- function(M, range=10) {
    rot_trans_2d(M,
        runif(1, 0, 360),
        runif(1, -range, range),
        runif(1, -range, range))
}

## ---------- simulation -----------

## Mean form
M <- rbind(
    L1=c(2, 0),
    L2=c(0, 2),
    L3=c(-2, 0),
    L4=c(-1, -5),
    L5=c(0, -6),
    L6=c(1, -5)
)
colnames(M) <- c("X", "Y")
M <- 10 * M

if (FALSE) {
set.seed(34)
plot(0, type="n", asp=1, xlim=c(-10, 10), ylim=c(-10, 10))
polygon(M)
for (i in 2:4)
    polygon(rnd_rot_trans_2d(M), border=i)
}

## SigmaK
S1 <- matrix(
  c("s1", NA,  NA, NA,  NA, NA,
    NA, "s1", NA, NA,  NA, NA,
    NA,  NA, "s1", NA,  NA, NA,
    NA,  NA,  NA, "s2", NA, NA,
    NA,  NA,  NA, NA, "s2", NA,
    NA,  NA,  NA, NA, NA, "s2"),
  nrow=6, ncol=6, byrow=TRUE)
dimnames(S1) <- list(rownames(M), rownames(M))
parm1 <- c("s1"=0.2, "s2"=4)
SigmaK <- EDMAinR:::.vec2mat(parm1, EDMAinR:::.mat2fac(S1))
dimnames(SigmaK) <- dimnames(S1)

## number of specimens
n <- 200

## simulate: natural space
sim0 <- edma_simulate_data(n=n, M, SigmaK)

## make an array for Procrustes
## rotate/shift
A <- as.array(sim0)
for (i in 1:n) {
    A[,,i] <- rnd_rot_trans_2d(A[,,i], range=100)
}
sim <- as.edma_data(A)

## quick check
plot_2d(sim)

## -------- natural space -----------

plot(as.matrix(sim0), type="n", asp=1, axes=FALSE, ann=FALSE)
#points(as.matrix(sim0), col="#00000044", pch=".")
for (i in 1:10) {
  polygon(sim0$data[[i]], border="#00000088", col=NA)
}
for (i in 1:nrow(M)) {
  polygon(EDMAinR:::.data_ellipse(t(as.array(sim0)[i,,])),
          col="#ff000044", border="#ff0000")
}
points(M, col=2, pch=3, cex=2)
#polygon(M, border=2, col=NA, lwd=2)

## --------- observed specimens ----------

plot(as.matrix(sim), type="n", asp=1, axes=FALSE, ann=FALSE)
#points(as.matrix(sim), col="#00000044", pch=".")
for (i in 1:10) {
  polygon(sim$data[[i]], border="#00000088", col=NA)
}
for (i in 1:nrow(M)) {
  polygon(EDMAinR:::.data_ellipse(t(as.array(sim0)[i,,])),
          col="#ff000044", border="#ff0000")
}
points(M, col=2, pch=3, cex=2)
polygon(M, border=2, col=NA)

## EDMA
fit <- edma_fit(sim)
e <- SigmaK_fit(fit, S1)
SigmaKstar(e)
SigmaK(e)
Meanform(e)

## Procrustes
p <- gpagen(as.array(sim))

## compare FM
## scaling estimate, need non-scaled Mean form estimate

## M hat
Me <- Meanform(e)
Mp <- p$consensus

## calculate unscaled residuals (RRR)
AA <- A
AA[] <- 0
RR <- AA
for (i in 1:n) {
  AA[,,i] <- p$coords[,,i] * p$Csize[i]
  RR[,,i] <- AA[,,i] - p$consensus * p$Csize[i]
}
AAA <- do.call(rbind, lapply(1:n, function(i) AA[,,i]))
RRR <- t(sapply(1:n, function(i) as.numeric(t(RR[,,i]))))

## plot
op <- par(mfrow=c(1,2))
plot(p, mean=FALSE)
plot(AAA, asp=1)
par(op)

## scaled (VSVs) and unscaled (VCV) covariances
VCVs <- cov(t(t(RRR) / p$Csize))
VCV <- cov(RRR)

## chack if we get the same as they
plot(as.numeric(VCVs), as.numeric(p$points.VCV))
abline(0,1)

## sigma hat
diag(SigmaK(e))
matrix(diag(VCV), ncol=2, byrow=TRUE)




D <- ncol(M)
K <- nrow(M)

system.time(v <- EDMAinR:::.edma_fit_np_old(stack(sim), n, K, D))
system.time(z <- EDMAinR:::.edma_fit_np(as.array(sim)))
v$M
z$M
z$M-v$M
max(abs(z$M-v$M))

v$SigmaKstar
z$SigmaKstar
z$SigmaKstar-v$SigmaKstar
max(abs(z$SigmaKstar-v$SigmaKstar))


library(magrittr)
library(shapes)
library(EDMAinR)

e <- edma_fit(sim) %>%
  SigmaK_fit(S1)

system.time(x1 <- .edma_simulate_data_old(n=n, M, SigmaK))
system.time(x2 <- .edma_simulate_data(n=n, M, SigmaK))
