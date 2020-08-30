library(EDMAinR)
library(geomorph)
library(shapes) # Morpho, Momocs

## -------- functions for rotating/translating/etc -----------

# R: rotation, t: translation, alpha: scaling, E: perturbation
#
# Model 1: alpha_i (M+E_i)R_i+t_i. (scaling after perturbation)
# Model 2: (alpha_i*M + E_i)R_i + t_i (scaling before perturbation)
#
# alpha_i vector to be provided by the user.
# We will generate it from some distribution for simulations.

## rotate M around its centroid (i.e. not origin)
rotate_2d <- function(M, deg=0, center=TRUE) {
    ## degree to radian
    deg2rad <- function(deg) deg * pi /180
    ## make a rotation matrix based on degrees
    Rmat2d <- function(deg) {
        rad <- deg2rad(deg)
        matrix(c(cos(rad), sin(rad), -sin(rad), cos(rad)), 2, 2)
    }
    tr <- if (center)
      colMeans(M) else rep(0, ncol(M))
    Mpr <- t(t(M) - tr)
    out <- Mpr %*% Rmat2d(deg)
    t(t(out) + tr)
}

## scale M around its centroid (i.e. not origin)
scale_2d <- function(M, scale=1, center=TRUE) {
    tr <- if (center)
      colMeans(M) else rep(0, ncol(M))
    Mpr <- t(t(M) - tr)
    out <- Mpr * scale
    t(t(out) + tr)
}

## translate M
translate_2d <- function(M, tx=0, ty=0) {
    t(t(M) + c(tx, ty))
}

## simulate a single specimen given M and SigmaK
simulate1_2d <- function(M, SigmaK, deg=0, tx=0, ty=0, scale=1, before=FALSE) {
    if (before)
        M <- scale_2d(M, scale)
    S <- EDMAinR:::.edma_simulate_data(1, M, SigmaK)$A[,,1]
    N <- S
    if (!before)
        S <- scale_2d(S, scale)
    S <- rotate_2d(S, deg)
    S <- translate_2d(S, tx, ty)
    dimnames(S) <- list(paste0("L", seq_len(nrow(S))), c("X", "Y"))
    dimnames(N) <- dimnames(S)
    attr(S, "natural") <- N
    S
}

simulate_2d <- function(n, M, SigmaK, deg=0, tx=0, ty=0, scale=1, before=FALSE) {
    nn <- seq_len(n)
    deg <- rep(deg, n)[nn]
    tx <- rep(tx, n)[nn]
    ty <- rep(ty, n)[nn]
    scale <- rep(scale, n)[nn]
    d <- lapply(nn, function(i) {
        simulate1_2d(M, SigmaK,
                     deg=deg[i],
                     tx=tx[i],
                     ty=ty[i],
                     scale=scale[i],
                     before=before)
    })
    names(d) <- paste0("S", nn)
    N <- d
    for (i in nn) {
        N[[i]] <- attr(d[[i]], "natural")
        attr(d[[i]], "natural") <- NULL
    }
    out <- list(name="Simulated data", data=d)
    class(out) <- c("edma_data_simul", "edma_data")
    nat <- out
    nat$data <- N
    out$natural <- as.array(nat)
    attr(out, "M") <- M
    attr(out, "SigmaK") <- SigmaK
    out
}

## plot the simulated objects
plot.edma_data_simul <- function(x, nmax=NULL, natural=FALSE,
chull=FALSE, ellipse=FALSE, ...) {
    A <- if (natural)
        x$natural else as.array(x)
    n <- dim(x)[3]
    if (is.null(nmax))
      nmax <- n
    AA <- rbind(x$M, do.call(rbind, lapply(1:n, function(i) A[,,i])))
    plot(AA, type="n", asp=1, axes=FALSE, ann=FALSE, ...)
    for (i in seq_len(nmax)) {
        polygon(A[,,i], border="#00000088", col=NA)
    }
    for (i in seq_len(dim(A)[1])) {
        AAA <- t(A[i,,])
        if (ellipse)
            polygon(EDMAinR:::.data_ellipse(AAA),
                col="#ff000044", border="#ff0000")
        if (chull)
            polygon(AAA[chull(AAA),],
                col="#ff000044", border="#ff0000")
    }
    invisible(x)
}

## ---------- simulation -----------

## Mean form
M1 <- rbind(
    L1=c(2, 0),
    L2=c(0, 2),
    L3=c(-2, 0),
    L4=c(-1, -5),
    L5=c(0, -6),
    L6=c(1, -5)
)
colnames(M1) <- c("X", "Y")
M1 <- 10 * M1

## SigmaK
S1 <- matrix(
  c("s1", NA,  NA, NA,  NA, NA,
    NA, "s1", NA, NA,  NA, NA,
    NA,  NA, "s1", NA,  NA, NA,
    NA,  NA,  NA, "s2", NA, NA,
    NA,  NA,  NA, NA, "s2", NA,
    NA,  NA,  NA, NA, NA, "s2"),
  nrow=6, ncol=6, byrow=TRUE)
dimnames(S1) <- list(rownames(M1), rownames(M1))
parm1 <- c("s1"=1, "s2"=5)
SigmaK1 <- make_Sigma(parm1, S1)

## simulation settings
set.seed(23)
n <- 100 # number of specimens
range <- 100
#alpha <- rlnorm(n, 0, 0.1) # set this 1 to remove scaling
before <- FALSE
alpha <- 1
deg <- runif(n, 0, 360)
tx <- runif(n, -range, range)
ty <- runif(n, -range, range)

## simulated object
sim <- simulate_2d(n, M1, SigmaK1,
    deg=deg, tx=tx, ty=ty, scale=alpha, before=before)
## has the following elements:
## $data:     EDMA data list after sclaing/perturbation
## $natural:  data as 3D array in natural space (before rotation/translation)

op <- par(mfrow=c(1,2))
plot(sim, nmax=5, natural=TRUE, chull=FALSE, ellipse=TRUE)
title("Natural space")
plot(sim, nmax=5)
polygon(M1, border=2, col=NA)
title("Observations")
par(op)

## here is how to get mean/var of the pairwise Eu distances
z <- EDMAinR:::.edma_fit_np(sim, less=FALSE)
str(z)
z$EuMean
z$EuVar

## EDMA
#sim <- edma_simulate_data(n, M1, SigmaK1)
fit <- edma_fit(sim)
e <- SigmaK_fit(fit, S1)
Meanform(e)
SigmaKstar(e)
SigmaK(e)

## Procrustes
A <- as.array(sim)
p <- gpagen(A)

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

## check if we get the same as they
plot(as.numeric(VCVs), as.numeric(p$points.VCV))
abline(0,1)

## sigma hat
diag(SigmaK(e))
matrix(diag(VCV), ncol=2, byrow=TRUE)


## using shapes

p <- procGPA(A, scale=FALSE)
plot(p$mshape,asp=1)
#plotshapes(A)




