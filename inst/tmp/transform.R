# Classification using composite likelihood ratio and computing the reliability

# Simulate the data
library(EDMAinR)

## data sets

## Crouzon
x1 <- read_xyz(system.file("extdata/crouzon/Crouzon_P0_Global_MUT.xyz", package="EDMAinR"))
x2 <- read_xyz(system.file("extdata/crouzon/Crouzon_P0_Global_NON-MUT.xyz", package="EDMAinR"))
all(landmarks(x1) == landmarks(x2))

## Apert
x1 <- read_xyz(system.file("extdata/apert/ALLAPT4.xyz", package="EDMAinR"))
x2 <- read_xyz(system.file("extdata/apert/ALLNORM4.xyz", package="EDMAinR"))
all(landmarks(x1) == landmarks(x2))

## Growth (another Crouzon data)
x1 <- read_xyz(system.file("extdata/growth/CZP0_mut_global.xyz", package="EDMAinR"))
x2 <- read_xyz(system.file("extdata/growth/CZP0_wt_global.xyz", package="EDMAinR"))
all(landmarks(x1) == landmarks(x2))


## Purple book: Sagittal synostosis calvarial data
x1 <- read_xyz(system.file("extdata/purplebook/craniosynostosiscalvarium.xyz", package="EDMAinR"))
x2 <- read_xyz(system.file("extdata/purplebook/normalcalvarium.xyz", package="EDMAinR"))
all(landmarks(x1) == landmarks(x2))

## Purple book: left hemi-mandible of trisomic mice
x1 <- read_xyz(system.file("extdata/purplebook/normalmandibles.xyz", package="EDMAinR"))
x2 <- read_xyz(system.file("extdata/purplebook/trisomicmandibles.xyz", package="EDMAinR"))
all(landmarks(x1) == landmarks(x2))

M <- cbind(c(1,2,3), c(0,1,0), c(0,0,0))

# R: rotation, t: translation, alpha: scaling, E: perturbation
#
# Model 1: alpha_i (M+E_i)R_i+t_i. (scaling after perturbation)
# Model 2: (alpha_i*M + E_i)R_i + t_i (scaling before perturbation)

transform_M <- function(M, tr=0, sc=1, rt=NA, center=TRUE) {
    D <- ncol(M)
    cs <- colMeans(M, na.rm=TRUE)
    if (center)
        M <- t(t(M) - cs)
    tr <- rep(tr, D)[seq_len(D)]
    sc <- rep(sc, D)[seq_len(D)]
    rt <- rep(rt, D)[seq_len(D)]
    rt <- rt * pi /180

    out <- cbind(M, 1)

    ## rotate
    if (!is.na(rt[1L])) {
        Xx <- diag(1, D+1L, D+1L)
        Xx[2L,2L] <- cos(rt[1L])
        Xx[3L,2L] <- sin(rt[1L])
        Xx[2L,3L] <- -sin(rt[1L])
        Xx[3L,3L] <- cos(rt[1L])
        out <- t((Xx %*% t(out)))
    }
    if (!is.na(rt[2L])) {
        Xy <- diag(1, D+1L, D+1L)
        Xy[1L,1L] <- cos(rt[2L])
        Xy[3L,3L] <- -sin(rt[2L])
        Xy[1L,1L] <- sin(rt[2L])
        Xy[3L,3L] <- cos(rt[2L])
        out <- t((Xy %*% t(out)))
    }
    if (!is.na(rt[3L])) {
        Xz <- diag(1, D+1L, D+1L)
        Xz[1L,1L] <- cos(rt[3L])
        Xz[2L,1L] <- sin(rt[3L])
        Xz[1L,2L] <- -sin(rt[3L])
        Xz[2L,2L] <- cos(rt[3L])
        out <- t((Xz %*% t(out)))
    }

    # scale
    Xs <- diag(c(sc, 1))
    out <- t((Xs %*% t(out)))

    ## translate
    Xt <- diag(D + 1)
    Xt[seq_len(D), D+1L] <- tr
    out <- t((Xt %*% t(out)))

    out <- out[,seq_len(D)]

    if (center)
        out <- t(t(out) + cs)
    out
}
plot(0, xlim=c(-5,5), ylim=c(-5,5), type="n")
polygon(M)
polygon(transform_M(M, sc=2, tr=-1, rt=c(45, NA, NA)), border=2)

