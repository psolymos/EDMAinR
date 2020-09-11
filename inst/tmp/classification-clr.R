# Classification using composite likelihood ratio and computing the reliability

# Simulate the data
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
n <- 10000 # number of specimens
range <- 100
# alpha <- rlnorm(n, 0, 0.1) # set this 1 to remove scaling
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


# Composite likelihood for a new observation

sim.new <- simulate_2d(1, M1, SigmaK1,
    deg=deg, tx=tx, ty=ty, scale=alpha, before=before)

eu_sim.new = as.vector(dist(sim.new))^2

# Now we compute the log-density for each component of the above vector. Each component is a non-central Chi-squared distribution with different parameters. We will compute the non-centrality parameter and scaling parameter for each first.

# Check if the distribution derivation is correct

plot(density(z$EuX[1,2,]))     #This is the empirical density of the observations.

# These are the parameters of the non-central chi-squared distribution.

delta_mat = as.matrix(dist(Meanform(e)))^2
phi_mat = (1/ncol(M1))*(z$EuMean - delta_mat)

# Generate data from the non-central chi-squre distribution and plot the density.

ncp12 = delta_mat[1,3]/phi_mat[1,3]
try = rchisq(10000,2,ncp=ncp12)

plot(density(try),xlim=range(try),ylim=c(0,0.015))
par(new=T)
try2 = z$EuX[1,3,]/phi_mat[1,3]
plot(density(try2),xlim=range(try),ylim=c(0,0.015))

# These match so the result about non-central chisquare distribution is correct.

# Now we can write the composite likelihood function.

scaled_obs = z$EuX   # Initiate the matrix

for (i in 1:n){
	scaled_obs[,,i] = z$EuX[,,i]/phi_mat
}

ncp_mat = delta_mat/phi_mat  # This is the matrix of non-centrality parameters. Only the upper diagonal is required.

CL1 = 0
for (i in 1:(nrow(M1)-1)){
	for (j in (i+1):nrow(M1))
	CL1 = CL1 + sum(dchisq(scaled_obs[i,j,],ncol(M1),ncp = ncp_mat[i,j],log=T))
}

# Now we have all the components for computing the CL based classification rule.

# Now generate data from two different populations (different means and/or variances) and see if we can do classification.

# Population 1 data

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
n1 <- 100 # number of specimens
range <- 100
# alpha <- rlnorm(n, 0, 0.1) # set this 1 to remove scaling
before <- FALSE
alpha <- 1
deg <- runif(n1, 0, 360)
tx <- runif(n1, -range, range)
ty <- runif(n1, -range, range)

## simulated object
sim1 <- simulate_2d(n1, M1, SigmaK1,
    deg=deg, tx=tx, ty=ty, scale=alpha, before=before)
## has the following elements:
## $data:     EDMA data list after sclaing/perturbation
## $natural:  data as 3D array in natural space (before rotation/translation)

## here is how to get mean/var of the pairwise Eu distances
z1 <- EDMAinR:::.edma_fit_np(sim1, less=FALSE)


## EDMA
#sim <- edma_simulate_data(n, M1, SigmaK1)
fit1 <- edma_fit(sim1)
e1 <- SigmaK_fit(fit1, S1)


# Estimate the parameters for the non-central chi-square

delta_mat1 = as.matrix(dist(Meanform(e1)))^2
phi_mat1 = (1/ncol(M1))*(z1$EuMean - delta_mat1)


# Data from population 2 (Change the mean form and sigma)


## Mean form
M2 <- rbind(
    L1=c(2, 0),
    L2=c(0, 2),
    L3=c(-2, 0),
    L4=c(-1, -5),
    L5=c(0, -6),
    L6=c(1, 5)
)
colnames(M1) <- c("X", "Y")
M2 <- 10 * M1

## SigmaK
S2 <- matrix(
  c("s1", NA,  NA, NA,  NA, NA,
    NA, "s1", NA, NA,  NA, NA,
    NA,  NA, "s1", NA,  NA, NA,
    NA,  NA,  NA, "s1", NA, NA,
    NA,  NA,  NA, NA, "s1", NA,
    NA,  NA,  NA, NA, NA, "s2"),
  nrow=6, ncol=6, byrow=TRUE)
dimnames(S2) <- list(rownames(M2), rownames(M2))
parm1 <- c("s1"=1, "s2"=5)
SigmaK2 <- make_Sigma(parm1, S2)

## simulation settings
set.seed(23)
n2 <- 200 # number of specimens
range <- 100
# alpha <- rlnorm(n, 0, 0.1) # set this 1 to remove scaling
before <- FALSE
alpha <- 1
deg <- runif(n1, 0, 360)
tx <- runif(n1, -range, range)
ty <- runif(n1, -range, range)

## simulated object
sim2 <- simulate_2d(n2, M2, SigmaK2,
    deg=deg, tx=tx, ty=ty, scale=alpha, before=before)
## has the following elements:
## $data:     EDMA data list after sclaing/perturbation
## $natural:  data as 3D array in natural space (before rotation/translation)

## here is how to get mean/var of the pairwise Eu distances
z2 <- EDMAinR:::.edma_fit_np(sim2, less=FALSE)

## EDMA
#sim <- edma_simulate_data(n, M1, SigmaK1)
fit2 <- edma_fit(sim2)
e2 <- SigmaK_fit(fit2, S2)

# Estimate the parameters for the non-central chi-square

delta_mat2 = as.matrix(dist(Meanform(e2)))^2
phi_mat2 = (1/ncol(M2))*(z2$EuMean - delta_mat2)


# Now generate one observation from population 1 or population 2 and see if we can classify it using CL ratio

sim.new <- simulate_2d(1, M2, SigmaK2,
    deg=deg, tx=tx, ty=ty, scale=alpha, before=before)

eu_sim.new = as.matrix(dist(sim.new))^2

# Scale according to population 1 parameters
scaled_new1 = eu_sim.new/phi_mat1

ncp_mat1 = delta_mat1/phi_mat1  # This is the matrix of non-centrality parameters for population 1.
CL1 = 0
for (i in 1:(nrow(M1)-1)){
	for (j in (i+1):nrow(M1))
	CL1 = CL1 + sum(dchisq(scaled_new1[i,j],ncol(M1),ncp = ncp_mat1[i,j],log=T))
}

# Scale according to population 1 parameters
scaled_new2 = eu_sim.new/phi_mat2

ncp_mat2 = delta_mat2/phi_mat2  # This is the matrix of non-centrality parameters for population 1.
CL2 = 0
for (i in 1:(nrow(M1)-1)){
	for (j in (i+1):nrow(M1))
	CL2 = CL2 + sum(dchisq(scaled_new2[i,j],ncol(M1),ncp = ncp_mat2[i,j],log=T))
}

CL2 - CL1



## Peter -------------

x <- sim.new # specimen to classify
fit <- fit1 # fit object
.compare_objects=EDMAinR:::.compare_objects
.Eu2=EDMAinR:::.Eu2
# need fit1 and fit2
fit1 <- edma_fit(sim1, B=10)
fit2 <- edma_fit(sim2, B=10)


## calculates phi and ncp matrices from mean form M and data array A
.get_mat <- function(M, A) {
    z <- .Eu2(A)
    delta_mat <- as.matrix(dist(M))^2
    phi_mat <- (1/ncol(M)) * (z$EuMean - delta_mat)
    ncp_mat <- delta_mat / phi_mat
    list(phi=phi_mat, ncp=ncp_mat)
}

## calculates composite likelihood for ith bootstrap sample
.get_cl <- function(x, fit, i=0) {
    eu <- as.matrix(dist(x)^2)
    if (i > 0) {
        if (is.null(fit$boot) || length(fit$boot) < i)
            stop("not enough bootstrap samples")
        M <- fit$boot[[i]]$M
        s <- attr(fit$boot, "samples")[,i]
        A <- as.array(as.edma_data(fit)[,,s])
    } else {
        M <- Meanform(fit)
        A <- as.array(as.edma_data(fit))
    }
    mat <- .get_mat(M, A)
    eu_scaled <- eu / mat$phi
    CL <- dchisq(
        x=as.numeric(as.dist(eu_scaled)),
        df=ncol(x),
        ncp=pmax(0, as.numeric(as.dist(mat$ncp))),
        log=TRUE)
    sum(CL)
}
## composite likelihood ratio
edma_clr <- function(x, fit1, fit2, boot=FALSE) {
    if (inherits(x, "edma_data")) {
        if (dim(x)[3] > 1L)
            stop("provide a single specimen only when x is an edma_data object")
        x <- x$data[[1L]]
    }
    .compare_objects(fit1, fit2)
    if (!identical(dimnames(fit1)[1:2], dimnames(x)))
        stop("specimen dimnames must match the EDMA fit objects")
    CL1 <- .get_cl(x, fit1)
    CL2 <- .get_cl(x, fit2)
    BOOT <- NULL
    B1 <- length(fit1$boot)
    B2 <- length(fit2$boot)
    if (boot && (B1 < 1L || B2 < 1L))
        warning("no bootstrap samples found: boot=TRUE ignored")
    if (boot && B1 > 0L && B2 > 0L) {
        CL1M <- sapply(seq_len(B1), function(i) .get_cl(x, fit1, i))
        CL2M <- sapply(seq_len(B2), function(i) .get_cl(x, fit2, i))
        CLRM <- t(outer(CL2M, CL1M, "-"))
        BOOT <- data.frame(
            complik1=as.numeric(CL1M),
            complik2=as.numeric(CL2M),
            complikr=as.numeric(CLRM))
    }
    ## CLR > 0 => x belongs to fit2
    ## CLR < 0 => x belongs to fit1
    out <- list(
        x=x,
        fit1=fit1,
        fit2=fit2,
        complik1=CL1,
        complik2=CL2,
        complikr=CL2 - CL1,
        boot=BOOT)
    class(out) <- "edma_clr"
    out
}



h <- edma_clr(x, fit1, fit2, boot=TRUE)

c(CLR=h$complikr, quantile(h$boot$complikr, c(0.025, 0.975)))

## -------------------


library(shapes)

data(schizophrenia.dat)
A <- schizophrenia.dat
x_con <- as.edma_data(A[,,1:14])
x_scz <- as.edma_data(A[,,15:28])

fit_con <- edma_fit(x_con)
fit_scz <- edma_fit(x_scz)

v <- numeric(28)
for (i in 1:14) {
    fit_con_less <- edma_fit(x_con[,,-i])
    fit_scz_less <- edma_fit(x_scz[,,-i])
    v[i] <- edma_clr(fit_con$data[[i]], fit_con_less, fit_scz)$complikr
    v[i+14] <- edma_clr(fit_scz$data[[i]], fit_con, fit_scz_less)$complikr
}

class_true <- c(rep(0, 14), rep(1, 14))
class_clr <- ifelse(v > 0, 1, 0)

cm <- table(class_true, class_clr)
1-sum(diag(cm)) / sum(cm)
cm

x1 <- as.edma_data(A[,,1:14])
x2 <- as.edma_data(A[,,15:28])

i <- 2
## loo function
B <- 0
level <- 0.95

loo <- function(x1, x2, B=0, level=0.95) {
    a <- c((1-level)/2, 1-(1-level)/2)
    x1 <- as.edma_data(x1)
    x2 <- as.edma_data(x2)
    n1 <- dim(x1)[3L]
    n2 <- dim(x2)[3L]
    x12 <- combine_data(x1, x2)
    Class <- data.frame(
        group=rep(1:2, c(n1, n2)),
        specimen=c(seq_len(n1), seq_len(n2)),
        id=seq_len(n1 + n2),
        complikr=NA, lower=NA, upper=NA)
    for (i in seq_len(n1 + n2)) {
        x <- x12$data[[i]]
        i1 <- which(Class$group == 1 & Class$id != i)
        i2 <- which(Class$group == 2 & Class$id != i)
        fit1 <- edma_fit(x12[,,i1], B=B)
        fit2 <- edma_fit(x12[,,i2], B=B)
        h <- edma_clr(x, fit1, fit2, boot=B > 0)
        Class$complikr[i] <- h$complikr
        if (B > 0 && !is.null(h$boot)) {
            q <- quantile(c(h$complikr, h$boot$complikr), a, na.rm=TRUE)
            Class$lower[i] <- q[1L]
            Class$upper[i] <- q[2L]
        }
    }
    Class$class <- ifelse(Class$complikr > 0, 2, 1)
    Class$signif <- !(Class$lower < 0 & Class$upper > 0)
    Class
}


## Bookstein's schizophrenia data
data(schizophrenia.dat)
x1 <- as.edma_data(schizophrenia.dat[,,1:14]) # control
x2 <- as.edma_data(schizophrenia.dat[,,15:28]) # schizo
l <- loo(x1, x2)
(cm <- table(l$group, l$class)) # confusion matrix
round(100*(1-sum(diag(cm)) / sum(cm)), 2) # error rate

## Gorillas
data(apes)
x1 <- as.edma_data(apes$x[,,apes$group == "gorf"]) # females
x2 <- as.edma_data(apes$x[,,apes$group == "gorm"]) # males
l <- loo(x1, x2)
(cm <- table(l$group, l$class)) # confusion matrix
round(100*(1-sum(diag(cm)) / sum(cm)), 2) # error rate

## Chimps
x1 <- as.edma_data(apes$x[,,apes$group == "panf"]) # females
x2 <- as.edma_data(apes$x[,,apes$group == "panm"]) # males
l <- loo(x1, x2)
(cm <- table(l$group, l$class)) # confusion matrix
round(100*(1-sum(diag(cm)) / sum(cm)), 2) # error rate

## Orangs
x1 <- as.edma_data(apes$x[,,apes$group == "pongof"]) # females
x2 <- as.edma_data(apes$x[,,apes$group == "pongom"]) # males
l <- loo(x1, x2)
(cm <- table(l$group, l$class)) # confusion matrix
round(100*(1-sum(diag(cm)) / sum(cm)), 2) # error rate

## Mice
data(mice)
x1 <- as.edma_data(mice$x[,,mice$group == "l"]) # large
x2 <- as.edma_data(mice$x[,,mice$group == "s"]) # small
l <- loo(x1, x2)
(cm <- table(l$group, l$class)) # confusion matrix
round(100*(1-sum(diag(cm)) / sum(cm)), 2) # error rate

