
source("R/read.R")
source("R/nonparametric.R")
file1 <- "inst/extdata/crouzon/Crouzon_P0_Global_MUT.xyz"
file2 <- "inst/extdata/crouzon/Crouzon_P0_Global_NON-MUT.xyz"
x1 <- read_xyz(file1)
x2 <- read_xyz(file2)
x1
x2
dim(x1)
dimnames(x1)
landmark_names(x1)
subset(x1, 1:10)
x1[1:10, 2:3, 1:5]
str(as.matrix(x1))
str(as.data.frame(x1))
str(as.array(x1))
str(stack(x1))

fit <- edma_fit(x1)
Meanform(fit)
SigmaKstar(fit)
head(stacked_dist(fit))
head(stacked_dist(fit, sort=TRUE, decreasing=TRUE))
head(stacked_dist(fit, sort=TRUE, decreasing=FALSE))

numerator <- edma_fit(x1, B=10)
denominator <- edma_fit(x2, B=10)
fd <- form_difference(numerator, denominator)
str(fd)
edma_test(numerator, denominator)
sfd <- stacked_form_difference(numerator, denominator, sort=TRUE)
str(sfd)

fit2 <- SigmaK_fit(fit)
str(fit2$optim)

## TODO:

## OK - edma_test: B+1 T values & P-value (value > Tobs / (B+1))
## OK - assess bootpstrat: !NULL, B1=B2
## OK - make stacked form diff with marginal CI

## OK - parametric fit for sig2*I
## OK - structural assessment
## - input structure

## - xlsx as output
## - 2D/3D plotting
## - Growth assessment (2x + 2x combo)


## parametric testing

# Generate the data and test the method

K <- 3 # number of landmarks
D <- 2 # dimension, 2 or 3

sig <- 0.75
rho <- 0
SigmaK <- sig^2*diag(1, K, K) + sig^2*rho*(1-diag(1, K, K))

M <- matrix(c(0,1,0,0,0,1), 3, 2)
M[,1] <- M[,1] - mean(M[,1])
M[,2] <- M[,2] - mean(M[,2])
M <- 10*M

n <- 100

Z <- matrix(nrow = n * K, ncol = D)
for (i in 1:n) {
    Z[((i - 1) * K + 1):(i * K), ] <- matrix(rnorm(K * D), nrow = K,
        ncol = D)
}
Cmat <- chol(SigmaK)
X <- matrix(nrow = n * K, ncol = D)
for (i in 1:n) {
    X[((i - 1) * K + 1):(i * K), ] <- crossprod(Cmat, Z[((i - 1) * K + 1):(i *
        K), ]) + M
}

I <- diag(1, K)
ones <- array(rep(1, K), c(1, K))
H <- I - (1/K) * crossprod(ones, ones)
SigmaKstar = H %*% SigmaK %*% H

#fit <- edma_fit.new(X, n, K, D)
fit <- .edma_fit_np(X, n, K, D)
o <- .SigmaK_fit(fit$SigmaKstar, fit$H)
o

# compare the non-parametric estimates
## sigma
c(true=sig, est=o$coefficients)
## M
cbind(true=dist(M), est=dist(fit$M))
## SigmaKstar
cbind(true=c(diag(SigmaKstar), SigmaKstar[upper.tri(SigmaKstar)]),
    nonpar=c(diag(fit$SigmaKstar), fit$SigmaKstar[upper.tri(fit$SigmaKstar)]))
## SigmaK
cbind(true=c(diag(SigmaK), SigmaK[upper.tri(SigmaK)]),
    param=c(diag(o$SigmaK), o$SigmaK[upper.tri(o$SigmaK)]))



## stuff from Subhash--------------------

## workhorse function to estimate centered mean form and SigmaKstar
## X: read.table(x) # data frame input, n tables (each K x D) stacked ???
## n: number of individuals (replicates)
## K: number of landmarks
## D: dimensions
edma_fit.new <- function(X, n, K, D) {
    ## step1: calculate Euclidean distance for X
    d <- function(x, y) {
        sum((x - y)^2)
    }
    EuX <- matrix(nrow = n * K, ncol = K)
    for (i in 1:n) {
        for (j in ((i - 1) * K + 1):(i * K)) {
            for (l in ((i - 1) * K + 1):(i * K)) {
                EuX[j, (l - (i - 1) * K)] <- d(X[j, ], X[l, ])
            }
        }
    }
    ## EuX step 2 calculate Eu(M)
    EuM <- matrix(nrow = K, ncol = K)
    for (i in 1:K) {
        for (j in 1:K) {
            sum1 <- 0
            for (l in 1:n) {
                sum1 <- sum1 + EuX[((l - 1) * K + i), j]
            }
            mean.squ.eudis <- sum1/n
            sum2 <- 0
            for (l in 1:n) {
                sum2 <- sum2 + (EuX[((l - 1) * K + i), j] - mean.squ.eudis)^2
            }
            var.squ.eudis <- sum2/n
            if (D == 2) {
                EuM[i, j] <- sqrt((mean.squ.eudis)^2 - var.squ.eudis)
            } else {
                EuM[i, j] <- sqrt((mean.squ.eudis)^2 - 1.5 * var.squ.eudis)
            }
        }
    }
    ## step 3 calculate B(M)
    I <- diag(1, K)
    ones <- array(rep(1, K), c(1, K))
    H <- I - (1/K) * crossprod(ones, ones)
    BM <- (-1/2) * H %*% EuM %*% H
    ## step 4 calculate eigenvalues and eigenvectors of B(M)
    EIG <- eigen(BM)
    ## step 5 estimate centred mean form M
    est.M <- matrix(nrow = K, ncol = D)
    for (i in 1:D) {
        est.M[, i] <- sqrt(EIG$values[i]) * EIG$vectors[, i]
    }
    ## step 6 estimate SigmaKstar=H*SigmaK*H
    BX <- matrix(nrow = n * K, ncol = K)
    for (i in 1:n) {
        BX[((i - 1) * K + 1):(i * K), ] <- (-1/2) * H %*% EuX[((i - 1) *
            K + 1):(i * K), ] %*% H
    }
    a <- 0
    for (i in 1:n) {
        a <- a + BX[((i - 1) * K + 1):(i * K), ]
    }
    est.SigmaKstar <- (a/n - BM)/D
out = list (est.M=est.M, est.SigmaKstar = est.SigmaKstar,H=H,K=K)
return(out)
}

# Create the parametric matrix as a function of parameters

V.fn = function(parms,K){
	v <- diag(1, K, K)
	sig <- exp(parms[1])
            rho <- tanh(parms[2])
            V <- sig^2 * v + sig^2 * rho * (1-v)
            return(V)
}

# Function to be minimized

fn.min = function(parms,K,est.SigmaKstar,H){

	out = max((est.SigmaKstar - (H %*% V.fn(parms,K) %*% t(H)))^2)
	return(10000*out)
}


# Generate the data and test the method

K <- 3 # number of landmarks
D <- 2 # dimension, 2 or 3

sig <- 0.5
rho <- 0.3
SigmaK <- sig^2*diag(1, K, K) + sig^2*rho*(1-diag(1, K, K))

M <- matrix(c(0,1,0,0,0,1), 3, 2)
M[,1] <- M[,1] - mean(M[,1])
M[,2] <- M[,2] - mean(M[,2])
M <- 10*M

n <- 10000
Z <- matrix(nrow = n * K, ncol = D)
for (i in 1:n) {
    Z[((i - 1) * K + 1):(i * K), ] <- matrix(rnorm(K * D), nrow = K,
        ncol = D)
}
C <- chol(SigmaK)
X <- matrix(nrow = n * K, ncol = D)
for (i in 1:n) {
    X[((i - 1) * K + 1):(i * K), ] <- crossprod(C, Z[((i - 1) * K + 1):(i *
        K), ]) + M
}

I <- diag(1, K)
    ones <- array(rep(1, K), c(1, K))
    H <- I - (1/K) * crossprod(ones, ones)
SigmaKstar = H %*% SigmaK %*% H

#fit <- edma_fit.new(X, n, K, D)
fit <- .edma_fit_np(X, n, K, D)

# compare the non-parametric estimates
cbind(true=dist(M), est=dist(fit$M))
cbind(true=c(diag(SigmaKstar), SigmaKstar[upper.tri(SigmaKstar)]),
    est=c(diag(fit$SigmaKstar), SigmaKstar[upper.tri(fit$SigmaKstar)]))

# Now do the minimization to get the parameters of the parametric covariance matrix


parms = c(log(sig+0.1),atanh(rho+0.1))
K = 3
est.SigmaKstar = fit$est.SigmaKstar
H = fit$H
fn.min(parms,K,est.SigmaKstar,H)

Sigma.par = optim(parms,fn.min,K=K,est.SigmaKstar=est.SigmaKstar,H=H)
Sigma.par

# Compare the estimated and true covariances.

parms.true = c(log(sig),atanh(rho))
cbind(V.fn(parms.true,K),V.fn(Sigma.par$par,K))
cbind(H %*% V.fn(parms.true,K) %*% H,H %*% V.fn(Sigma.par$par,K) %*% H)






SigmaK_fit <- function(fit,
type=c("sig", "sig_rho"),
method = "Nelder-Mead", control = list(), hessian = FALSE)
{
    type <- match.arg(type)
    v <- diag(1, fit$K, fit$K)
    if (type == "sig") {
        init <- log(sqrt(mean(diag(fit$est.SigmaKstar))))
        fun <- function(par) {
            sig <- exp(par[1L])
            V <- sig^2 * v
            sum((fit$SigmaK - (fit$H %*% V %*% t(fit$H)))^2)
        }
        o <- suppressWarnings({
            optim(init, fun, method=method, control=control, hessian=hessian)
        })
        o$coefficients <- c(sigma=exp(o$par))
        o$SigmaKhat <- exp(o$par)^2 * v
    }
    if (type == "sig_rho") {
        init <- c(log(sqrt(mean(diag(fit$est.SigmaKstar)))),
                  atanh(mean(fit$SigmaK[lower.tri(fit$est.SigmaKstar)])))
        fun <- function(par) {
            sig <- exp(par[1L])
            rho <- tanh(par[2L])
            V <- sig^2 * v + sig^2 * rho * (1-v)
            sum((fit$SigmaK - (fit$H %*% V %*% t(fit$H)))^2)
        }
        o <- optim(init, fun, method=method, control=control, hessian=hessian)
        o$coefficients <- c(sigma=exp(o$par[1L]),
                            rho=tanh(o$par[2L]))
        o$SigmaKhat <- exp(o$par[1L])^2 * v +
            exp(o$par[1L])^2 * tanh(o$par[2L]) * (1-v)
    }
    fit$optim <- o
    fit
}

SigmaK_fit(fit)

## Cheverud et al. 1997
Mandible_data_meanform <- matrix(c(
    -2.5642147, 2.9238803,
    -4.3598175, 2.1062128,
    -5.0189366, 2.0386707,
    -5.7237441, 0.6476102,
    -3.4644795, -1.1693193,
    -4.5855937, -2.9159961,
    -3.3321319, -3.2438449,
    0.1369185, -1.7688938,
    2.4780984, -2.1174301,
    4.6209991, -1.6495556,
    6.0020391, -0.5078821,
    6.5324078, 1.7178527,
    4.2879408, 0.3604548,
    3.0344353, 1.1350108,
    1.2218793, 1.1121225,
    0.7341997, 1.3311070
), ncol=2, byrow=TRUE)
#1 Coronoid tip
#2 Anterior condylar facet
#3 Posterior condylar facet
#4 Inferior condylar facet
#5 Deepest point of mandibular notch
#6 Posterior angular process
#7 Inferior angular process
#8 Anterior angular process
#9 Posterior inferior corpus
#10 Anterior inferior corpus
#11 Inferior incisor alveolus
#12 Superior incisor alveolus
#13 Deepest point of incisive notch
#14 Anterior molar alveolus
#15 Posterior molar alveolus
#16 Coronoid base
