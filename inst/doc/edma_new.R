
source("R/read.R")
source("R/nonparametric.R")
file <- "inst/extdata/crouzon/Crouzon_P0_Global_MUT.xyz"
x <- read_xyz(file)
x
dim(x)
dimnames(x)
subset(x, 1:10)
x[1:10, 2:3, 1:5]
str(as.matrix(x))
str(as.data.frame(x))
str(as.array(x))
str(X <- stack(x, TRUE))


fit <- edma_fit(x)
object <- edma_fit(x, B=10)
Meanform(fit)
SigmaKstar(fit)

print.edma_fit_np <- function(x, ...) {
    cat("EDMA nonparametric fit: ", x$name, "\n",
        ncol(x$data[[1L]]), " dimensions, ",
        nrow(x$data[[1L]]), " landmarks, ",
        length(x$data), " replicates, ",
        if (length(x$boot))
            paste(length(x$boot) + 1L, "bootstrap runs") else "no bootstrap",
        sep="")
    invisible(x)
}

Meanform <- function (object, ...) UseMethod("Meanform")
Meanform.edma_fit <- function (object, ...) object$M

SigmaKstar <- function (object, ...) UseMethod("SigmaKstar")
SigmaKstar.edma_fit_np <- function (object, ...) object$SigmaKstar

as.dist.edma_fit <- function(m, diag = FALSE, upper = FALSE) {
    out <- dist(Meanform(m), diag=diag, upper=upper)
    class(out) <- c(class(out), "edma_dist")
    out
}

stack.dist <- function(x) {
    id <- as.matrix(x)
    id[] <- 0
    id[lower.tri(id)] <- 1
    rm <- row(id)
    cm <- col(id)
    rm <- array(rm)[array(id) == 1]
    cm <- array(cm)[array(id) == 1]
    d <- as.vector(x)
    out <- data.frame(row=rm, col=cm, dist=d)
    out$row <- as.factor(out$row)
    out$col <- as.factor(out$col)
    levels(out$row) <- rownames(id)[-1]
    levels(out$col) <- colnames(id)[-ncol(id)]
    out
}

stacked_dist <- function (object, ...) UseMethod("stacked_dist")
stacked_dist.edma_fit <- function (object, sort=FALSE, ...) {
    out <- stack(as.dist(object, diag = FALSE, upper = FALSE))
    if (sort)
        out <- out[order(out$dist, ...),]
    out
}
head(stacked_dist(fit))
head(stacked_dist(fit, sort=TRUE, decreasing=TRUE))
head(stacked_dist(fit, sort=TRUE, decreasing=FALSE))

form_difference <- function (numerator, denominator, ...)
    UseMethod("form_difference")

file1 <- "inst/extdata/crouzon/Crouzon_P0_Global_MUT.xyz"
file2 <- "inst/extdata/crouzon/Crouzon_P0_Global_NON-MUT.xyz"
x1 <- read_xyz(file1)
x2 <- read_xyz(file2)
numerator <- edma_fit(x1)
denominator <- edma_fit(x2)

form_difference.edma_fit <- function (numerator, denominator, ...) {
    f <- function(a, b) as.dist(a) / as.dist(b)
    r <- f(numerator, denominator)
    attr(r, "method") <- "euclidean_distance_ratio"
    attr(r, "call") <- NULL
    attr(r, "T") <- max(r) / min(r)
    r
}

## TODO:

## - edma_test: B+1 T values & P-value (value > Tobs / (B+1))
## - assess bootpstrat: !NULL, B1=B2
## - make stacked form diff with marginal CI

## - parametric fit for sig2*I
## - structural assessment
## - input structure

## - xlsx as output
## - 2D/3D plotting
## - Growth assessment (2x + 2x combo)


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
