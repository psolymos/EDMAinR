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
fit <- edma_fit.new(X, n, K, D)

# compare the non-parametric estimates
cbind(dist(M),dist(fit$est.M))
cbind(SigmaKstar,fit$est.SigmaKstar)

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