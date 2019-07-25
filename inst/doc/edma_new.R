
source("R/read.R")
source("R/nonparametric.R")
file1 <- "inst/extdata/crouzon/Crouzon_P0_Global_MUT.xyz"
file2 <- "inst/extdata/crouzon/Crouzon_P0_Global_NON-MUT.xyz"
x1 <- read_xyz(file1)
x2 <- read_xyz(file2)
dim(x1)
dimnames(x1)
subset(x1, 1:10)
x1[1:10, 2:3, 1:5]
str(as.matrix(x1))
str(as.data.frame(x1))
str(as.array(x1))
str(stack(x1, TRUE))


fit <- edma_fit(x1)
object <- edma_fit(x1, B=10)
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
    d <- as.dist(object, diag = FALSE, upper = FALSE)
    out <- stack(d)
    attr(out, "method") <- attr(d, "method")
    attr(out, "Tval") <- attr(d, "Tval")
    if (sort)
        out <- out[order(out$dist, ...),]
    out
}
head(stacked_dist(fit))
head(stacked_dist(fit, sort=TRUE, decreasing=TRUE))
head(stacked_dist(fit, sort=TRUE, decreasing=FALSE))

## a and b are edma_fit objects
.compare_objects <- function (a, b, ...) {
    if (nrow(a$data[[1L]]) != nrow(b$data[[1L]]))
        stop("number of landmarks must be identical")
    if (ncol(a$data[[1L]]) != ncol(b$data[[1L]]))
        stop("number of dimensions must be identical")
    if (length(a$boot) != length(b$boot))
        stop("number of bootstrap runs  must be identical")
    if (!all(rownames(a$data[[1L]]) == rownames(b$data[[1L]])))
        stop("landmark names and ordering must be identical")
    if (!all(colnames(a$data[[1L]]) == colnames(b$data[[1L]])))
        stop("dimension names and ordering must be identical")
    invisible(NULL)
}

## inputs are meanform matrices
.formdiff <- function(M1, M2) {
    r <- dist(M1) / dist(M2)
    attr(r, "method") <- "euclidean_distance_ratio"
    attr(r, "call") <- NULL
    attr(r, "Tval") <- max(r) / min(r)
    r
}
form_difference <- function (numerator, denominator, ...)
    UseMethod("form_difference")
form_difference.edma_fit <- function (numerator, denominator, ...) {
    .compare_objects(numerator, denominator)
    .formdiff(Meanform(numerator), Meanform(denominator))
}

numerator <- edma_fit(x1, B=10)
denominator <- edma_fit(x2, B=10)

edma_test <- function (numerator, denominator) {
    .compare_objects(numerator, denominator)
    DNAME <- paste(numerator$name, denominator$name, sep = ", ")
    METHOD <- "Bootstrap based EDMA T-test"
    B <- length(numerator$boot)
    Tval <- attr(form_difference(numerator, denominator), "Tval")
    Tvals <- c(Tval, sapply(seq_len(B), function(i) {
        attr(.formdiff(numerator$boot[[i]]$M, denominator$boot[[i]]$M), "Tval")
    }))
    PVAL = sum(Tvals > Tval) / (B + 1)
    PARAMETER <- B + 1L
    names(Tval) <- "T-value"
    names(PARAMETER) <- "B"
    structure(list(statistic = Tval, parameter = PARAMETER,
        p.value = PVAL, method = METHOD, data.name = DNAME,
        Tvals=Tvals),
        class = "htest")
}
stacked_form_difference <- function (numerator, denominator,
sort=FALSE, level=0.95) {
    .compare_objects(numerator, denominator)
    d <- form_difference(numerator, denominator)
    out <- stack(d)
    B <- length(numerator$boot)
    fd <- sapply(seq_len(B), function(i) {
        stack(.formdiff(numerator$boot[[i]]$M, denominator$boot[[i]]$M))$dist
    })
    a <- c((1-level)/2, 1-(1-level)/2)
    q <- t(apply(cbind(out$dist, fd), 1, quantile, a))
    attr(out, "method") <- attr(d, "method")
    attr(out, "Tval") <- attr(d, "Tval")
    attr(out, "level") <- level
    out$lower <- q[,1L]
    out$upper <- q[,2L]
    out$inside <- out$dist <= out$upper & out$dist >= out$lower
    if (sort)
        out <- out[order(out$dist, ...),]
    out
}
## TODO:

## OK - edma_test: B+1 T values & P-value (value > Tobs / (B+1))
## OK - assess bootpstrat: !NULL, B1=B2
## OK - make stacked form diff with marginal CI

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
