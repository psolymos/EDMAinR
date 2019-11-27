## take design matrix and returns a factor for parameter matching
## NA --> 0
## other unique values will be factor levels
.mat2fac <- function(m) {
    if (dim(m)[1L] != dim(m)[1L])
        stop("matrix must be square matrix")
    if (any(is.na(diag(m))))
        stop("values in the diagonal must not be NA")
    K <- nrow(m)
    m[upper.tri(m)] <- m[lower.tri(m)]
    m <- factor(as.character(m))
    attr(m, "K") <- K
    m
}

## takes a named param vector and reconstructs matrix
.vec2mat <- function(parm, fac) {
    x <- matrix(0, attr(fac, "K"), attr(fac, "K"))
    fac <- droplevels(fac)
    if (!all(sort(names(parm)) == sort(levels(fac))))
        stop("names of parm and levels of fact must match")
    for (i in names(parm))
        x[!is.na(fac) & fac == i] <- parm[i]
    x
}

## simulate stuff -- this is already there as edma_simulate_data !!!
.edma_simulate <- function(n, M, SigmaK) {
    D <- ncol(M)
    K <- nrow(M)
    Z <- matrix(nrow = n * K, ncol = D)
    for (i in 1:n) {
        Z[((i - 1) * K + 1):(i * K), ] <-
            matrix(rnorm(K * D), nrow = K, ncol = D)
    }
    C <- chol(SigmaK)
    X <- matrix(nrow = n * K, ncol = D)
    for (i in 1:n) {
        X[((i - 1) * K + 1):(i * K), ] <-
            crossprod(C, Z[((i - 1) * K + 1):(i * K), ]) + M
    }

    I <- diag(1, K)
        ones <- array(rep(1, K), c(1, K))
        H <- I - (1/K) * crossprod(ones, ones)
    SigmaKstar = H %*% SigmaK %*% H
    list(
        M=M, SigmaK=SigmaK, SigmaKstar=SigmaKstar, H=H,
        X=X, D=D, K=K, n=n)
}

m <- matrix(c(
    "a", NA, NA, NA,
    "c", "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, "d", "b"
), 4, 4, byrow=TRUE)
parm <- c(a=1, b=2, c=3, d=4)
fac <- .mat2fac(m)
S <- .vec2mat(parm, fac)
.estimable_SigmaK(S)

# Generate the data and test the method


m <- matrix(c(
    "a", NA, NA, NA,
    "c", "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, "d", "b"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.3, c=0.075, d=0.09)

## this works
m <- matrix(c(
    "a", NA, NA, NA,
    NA, "a", NA, NA,
    NA,  NA, "a", NA,
    NA,  NA, NA, "a"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25)

## this works
m <- matrix(c(
    "a", NA, NA, NA,
    NA, "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, NA, "b"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.35)

# this worked
m <- matrix(c(
    "a", NA, NA, NA,
    NA, "b", NA, NA,
    NA,  NA, "c", NA,
    NA,  NA, NA, "d"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.35, c=0.1, d=0.05)

## this works
m <- matrix(c(
    "a", NA, NA, NA,
    "b", "a", NA, NA,
    NA,  NA, "a", NA,
    NA,  NA, "b", "a"
), 4, 4, byrow=TRUE)
parm <- c(a=0.15, b=0.05)

## this does not works
m <- matrix(c(
    "a", NA, NA, NA,
    "b", "a", NA, NA,
    "b", "b", "a", NA,
    "b", "b", "b", "a"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.07)

.SigmaK_fit <- function(SigmaKstar, H, pattern, init,
method = "Nelder-Mead", control = list(), hessian = FALSE) {
    ## test that SigmaKstar is of the right (K-1) rank
    K <- nrow(SigmaKstar)
    fac <- .mat2fac(pattern)
    p <- nlevels(droplevels(fac))
    if (p > K*(K-1)/2)
        stop(
            sprintf(
                "pattern with %s unknowns implies non estimable SigmaK",
            p))
    if (missing(init))
        init <- structure(numeric(p), names=levels(fac))
    num_max <- .Machine$double.xmax^(1/3)
    ## we might need constraints here, i.e. >0 diag values
    fun <- function(parms){
        SigmaK <- .vec2mat(parms, fac)
        if (any(diag(SigmaK) <= 0))
            return(num_max)
        10^4 * max((SigmaKstar - (H %*% SigmaK %*% H))^2)
    }
    if (!is.null(control$fnscale) && control$fnscale < 0)
        stop("control$fnscale can not be negative")
    o <- suppressWarnings({
        optim(init, fun, method=method, control=control, hessian=hessian)
    })
    o$SigmaK <- .vec2mat(o$par, fac)
    o
}

M <- structure(c(-2.5, 7.5, -2.5, -2.5, -7.5, 2.5, 2.5, 4.5),
    .Dim = c(4L, 2L))
SigmaK <- .vec2mat(parm, .mat2fac(m))
sim <- .edma_simulate(n=1000, M, SigmaK)

fit <- EDMAinR:::.edma_fit_np(sim$X, sim$n, sim$K, sim$D)
o <- .SigmaK_fit(fit$SigmaKstar, fit$H, m, init=c(a=runif(0, 1), b=0))
o

summary(t(replicate(10, unlist(.SigmaK_fit(fit$SigmaKstar, fit$H, m, init=c(a=runif(1), b=0))[1:2]))))

cbind(true=parm, est=o$par)

