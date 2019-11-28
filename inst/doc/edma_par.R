## take design matrix and returns a factor for parameter matching
## NA --> 0
## other unique values will be factor levels
.mat2fac <- function(m) {
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

.SigmaK_fit <- function(SigmaKstar, H, pattern, init,
method = "Nelder-Mead", control = list(), hessian = FALSE) {
    ## test that SigmaKstar is of the right (K-1) rank
    r <- qr(SigmaKstar)$rank
    K <- nrow(SigmaKstar)
    if (r != K-1)
        warning(sprintf("rank of SigmaKstar was %s instead of %s", r, K-1))
    if (any(is.na(diag(pattern))))
        stop("pattern matrix must have parameters in the diagonal")
    if (dim(pattern)[1L] != dim(pattern)[1L])
        stop("pattern matrix must be square matrix")
    pattern1 <- pattern
    pattern1[] <- NA
    diag(pattern1) <- diag(pattern)
    pattern0 <- pattern
    diag(pattern0) <- NA
    fac <- .mat2fac(pattern) # all parms
    lev1 <- levels(.mat2fac(pattern1)) # parms in diag
    lev0 <- levels(.mat2fac(pattern0)) # parms off-diag
    if (length(intersect(lev1, lev0)) > 0)
        stop("diagonal and off-diagonal parameters must not overlap")
    p <- nlevels(droplevels(fac))
    if (p > K*(K-1)/2)
        stop(
            sprintf(
                "pattern with %s unknowns implies non estimable SigmaK",
            p))
    ## make sure diags are >0
    ## generate random starting values using runif()
    init0 <- structure(numeric(p), names=levels(fac))
    if (missing(init)) {
        init0[lev1] <- runif(length(lev1))
    } else {
        init <- init[names(init) %in% names(init0)]
        init0[names(init)] <- init
    }
    if (any(init0[lev1] <= 0))
        stop("inits for diagonal elements must be > 0")
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
        optim(init0, fun, method=method, control=control, hessian=hessian)
    })
    o$init <- init0
    o$SigmaK <- .vec2mat(o$par, fac)
    o
}


## test pattern matrix manipulation
m <- matrix(c(
    "a", NA, NA, NA,
    "c", "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, "d", "b"
), 4, 4, byrow=TRUE)
parm <- c(a=1, b=2, c=3, d=4)
(fac <- .mat2fac(m))
(S <- .vec2mat(parm, fac))

# Generate the data and test the method

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

# this workes, but it is harder to find
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

m <- matrix(c(
    "a", NA, NA, NA,
    "c", "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, "d", "b"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.3, c=0.075, d=0.09)


f <- function(fit)
    unlist(.SigmaK_fit(fit$SigmaKstar, fit$H, m)[1:2])
M <- structure(c(-2.5, 7.5, -2.5, -2.5, -7.5, 2.5, 2.5, 4.5),
    .Dim = c(4L, 2L))
SigmaK <- .vec2mat(parm, .mat2fac(m))

# barebones version
sim <- EDMAinR:::.edma_simulate_data(n=1000, M, SigmaK)
fit <- EDMAinR:::.edma_fit_np(sim$X, sim$n, sim$K, sim$D)
str(o <- .SigmaK_fit(fit$SigmaKstar, fit$H, m))
cbind(true=parm, est=o$par)
summary(t(replicate(10, f(fit))))

# nicely formatted version
sim2 <- EDMAinR::edma_simulate_data(n=1000, M, SigmaK)
fit2 <- EDMAinR::edma_fit(sim2)
str(o <- .SigmaK_fit(fit2$SigmaKstar, fit2$H, m))
cbind(true=parm, est=o$par)
summary(t(replicate(10, f(fit))))

