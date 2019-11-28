library(EDMAinR)

## test pattern matrix manipulation
m <- matrix(c(
    "a", NA, NA, NA,
    "c", "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, "d", "b"
), 4, 4, byrow=TRUE)
parm <- c(a=1, b=2, c=3, d=4)
(fac <- EDMAinR:::.mat2fac(m))
(S <- EDMAinR:::.vec2mat(parm, fac))

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

## this does not work
m <- matrix(c(
    "a", NA, NA, NA,
    "b", "a", NA, NA,
    "b", "b", "a", NA,
    "b", "b", "b", "a"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.07)

## this does not work
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
    "c", "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, NA, "b"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.3, c=0.075)

f <- function(fit)
    unlist(.SigmaK_fit(fit$SigmaKstar, fit$H, m)[1:2])
M <- structure(c(-2.5, 7.5, -2.5, -2.5, -7.5, 2.5, 2.5, 4.5),
    .Dim = c(4L, 2L))
SigmaK <- EDMAinR:::.vec2mat(parm, EDMAinR:::.mat2fac(m))

if (FALSE) {
# barebones version
#sim <- EDMAinR:::.edma_simulate_data(n=1000, M, SigmaK)
sim <- .edma_simulate_data(n=1000, M, SigmaK)
fit <- EDMAinR:::.edma_fit_np(sim$X, sim$n, sim$K, sim$D)
str(o <- .SigmaK_fit(fit$SigmaKstar, fit$H, m))
cbind(true=parm, est=o$par)
summary(t(replicate(10, f(fit))))
}

# nicely formatted version
sim <- edma_simulate_data(n=1000, M, SigmaK)
fit <- edma_fit(sim, B=10)
o <- SigmaK_fit(fit, m)
cbind(true=parm, est=o$results$par)
sensitivity(o)

## Liangyuan's estimator

## pattern is a 0/1 matrix: 1 indicating unknowns
.full_p_fun <- function(object, pattern) {
    est.M <- Meanform(object)
    est.SigmaKstar <- SigmaKstar(object)
    if (!all(dim(est.SigmaKstar) == dim(pattern)))
        stop("Dimension mismatch for pattern.")
    K <- nrow(est.M)
    D <- ncol(est.M)
    ## Use maple to solve Y such that YH=L
    fcol <- array(c(rep(-1, K - 2), -2), c(K - 1, 1))
    identity <- diag(1, K - 2, K - 1)
    augment <- array(c(rep(-1, K - 2), 0), c(1, K - 1))
    Y <- cbind(fcol, rbind(identity, augment))
    ## from est.SigmaKstar to est.SigmaK~
    est.SigmaKtilde <- Y %*% est.SigmaKstar %*% t(Y)
    L <- cbind(c(rep(-1, K - 1)), diag(1, K - 1))
    ## L%*%SigmaK%*%t(L) # model for diagonal matrix
    cvec <- as.vector(est.SigmaKtilde)
    A <- matrix(nrow = (K - 1)^2, ncol = K)
    iden <- diag(1, K - 1)
    for (i in 1:(K - 1)) {
        A[(i - 1) * K + 1, ] <- c(1, iden[i, ])
    }
    for (i in 1:(K - 2)) {
        A[((i - 1) * K + 2):(i * K), ] <- matrix(rep(c(1, rep(0, K - 1)),
            K - 1), nrow = K - 1, byrow = T)
    }
    ## A
    est.vect <- round((solve(t(A) %*% A)) %*% t(A) %*% cvec, 3)
    est.SigmaK <- diag(c(est.vect))
    ## general case
    #c <- as.vector(est.SigmaKtilde)
    #b.ori <- as.vector(SigmaK)
    m <- 0
    b.dis <- rep(0, K * (K + 1)/2)
    for (i in 1:K) {
        for (j in i:K) {
            m <- m + 1
            b.dis[m] <- pattern[i, j]
        }
    }
    d <- 0
    for (i in 1:(K * (K + 1)/2)) {
        if (b.dis[i] != 0)
            d <- c(d, b.dis[i])
    }
    b <- d[-1]
    c.dis <- rep(0, K * (K - 1)/2)
    q <- 0
    for (i in 1:(K - 1)) {
        for (j in i:(K - 1)) {
            q <- q + 1
            c.dis[q] <- est.SigmaKtilde[i, j]
        }
    }
    iden <- diag(1, K - 1)
    A <- matrix(nrow = (K - 1)^2, ncol = K^2)
    A[, K + 1] <- 0
    for (i in 1:(K - 1)) {
        A[, 1] <- 1
        A[((i - 1) * (K - 1) + 1):(i * (K - 1)), 2:K] <- diag(-1, K - 1)
        A[((i - 1) * (K - 1) + 1):(i * (K - 1)), i + 1] <- -1
        A[(i - 1) * K + 1, i + 1] <- -2
        A[((i - 1) * (K - 1) + 1):(i * (K - 1)), (i * K + 2):((i + 1) *
            K)] <- iden
        A[((i - 1) * (K - 1) + 1):(i * (K - 1)), c(-(1:(K + 1)), -((i *
            K + 2):((i + 1) * K)))] <- 0
    }
    pcol <- rep(0, K * (K - 1)/2)
    t <- 1
    for (i in 1:(K - 1)) {
        for (j in 1:i) {
            pcol[t] <- i * K + j
            t <- t + 1
        }
    }
    B <- matrix(nrow = (K - 1)^2, ncol = K * (K + 1)/2)
    pcol.new <- c(pcol, 0)
    i <- 1
    ColB <- 1
    for (ColA in 1:K^2) {
        if (ColA == pcol.new[i])
            i <- i + 1 else {
            B[, ColB] <- A[, ColA]
            ColB <- ColB + 1
        }
    }
    prow <- rep(0, (K - 2) * (K - 1)/2)
    t <- 1
    for (i in 1:(K - 2)) {
        for (j in 1:i) {
            prow[t] <- i * (K - 1) + j
            t <- t + 1
        }
    }
    Cmat <- matrix(nrow = K * (K - 1)/2, ncol = K * (K + 1)/2)
    prow.new <- c(prow, 0)
    i <- 1
    RowC <- 1
    for (RowB in 1:(K - 1)^2) {
        if (RowB == prow.new[i])
            i <- i + 1 else {
            Cmat[RowC, ] <- B[RowB, ]
            RowC <- RowC + 1
        }
    }
    Q <- 0
    for (i in 1:(K * (K + 1)/2)) {
        if (b.dis[i] != 0)
            Q <- cbind(Q, Cmat[, i])
    }
    A.aug <- Q[, -1]
    Kstar <- length(b)
    ## this part here doesnt make sense: solve expects a to be square matrix
#    if (nrow(A.aug) != ncol(A.aug)) {
#        est.vect <- solve(A.aug, c.dis)
#    } else if (ncol(A.aug) == Kstar) {
#        est.vect <- (solve(t(A.aug) %*% A.aug)) %*% t(A.aug) %*% c.dis
#    } else stop("SigmaK is not identifiable\n")
    est.vect <- (solve(t(A.aug) %*% A.aug)) %*% t(A.aug) %*% c.dis
    ## general case est.SigmaK
    est.SigmaK <- matrix(nrow = K, ncol = K)
    t <- 1
    for (i in 1:K) {
        for (j in i:K) {
            if (pattern[i, j] == 0) {
                est.SigmaK[i, j] <- 0
                est.SigmaK[j, i] <- 0
            } else {
                est.SigmaK[i, j] <- est.vect[t]
                est.SigmaK[j, i] <- est.vect[t]
                t <- t + 1
            }
        }
    }
    dimnames(est.SigmaK) <- dimnames(est.SigmaKstar)
    est.SigmaK
}

