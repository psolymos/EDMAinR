other_fun <- function(z) {
    ## Use maple to solve Y such that YH=L
    fcol <- array(c(rep(-1, K - 2), -2), c(K - 1, 1))
    identity <- diag(1, K - 2, K - 1)
    augment <- array(c(rep(-1, K - 2), 0), c(1, K - 1))
    Y <- cbind(fcol, rbind(identity, augment))
    ## from est.SigmaKstar to est.SigmaK~
    est.SigmaKtilde <- Y %*% est.SigmaKstar %*% t(Y)
    L <- cbind(c(rep(-1, K - 1)), diag(1, K - 1))
    ## L%*%SigmaK%*%t(L) # model for diagonal matrix
    c <- as.vector(est.SigmaKtilde)
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
    est.vect <- round((solve(t(A) %*% A)) %*% t(A) %*% c, 3)
    est.SigmaK <- diag(c(est.vect))
    ## general case
    c <- as.vector(est.SigmaKtilde)
    b.ori <- as.vector(SigmaK)
    m <- 0
    b.dis <- rep(0, K * (K + 1)/2)
    for (i in 1:K) {
        for (j in i:K) {
            m <- m + 1
            b.dis[m] <- SigmaK[i, j]
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
    if (nrow(A.aug) != ncol(A.aug)) {
        est.vect <- solve(A.aug, c.dis) # this is where rho has issues
    } else if (ncol(A.aug) == Kstar) {
        est.vect <- (solve(t(A.aug) %*% A.aug)) %*% t(A.aug) %*% c.dis
    } else stop(" SigmaK is not identifiable\n")
    ## general case est.SigmaK
    est.SigmaK <- matrix(nrow = K, ncol = K)
    t <- 1
    for (i in 1:K) {
        for (j in i:K) {
            if (SigmaK[i, j] == 0) {
                est.SigmaK[i, j] <- 0
                est.SigmaK[j, i] <- 0
            } else {
                est.SigmaK[i, j] <- est.vect[t]
                est.SigmaK[j, i] <- est.vect[t]
                t <- t + 1
            }
        }
    }
    list(M=est.M, SigmaK=est.SigmaK, H=H, D=D, K=K)
}

SigmaK_fit <- function(fit,
type=c("sig", "sig_rho"),
method = "Nelder-Mead", control = list(), hessian = FALSE)
{
    type <- match.arg(type)
    v <- diag(1, fit$K, fit$K)
    if (type == "sig") {
        init <- log(sqrt(mean(diag(fit$SigmaK))))
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
        init <- c(log(sqrt(mean(diag(fit$SigmaK)))),
                  atanh(mean(fit$SigmaK[lower.tri(fit$SigmaK)])))
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
