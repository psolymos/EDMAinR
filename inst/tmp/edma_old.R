#devtools::install_github("psolymos/EDMAinR")
source("R/edma_fit.R")

K <- 3 # number of landmarks
D <- 2 # dimension, 2 or 3

sig <- 0.2
rho <- 0
SigmaK <- sig^2*diag(1, K, K) + sig^2*rho*(1-diag(1, K, K))

M <- matrix(c(0,1,0,0,0,1), 3, 2)
M[,1] <- M[,1] - mean(M[,1])
M[,2] <- M[,2] - mean(M[,2])
M <- 10*M

n <- 1000
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

(fit <- edma_fit(X, n, K, D))
SigmaK_fit(fit, "sig")$optim
SigmaK_fit(fit, "sig_rho")$optim

