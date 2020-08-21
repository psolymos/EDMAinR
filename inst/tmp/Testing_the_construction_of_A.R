# Testing the program

K = 6

H = diag(rep(1,K),K,K) - matrix(rep(1/K,K^2),K,K)

L = cbind(rep(-1,K-1),diag(1,K-1))

L %*% H - L     # Testing if LH = L.

# Construct a true Sigma_K with appropriate number of 0s. This should be positive definite, symmetric, square matrix.

G = matrix(rnorm(K^2),c(K,K))
Sigma_K = t(G) %*% G

# Equate some of the covariances to 0. This should not affect the positive definiteness.

Sigma_K[1,2] = 0
Sigma_K[2,1] = 0

Sigma_K[1,3] = 0
Sigma_K[3,1] = 0

Sigma_K[1,4] = 0
Sigma_K[4,1] = 0

# Sigma_K[2,3] = 0
# Sigma_K[3,2] = 0

Sigma_K[2,4] = 0
Sigma_K[4,2] = 0

Sigma_K[3,4] = 0
Sigma_K[4,3] = 0

if (K > 4) {
    Sigma_K[5,1:4] <- Sigma_K[1:4, 5] <- 0
    Sigma_K[6,1:4] <- Sigma_K[1:4, 6] <- 0
    Sigma_K[6,5] <- Sigma_K[5,6] <- 0
}

pattern <- ifelse(Sigma_K != 0, Sigma_K, NA)

Sigma_K
IDs <- ifelse(Sigma_K != 0, matrix(seq_len(length(Sigma_K)), K, K), NA)
IDs <- IDs[lower.tri(IDs, diag=TRUE)]


# Sigma_K_curl

Sigma_K_curl = L %*% Sigma_K %*% t(L)

Sigma_K_curl_vec = Sigma_K_curl[lower.tri(Sigma_K_curl,diag=T)]

# Now use the constructed A to reconstruct the true Sigma_K.


L = cbind(rep(-1,K-1),diag(1,K-1))

A.0 = kronecker(L,L)

# Covariance matrix Simga_K_curl is symmetric. There are duplicate entries corresponding to the covariances. We need to remove the corresponding rows to reduce the dimension of the system of equations. So the row dimension corresponds to the number of quantities available from the estimator of the singular matrix.

Index.matrix0 = matrix(seq(1:(K-1)^2),c(K-1,K-1))
Index.vec0 = 0
for (i in 1:K-1){
	for (j in 1:K-1){
		if (i < j){Index.vec0 = c(Index.vec0,Index.matrix0[i,j])}
	}
}

Index.vec0 = Index.vec0[-1]

A.1 = A.0[-Index.vec0,]

# Now we use the symmetry in Sigma_K to reduce the column dimension. We need to add the corresponding columns. For example, [1,2] and [2,1] are identical so the corresponding columns should be added and one of them need to removed.

Index.matrix1 = matrix(seq(1:K^2),c(K,K))
Index.mat1 = c(0,0)
for (i in 1:K){
	for (j in 1:K){
		if (i < j){Index.mat1 = rbind(Index.mat1,c(Index.matrix1[i,j],Index.matrix1[j,i]))}
	}
}

Index.mat1 = Index.mat1[-1,]

A.2 = A.1
A.2[,Index.mat1[,2]] = A.1[,Index.mat1[,2]] + A.1[,Index.mat1[,1]]
A.2 = A.2[,-Index.mat1[,1]]

# Now remove the columns that correspond to the structural 0s in the model covariance matrix. The indices for the model covariance matrix should go from 1 to K(K+1)/2. This should equal the number of columns in A.2

check = ncol(A.2) - K*(K+1)/2
if (check == 0){cat("Good to go")}

# Find the indices of the non-zero entries in the model Sigma_K in vectorized form. Order is important here.

which(Sigma_K[lower.tri(Sigma_K,diag=TRUE)] != 0)
which(!is.na(pattern[lower.tri(pattern,diag=TRUE)]))
Index.nonzero = which(Sigma_K[lower.tri(Sigma_K,diag=T)] != 0)   # This changes as per the model for covariance. You will need to change the non-equality to 0 to non-equality to "NA".

A = A.2[,Index.nonzero]

Sigma_hat = solve(t(A) %*% A) %*% t(A) %*% Sigma_K_curl_vec

# This Sigma_hat should have the same entries as Sigma. Need to check the order.

Sigma_K[IDs[!is.na(IDs)]]
drop(Sigma_hat)

SigmaKhat <- matrix(0, K, K)
SigmaKhat[IDs[!is.na(IDs)]] <- Sigma_hat
SigmaKhat <- t(SigmaKhat)
SigmaKhat[IDs[!is.na(IDs)]] <- Sigma_hat

Sigma_K
SigmaKhat

max(abs(SigmaKhat - Sigma_K))



## --------------

.make_A <- function(K, pattern) {
    if (K < 3)
        stop("K must be at least 3")
    K <- as.integer(K)
    L <- cbind(rep(-1, K-1), diag(1, K-1, K-1))
    A.0 <- kronecker(L, L)
    Index.matrix0 <- matrix(seq(1:(K-1)^2),c(K-1,K-1))
    Index.vec0 <- 0
    for (i in 1:K-1){
        for (j in 1:K-1){
            if (i < j) {
                Index.vec0 <- c(Index.vec0, Index.matrix0[i,j])
            }
        }
    }
    Index.vec0 <- Index.vec0[-1]
    A.1 <- A.0[-Index.vec0,]
    Index.matrix1 <- matrix(seq(1:K^2),c(K,K))
    Index.mat1 <- c(0, 0)
    for (i in 1:K){
        for (j in 1:K){
            if (i < j) {
                Index.mat1 <- rbind(
                    Index.mat1,
                    c(Index.matrix1[i,j], Index.matrix1[j,i]))
            }
        }
    }
    Index.mat1 <- Index.mat1[-1,]
    A.2 <- A.1
    A.2[,Index.mat1[,2]] <- A.1[,Index.mat1[,2]] + A.1[,Index.mat1[,1]]
    A.2 <- A.2[,-Index.mat1[,1]]
    if (ncol(A.2) != K*(K+1)/2)
        stop("Something went wrong with constructing A.")
    Index.nonzero <- which(!is.na(pattern[lower.tri(pattern,diag=TRUE)]))
    A <- A.2[,Index.nonzero]
    A
}

.SigmaK_fit_full <- function(SigmaKstar, pattern) {
    K <- nrow(SigmaKs)
    L <-  cbind(rep(-1, K-1), diag(1, K-1, K-1))
    A <- .make_A(K, pattern)

    IDs <- ifelse(!is.na(pattern), matrix(seq_len(length(pattern)), K, K), NA)
    IDs <- IDs[lower.tri(IDs, diag=TRUE)]
    IDs <- IDs[!is.na(IDs)]

    SigmaKc <- L %*% SigmaKstar %*% t(L)
    SigmaKcvec = SigmaKc[lower.tri(SigmaKc, diag=TRUE)]
    SigmaKhat1 = solve(t(A) %*% A) %*% t(A) %*% SigmaKcvec

    SigmaKhat <- matrix(0, K, K)
    SigmaKhat[IDs] <- SigmaKhat1
    SigmaKhat <- t(SigmaKhat)
    SigmaKhat[IDs] <- SigmaKhat1
    dimnames(SigmaKhat) <- dimnames(SigmaKstar)
    SigmaKhat
}

.SigmaK_fit_new <- function(SigmaKstar, pattern, init,
method = "Nelder-Mead", control = list()) {
    K <- nrow(SigmaKstar)
    IDs <- ifelse(!is.na(pattern), matrix(seq_len(length(pattern)), K, K), NA)
    IDs <- IDs[lower.tri(IDs, diag=TRUE)]
    IDs <- IDs[!is.na(IDs)]
    IDd <- which(diag(1, K, K) > 0)
    isDiag <- IDs %in% IDd
    names(IDs) <- pattern[IDs]
    uni <- !duplicated(names(IDs))
    attr(IDs, "diag") <- isDiag
    attr(IDs, "unique") <- uni
    isDiagUni <- isDiag[uni]

    SigmaKfull <- .SigmaK_fit_full(SigmaKstar, pattern)
    parms_full <- SigmaKfull[IDs]
    parms_uni <- parms_full[uni]
    for (i in names(IDs)[uni])
        parms_uni[names(IDs)[uni] == i] <- mean(parms_full[names(IDs) == i])
    if (missing(init))
        init <- parms_uni
    if (length(init) != length(parms_uni))
        stop(sprintf("init length must be %s", length(parms_uni)))
    num_max <- .Machine$double.xmax^(1/3)
    fun <- function(parms) {
        if (any(parms[isDiagUni] <= 0))
            return(num_max)
        sum((parms - parms_uni)^2)
    }
    if (!is.null(control$fnscale) && control$fnscale < 0)
        stop("control$fnscale can not be negative")
    o <- suppressWarnings({
        optim(init, fun, method=method, control=control, hessian=FALSE)
    })

    parms_back <- parms_full
    for (i in names(IDs)[uni])
        parms_back[names(IDs) == i] <- o$par[names(IDs)[uni] == i]

    SigmaKhat <- matrix(0, K, K)
    SigmaKhat[IDs] <- parms_back
    SigmaKhat <- t(SigmaKhat)
    SigmaKhat[IDs] <- parms_back
    dimnames(SigmaKhat) <- dimnames(SigmaKfull)

    o$init <- init0
    o$SigmaK <- SigmaKhat
    o$SigmaKfull <- SigmaKfull
    o$method <- method
    o$control <- control
    o$id <- IDs
    o
}


## SigmaK
S2 <- matrix(
  c("s1", "r1",  NA, NA,  NA, NA,
    "r1", "s2", NA, NA,  NA, NA,
    NA,  NA, "s1", NA,  NA, NA,
    NA,  NA,  NA, "s3", "r2", NA,
    NA,  NA,  NA, "r2", "s3", NA,
    NA,  NA,  NA, NA, NA, "s3"),
  nrow=6, ncol=6, byrow=TRUE)
dimnames(S2) <- list(rownames(M), rownames(M))
parm2 <- c("s1"=0.2, "s2"=4, "s3"=4,"r1"=0.22, "r2"=0.33)
SigmaK2 <- EDMAinR:::.vec2mat(parm2, EDMAinR:::.mat2fac(S2))
dimnames(SigmaK2) <- dimnames(S2)

n <- 1000
sim0 <- edma_simulate_data(n=n, M, SigmaK2)
fit <- edma_fit(sim0)

SigmaKhat1 <- EDMAinR:::.SigmaK_fit_full(SigmaKstar(fit), S2)
old <- EDMAinR:::.SigmaK_fit_old(SigmaKstar(fit), fit$H, S2)
est <- SigmaK_fit(fit, S2)
SigmaKhat3 <- est$SigmaK
id <- est$results$id

SigmaK2
round(SigmaKhat1, 4)
round(SigmaKhat3, 4)
max(abs(SigmaKhat1 - SigmaK2))
max(abs(SigmaKhat3 - SigmaK2))

df <- data.frame(id=id, parm=names(id),
    true=SigmaK2[id],
    full=SigmaKhat1[id],
    const=SigmaKhat3[id])
df[order(df$parm),]
