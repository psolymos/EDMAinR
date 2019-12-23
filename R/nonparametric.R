## nonparametric estimate of centered mean form and SigmaKstar
## X: data frame input, n tables (each K x D) stacked
## n: number of individuals (replicates)
## K: number of landmarks
## D: number of dimensions (2-3)
.edma_fit_np <- function(X, n, K, D) {
    ## calculate Euclidean distance for X
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
    ## calculate Eu(M)
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
    ## calculate B(M)
    I <- diag(1, K)
    ones <- array(rep(1, K), c(1, K))
    H <- I - (1/K) * crossprod(ones, ones)
    BM <- (-1/2) * H %*% EuM %*% H
    ## calculate eigenvalues and eigenvectors of B(M)
    EIG <- eigen(BM)
    ## estimate centred mean form M
    est.M <- matrix(nrow = K, ncol = D)
    for (i in 1:D) {
        est.M[, i] <- sqrt(EIG$values[i]) * EIG$vectors[, i]
    }
    ## estimate SigmaKstar=H*SigmaK*H
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

    ## check rank of SigmaKstar
    r <- qr(est.SigmaKstar)$rank
    if (r != K-1)
        warning(sprintf("rank of SigmaKstar was %s instead of %s", r, K-1))

    ## meanform eigenvalues: first D values > 0, rest =0
    mds <- cmdscale(dist(est.M), k=D, eig=TRUE)
    if (!all(mds$eig[seq_len(D)] > 0))
        warning(sprintf("%s of the first %s eigenvalues <= 0",
            sum(mds$eig[seq_len(D)] > 0), D))
    if (sum(mds$eig[seq_len(D)]) / sum(mds$eig) < 0.99)
        warning(sprintf("the first %s eigenvalues summed to %s",
            sum(mds$eig[seq_len(D)]), D))

    out <- list(
        M=est.M,
        SigmaKstar=est.SigmaKstar,
        H=H)
    out
}

## x is edma_data object
## use same pbapply setup as in opticut with cl arg etc
## bootstrap here is used for M uncertainty etc and not for T-test
edma_fit <- function(x, B=0) {
    if (!inherits(x, "edma_data"))
        stop("x must be of class edma_data")
    if (B < 0)
        stop("B must be non-negative")
    B <- as.integer(B)
    DIM <- dim(x)
    fit <- .edma_fit_np(stack(x), DIM[3L], DIM[1L], DIM[2L])
    dimnames(fit$M) <- dimnames(x)[1:2]
    dimnames(fit$SigmaKstar) <- dimnames(x)[c(1,1)]
    if (B > 0) {
        boot <- pbapply::pblapply(seq_len(B), function(i) {
            j <- sample(DIM[3L], replace=TRUE)
            z <- subset(x, j)
            .edma_fit_np(stack(z), DIM[3L], DIM[1L], DIM[2L])
        })
    } else {
        boot <- NULL
    }
    out <- c(x, fit)
    out$call <- match.call()
    out$boot <- boot
    class(out) <- c("edma_fit_np", "edma_fit", "edma_data")
    out
}

print.edma_fit_np <- function(x, ...) {
    cat("EDMA nonparametric fit: ", x$name, "\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        ncol(x$data[[1L]]), " dimensions, ",
        nrow(x$data[[1L]]), " landmarks, ",
        length(x$data), " replicates, ",
        if (length(x$boot))
            paste(length(x$boot), "bootstrap runs") else "no bootstrap",
        sep="")
    invisible(x)
}

Meanform <- function (object, ...) UseMethod("Meanform")
Meanform.edma_fit <- function (object, ...) object$M

SigmaKstar <- function (object, ...) UseMethod("SigmaKstar")
SigmaKstar.edma_fit_np <- function (object, ...) object[["SigmaKstar"]]

as.dist.edma_fit <- function(m, diag = FALSE, upper = FALSE) {
    out <- dist(Meanform(m), diag=diag, upper=upper)
    class(out) <- c("edma_dist", class(out))
    out
}

stack.dist <- function(x, ...) {
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

## form matrix, stacked
get_fm <- function (object, ...) UseMethod("get_fm")
get_fm.edma_fit <- function (object, sort=FALSE, ...) {
    d <- as.dist(object, diag = FALSE, upper = FALSE)
    out <- stack(d)
    if (sort)
        out <- out[order(out$dist, ...),]
    class(out) <- c("fm", class(out))
    out
}

.get_data <- function(object) {
    object <- list(name=object$name, data=object$data)
    class(object) <- "edma_data"
    object
}

landmarks.edma_fit <- function(x, ...)
    landmarks(.get_data(x))
