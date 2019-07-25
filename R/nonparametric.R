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
    out <- list(
        M=est.M,
        SigmaKstar=est.SigmaKstar,
        H=H)
    out
}

## x is edma_data object
## use same pbapply setup as in opticut with cl arg etc
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
    class(out) <- c("edma_fit", "edma_fit_np")
    out
}

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
sort=FALSE, level=0.95, ...) {
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
    attr(out, "T-test") <- edma_test(numerator, denominator)
    attr(out, "level") <- level
    out$lower <- q[,1L]
    out$upper <- q[,2L]
    out$inside <- out$dist <= out$upper & out$dist >= out$lower
    if (sort)
        out <- out[order(out$dist, ...),]
    out
}

