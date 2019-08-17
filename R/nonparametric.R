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
    class(out) <- c("edma_fit_np", "edma_fit")
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
edma_fm <- function (object, ...) UseMethod("edma_fm")
edma_fm.edma_fit <- function (object, sort=FALSE, ...) {
    d <- as.dist(object, diag = FALSE, upper = FALSE)
    out <- stack(d)
    attr(out, "method") <- attr(d, "method")
    attr(out, "Tval") <- attr(d, "Tval")
    if (sort)
        out <- out[order(out$dist, ...),]
    class(out) <- c("edma_fm", class(out))
    out
}

## a and b are edma_fit objects
.compare_objects <- function (a, b, text="", ...) {
    if (!inherits(b, "edma_fit") || !inherits(a, "edma_fit"))
        stop("input must be edma_fit object")
    if (nrow(a$data[[1L]]) != nrow(b$data[[1L]]))
        stop(paste(text, "number of landmarks must be identical"))
    if (ncol(a$data[[1L]]) != ncol(b$data[[1L]]))
        stop(paste(text, "number of dimensions must be identical"))
    if (length(a$boot) != length(b$boot))
        stop(paste(text, "number of bootstrap runs  must be identical"))
    if (!all(rownames(a$data[[1L]]) == rownames(b$data[[1L]])))
        stop(paste(text, "landmark names and ordering must be identical"))
    if (!all(colnames(a$data[[1L]]) == colnames(b$data[[1L]])))
        stop(paste(text, "dimension names and ordering must be identical"))
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

## inputs are edma_fit objects
## this returns the distance matrix, not stacked
formdiff <- function (numerator, denominator, ...) {
    .compare_objects(numerator, denominator)
    .formdiff(Meanform(numerator), Meanform(denominator))
}

## T-test for 2 edma_fit objects
edma_test <- function (numerator, denominator) {
    .compare_objects(numerator, denominator)
    DNAME <- paste(deparse(substitute(numerator)),
        deparse(substitute(denominator)), sep = " / "))
    METHOD <- "Bootstrap based EDMA T-test"
    B <- length(numerator$boot)
    Tval <- attr(formdiff(numerator, denominator), "Tval")
    Tvals <- c(Tval, sapply(seq_len(B), function(i) {
        attr(.formdiff(numerator$boot[[i]]$M, denominator$boot[[i]]$M), "Tval")
    }))
    PVAL = sum(Tvals > Tval) / (B + 1)
    PARAMETER <- B + 1L
    names(Tval) <- "T-value"
    names(PARAMETER) <- "B+1"
    out <- list(statistic = Tval, parameter = PARAMETER,
        p.value = PVAL, method = METHOD, data.name = DNAME,
        Tvals=Tvals)
    class(out) <- c("edma_test", "htest")
    out
}

plot.edma_test <- function(x, ...) {
    hist(x$Tvals, xlab="T-values", main=x$data.name, ...)
    abline(v=x$statistic, col=2, lwd=2)
}

.sdm <- function (numerator, denominator, sort=FALSE, level=0.95, ...) {
    d <- formdiff(numerator, denominator)
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
    out$inside <- 1 <= out$upper & 1 >= out$lower
    if (sort)
        out <- out[order(out$dist, ...),]
    out
}

## stacked form difference matrix
## inputs are edma_fit objects
edma_fdm <- function (b, a, sort=FALSE, level=0.95, ...) {
    .compare_objects(b, a, "b vs a:")
    out <- .sdm(b, a, sort, level, ...)
    attr(out, "numerator") <- deparse(substitute(b))
    attr(out, "denominator") <- deparse(substitute(q))
    class(out) <- c("edma_fdm", class(out))
    out
}

## stacked growth matrix
## inputs are edma_fit objects
edma_gm <- function (a1, a2, sort=FALSE, level=0.95, ...) {
    .compare_objects(a1, a2, "a1 vs a2:")
    out <- .sdm(a2, a1, sort, level, ...)
    attr(out, "age1") <- deparse(substitute(a1))
    attr(out, "age2") <- deparse(substitute(a2))
    class(out) <- c("edma_gm", class(out))
    out
}

## stacked growth difference matrix
## inputs are edma_fit objects
edma_gdm <- function (a1, a2, b1, b2,
sort=FALSE, level=0.95, ...) {
    .compare_objects(a1, a2, "a1 vs a2:")
    .compare_objects(b1, b2, "b1 vs b2:")
    .compare_objects(a1, b1, "a1 vs b1:")
    .compare_objects(a2, b2, "a2 vs b2:")
    gma <- formdiff(numerator=a2, denominator=a1)
    gmb <- formdiff(numerator=b2, denominator=b1)
    out <- stack(gmb / gma)
    B <- length(a1$boot)
    fd <- sapply(seq_len(B), function(i) {
        gma <- .formdiff(a2$boot[[i]]$M, a1$boot[[i]]$M)
        gmb <- .formdiff(b2$boot[[i]]$M, b1$boot[[i]]$M)
        stack(gmb / gma)$dist
    })
    a <- c((1-level)/2, 1-(1-level)/2)
    q <- t(apply(cbind(out$dist, fd), 1, quantile, a))
    attr(out, "method") <- "growth_difference"
    attr(out, "level") <- level
    out$lower <- q[,1L]
    out$upper <- q[,2L]
    out$inside <- 1 <= out$upper & 1 >= out$lower
    if (sort)
        out <- out[order(out$dist, ...),]
    attr(out, "a1") <- deparse(substitute(a1))
    attr(out, "a2") <- deparse(substitute(a2))
    attr(out, "b1") <- deparse(substitute(b1))
    attr(out, "b2") <- deparse(substitute(b2))
    class(out) <- c("edma_gdm", class(out))
    out
}

.sdmplot <- function(x, xlab="", ylab="",
    bottom=1.5, xcex=0.5, xshow=TRUE, ...) {
    x <- x[order(x$dist),]
    k <- nrow(x)
    xv <- seq_len(k)
    r <- range(x$dist, x$lower, x$upper)
    op <- par(srt=90, xpd = TRUE, mar=par()$mar*c(bottom, 1, 1, 1))
    on.exit(par(op), add=TRUE)
    plot(xv, x$dist, ylim=r, type="n",
        xlab=xlab, ylab=ylab, axes=FALSE)
    polygon(c(xv, rev(xv)), c(x$lower, rev(x$upper)), col="lightblue", border=NA)
    abline(h=1, col="grey")
    for (i in seq_len(k))
        lines(xv[i]+c(-0.5, 0.5), c(1, 1),
            col=if (1 > x$upper[i] || 1 < x$lower[i]) "red" else "grey")
    lines(xv, x$dist, col="blue")
    axis(2)
    if (xshow) {
        lab <- paste0(as.character(x$row), "-", as.character(x$col))
        text(xv, min(r) - 0.02 * diff(r), lab, adj=c(1, 0.5), cex=xcex)
    }
    invisible(x)
}
plot.edma_fdm <- function(x, ...) {
    .sdmplot(x, ylab="FDM Ratio", ...)
}
plot.edma_gdm <- function(x, ...) {
    .sdmplot(x, ylab="GDM Ratio", ...)
}
