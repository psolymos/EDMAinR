## nonparametric estimate of centered mean form and SigmaKstar
## X: data frame input, n tables (each K x D) stacked
## n: number of individuals (replicates)
## K: number of landmarks
## D: number of dimensions (2-3)
## this _old function is here for comparison (~250% slower)
.edma_fit_np_old <- function(X, n, K, D) {
    ## calculate squared Euclidean distance for X
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
        warning(sprintf("Rank of SigmaKstar was %s instead of %s.", r, K-1))

    ## meanform eigenvalues: first D values > 0, rest =0
    mds <- cmdscale(dist(est.M), k=D, eig=TRUE)
    if (!all(mds$eig[seq_len(D)] > 0))
        warning(sprintf("%s of the first %s eigenvalues <= 0.",
            sum(mds$eig[seq_len(D)] > 0), D))
    if (sum(mds$eig[seq_len(D)]) / sum(mds$eig) < 0.99)
        warning(sprintf("The first %s eigenvalues summed to %s.",
            sum(mds$eig[seq_len(D)]), D))

    out <- list(
        M=est.M,
        SigmaKstar=est.SigmaKstar)
    out
}

## nonparametric estimate of centered mean form and SigmaKstar
## A: array input, K x D x n
## this new function makes use of vectorized base functions
## to gain speed and precistion
## also returns some intermediate results for future use cases
.Eu2 <- function(A) {
    K <- dim(A)[1L]
    D <- dim(A)[2L]
    n <- dim(A)[3L]
    ## calculate squared Euclidean distances for X
    EuX <- array(0, c(K, K, n))
    dimnames(EuX) <- dimnames(A)[c(1L,1L,3L)]
    for (i in seq_len(n)) {
        EuX[,,i] <- as.matrix(dist(A[,,i])^2)
    }
    ## calculate Eu(M)
    EuM <- matrix(nrow = K, ncol = K)
    dimnames(EuM) <- dimnames(A)[c(1L,1L)]
    Mean <- Var <- EuM
    for (i in seq_len(K)) {
        for (j in 1:K) {
            Mean[i,j] <- mean(EuX[i,j,])
            Var[i,j] <- sd(EuX[i,j,])^2
            if (D == 2L) {
                EuM[i,j] <- sqrt((Mean[i,j])^2 - Var[i,j])
            } else {
                EuM[i,j] <- sqrt((Mean[i,j])^2 - 1.5 * Var[i,j])
            }
        }
    }
    list(EuX=EuX, EuM=EuM, EuMean=Mean, EuVar=Var)
}
.edma_fit_np <- function(A, less=TRUE) {
    K <- dim(A)[1L]
    D <- dim(A)[2L]
    n <- dim(A)[3L]
    if (K <= D)
        stop("K should be larger than D for EDMA to work")

    ## calculate arrays using .Eu2()
    Eu2 <- .Eu2(A)
    EuX <- Eu2$EuX
    EuM <- Eu2$EuM
    #Mean <- Eu2$EuMean=Mean
    #Var <- Eu2$EuVar

    ## calculate B(M)
    I <- diag(1, K)
    ones <- array(rep(1, K), c(1, K))
    H <- I - (1/K) * crossprod(ones, ones)
    BM <- (-1/2) * H %*% EuM %*% H
    ## calculate eigenvalues and eigenvectors of B(M)
    EIG <- eigen(BM)
    ## estimate centred mean form M
    M_hat <- matrix(nrow = K, ncol = D)
    for (i in seq_len(D)) {
        M_hat[, i] <- sqrt(EIG$values[i]) * EIG$vectors[, i]
    }
    dimnames(M_hat) <- dimnames(A)[1:2]
    ## estimate SigmaKstar=H*SigmaK*H
    a <- 0
    for (i in seq_len(n)) {
        a <- a + (-1/2) * H %*% EuX[,,i] %*% H
    }
    SigmaKstar_hat <- (a/n - BM)/D
    dimnames(SigmaKstar_hat) <- dimnames(EuM)

    ## check rank of SigmaKstar
    r <- qr(SigmaKstar_hat)$rank
    if (r != K-1)
        warning(sprintf("Rank of SigmaKstar was %s instead of %s.", r, K-1))

    ## meanform eigenvalues: first D values > 0, rest =0
    mds <- cmdscale(dist(M_hat), k=D, eig=TRUE)
    if (!all(mds$eig[seq_len(D)] > 0))
        warning(sprintf("%s of the first %s eigenvalues <= 0.",
            sum(mds$eig[seq_len(D)] > 0), D))
    if (sum(mds$eig[seq_len(D)]) / sum(mds$eig) < 0.99)
        warning(sprintf("The first %s eigenvalues summed to %s.",
            sum(mds$eig[seq_len(D)]), D))
    ## usually we don't want to accumulate everything
    out <- list(
        M=M_hat,
        SigmaKstar=SigmaKstar_hat)
    if (!less)
        out <- c(out, Eu2)
    out
}

## x is edma_data object
## use same pbapply setup as in opticut with cl arg etc
## bootstrap here is used for M uncertainty etc and not for T-test
edma_fit <- function(x, B=0, ncores=getOption("Ncpus", 1L)) {
    if (!inherits(x, "edma_data"))
        stop("x must be of class edma_data.")
    if (B < 0)
        stop("B must not be negative.")
    B <- as.integer(B)
    DIM <- dim(x)
    fit <- .edma_fit_np(as.array(x))
    dimnames(fit$M) <- dimnames(x)[1:2]
    dimnames(fit$SigmaKstar) <- dimnames(x)[c(1,1)]
    if (B > 0) {
        #fun <- .edma_fit_np
        BB <- replicate(B, sample(DIM[3L], replace=TRUE))
        A <- as.array(x)
        ncores <- as.integer(ncores)
        if (ncores > 1L) {
            ncores <- min(ncores, parallel::detectCores(TRUE), na.rm=TRUE)
            cl <- parallel::makeCluster(ncores)
            on.exit(parallel::stopCluster(cl))
            #parallel::clusterExport(cl, c("fun"))
            #parallel::clusterEvalQ(cl, library(EDMAinR))
        } else {
            cl <- NULL
        }
        boot <- pbapply::pbapply(BB, 2L, function(j, A) {
            EDMAinR::.edma_fit_np(A[,,j])
        }, A=A, cl=cl)
        attr(boot, "samples") <- BB
    } else {
        boot <- NULL
    }
    out <- c(x, fit)
    out$call <- match.call()
    out$boot <- boot
    class(out) <- c("edma_fit_np", "edma_fit", "edma_data")
    out
}

## this fits GPA or weighted GPA to the 3d array object
.gpa_fit <- function(A, ..., weighted=FALSE, gpa_results=FALSE) {
    SigmaKstar_hat <- matrix(NA, dim(A)[1L], dim(A)[1L])
    dimnames(SigmaKstar_hat) <- list(dimnames(A)[[1L]], dimnames(A)[[1L]])
    m <- if (weighted)
        shapes::procWGPA(A, ...) else shapes::procGPA(A, ...)
    m$weighted <- weighted
    M_hat <- as.matrix(m$mshape)
    dimnames(M_hat) <- dimnames(A)[1:2]
    out <- list(M=M_hat,
        SigmaKstar=SigmaKstar_hat)
    if (gpa_results)
        out$gpa_results <- m
    out
}
gpa_fit <- function(x, B=0, ncores=getOption("Ncpus", 1L), weighted=FALSE, ...) {
    if (!requireNamespace("shapes"))
        stop("The shapes package is required for GPA, use install.packages('shapes')")
    if (!inherits(x, "edma_data"))
        stop("x must be of class edma_data.")
    if (B < 0)
        stop("B must not be negative.")
    B <- as.integer(B)
    DIM <- dim(x)
    fit <- .gpa_fit(as.array(x), ..., weighted=weighted)
    dimnames(fit$M) <- dimnames(x)[1:2]
    dimnames(fit$SigmaKstar) <- dimnames(x)[c(1,1)]
    if (B > 0) {
        BB <- replicate(B, sample(DIM[3L], replace=TRUE))
        A <- as.array(x)
        ncores <- as.integer(ncores)
        if (ncores > 1L) {
            ncores <- min(ncores, parallel::detectCores(TRUE), na.rm=TRUE)
            cl <- parallel::makeCluster(ncores)
            on.exit(parallel::stopCluster(cl))
        } else {
            cl <- NULL
        }
        boot <- pbapply::pbapply(BB, 2L, function(j, A, ..., weighted) {
            EDMAinR::.gpa_fit(A[,,j], ..., weighted=weighted)
        }, A=A, cl=cl, ..., weighted=weighted)
        attr(boot, "samples") <- BB
    } else {
        boot <- NULL
    }
    out <- c(x, fit)
    out$call <- match.call()
    out$boot <- boot
    out$weighted <- weighted
    class(out) <- c("gpa_fit", "edma_fit_np", "edma_fit", "edma_data")
    out
}


## print methods for np fit object
print.edma_fit_np <- function(x, truncate=40, ...) {
    cat("EDMA nonparametric fit: ",
        .shorten_name(x$name, truncate), "\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        nrow(x$data[[1L]]), " landmarks, ",
        ncol(x$data[[1L]]), " dimensions, ",
        length(x$data), " specimens, ",
        if (length(x$boot))
            paste(length(x$boot), "bootstrap runs") else "no bootstrap",
        sep="")
    invisible(x)
}

## print methods for GPA fit object
print.gpa_fit <- function(x, truncate=40, ...) {
    cat("GPA fit (mean form only): ",
        .shorten_name(x$name, truncate), "\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        nrow(x$data[[1L]]), " landmarks, ",
        ncol(x$data[[1L]]), " dimensions, ",
        length(x$data), " specimens, ",
        if (length(x$boot))
            paste(length(x$boot), "bootstrap runs") else "no bootstrap",
        sep="")
    invisible(x)
}

## mean form
Meanform <- function (object, ...) UseMethod("Meanform")
Meanform.edma_fit <- function (object, ...) object$M

## SigmaKstar
SigmaKstar <- function (object, ...) UseMethod("SigmaKstar")
SigmaKstar.edma_fit <- function (object, ...) object[["SigmaKstar"]]

## KxK distance matrix based on mean form
as.dist.edma_fit <- function(m, diag = FALSE, upper = FALSE) {
    out <- dist(Meanform(m), diag=diag, upper=upper)
    class(out) <- c("edma_dist", class(out))
    out
}

## stacked distances based on an input object of class dist
## this is a generic method
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

## CI based on the 1 input object boot sample
confint.edma_fit <- function (object, parm, level=0.95, ...) {
    d <- stack(as.dist(object))
    if (missing(parm))
        parm <- seq_len(nrow(d))
    a <- c((1-level)/2, 1-(1-level)/2)
    b <- if (is.null(object$boot)) {
        data.matrix(d$dist)
    } else {
        cbind(d$dist,
            sapply(object$boot, function(z) stack(dist(z$M))$dist))
    }
    out <- t(apply(b, 1, quantile, a))
    if (is.null(object$boot))
        out[] <- NA
    rownames(out) <- paste0(as.character(d$row), "-", as.character(d$col))
    out[parm,,drop=FALSE]
}

## form matrix, stacked distances from mean form
## this is the intended EDMA interface
get_fm <- function (object, ...) UseMethod("get_fm")
get_fm.edma_fit <- function (object, sort=FALSE,
level = 0.95, ...) {
    d <- as.dist(object, diag = FALSE, upper = FALSE)
    out <- stack(d)
    ci <- confint(object, level=level)
    out$lower <- ci[,1L]
    out$upper <- ci[,2L]
    if (sort)
        out <- out[order(out$dist, ...),]
    class(out) <- c("fm", class(out))
    attr(out, "level") <- level
    out
}

## reduces a fit object to a data object
.get_data <- function(object) {
    object <- list(name=object$name, data=object$data)
    class(object) <- "edma_data"
    object
}

as.edma_data.edma_fit <- function(x, ...) .get_data(x)

