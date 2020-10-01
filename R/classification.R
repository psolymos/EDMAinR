## calculates phi and ncp matrices from mean form M and data array A
.get_mat <- function(M, A) {
    z <- .Eu2(A)
    delta_mat <- as.matrix(dist(M))^2
    phi_mat <- (1/ncol(M)) * (z$EuMean - delta_mat)
    ncp_mat <- delta_mat / phi_mat
    list(phi=phi_mat, ncp=ncp_mat)
}

## pairs of pairs of distances
.double_pairs <- function(A) {
    n <- dim(A)[3L]
    EuX <- .Eu2(A)$EuX
    tmp <- EuX[,,1]
    tmp <- tmp[lower.tri(tmp)]

    ind <- seq(1:length(tmp))
    index.mat <- expand.grid(ind,ind)
    index.mat <- subset(index.mat, index.mat[,1] < index.mat[,2])

    tmp2 <- cbind(tmp[index.mat[,1]],tmp[index.mat[,2]])

    tmp3 <- array(0, c(nrow(tmp2),2,n))
    tmp3[,,1] <- tmp2

    for (i in seq_len(n)[-1]){
        tmp <- EuX[,,i]
        tmp <- tmp[lower.tri(tmp)]
        tmp2 <- cbind(tmp[index.mat[,1]],tmp[index.mat[,2]])
        tmp3[,,i] <- tmp2
    }
    tmp3
}

## means and covariances for bivariate normal
.double_mean_cov <- function(x) {
    m <- dim(x)[1L]
    mean_mat <- matrix(0, 2L, m)
    cov_array <- array(0, c(2L, 2L, m))
    for (i in seq_len(m)) {
      mean_mat[,i] <- rowMeans(x[i,,])
      cov_array[,,i] <- cov(t(x[i,,]))
    }
    list(means=mean_mat, covariances=cov_array)
}

## calculates composite likelihood for ith bootstrap sample
.get_cl <- function(x, fit, i=0, method=c("cip", "chisq", "bnorm", "norm")) {
    method <- match.arg(method)
    if (i > 0) {
        if (is.null(fit$boot) || length(fit$boot) < i)
            stop("not enough bootstrap samples")
        M <- fit$boot[[i]]$M
        s <- attr(fit$boot, "samples")[,i]
        A <- as.array(as.edma_data(fit)[,,s])
    } else {
        M <- Meanform(fit)
        A <- as.array(as.edma_data(fit))
    }
    K <- dim(A)[1L]
    n <- dim(A)[3L]
    if (method == "chisq") {
        mat <- .get_mat(M, A)
        eu <- as.matrix(dist(x)^2)
        eu_scaled <- eu / mat$phi
        CL <- sum(dchisq(
            x=as.numeric(as.dist(eu_scaled)),
            df=ncol(x),
            ncp=pmax(0, as.numeric(as.dist(mat$ncp))),
            log=TRUE))
    }
    if (method == "bnorm") {
        eu <- .double_pairs(array(x, c(dim(x), 1L)))[,,1L]
        mc <- .double_mean_cov(.double_pairs(A))
        CL <- sum(sapply(seq_len(nrow(eu)), function(i) {
            mvtnorm::dmvnorm(
                x=eu[i,],
                mean=mc$means[,i],
                sigma=mc$covariances[,,i],
                log=TRUE)
        }))
    }
    if (method == "norm") {
        dm1 <- as.numeric(dist(x))^2
        dm <- sapply(seq_len(n), function(i) as.numeric(dist(A[,,i]))^2)
        Means <- rowMeans(dm)
        Vars <- apply(dm, 1, sd)^2
        CL <- sum(dnorm(dm1,
            mean=rowMeans(dm),
            sd=apply(dm, 1, sd),
            log=TRUE))
    }
    if (method == "cip") {
        CM <- apply(A, c(2, 3), scale, scale=FALSE)
        j <- !upper.tri(matrix(0, K, K))
        CIP <- sapply(seq_len(n), function(i) {
            (CM[,,i] %*% t(CM[,,i]))[j]
        })
        Ex <- apply(CIP, 1L, mean)
        Va <- apply(CIP, 1L, var)
        Cnew <- scale(x, scale=FALSE)
        CIPnew <- (Cnew %*% t(Cnew))[j]
        d <- CIPnew - Ex
        # taking the negative here to align with other CL types
        CL <- -as.numeric(t(d) %*% diag(Va^-1) %*% d + sum(log(Va)))
    }
    CL
}

## composite likelihood ratio and CIP
edma_class <- function(x, fit1, fit2, boot=FALSE,
method=c("cip", "chisq", "bnorm", "norm")) {
    if (inherits(x, "edma_data")) {
        if (dim(x)[3] > 1L)
            stop("provide a single specimen only when x is an edma_data object")
        x <- x$data[[1L]]
    }
    .compare_objects(fit1, fit2)
    if (!identical(dimnames(fit1)[1:2], dimnames(x)))
        stop("specimen dimnames must match the EDMA fit objects")
    CL1 <- .get_cl(x, fit1, method=method)
    CL2 <- .get_cl(x, fit2, method=method)
    BOOT <- NULL
    B1 <- length(fit1$boot)
    B2 <- length(fit2$boot)
    if (boot && (B1 < 1L || B2 < 1L))
        warning("no bootstrap samples found: boot=TRUE ignored")
    if (boot && B1 > 0L && B2 > 0L) {
        CL1M <- sapply(seq_len(B1), function(i) .get_cl(x, fit1, i, method=method))
        CL2M <- sapply(seq_len(B2), function(i) .get_cl(x, fit2, i, method=method))
        CLRM <- t(outer(CL2M, CL1M, "-"))
        BOOT <- data.frame(
            complik1=as.numeric(CL1M),
            complik2=as.numeric(CL2M),
            complikr=as.numeric(CLRM))
    }
    ## CLR > 0 => x belongs to fit2
    ## CLR < 0 => x belongs to fit1
    out <- list(
        x=x,
        fit1=fit1,
        fit2=fit2,
        complik1=CL1,
        complik2=CL2,
        complikr=CL2 - CL1,
        boot=BOOT)
    class(out) <- "edma_class"
    out
}

## leave-one-out cross validation
loo <- function(fit1, fit2, B=0, level=0.95,
method=c("cip", "chisq", "bnorm", "norm")) {
    x1 <- as.edma_data(fit1)
    x2 <- as.edma_data(fit2)
    n1 <- dim(x1)[3L]
    n2 <- dim(x2)[3L]
    x12 <- combine_data(x1, x2)
    a <- c((1-level)/2, 1-(1-level)/2)
    Class <- data.frame(
        group=rep(1:2, c(n1, n2)),
        specimen=c(seq_len(n1), seq_len(n2)),
        id=seq_len(n1 + n2),
        complikr=NA, lower=NA, upper=NA)
    for (i in seq_len(n1 + n2)) {
        x <- x12$data[[i]]
        i1 <- which(Class$group == 1 & Class$id != i)
        i2 <- which(Class$group == 2 & Class$id != i)
        fit1 <- edma_fit(x12[,,i1], B=B)
        fit2 <- edma_fit(x12[,,i2], B=B)
        h <- edma_class(x, fit1, fit2, boot=B > 0, method=method)
        Class$complikr[i] <- h$complikr
        if (B > 0 && !is.null(h$boot)) {
            q <- quantile(c(h$complikr, h$boot$complikr), a, na.rm=TRUE)
            Class$lower[i] <- q[1L]
            Class$upper[i] <- q[2L]
        }
    }
    Class$class <- ifelse(Class$complikr > 0, 2, 1)
    Class$signif <- !(Class$lower < 0 & Class$upper > 0)
    cm <- table(Class$group, Class$class)
    attr(Class, "accuracy") <- sum(diag(cm)) / sum(cm)
    class(Class) <- c("edma_loo", class(Class))
    Class
}
