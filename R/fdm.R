## a and b are edma_fit objects
.compare_data <- function (a, b, text="", ...) {
    if (nrow(a$data[[1L]]) != nrow(b$data[[1L]]))
        stop(paste(text, "number of landmarks must be identical"))
    if (ncol(a$data[[1L]]) != ncol(b$data[[1L]]))
        stop(paste(text, "number of dimensions must be identical"))
    if (!all(rownames(a$data[[1L]]) == rownames(b$data[[1L]])))
        stop(paste(text, "landmark names and ordering must be identical"))
    if (!all(colnames(a$data[[1L]]) == colnames(b$data[[1L]])))
        stop(paste(text, "dimension names and ordering must be identical"))
    invisible(NULL)
}
.compare_objects <- function (a, b, text="", ...) {
    if (!inherits(b, "edma_fit") || !inherits(a, "edma_fit"))
        stop("input must be edma_fit object")
    .compare_data(a, b, text, ...)
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
## random mixed bootstrap index based on 2 objects, singe draw
.Test_index1 <- function(n1, n2, sample=TRUE, replace=TRUE, ...) {
    ii <- c(-seq_len(n1), seq_len(n2))
    if (sample)
        ii <- sample(ii, length(ii), replace=replace, ...)
    ii
}
## random mixed bootstrap index based on 2 objects, multiple
.Ttest_index <- function(numerator, denominator, B=0) {
    n1 <- length(numerator$data)
    n2 <- length(denominator$data)
    v0 <- .Test_index1(n1, n2, sample=FALSE)
    out <- if (B > 0) {
        cbind(v0, replicate(B, .Test_index1(n1, n2, sample=TRUE)))
    } else {
        data.matrix(v0)
    }
    unname(out)
}
## make a fit object like input data
.get_data <- function(object) {
    object <- list(name=object$name, data=object$data)
    class(object) <- "edma_data"
    object
}
## form diff matrix based on mixed bootstrap, single draw
.Ttest_fit1 <- function(i, d1, d2) {
    n1 <- length(d1$data)
    n2 <- length(d2$data)
    i1 <- i[seq_len(n1)]
    i2 <- i[seq_len(n2)+n1]
    d1o <- c(d1$data[-i1[i1 < 0]], d2$data[i1[i1 > 0]])
    d2o <- c(d1$data[-i2[i2 < 0]], d2$data[i2[i2 > 0]])
    d1$data <- d1o
    d2$data <- d2o
    f1 <- edma_fit(d1)
    f2 <- edma_fit(d2)
    .formdiff(Meanform(f1), Meanform(f2))
}
## form diff matrix based on mixed bootstrap, multiple runs
.Ttest_fit <- function(numerator, denominator, B=0) {
    d1 <- .get_data(numerator)
    d2 <- .get_data(denominator)
    ii <- .Ttest_index(numerator, denominator, B=B)
    fd <- pbapply::pblapply(seq_len(B+1L),
        function(i) .Ttest_fit1(ii[,i], d1, d2))
    Tval <- sapply(fd, function(z) attr(z, "Tval"))
    out <- do.call(cbind, fd)
    attr(out, "Tval") <- Tval
    out
}

## form difference matrix with mixed bootstrap
edma_fdm <- function(numerator, denominator, B=0) {
    .compare_objects(numerator, denominator)
    fd <- stack(formdiff(numerator, denominator))
    b <- .Ttest_fit(numerator, denominator, B=B)
    out <- list(
        call=match.call(),
        numerator=numerator,
        denominator=denominator,
        B=B,
        dm=fd,
        boot=b)
    class(out) <- c("edma_fdm", class(out))
    out
}
.print_edma_fdm <- function(x, title="EDMA", truncate=20, ...) {
    H <- T_test(x)
    cat(title, "\n",
        #"Numerator: ", .shorten_name(x$numerator$name, truncate), "\n",
        #"Denominator: ", .shorten_name(x$denominator$name, truncate), "\n",
        if (x$B)
            paste(x$B, "mixed bootstrap runs") else "no mixed bootstrap",
        ", T=", round(H$statistic, 4), ", p=", round(H$p.value, 4),
        sep="")
    invisible(x)
}
print.edma_fdm <- function(x, ...) {
    .print_edma_fdm(x, "EDMA form difference matrix", ...)
}


## this pulls out the stacked form difference matrix
get_fdm <- function (object, ...) UseMethod("get_fdm")
get_fdm.edma_fdm <- function (object, sort=FALSE,
level = 0.95, ...) {
    out <- object$dm
    ci <- confint(object, level=level)
    out$lower <- ci[,1L]
    out$upper <- ci[,2L]
    if (sort)
        out <- out[order(out$dist, ...),]
    class(out) <- c("fdm", class(out))
    attr(out, "level") <- level
    out
}

## T-test for 2 edma_fit objects
T_test <- function (object, ...) UseMethod("T_test")
T_test.edma_fdm <- function (object, ...) {
    DNAME <- paste(object$call$numerator, object$call$denominator, sep = " / ")
    METHOD <- "Bootstrap based EDMA T-test"
    Tval <- attr(object$boot, "Tval")
    PVAL <- sum(Tval[1L] <= Tval) / length(Tval)
    PARAMETER <- length(Tval) - 1L
    names(Tval) <- "T-value"
    names(PARAMETER) <- "B"
    out <- list(statistic = Tval[1L], parameter = PARAMETER,
        p.value = PVAL, method = METHOD, data.name = DNAME,
        Tvals=Tval)
    class(out) <- c("edma_test", "htest")
    out
}

## CI based on the 2x input object boot sample
confint.edma_fdm <- function (object, parm, level=0.95, ...) {
    d <- stack(formdiff(object$numerator, object$denominator))
    if (missing(parm))
        parm <- seq_len(nrow(d))
    if (is.null(object$numerator$boot) || is.null(object$denominator$boot)) {
        fd <- NULL
    } else {
        B <- min(length(object$numerator$boot), length(object$denominator$boot))
        fd <- sapply(seq_len(B), function(i) {
            stack(.formdiff(object$numerator$boot[[i]]$M,
                            object$denominator$boot[[i]]$M))$dist
        })
    }
    a <- c((1-level)/2, 1-(1-level)/2)
    out <- t(apply(cbind(d$dist, fd), 1, quantile, a))
    if (is.null(fd))
        out[] <- NA
    rownames(out) <- rownames(d)
    out[parm,,drop=FALSE]
}

.sdmplot_ci <- function(x, xlab="", ylab="",
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
    #abline(h=1, col="grey")
    lines(xv, rep(1, k), col="grey")
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

## plot global/local test based on FDM
.plot_edma_fdm <- function(x, type=c("global", "local"),
                           ylab="Distance Ratio", ...) {
    type <- match.arg(type)
    if (type == "global") {
        z <- T_test(x)
        hist(z$Tvals, xlab="T-values", main=z$data.name, ...)
        abline(v=z$statistic, col=2, lwd=2)
    } else {
        z <- get_fdm(x)
        .sdmplot_ci(z, ylab=ylab, ...)
    }
    invisible(x)
}
plot.edma_fdm <- function(x, type=c("global", "local"), ylab, ...) {
    if (missing(ylab))
        ylab <- "FDM Ratio"
    .plot_edma_fdm(x, type, ylab, ...)
}

## stacked growth matrix
## inputs are edma_fit objects
edma_gm <- function (a1, a2, B=0) {
    .compare_objects(a1, a2, "a1 vs a2:")
    out <- edma_fdm(numerator=a2, denominator=a1, B=B)
    class(out) <- c("edma_gm", class(out))
    out$call <- match.call()
    out
}
print.edma_gm <- function(x, ...) {
    .print_edma_fdm(x, "EDMA growth matrix", ...)
}
plot.edma_gm <- function(x, type=c("global", "local"), ylab, ...) {
    if (missing(ylab))
        ylab <- "GM Ratio"
    .plot_edma_fdm(x, type, ylab, ...)
}

## stacked growth difference matrix
## inputs are edma_fit objects
edma_gdm <- function (a1, a2, b1, b2, B=0) {
    .compare_objects(a1, a2, "a1 vs a2:")
    .compare_objects(b1, b2, "b1 vs b2:")
    .compare_objects(a1, b1, "a1 vs b1:")
    .compare_objects(a2, b2, "a2 vs b2:")
    gma <- edma_fdm(numerator=a2, denominator=a1, B=B)
    gmb <- edma_fdm(numerator=b2, denominator=b1, B=B)
    out <- list(
        call=match.call(),
        a1=a1, a2=a2, b1=b1, b2=b2,
        B=B,
        gdm=gma$dm,
        boot=gma$boot)
    out$gdm$dist <- gmb$dm$dist / gma$dm$dist
    out$boot <- gmb$boot / gma$boot
    attr(out$boot, "Tval") <- apply(out$boot, 2, max) / apply(out$boot, 2, min)
    class(out) <- c("edma_gdm", class(out))
    out
}

plot.edma_gdm <- function(x, ...) {
    .sdmplot_ci(x, ylab="GDM Ratio", ...)
}


