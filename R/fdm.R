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

## x: object w bootstrap
## ref: reference object (w or w/o bootstrap)
.Ttest_fit <- function(x, ref, boot=TRUE) {
    if (!boot)
        x$boot <- NULL
    FMref <- as.numeric(as.dist(ref))
    FM <- if (is.null(x$boot)) {
        data.matrix(as.numeric(as.dist(x)))
    } else {
        cbind(as.numeric(as.dist(x)),
            sapply(x$boot, function(z) as.numeric(dist(z$M))))
    }
    FDM <- FM/FMref
    Tval <- apply(FDM, 2, function(z) max(z)/min(z))
    attr(FDM, "Tval") <- Tval
    FDM
}

## form difference matrix with mixed bootstrap
edma_fdm <- function(numerator, denominator, ref_denom=TRUE) {
    .compare_objects(numerator, denominator)
    fd <- stack(formdiff(numerator, denominator))
    b <- if (ref_denom) {
        .Ttest_fit(numerator, denominator)
    } else {
        .Ttest_fit(denominator, numerator)
    }
    out <- list(
        call=match.call(),
        numerator=numerator,
        denominator=denominator,
        ref_denom=ref_denom,
        dm=fd,
        B=ncol(b)-1,
        boot=b)
    class(out) <- c("edma_fdm", class(out))
    out
}
.print_edma_fdm <- function(x, title="EDMA", truncate=20, ...) {
    H <- T_test(x)
    cat(title, "\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        if (x$B)
            paste(x$B, "bootstrap runs") else "no bootstrap",
        " (ref: ", if (x$ref_denom) "denominator" else "numerator", ")",
        "\nT=", round(H$statistic, 4),
        if (x$B)
            paste0(", p=", round(H$p.value, 4)) else "",
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
    a <- c((1-level)/2, 1-(1-level)/2)
    out <- t(apply(object$boot, 1, quantile, a))
    if (object$B < 1)
        out[] <- NA
    rownames(out) <- rownames(d)
    out[parm,,drop=FALSE]
}

landmarks.edma_fdm <- function(x, ...)
    landmarks(x$numerator)

.sdmplot_ci <- function(x, xlab="", ylab="",
    bottom=1.5, xcex=0.5, xshow=TRUE, ...) {
    x <- x[order(x$dist),]
    k <- nrow(x)
    xv <- seq_len(k)
    r <- range(x$dist, x$lower, x$upper, 1)
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

plot.edma_gdm <- function(x, ...) {
    .sdmplot_ci(x, ylab="GDM Ratio", ...)
}

## influential landmarks
## no bootstrap for quick=FALSE
.influence <- function(i, object, quick=TRUE, level=0.95) {
    ls <- landmarks(object)
    names(ls) <- ls
    i <- ls[i]
    lsd <- ls[!(ls %in% i)]
    ci <- c(NA, NA)
    if (quick) {
        a <- c((1-level)/2, 1-(1-level)/2)
        keep <- !(object$dm$row %in% i) & !(object$dm$col %in% i)
        FDMkeep <- object$boot[keep,,drop=FALSE]
        Tvals <- apply(FDMkeep, 2, function(z) max(z)/min(z))
        Tval <- Tvals[1L]
        if (object$B > 0)
            ci <- unname(quantile(Tvals, a))
    } else {
        if (object$ref_denom) {
            x <- edma_fit(.get_data(object$numerator)[lsd,,])
            ref <- edma_fit(.get_data(object$denominator)[lsd,,])
        } else {
            x <- edma_fit(.get_data(object$denominator)[lsd,,])
            ref <- edma_fit(.get_data(object$numerator)[lsd,,])
        }
        Tval <- attr(.Ttest_fit(x, ref, boot=FALSE), "Tval")
    }
    c(Tval, ci)
}

get_influence <- function (object, ...) UseMethod("get_influence")
get_influence.edma_fdm <- function (object,
quick=TRUE, level=0.95, ...) {
    ls <- landmarks(object)
    Tval <- T_test(object)$statistic
    Tvals <- if (quick) {
        t(sapply(ls, .influence,
            object=object, quick=quick, level=level))
    } else {
        t(pbapply::pbsapply(ls, .influence,
            object=object, quick=quick, level=level))
    }
    colnames(Tvals) <- c("Tdrop", "lower", "upper")
    out <- data.frame(landmark=ls, Tvals)
    rownames(out) <- NULL
    attr(out, "Tval") <- Tval
    attr(out, "level") <- level
    attr(out, "quick") <- quick
    class(out) <- c("edma_influence", class(out))
    out
}

.sdmplot_ci2 <- function(x, xlab="", ylab="",
    bottom=1.5, xcex=0.5, xshow=TRUE, ...) {
    x <- x[order(x$Tdrop),]
    k <- nrow(x)
    xv <- seq_len(k)
    r <- range(x$Tdrop, x$lower, x$upper, 1, na.rm=TRUE)
    op <- par(srt=90, xpd = TRUE, mar=par()$mar*c(bottom, 1, 1, 1))
    on.exit(par(op), add=TRUE)
    plot(xv, x$Tdrop, ylim=r, type="n",
        xlab=xlab, ylab=ylab, axes=FALSE)
    polygon(c(xv, rev(xv)), c(x$lower, rev(x$upper)), col="lightblue", border=NA)
    lines(xv, rep(1, k), col="grey")
    Tv <- attr(x, "Tval")
    for (i in seq_len(k))
        if (!is.na(x$upper[i]) && !is.na(x$lower[i]))
        lines(xv[i]+c(-0.5, 0.5), c(Tv, Tv),
            col=if (Tv > x$upper[i] ||
                    Tv < x$lower[i]) "red" else "grey")
    lines(xv, x$Tdrop, col="blue")
    axis(2)
    if (xshow) {
        lab <- x$landmark
        text(xv, min(r) - 0.02 * diff(r), lab, adj=c(1, 0.5), cex=xcex)
    }
    invisible(x)
}
plot.edma_influence <- function(x, ...) {
    .sdmplot_ci2(x, ylab="T-value", ...)
    invisible(x)
}


## stacked growth matrix
## inputs are edma_fit objects
edma_gm <- function (a1, a2, ...) {
    .compare_objects(a1, a2, "a1 vs a2:")
    out <- edma_fdm(numerator=a2, denominator=a1, ...)
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
edma_gdm <- function (a1, a2, b1, b2, ...) {
    .compare_objects(a1, a2, "a1 vs a2:")
    .compare_objects(b1, b2, "b1 vs b2:")
    .compare_objects(a1, b1, "a1 vs b1:")
    .compare_objects(a2, b2, "a2 vs b2:")
    # ref is a1 an b1
    gma <- edma_fdm(numerator=a2, denominator=a1, ...)
    gmb <- edma_fdm(numerator=b2, denominator=b1, ...)
    B <- min(gma$B, gmb$B)
    gdm <- gma$dm
    gdm$dist <- gmb$dm$dist / gma$dm$dist
    i <- seq_len(B+1L)
    b <- gmb$boot[,i,drop=FALSE] / gma$boot[,i,drop=FALSE]

    out <- list(
        call=match.call(),
        a1=a1, a2=a2, b1=b1, b2=b2,
        B=B,
        ref_denom=gma$ref_denom,
        dm=gdm,
        boot=b)
    attr(out$boot, "Tval") <- apply(out$boot, 2, max) / apply(out$boot, 2, min)
#    class(out) <- c("edma_gdm", class(gma))
    class(out) <- c("edma_gdm", class(out))
    out
}

print.edma_gdm <- function(x, ...) {
    .print_edma_fdm(x, "EDMA growth difference matrix", ...)
}


