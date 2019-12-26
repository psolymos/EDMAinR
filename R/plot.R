## visualization generics

plot_2d <- function (x,...) UseMethod("plot_2d")
plot_3d <- function (x,...) UseMethod("plot_3d")
plot_ord <- function (x,...) UseMethod("plot_ord")
plot_clust <- function (x,...) UseMethod("plot_clust")
plot_Ttest <- function (x,...) UseMethod("plot_Ttest")
plot_ci <- function (x,...) UseMethod("plot_ci")

## diagnostics plot for xyz objects

.plot_edma_data <- function(x, which=NULL,
col_chull="#44444444", col_spec=2, ...) {
    n <- length(x$data)
    K <- nrow(x$data[[1]])
    fm <- lapply(x$data, dist)
    fma <- fm[[1]]
    ii <- seq_len(n)
    if (is.null(which)) {
        fma[] <- rowMeans(sapply(fm, function(z) as.numeric(z)))
    } else {
        which <- as.integer(which[[1]])
        fma[] <- rowMeans(sapply(fm[-which], function(z) as.numeric(z)))
        ii <- ii[-which]
    }
    pca <- cmdscale(fma, k=2)
    pci <- lapply(seq_len(n), function(i) {
        pci <- cmdscale(fm[[i]], k=2)
        tmp1 <- apply(abs(pca-pci), 2, max)
        tmp2 <- apply(abs(pca+pci), 2, max)
        if (tmp1[1] > tmp2[1]) {
            pci[,1] <- -pci[,1]
        }
        if (tmp1[2] > tmp2[2]) {
            pci[,2] <- -pci[,2]
        }
        pci
    })
    plot(pca, pch=3, axes=FALSE, ann=FALSE, ...)
    for (j in seq_len(K)) {
        ll <- t(sapply(pci[ii], function(z) z[j,]))
        polygon(ll[chull(ll),], col=col_chull, border=NA)
    }
    if (!is.null(which)) {
        segments(pca[,1], pca[,2],
            pci[[which]][,1], pci[[which]][,2], col=col_spec)
        points(pci[[which]][,1], pci[[which]][,2],
            pch=19, cex=0.5, col=col_spec)
    }
    invisible(x)
}

plot_2d.edma_data <- function(x, which=NULL,
ask = dev.interactive(),
col_chull="#44444444", col_spec=2, ...) {
    n <- length(x$data)
    if (!is.null(which)) {
        if (which < 1 || which > n)
            stop(sprintf("which must be <= %s", n))
        .plot_edma_data(x, which, col_chull, col_spec, ...)
        title(main=paste("Specimen", which))
    } else {
        if (ask) {
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))
        }
        .plot_edma_data(x, NULL, col_chull, col_spec, ...)
        title(main="All specimens")
        for (i in seq_len(n)) {
            .plot_edma_data(x, i, col_chull, col_spec, ...)
            title(main=paste("Specimen", i))
        }
    }
    invisible(x)
}

## ordination and cluster plots

.plot_specimens_ord1 <- function (x, cex.label=0.6, plot=TRUE, ...) {
    d <- as.dist(.get_data(x))
    mds <- cmdscale(sqrt(d), k=2, add=TRUE)
    if (plot) {
        xlim <- range(mds$points[,1L])
        xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
        ylim <- range(mds$points[,2L])
        ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
        plot(mds$points, type="p",
            ylim=ylim, xlim=xlim,
            xlab="Axis 1", ylab="Axis 2", main="Ordination", ...)
        text(mds$points, labels=specimens(x), cex=cex.label, pos=1)
    }
    invisible(mds)
}

.plot_specimens_clust1 <- function (x, cex.label=0.6, plot=TRUE, ...) {
    d <- as.dist(.get_data(x))
    h <- hclust(d, method="ward.D2")
    if (plot) {
        plot(as.phylo(h), font=1,
            main="Dendrogram", cex=cex.label, ...)
    }
    invisible(h)
}

.plot_specimens_ord2 <- function (x, cex.label=0.6, plot=TRUE, ...) {
    xx <- .combine_data(
        .get_data(x$numerator),
        .get_data(x$denominator))
    g <- c(2, 4)[xx$groups]
    d <- as.dist(xx)
    mds <- cmdscale(sqrt(d), k=2, add=TRUE)
    if (plot) {
        xlim <- range(mds$points[,1L])
        xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
        ylim <- range(mds$points[,2L])
        ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
        plot(mds$points, type="p",
            ylim=ylim, xlim=xlim,
            xlab="Axis 1", ylab="Axis 2", main="Ordination",
            col=g, ...)
        text(mds$points, labels=specimens(xx), cex=cex.label,
            col=g, pos=1)
    }
    invisible(list(mds=mds, data=xx))
}

.plot_specimens_clust2 <- function (x, cex.label=0.6, plot=TRUE, ...) {
    xx <- .combine_data(
        .get_data(x$numerator),
        .get_data(x$denominator))
    g <- c(2, 4)[xx$groups]
    d <- as.dist(xx)
    h <- hclust(d, method="ward.D2")
    if (plot) {
        plot(as.phylo(h), font=1,
            main="Dendrogram", tip.color=g, cex=cex.label, ...)
    }
    invisible(list(hclust=hclust, data=xx))
}

.plot_specimens_ord4 <- function (x, cex.label=0.6, plot=TRUE, ...) {
    xx <- .combine_data4(
        .get_data(x$a1),
        .get_data(x$a2),
        .get_data(x$b1),
        .get_data(x$b2))
    g <- xx$groups
    d <- as.dist(xx)
    mds <- cmdscale(sqrt(d), k=2, add=TRUE)
    if (plot) {
        xlim <- range(mds$points[,1L])
        xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
        ylim <- range(mds$points[,2L])
        ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
        plot(mds$points, type="p",
            ylim=ylim, xlim=xlim,
            xlab="Axis 1", ylab="Axis 2", main="Ordination",
            col=g, ...)
        text(mds$points, labels=specimens(xx), cex=cex.label,
            col=g, pos=1)
    }
    invisible(list(mds=mds, data=xx))
}

.plot_specimens_clust4 <- function (x, cex.label=0.6, plot=TRUE, ...) {
    xx <- .combine_data4(
        .get_data(x$a1),
        .get_data(x$a2),
        .get_data(x$b1),
        .get_data(x$b2))
    g <- xx$groups
    d <- as.dist(xx)
    h <- hclust(d, method="ward.D2")
    if (plot) {
        plot(as.phylo(h), font=1,
            main="Dendrogram", tip.color=g, cex=cex.label, ...)
    }
    invisible(list(hclust=hclust, data=xx))
}

## gm inherits from fdm
plot_ord.edma_data <- function(x, ...) .plot_specimens_ord1(x, ...)
plot_ord.edma_fit <- function(x, ...) .plot_specimens_ord1(x, ...)
plot_ord.edma_fdm <- function(x, ...) .plot_specimens_ord2(x, ...)
plot_ord.edma_gdm <- function(x, ...) .plot_specimens_ord4(x, ...)
plot_clust.edma_data <- function(x, ...) .plot_specimens_clust1(x, ...)
plot_clust.edma_fit <- function(x, ...) .plot_specimens_clust1(x, ...)
plot_clust.edma_fdm <- function(x, ...) .plot_specimens_clust2(x, ...)
plot_clust.edma_gdm <- function(x, ...) .plot_specimens_clust4(x, ...)

## global T-test plots

plot_Ttest.edma_dm <- function(x, ...) {
    z <- .T_test(x)
    hist(z$Tvals, xlab="T-values", ...)
    abline(v=z$statistic, col=2, lwd=2)
    invisible(x)
}

## plot for influence

.plot_influence <- function(x, xlab="", ylab="",
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
    .plot_influence(x, ylab="T-value", ...)
    invisible(x)
}



## CI plot for FDM, GM, GDM

.plot_ci <- function(x, xlab="", ylab="",
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

plot_ci.edma_fdm <- function(x, ...)
    .plot_ci(get_fdm(x), ylab="FDM Ratio", ...)
plot_ci.edma_gm <- function(x, ...)
    .plot_ci(get_gm(x), ylab="GM Ratio", ...)
plot_ci.edma_gdm <- function(x, ...)
    .plot_ci(get_gdm(x), ylab="GDM Ratio", ...)

## 2D and 3D plots

.plot_d_data <- function(proto, d3=TRUE, ...) {
    xyz <- Meanform(proto)
    if (!d3) {
        xyz[,1:2] <- cmdscale(dist(xyz), k=2, add=TRUE)$points
        xyz[,3L] <- 0
    }
    V <- sqrt(diag(SigmaKstar(proto)))
    V <- 0.2 + 0.8 * V/max(V)
    if (d3) {
        plot3d(xyz[,1L], xyz[,2L], xyz[,3L],
            type="s",
            ann=FALSE, axes=FALSE,
            xlab="", ylab="", zlab="",
            radius=V*diff(range(xyz))/50)
    } else {
        plot(xyz[,1:2], type="p", axes=FALSE, ann=FALSE,
            cex=V, ...)
    }
    invisible(proto)
}
plot_2d.edma_fit <- function(x, ...) .plot_d_data(x, d3=FALSE, ...)
plot_3d.edma_fit <- function(x, ...) .plot_d_data(x, d3=TRUE, ...)

.plot_d_dm <- function(x, d3=TRUE, pal=NULL, ...) {
    if (inherits(x, "edma_fdm")) {
        proto <- if (x$ref_denom)
            x$denominator else x$numerator
    }
    if (inherits(x, "edma_gdm")) {
        proto <- if (x$ref_denom)
            x$a1 else x$a2
    }
    xyz <- Meanform(proto)
    if (!d3) {
        xyz[,1:2] <- cmdscale(dist(xyz), k=2, add=TRUE)$points
        xyz[,3L] <- 0
    }
    i <- get_influence(x)
    f <- x$dm
    ci <- confint(x)
    f$lower <- ci[,1L]
    f$upper <- ci[,2L]
    f <- f[order(abs(log(f$dist))),]
    iSig <- !(i$lower < attr(i, "Tval") & i$upper > attr(i, "Tval"))
    names(iSig) <- rownames(xyz)
    fSig <- !(f$lower < 1 & f$upper > 1)
    Max <- max(1/min(1, f$dist), max(f$dist))
    v <- (Max-1) * c(0, 0.1, 0.25, 0.5, 1) + 1
    v <- c(0, rev(1/v[-1]), v[-1], Inf)
    f$cut <- cut(f$dist, v, include.lowest=TRUE, labels=FALSE)
    if (is.null(pal))
        pal <- colorRampPalette(
            rev(c('#d7191c','#fdae61','#eeeeee','#abd9e9','#2c7bb6'))
            )(length(v)-1)
    if (d3) {
        plot3d(xyz[,1L], xyz[,2L], xyz[,3L],
            type="s",
            ann=FALSE, axes=FALSE,
            xlab="", ylab="", zlab="",
            col=c("grey", "red")[iSig+1], radius=0.1)
        for (j in which(fSig & f$cut != 5)) {
            xyz1 <- rbind(xyz[as.character(f$row[j]),],
                xyz[as.character(f$col[j]),])
            lines3d(xyz1[,1L], xyz1[,2L], xyz1[,3L],
                col=as.character(pal[f$cut[j]]),
                lwd= 2)
        }

    } else {
        plot(xyz[,1:2], type="n", axes=FALSE, ann=FALSE)
        for (j in which(fSig)) {
            xy1 <- xyz[as.character(f$row[j]),1:2]
            xy2 <- xyz[as.character(f$col[j]),1:2]
            lines(rbind(xy1, xy2),
                col=paste0(pal[f$cut[j]], "ff"),
                lwd=if (f$cut[j] == 5) 0.5 else 2)
        }
        points(xyz[iSig,1:2], pch=19, col=2)
        points(xyz[,1:2])
    }
    invisible(x)
}

plot_2d.edma_dm <- function(x, ...) .plot_d_dm(x, d3=FALSE, ...)
plot_3d.edma_dm <- function(x, ...) .plot_d_dm(x, d3=TRUE, ...)


## default plot methods

plot.edma_data <- plot_2d.edma_data
plot.edma_fit <- plot_2d.edma_fit
plot.edma_dm <- plot_2d.edma_dm

