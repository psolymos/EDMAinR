plot_specimens <- function (x, ...)
    UseMethod("plot_specimens")
plot_specimens.default <- function (x,
method=c("cluster", "ordination"), cex.label=0.6, ...) {
    method <- match.arg(method)
    d <- as.dist(.get_data(x))
    if (method == "cluster") {
        h <- hclust(d, method="ward.D2")
        plot(as.phylo(h), font=1,
            main="Dendrogram", cex=cex.label, ...)
    }
    if (method == "ordination") {
        mds <- cmdscale(sqrt(d), k=2, add=TRUE)
        xlim <- range(mds$points[,1L])
        xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
        ylim <- range(mds$points[,2L])
        ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
        plot(mds$points, type="p",
            ylim=ylim, xlim=xlim,
            xlab="Axis 1", ylab="Axis 2", main="Ordination", ...)
        text(mds$points, labels=specimens(x), cex=cex.label, pos=1)
    }
    invisible(x)
}

plot_specimens.edma_fdm <- function (x,
method=c("cluster", "ordination"), cex.label=0.6, ...) {
    method <- match.arg(method)
    xx <- .combine_data(
        .get_data(x$numerator),
        .get_data(x$denominator))
    g <- c(2, 4)[xx$groups]
    d <- as.dist(xx)
    if (method == "cluster") {
        h <- hclust(d, method="ward.D2")
        plot(as.phylo(h), font=1,
            main="Dendrogram", tip.color=g, cex=cex.label, ...)
    }
    if (method == "ordination") {
        mds <- cmdscale(sqrt(d), k=2, add=TRUE)
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
    invisible(xx)
}
