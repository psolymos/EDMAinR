## visualization generics

plot_2d <- function (x,...) UseMethod("plot_2d")
plot_3d <- function (x,...) UseMethod("plot_3d")
plot_ord <- function (x,...) UseMethod("plot_ord")
plot_clust <- function (x,...) UseMethod("plot_clust")
plot_test <- function (x,...) UseMethod("plot_test")
plot_ci <- function (x,...) UseMethod("plot_ci")

## color options

edma_colors <- function(n,
type=c("diverging", "sequential", "qualitative"), alpha=1, rev=FALSE) {
    if (n < 1L)
        stop("Number of colors must be > 0.")
    type <- match.arg(type)
    op <- getOption("edma_options")
    ## avoid surprises
    opq <- op$qualitative
    if (is.null(opq) || any(is.na(opq)))
        opq <- "Set 2"
    opd <- op$diverging
    if (is.null(opd) || any(is.na(opd)))
        opd <- "Blue-Red"
    ## qualitative
    if (type == "qualitative") {
        pals <- hcl.pals("qualitative")
        if (length(opq) == 1L && opq %in% pals) {
            col <- hcl.colors(n, opq)
        } else {
            cols <- opq
            if (n > length(cols))
                cols <- rep(cols, n)
            col <- cols[seq_len(n)]
        }
    } else {
        palsd <- hcl.pals("diverging")
        palss <- hcl.pals("sequential")
        if (type == "diverging") {
            if (length(opd) == 1L && opd %in% palsd) {
                col <- hcl.colors(n, opd)
            } else {
                col <- colorRampPalette(opd)(n)
            }
        }
        if (type == "sequential") {
            if (length(opd) == 1L && opd %in% palss) {
                col <- hcl.colors(n, opd)
            } else {
                if (length(opd) == 1L && opd %in% palsd) {
                    m <- n+(n-1L)
                    col <- hcl.colors(m, opd)[n:m]
                } else {
                    col <- colorRampPalette(opd)(n)
                }
            }
        }
    }
    ## alpha
    if (is.null(alpha) || any(is.na(alpha)))
        alpha <- 1
    alpha <- max(0, min(1, alpha[1L]))
    rgba <- col2rgb(col, alpha=TRUE)
    rgba["alpha",] <- round(alpha * 255)
    col <- rgb(rgba["red",], rgba["green",], rgba["blue",],
        rgba["alpha",], maxColorValue=255)
    ## reverse
    if (rev)
        rev(col) else col
}

plot_edma_colors <- function(n=9, maxq=9) {
    op <- par(mar=c(0,0,0,0)+0.1)
    on.exit(par(op))
    cold <- edma_colors(n, "diverging")
    cols <- edma_colors(n, "sequential")
    nq <- min(n, maxq)
    colq <- edma_colors(nq, "qualitative")
    plot(0, ann=FALSE, axes=FALSE, ylim=c(3,0), xlim=c(0, n), type="n")
    for (i in seq_len(n)) {
        polygon(c(i-1, i-1, i, i), c(0.1, 0.9, 0.9, 0.1),
                border=cold[i], col=cold[i])
        polygon(c(i-1, i-1, i, i), c(1.1, 1.9, 1.9, 1.1),
                border=cols[i], col=cols[i])
    }
    u <- n/nq
    for (i in seq_len(nq)) {
        polygon(c((i-1)*u, (i-1)*u, i*u, i*u),
                c(2.1, 2.9, 2.9, 2.1),
                border=colq[i], col=colq[i])
    }
    text(n/2, 0.5, "diverging")
    text(n/2, 1.5, "sequential")
    text(n/2, 2.5, "qualitative")
    invisible(list(diverging=cold, sequential=cols, qualitative=colq))
}

## diagnostics plot for xyz objects

.data_ellipse <- function (x, y, level=0.95, segments = 51, ...) {
    if (missing(y)) {
      y <- x[,2L]
      x <- x[,1L]
    }
    dfn <- 2
    dfd <- length(x) - 1
    v <- cov.wt(cbind(x, y))
    shape <- v$cov
    center <- v$center
    radius <- sqrt(dfn * qf(level, dfn, dfd))
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    Q <- chol(shape, pivot = TRUE)
    order <- order(attr(Q, "pivot"))
    ell <- t(center + radius * t(unit.circle %*% Q[, order]))
    colnames(ell) <- c("x", "y")
    ell
}

## issue when only one specimen is present:
## rowMeans(sapply(fm[-which], ...) gives error
.plot_edma_data <- function(x, which=NULL,
col_chull=NA, col_spec=2, hull=TRUE, level=0.95, segments=51,
xlim=NULL, ylim=NULL, labels=FALSE, label_cex=1, label_off=c(0,0), ...) {
    c3 <- edma_colors(3, "diverging")
    if (is.na(col_chull))
        col_chull <- c3[1L]
    if (is.na(col_spec))
        col_spec <- c3[3L]
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
    colnames(pca) <- c("X", "Y")
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
    lims <- apply(do.call(rbind, pci), 2, range)
    if (is.null(xlim))
      xlim <- lims[,1L]
    if (is.null(ylim))
      ylim <- lims[,2L]
    op <- par(xpd=NA)
    on.exit(par(op))
    plot(pca, pch=3, axes=FALSE, ann=FALSE, xlim=xlim, ylim=ylim, asp=1, ...)
    if (labels) {
      text(pca[,1L]+label_off[1L], pca[,2L]+label_off[2L],
        labels=landmarks(x), cex=label_cex)
    }
    for (j in seq_len(K)) {
        ll <- t(sapply(pci[ii], function(z) z[j,]))
        if (hull) {
            polygon(ll[chull(ll),],
                col=paste0(substr(col_chull, 1, 7), "44"), border=col_chull)
        } else {
            polygon(.data_ellipse(ll, level=level, segments=segments),
                col=paste0(substr(col_chull, 1, 7), "44"), border=col_chull)
        }
    }
    if (!is.null(which)) {
        segments(pca[,1], pca[,2],
            pci[[which]][,1], pci[[which]][,2], col=col_spec)
        points(pci[[which]][,1], pci[[which]][,2],
            pch=19, cex=0.5, col=col_spec)
    }
    attr(x, "coordinates") <- pca
    invisible(x)
}

plot.edma_data <- function(x, which=NULL,
ask = dev.interactive(), ...) {
    n <- length(x$data)
    if (is.na(ask)) {
        x <- .plot_edma_data(x, NULL,  ...)
    } else {
        if (!is.null(which)) {
            if (which < 1 || which > n)
                stop(sprintf("which must be <= %s", n))
            x <- .plot_edma_data(x, which, ...)
            title(main=paste("Specimen", which))
        } else {
            if (ask) {
                oask <- devAskNewPage(TRUE)
                on.exit(devAskNewPage(oask))
            }
            x <- .plot_edma_data(x, NULL,  ...)
            title(main="All specimens")
            for (i in seq_len(n)) {
                .plot_edma_data(x, i, ...)
                title(main=paste("Specimen", specimens(x)[i]))
            }
        }
    }
    invisible(x)
}
plot_2d.edma_data <- function(x, which=NULL, ...)
    plot.edma_data(x, which=which, ask=NA, ...)

## ordination and cluster plots

.plot_specimens_ord1 <- function (x, cex.label=0.6, col=1, plot=TRUE, ...) {
    d <- as.dist(.get_data(x))
    mds <- cmdscale(sqrt(d), k=2, add=TRUE)
    if (plot) {
        xlim <- range(mds$points[,1L])
        xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
        ylim <- range(mds$points[,2L])
        ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
        op <- par(xpd=NA)
        on.exit(par(op))
        plot(mds$points, type="p",
            ylim=ylim, xlim=xlim,
            xlab="Axis 1", ylab="Axis 2", main="Ordination", col=col, ...)
        text(mds$points, labels=specimens(x), cex=cex.label, pos=1, col=col)
    }
    colnames(mds$points) <- c("Axis 1", "Axis 2")
    invisible(mds)
}

.plot_specimens_clust1 <- function (x, cex.label=0.6,
    plot=TRUE, method="ward.D2", ...) {
    d <- as.dist(.get_data(x))
    h <- hclust(d, method=method)
    if (plot) {
        plot(as.phylo(h), font=1,
            main="Dendrogram", cex=cex.label, ...)
    }
    invisible(h)
}

.plot_specimens_ord2 <- function (x, cex.label=0.6, plot=TRUE,
                                  legend="topleft", ...) {
    xx <- combine_data(
        .get_data(x$numerator),
        .get_data(x$denominator))
    c2 <- edma_colors(2, "qualitative")
    g <- c2[xx$groups]
    d <- as.dist(xx)
    mds <- cmdscale(sqrt(d), k=2, add=TRUE)
    if (plot) {
        xlim <- range(mds$points[,1L])
        xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
        ylim <- range(mds$points[,2L])
        ylim <- ylim + c(-0.1, 0.1) * diff(ylim)
        op <- par(xpd=NA)
        on.exit(par(op))
        plot(mds$points, type="p",
            ylim=ylim, xlim=xlim,
            xlab="Axis 1", ylab="Axis 2", main="Ordination",
            col=g, ...)
        text(mds$points, labels=specimens(xx), cex=cex.label,
            col=g, pos=1)
        legend(legend, legend=c("Numerator", "Denominator"),
              pch=1, col=c2, bty="n", horiz=TRUE)
    }
    invisible(list(mds=mds, data=xx))
}

.plot_specimens_clust2 <- function (x, cex.label=0.6,
    plot=TRUE, method="ward.D2", legend="topleft", ...) {
    xx <- combine_data(
        .get_data(x$numerator),
        .get_data(x$denominator))
    c2 <- edma_colors(2, "qualitative")
    g <- c2[xx$groups]
    d <- as.dist(xx)
    h <- hclust(d, method=method)
    if (plot) {
        plot(as.phylo(h), font=1,
            main="Dendrogram", tip.color=g, cex=cex.label, ...)
        legend(legend, legend=c("Numerator", "Denominator"),
              lty=1, col=c2, bty="n", horiz=TRUE)
    }
    invisible(list(hclust=hclust, data=xx))
}

.plot_specimens_ord4 <- function (x, cex.label=0.6, plot=TRUE, ...) {
    xx <- combine_data4(
        .get_data(x$a1),
        .get_data(x$a2),
        .get_data(x$b1),
        .get_data(x$b2))
    c4 <- edma_colors(4, "qualitative")
    g <- c4[xx$groups]
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

.plot_specimens_clust4 <- function (x, cex.label=0.6,
    plot=TRUE, method="ward.D2", ...) {
    xx <- combine_data4(
        .get_data(x$a1),
        .get_data(x$a2),
        .get_data(x$b1),
        .get_data(x$b2))
    c4 <- edma_colors(4, "qualitative")
    g <- c4[xx$groups]
    d <- as.dist(xx)
    h <- hclust(d, method=method)
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

## global Tobs-test plots

plot_test.edma_dm <- function(x, ...) {
    c5 <- edma_colors(5, "diverging")
    z <- .global_test(x)
    hist(z$Tvals, xlab="Tobs values", main="",
        col=c5[2L], border=c5[1L],
        xlim=range(c(z$statistic, z$Tvals)), ...)
    abline(v=z$statistic, lwd=2, col=c5[5L])
    invisible(x)
}

## plot for influence

.plot_influence <- function(x, xlab="", ylab="",
    bottom=1.5, xcex=0.5, xshow=TRUE, ...) {
    WHAT <- paste0(attr(x, "statistic"), "val")
    DROP <- paste0(attr(x, "statistic"), "drop")
    x <- x[order(x[[DROP]]),]
    k <- nrow(x)
    xv <- seq_len(k)
    r <- range(x[[DROP]], x$lower, x$upper, attr(x, "null"), na.rm=TRUE)
    op <- par(srt=90, xpd = TRUE, mar=par()$mar*c(bottom, 1, 1, 1))
    on.exit(par(op), add=TRUE)
    c5 <- edma_colors(5, "diverging")
    plot(xv, x[[DROP]], ylim=r, type="n",
        xlab=xlab, ylab=ylab, axes=FALSE)
    polygon(c(xv, rev(xv)), c(x$lower, rev(x$upper)), border=NA,
        col=c5[2L])
    lines(xv, rep(attr(x, "null"), k), col=c5[3L])
    Tv <- attr(x, WHAT)
    for (i in seq_len(k))
        if (!is.na(x$upper[i]) && !is.na(x$lower[i]))
        lines(xv[i]+c(-0.5, 0.5), c(Tv, Tv),
            col=if (Tv > x$upper[i] ||
                    Tv < x$lower[i]) c5[5L] else c5[3L])
    lines(xv, x[[DROP]], col=c5[1L])
    axis(2)
    if (xshow) {
        lab <- x$landmark
        text(xv, min(r) - 0.02 * diff(r), lab, adj=c(1, 0.5), cex=xcex)
    }
    invisible(x)
}
plot.edma_influence <- function(x, ...) {
    ylab <- paste0(attr(x, "statistic"), "-value")
    .plot_influence(x, ylab=ylab, ...)
    invisible(x)
}


## CI plot for FDM, GM, GDM

.plot_ci <- function(x, xlab="", ylab="",
    bottom=1.5, xcex=0.5, xshow=TRUE, what="dist", ...) {
    x <- x[order(x[[what]]),]
    k <- nrow(x)
    xv <- seq_len(k)
    r <- range(x[[what]], x$lower, x$upper, 1, na.rm=TRUE)
    op <- par(srt=90, xpd = TRUE, mar=par()$mar*c(bottom, 1, 1, 1))
    on.exit(par(op), add=TRUE)
    c5 <- edma_colors(5, "diverging")
    plot(xv, x[[what]], ylim=r, type="n",
        xlab=xlab, ylab=ylab, axes=FALSE)
    polygon(c(xv, rev(xv)), c(x$lower, rev(x$upper)), border=NA,
        col=c5[2L])
    lines(xv, rep(1, k), col=c5[3L])
    for (i in seq_len(k))
        if (!is.na(x$upper[i]) && !is.na(x$lower[i]))
            lines(xv[i]+c(-0.5, 0.5), c(1, 1),
                col=if (1 > x$upper[i] || 1 < x$lower[i])
                    c5[5L] else c5[3L])
    lines(xv, x[[what]], col=c5[1L])
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
plot_ci.edma_sdm <- function(x, ...)
    .plot_ci(get_sdm(x), ylab="Shape Difference", what="sdm", ...)

## 2D and 3D plots

.plot_d_data <- function(proto, d3=TRUE,
cex=NULL, pch=19, col=NULL, alpha=0.8,
labels=FALSE, label_cex=1, ...) {
    xyz <- Meanform(proto)
    if (!d3) {
        if (ncol(xyz) > 2L) {
            xyz[,1:2] <- cmdscale(dist(xyz), k=2, add=TRUE)$points
            xyz <- xyz[,1:2]
        }
    }
    V <- diag(SigmaKstar(proto))
    if (any(V < 0)) # diag can have negative values
      V <- V + abs(min(V))
    V <- sqrt(V)
    V <- V - min(V)
    V <- V / max(V)
    V <- 0.25 + 1.75 * V
    if (!is.null(cex)) {
        if (length(cex) > 1L) {
            V <- cex
        } else {
            V <- V * cex
        }
    }
    if (is.null(col)) {
        pal <- edma_colors(5, "sequential", alpha=alpha)
        coli <- cut(V, length(pal)-1L, labels=FALSE, include.lowest=TRUE)
        if (length(unique(coli)) < 2L)
          coli[] <- 1L
        col <- pal[coli]
    }
    if (d3) {
        requireNamespace("rgl")
        rgl::plot3d(xyz[,1L], xyz[,2L], xyz[,3L],
            type="s",
            ann=FALSE, axes=FALSE,
            xlab="", ylab="", zlab="",
            radius=V*diff(range(xyz))/50, col=col, ...)
        if (labels)
            rgl::text3d(xyz, texts=rownames(xyz), pos=1, cex=label_cex)
        rgl::aspect3d("iso")
    } else {
        plot(xyz[,1:2], type="p", axes=FALSE, ann=FALSE,
            cex=V, pch=pch, col=col, asp=1, ...)
    }
    invisible(xyz)
}
plot_2d.edma_fit <- function(x, ...) .plot_d_data(x, d3=FALSE, ...)
plot_3d.edma_fit <- function(x, ...) .plot_d_data(x, d3=TRUE, ...)

## midpoints for pairwise distances
## proto is an edma_fit object
.midpoints <- function(proto) {
  xyz <- Meanform(proto)
  M <- Meanform(proto)
  d <- get_fm(proto)[,c("dist", "row", "col")]
  d$row <- as.character(d$row)
  d$col <- as.character(d$col)
  Mi <- M[match(d$row, rownames(M)),]
  Mj <- M[match(d$col, rownames(M)),]
  0.5 * (Mi + Mj)
}

.plot_d_dm <- function(x, d3=TRUE, pal=NULL, pch=19,
cex=1, alpha=0.8, signif_only=FALSE, q=0.1,
labels=FALSE, label_cex=1, proto=NULL, level = 0.95, ...) {
    if (is.null(proto)) {
        if (inherits(x, "edma_fdm")) {
            proto <- if (x$ref_denom)
                x$denominator else x$numerator
        }
        if (inherits(x, "edma_gdm")) {
            proto <- if (x$ref_denom)
                x$a1 else x$a2
        }
        xyz <- Meanform(proto)
    } else {
        xyz <- proto
    }

    if (!d3) {
        if (ncol(xyz) > 2L) {
            xyz[,1:2] <- cmdscale(dist(xyz), k=2, add=TRUE)$points
            xyz <- xyz[,1:2]
        }
    }
    ## landmarks
    i <- get_influence(x)
    td <- attr(i, "Tval") - i$Tdrop
    td <- td/max(td)
    icol <- cut(td, seq(0, 1, 0.25), labels=FALSE) + 1L
    icol[is.na(icol)] <- 1L

    ## pairwise distances
    #f <- get_fdm(x, level=level, sort=TRUE)
    f <- x$dm
    ci <- confint(x, level=level)
    f$lower <- ci[,1L]
    f$upper <- ci[,2L]
    f <- f[order(f$dist),]
    ## pairwise distances: significance
    f$sign <- 0
    f$sign[f$lower < 1 & f$upper < 1] <- -1
    f$sign[f$lower > 1 & f$upper > 1] <- 1
    ## pairwise distances: colors
    f$lcol <- 1 - ifelse(f$dist < 1, f$dist, 1/f$dist)
    f$lcol[f$dist < 1] <- -f$lcol[f$dist < 1]
    f$lcol[f$dist > 1] <- f$lcol[f$dist > 1] / max(f$lcol[f$dist > 1])
    f$lcol[f$dist < 1] <- f$lcol[f$dist < 1] / abs(min(f$lcol[f$dist < 1]))
    f$lcol <- round(5 * f$lcol) + 6
    f$lcol[f$lcol <= 1] <- 1
    f$lcol[f$lcol >= 11] <- 11
    ## quantiles
    f$quant <- (seq_len(nrow(f))-0.5)/nrow(f)

    s5 <- edma_colors(5, "sequential", alpha=alpha)
    c5 <- edma_colors(5, "diverging", alpha=alpha)
    pal <- edma_colors(11, "diverging", alpha=alpha)

    lo <- f[f$dist < 1,,drop=FALSE]
    hi <- f[f$dist >= 1,,drop=FALSE]
    lo <- lo[lo$quant <= q/2,,drop=FALSE]
    hi <- hi[hi$quant >= 1-q/2,,drop=FALSE]
    if (signif_only) {
      lo <- lo[lo$sign < 0,,drop=FALSE]
      hi <- hi[hi$sign > 0,,drop=FALSE]
    }

    if (d3) {
        requireNamespace("rgl")
        rgl::plot3d(xyz[,1L], xyz[,2L], xyz[,3L],
            type="s",
            ann=FALSE, axes=FALSE,
            xlab="", ylab="", zlab="",
            col=s5[icol], radius=0.1*cex)
        for (j in seq_len(nrow(lo))) {
            xyz1 <- rbind(xyz[as.character(lo$row[j]),],
                xyz[as.character(lo$col[j]),])
            rgl::lines3d(xyz1[,1L], xyz1[,2L], xyz1[,3L],
                col=pal[lo$lcol[j]],
                lwd=2)
        }
        for (j in seq_len(nrow(hi))) {
            xyz1 <- rbind(xyz[as.character(hi$row[j]),],
                xyz[as.character(hi$col[j]),])
            rgl::lines3d(xyz1[,1L], xyz1[,2L], xyz1[,3L],
                col=pal[hi$lcol[j]],
                lwd=2)
        }
        if (labels)
            rgl::text3d(xyz, texts=rownames(xyz), pos=1, cex=label_cex)
        rgl::aspect3d("iso")
    } else {
        op <- par(xpd=NA)
        on.exit(par(op))
        plot(xyz[,1:2], type="n", axes=FALSE, ann=FALSE, asp=1)
        segments(
            x0=xyz[as.character(lo$row),1],
            y0=xyz[as.character(lo$row),2],
            x1=xyz[as.character(lo$col),1],
            y1=xyz[as.character(lo$col),2],
            col=pal[lo$lcol],
            lwd=1)
        segments(
            x0=xyz[as.character(hi$row),1],
            y0=xyz[as.character(hi$row),2],
            x1=xyz[as.character(hi$col),1],
            y1=xyz[as.character(hi$col),2],
            col=pal[hi$lcol],
            lwd=1)
        points(xyz[,1:2], pch=pch, col=s5[icol], cex=cex)
        if (labels)
            text(xyz[,1:2], labels=rownames(xyz), pos=1, cex=label_cex)
    }
    invisible(xyz)
}


## this version messes up the colors
## q is the % shown, 0.1=10% --> 5% top & 5% bottom
.plot_d_dm_old <- function(x, d3=TRUE, pal=NULL, pch=19,
cex=1, alpha=0.8, signif_only=TRUE, q=0.1,
midpoints=FALSE, breaks=NULL,
labels=FALSE, label_cex=1, proto=NULL, ...) {
    if (is.null(proto)) {
        if (inherits(x, "edma_fdm")) {
            proto <- if (x$ref_denom)
                x$denominator else x$numerator
        }
        if (inherits(x, "edma_gdm")) {
            proto <- if (x$ref_denom)
                x$a1 else x$a2
        }
        xyz <- if (midpoints)
          .midpoints(proto) else Meanform(proto)
    } else {
        xyz <- proto
    }
    all <- !signif_only && q >= 1

    if (!d3) {
        if (ncol(xyz) > 2L) {
            xyz[,1:2] <- cmdscale(dist(xyz), k=2, add=TRUE)$points
            xyz <- xyz[,1:2]
        }
    }
    ## landmarks
    i <- get_influence(x)
    #iSig <- !(i$lower < attr(i, "Tval") & i$upper > attr(i, "Tval"))
    #names(iSig) <- rownames(xyz)
    td <- attr(i, "Tval") - i$Tdrop
    td <- td/max(td)
    icol <- cut(td, seq(0, 1, 0.25), labels=FALSE) + 1L
    icol[is.na(icol)] <- 1L

    ## pairwise distances
    f <- x$dm
    ci <- confint(x)
    f$lower <- ci[,1L]
    f$upper <- ci[,2L]
    o <- order(abs(log(f$dist)))
    f <- f[o,]
    if (midpoints)
        xyz <- xyz[o,]
    fSig <- f$dist < f$lower | f$dist > f$upper
    Max <- max(1/min(1, f$dist), max(f$dist))
    if (midpoints) {
      Sv <- if (is.null(breaks))
        c(0, 0.4, 0.7, 0.9, 0.95, 1) else breaks
      v1 <- quantile(f$dist[f$dist < 1], rev(1-Sv))
      v2 <- quantile(f$dist[f$dist > 1], Sv)
      v1[1L] <- 0
      v2[length(v2)] <- Inf
      v <- c(v1[-length(v1)], 1, v2[-1L])
    } else {
      v <- (Max-1) * c(0, 0.1, 0.25, 0.5, 1) + 1
      v <- c(0, rev(1/v[-1]), v[-1], Inf)
    }
    f$cut <- cut(f$dist, v, include.lowest=TRUE, labels=FALSE)

    s5 <- edma_colors(5, "sequential", alpha=alpha)
    c5 <- edma_colors(5, "diverging", alpha=alpha)
    if (is.null(pal))
        pal <- edma_colors(max(1L, length(v)-1L), "diverging", alpha=alpha)
    if (!signif_only && !all) {
        f$cut <- 0
        iii <- seq_len(nrow(f))
        iii1 <- round(nrow(f)*q/2)
        iii2 <- nrow(f) - iii1
        f$cut[iii <= iii1 & f$dist < 1] <- 1
        f$cut[iii >= iii2 & f$dist > 1] <- 2
        pal <- edma_colors(2, "diverging", alpha=alpha)
    }
    if (d3) {
        requireNamespace("rgl")
        if (midpoints) {
          rgl::plot3d(xyz[,1L], xyz[,2L], xyz[,3L],
              type="s",
              ann=FALSE, axes=FALSE,
              xlab="", ylab="", zlab="",
              col=pal[f$cut], radius=0.1*cex)
        } else {
          rgl::plot3d(xyz[,1L], xyz[,2L], xyz[,3L],
              type="s",
              ann=FALSE, axes=FALSE,
              xlab="", ylab="", zlab="",
              col=s5[icol], radius=0.1*cex)
          if (all) {
              for (j in 1:nrow(f)) {
                  xyz1 <- rbind(xyz[as.character(f$row[j]),],
                      xyz[as.character(f$col[j]),])
                  rgl::lines3d(xyz1[,1L], xyz1[,2L], xyz1[,3L],
                      col=pal[f$cut[j]],
                      lwd=2)
              }
          } else {
            if (signif_only) {
              for (j in which(fSig & f$cut != 5)) {
                  xyz1 <- rbind(xyz[as.character(f$row[j]),],
                      xyz[as.character(f$col[j]),])
                  rgl::lines3d(xyz1[,1L], xyz1[,2L], xyz1[,3L],
                      col=pal[f$cut[j]],
                      lwd=2)
              }
            } else {
              for (j in which(f$cut > 0)) {
                  xyz1 <- rbind(xyz[as.character(f$row[j]),],
                      xyz[as.character(f$col[j]),])
                  rgl::lines3d(xyz1[,1L], xyz1[,2L], xyz1[,3L],
                      col=pal[f$cut[j]],
                      lwd=2)
              }
            }
          }
          if (labels)
            rgl::text3d(xyz, texts=rownames(xyz), pos=1, cex=label_cex)
        }
        rgl::aspect3d("iso")
    } else {
      op <- par(xpd=NA)
      on.exit(par(op))
      if (midpoints) {
        plot(xyz[,1:2], axes=FALSE, ann=FALSE,
          col=pal[f$cut], pch=pch, cex=cex, asp=1)
      } else {
        plot(xyz[,1:2], type="n", axes=FALSE, ann=FALSE, asp=1)
        if (all) {
            segments(
                x0=xyz[as.character(f$row),1],
                y0=xyz[as.character(f$row),2],
                x1=xyz[as.character(f$col),1],
                y1=xyz[as.character(f$col),2],
                col=pal[f$cut],
                lwd=1)
        } else {
          if (signif_only) {
            for (j in which(fSig)) {
                xy1 <- xyz[as.character(f$row[j]),1:2]
                xy2 <- xyz[as.character(f$col[j]),1:2]
                lines(rbind(xy1, xy2),
                    col=pal[f$cut[j]],
                    lwd=if (f$cut[j] == 5) 0.5 else 2)
            }
          } else {
            for (j in which(f$cut > 0)) {
                xy1 <- xyz[as.character(f$row[j]),1:2]
                xy2 <- xyz[as.character(f$col[j]),1:2]
                lines(rbind(xy1, xy2),
                    col=pal[f$cut[j]],
                    lwd=if (f$cut[j] == 5) 0.5 else 2)
            }
          }
        }
        points(xyz[,1:2], pch=pch, col=s5[icol], cex=cex)
        if (labels)
          text(xyz[,1:2], labels=rownames(xyz), pos=1, cex=label_cex)
      }
    }
    invisible(xyz)
}

plot_2d.edma_dm <- function(x, ...) .plot_d_dm(x, d3=FALSE, ...)
plot_3d.edma_dm <- function(x, ...) .plot_d_dm(x, d3=TRUE, ...)


## default plot methods

#plot.edma_data <- plot_2d.edma_data
plot.edma_fit <- plot_2d.edma_fit
plot.edma_dm <- plot_2d.edma_dm


## Global plots for shape difference

plot_Ztest <- function (x,...) UseMethod("plot_Ztest")

plot_Ztest.edma_sdm <- function(x, statistic=c("Z", "C"), level = 0.95, ...) {
    statistic <- match.arg(statistic)
    c5 <- edma_colors(5, "diverging")
    if (statistic == "Z") {
        xlab <- "Z-values"
        z <- x$boot$Zval
    } else {
        xlab <- "C-values"
        z <- x$boot$Cval
    }
    a <- c((1-level)/2, 1-(1-level)/2)
    ci <- quantile(z, a)
    hist(z, xlab=xlab, main="",
        col=c5[2L], border=c5[1L], ...)
    abline(v=0, lwd=2, col=c5[5L])
    abline(v=ci, lwd=2, col=c5[1L], lty=2)
    invisible(x)
}


