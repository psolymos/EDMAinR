## read in xyz landmark data
## file: file name
## ...: passed to read.table for things like dec
## Structure is the following:
## Header
## XYZ
## 42L 3 9    (which means 42 landmarks, 3 dimensions, 9 specimens)
## amsph bas ethma ethmp intpar laalf lalp lflac lfppm liohd lnasapl loci lpalf lpfl lplpp lpmx lpns lpsq lpto lptyp lsqu lsyn lzyt raalf ralp rflac rfppm riohd rmaxi rnasapl roci rpalf rpfl rplpp rpmx rpns rpsq rpto rptyp rsqu rsyn rzyt
## 3.27773899999999995813 3.08274600000000020827 3.82010100000000019094
## 1.69021700000000008046 3.48286000000000006693 2.04403800000000002157
## ...I truncated the XYZ data itself here...
## a blank line, then the date on of scans for each specimen
##
## CZEM_087 Scanned on
## CZEM_094 Scanned on
## CZEM_097 Scanned on
## etc, this part is particularly useful because then you actually know the specimen numbers for the data within the file. Space is the delimiter in the XYZ files.
read_xyz <- function(file, split_spec_names=TRUE, ...) {
    h <- readLines(file, n=4L)
    NAME <- h[1L]
    DIMS <- character(nchar(h[2L]))
    for (i in seq_along(DIMS))
        DIMS[i] <- substr(h[2L], i, i)
    tmp <- gsub("\\D+","", strsplit(h[3L], " ")[[1L]])
    K <- as.integer(tmp[1L])
    D <- as.integer(tmp[2L])
    DIMS <- DIMS[seq_len(D)]
    n <- as.integer(tmp[3L])
    LABELS <- make.names(strsplit(h[4L], " ")[[1L]])
    if (length(LABELS) < 1L) { # missing
        LABELS <- paste0("L", seq_len(K))
    } else {
        if (length(LABELS) != K) {
            warning("Labels are ignored, length did not match data.")
        }
    }
    X <- read.table(file, header=FALSE, sep=" ", skip=4L, nrows=n*K, ...)
    SPECIMENS <- try(read.table(file, header=FALSE, sep="\t", skip=4L+n*K,
        stringsAsFactors=FALSE)[[1L]], silent=TRUE)
    if (inherits(SPECIMENS, "try-error"))
        SPECIMENS <- paste0("S", seq_len(n))
    if (split_spec_names)
        SPECIMENS <- sapply(strsplit(SPECIMENS, " "), "[[", 1L)
    if (length(SPECIMENS) != n)
        stop("length of specimen names must match n")
    if (any(duplicated(SPECIMENS)))
        stop("specimen names must be unique")
    if (any(duplicated(LABELS)))
        stop("landmark names must be unique")
    DATA <- list()
    for (i in seq_len(n)) {
        DATA[[SPECIMENS[i]]] <- as.matrix(X[((i-1)*K+1):(i*K),,drop=FALSE])
        dimnames(DATA[[i]]) <- list(LABELS, DIMS)
    }
    out <- list(
        name=NAME,
        data=DATA
    )
    class(out) <- c("edma_data")
    out
}

## this turns the data into X expected by fitting functions
## is centering useful here?
stack.edma_data <- function(x, ...) {
    d <- x$data
    out <- do.call(rbind, d)
    rownames(out) <- paste0(rep(specimens(x), each=nrow(x$data[[1L]])),
        "_", rownames(x$data[[1L]]))
    out
}

.shorten_name <- function(x, truncate=20) {
    n <- nchar(x[1L])
    x <- trimws(substr(x[1L], 1, min(truncate, n)))
    if (n > truncate)
        x <- paste0(x, "...")
    x
}

## print function
print.edma_data <- function(x, truncate=20, ...) {
    cat("EDMA data: ", .shorten_name(x$name, truncate), "\n",
        ncol(x$data[[1L]]), " dimensions, ",
        nrow(x$data[[1L]]), " landmarks, ",
        length(x$data), " specimens", sep="")
    invisible(x)
}

## subset specimens in the the data list
subset.edma_data <- function(x, subset, ...) {
    if (missing(subset))
        subset <- seq_along(x$data)
    x$data <- x$data[subset]
    x
}

## subset landmarks, dimensions, specimens
## [landmarks, dims, specimens]
`[.edma_data` <- function (x, i, j, k) {
    if (missing(i))
        i <- seq_len(nrow(x$data[[1L]]))
    if (missing(j))
        j <- seq_len(ncol(x$data[[1L]]))
    if (missing(k))
        k <- seq_along(x$data)
    x <- subset(x, k)
    for (h in names(x$data))
        x$data[[h]] <- x$data[[h]][i,j,drop=FALSE]
    x
}

## retrieve dimensions and dimnames
`dim.edma_data` <- function(x) {
    c(nrow(x$data[[1L]]), ncol(x$data[[1L]]), length(x$data))
}
`dimnames.edma_data` <- function(x) {
    c(dimnames(x$data[[1L]]), list(names(x$data)))
}
landmarks <- function (x, ...) UseMethod("landmarks")
landmarks.edma_data <- function(x, ...) dimnames(x)[[1L]]
specimens <- function (x, ...) UseMethod("specimens")
specimens.edma_data <- function(x, ...) names(x$data)


## coercion methods
as.matrix.edma_data <- function(x, ...) stack(x, ...)
as.data.frame.edma_data <- function(x, ...) as.data.frame(stack(x, ...))
as.array.edma_data <- function (x, ...) {
    out <- array(0, dim(x), dimnames(x))
    for (i in seq_along(x$data))
        out[,,i] <- x$data[[i]]
    out
}

.edma_simulate_data <- function(n, M, SigmaK, H=NULL) {
    K <- nrow(M)
    D <- ncol(M)
    if (D > 3 || D < 2)
        stop("Mean form must have 2 or 3 dimensions.")
    Z <- matrix(nrow = n * K, ncol = D)
    for (i in 1:n) {
        Z[((i - 1) * K + 1):(i * K), ] <- matrix(rnorm(K * D), nrow = K,
            ncol = D)
    }
    Cmat <- chol(SigmaK)
    X <- matrix(nrow = n * K, ncol = D)
    for (i in 1:n) {
        X[((i - 1) * K + 1):(i * K), ] <- crossprod(Cmat,
            Z[((i - 1) * K + 1):(i * K), ]) + M
    }
    if (is.null(H)) {
        ones <- array(rep(1, K), c(1, K))
        H <- diag(1, K) - (1/K) * crossprod(ones, ones)
    }
    SigmaKstar = H %*% SigmaK %*% H
    list(
        M=M, SigmaK=SigmaK, SigmaKstar=SigmaKstar, H=H,
        X=X, D=D, K=K, n=n)
}
edma_simulate_data <- function(n, M, SigmaK, H=NULL) {
    z <- .edma_simulate_data(n, M, SigmaK, H)
    DATA <- list()
    for (i in seq_len(z$n)) {
        DATA[[paste0("S", i)]] <- as.matrix(z$X[((i-1)*z$K+1):(i*z$K),,
                                                drop=FALSE])
        dimnames(DATA[[i]]) <- list(
            paste0("L", seq_len(z$K)),
            c("X", "Y", "Z")[seq_len(z$D)])
    }
    out <- list(
        name="Simulated landmark data",
        data=DATA
    )
    class(out) <- c("edma_data")
    attr(out, "simulation_settings") <- list(M=z$M, SigmaK=z$SigmaK,
        SigmaKstar=z$SigmaKstar, H=z$H)
    out
}

## diagnostics

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

plot.edma_data <- function(x, which=NULL,
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

.combine_data <- function(a, b, ga="G1", gb="G2") {
    ls <- intersect(landmarks(a), landmarks(b))
    if (!all(ls %in% union(landmarks(a), landmarks(b))))
        stop("landmarks must be the same in the two objects")
    ds <- intersect(dimnames(a)[[2L]], dimnames(b)[[2L]])
    if (!all(ds %in% union(dimnames(a)[[2L]], dimnames(b)[[2L]])))
        stop("dimensions must be the same in the two objects")
    names(a$data) <- paste0(ga, "_", names(a$data))
    names(b$data) <- paste0(gb, "_", names(b$data))
    a <- a[ls,ds,]
    b <- b[ls,ds,]
    ab <- a
    na <- length(a$data)
    nb <- length(b$data)
    ab$name <- "data with 2 groups"
    ab$data <- c(a$data, b$data)
    ab$groups <- rep(1:2, c(na, nb))
    ab
}

as.dist.edma_data <- function(m, diag = FALSE, upper = FALSE) {
    n <- dim(m)[3L]
    nam <- dimnames(m)[[3L]]
    mat <- matrix(0, n, n)
    dimnames(mat) <- list(nam, nam)
    for (i in seq_len(n)) {
        for (j in seq_len(i)) {
            di <- dist(m$data[[i]])
            dj <- dist(m$data[[j]])
            fdm <- di / dj
            mat[i,j] <- mat[j,i] <- log(max(fdm) / min(fdm))
        }
    }
    out <- as.dist(mat, diag=diag, upper=upper)
    class(out) <- c("edma_logT", class(out))
    out
}

