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
read_xyz <- function(file, ...) {
    h <- readLines(file, n=4L)
    NAME <- h[1L]
    DIMS <- character(nchar(h[2L]))
    for (i in seq_along(DIMS))
        DIMS[i] <- substr(h[2L], i, i)
    tmp <- gsub("\\D+","", strsplit(h[3L], " ")[[1L]])
    K <- as.integer(tmp[1L])
    D <- as.integer(tmp[2L])
    if (K <= D)
        warning("K <= D: EDMA will not work.")
    if (K < 2L)
        warning("K < 2: EDMA will not work.")
    if (D > 3L || D < 2L)
        warning("D should be 2 or 3: EDMA will not work.")
    DIMS <- DIMS[seq_len(D)]
    n <- as.integer(tmp[3L])
    ## don't standardize names: no known risk at the moment
    #LABELS <- make.names(strsplit(h[4L], " ")[[1L]])
    LABELS <- strsplit(h[4L], " ")[[1L]]
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
    SPECIMENS <- sapply(strsplit(SPECIMENS, " "), "[[", 1L)
    if (length(SPECIMENS) != n)
        stop("Length of specimen names must match n.")
    if (any(duplicated(SPECIMENS)))
        stop("Specimen names must be unique.")
    if (any(duplicated(LABELS)))
        stop("Landmark names must be unique.")
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

write_xyz <- function(x, file) {
    if (!inherits(x, "edma_data"))
        stop("x must be an edma_data object.")
    text <- c(
        x$name,
        "XYZ",
        paste0(dim(x)[1L], "L ", dim(x)[2L], " ", dim(x)[3L]),
        paste(landmarks(x), collapse=" "),
        apply(as.matrix(x), 1, paste, collapse=" "),
        "",
        specimens(x)
    )
    if (missing(file))
        file <- stdout()
    writeLines(text, con = file, sep = "\n")
}


## this turns the data into X expected by fitting functions
stack.edma_data <- function(x, ...) {
    d <- x$data
    out <- do.call(rbind, d)
    rownames(out) <- paste0(rep(specimens(x), each=nrow(x$data[[1L]])),
        "_", rownames(x$data[[1L]]))
    out
}

## shortens the xyz object title, it can be extra long
.shorten_name <- function(x, truncate=40) {
    n <- nchar(x[1L])
    x <- trimws(substr(x[1L], 1, min(truncate, n)))
    if (n > truncate)
        x <- paste0(x, "...")
    x
}

## print function
print.edma_data <- function(x, truncate=40, ...) {
    cat("EDMA data: ", .shorten_name(x$name, truncate), "\n",
        nrow(x$data[[1L]]), " landmarks, ",
        ncol(x$data[[1L]]), " dimensions, ",
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
landmarks.edma_data <- function(x, ...) rownames(x$data[[1L]])
dimensions <- function (x, ...) UseMethod("dimensions")
dimensions.edma_data <- function(x, ...) colnames(x$data[[1L]])
specimens <- function (x, ...) UseMethod("specimens")
specimens.edma_data <- function(x, ...) names(x$data)

`dimnames<-.edma_data` <- function(x, value) {
    names(x$data) <- value[[3L]]
    for (i in seq_along(x$data))
        dimnames(x$data[[i]]) <- value[1:2]
    x
}
`landmarks<-` <- function (x, value) UseMethod("landmarks<-")
`dimensions<-` <- function (x, value) UseMethod("dimensions<-")
`specimens<-` <- function (x, value) UseMethod("specimens<-")
`landmarks<-.edma_data` <- function(x, value) {
    dn <- dimnames(x)
    dn[[1L]] <- value
    dimnames(x) <- dn
    x
}
`dimensions<-.edma_data` <- function(x, value) {
    dn <- dimnames(x)
    dn[[2L]] <- value
    dimnames(x) <- dn
    x
}
`specimens<-.edma_data` <- function(x, value) {
    dn <- dimnames(x)
    dn[[3L]] <- value
    dimnames(x) <- dn
    x
}


## coercion methods
as.matrix.edma_data <- function(x, ...) stack(x, ...)
as.data.frame.edma_data <- function(x, ...) as.data.frame(stack(x, ...))
as.array.edma_data <- function (x, ...) {
    out <- array(0, dim(x), dimnames(x))
    for (i in seq_along(x$data))
        out[,,i] <- x$data[[i]]
    out
}

## simulate data xyz sets
.edma_simulate_data_old <- function(n, M, SigmaK, H=NULL) {
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
    SigmaKstar <- H %*% SigmaK %*% H
    list(
        M=M, SigmaK=SigmaK, SigmaKstar=SigmaKstar, H=H,
        X=X, D=D, K=K, n=n)
}
.edma_simulate_data <- function(n, M, SigmaK) {
    K <- nrow(M)
    if (K < 2L)
        stop("Mean form must have at leats 2 landmarks (rows).")
    D <- ncol(M)
    if (D > 3L || D < 2L)
        stop("Mean form must have 2 or 3 dimensions (columns).")
    if (K <= D)
        warning("K <= D: EDMA will not work.")
    Z <- array(rnorm(K * D * n), c(K, D, n))
    Cmat <- chol(SigmaK)
    A <- array(0, c(K, D, n))
    for (i in seq_len(n)) {
        A[,,i] <- crossprod(Cmat, Z[,,i]) + M
    }
    ones <- array(rep(1, K), c(1, K))
    H <- diag(1, K) - (1/K) * crossprod(ones, ones)
    SigmaKstar <- H %*% SigmaK %*% H
    list(
        M=unname(M),
        SigmaK=unname(SigmaK),
        SigmaKstar=unname(SigmaKstar),
        A=unname(A),
        D=D, K=K, n=n)
}
edma_simulate_data <- function(n, M, SigmaK) {
    z <- .edma_simulate_data(n, M, SigmaK)
    DATA <- list()
    LM <- paste0("L", seq_len(z$K))
    for (i in seq_len(z$n)) {
        DATA[[paste0("S", i)]] <- matrix(z$A[,,i], z$K, z$D)
        dimnames(DATA[[i]]) <- list(LM, c("X", "Y", "Z")[seq_len(z$D)])
    }
    out <- list(
        name="Simulated landmark data",
        data=DATA
    )
    class(out) <- c("edma_data")
    dimnames(z$M) <- dimnames(DATA[[1L]])
    dimnames(z$SigmaK) <- dimnames(z$SigmaKstar) <- list(LM, LM)
    attr(out, "simulation_settings") <- list(
        M=z$M,
        SigmaK=z$SigmaK,
        SigmaKstar=z$SigmaKstar)
    out
}

## combines 2 xyz data sets into 1 object
## for plotting etc.
combine_data <- function(a, b, ga="G1", gb="G2") {
    ls <- intersect(landmarks(a), landmarks(b))
    if (!all(ls %in% union(landmarks(a), landmarks(b))))
        stop("Landmarks must be the same in the 2 objects.")
    ds <- intersect(dimnames(a)[[2L]], dimnames(b)[[2L]])
    if (!all(ds %in% union(dimnames(a)[[2L]], dimnames(b)[[2L]])))
        stop("Dimensions must be the same in the 2 objects.")
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
## combines 4 xyz data sets into 1 object
## for plotting etc.
combine_data4 <- function(a1, a2, b1, b2,
ga1="A1", ga2="A2", gb1="B1", gb2="B2") {
    ls <- intersect(landmarks(a1), landmarks(a2))
    if (!all(ls %in% union(landmarks(a1), landmarks(a2))))
        stop("Landmarks must be the same in the 4 objects.")
    ls <- intersect(landmarks(b1), landmarks(b2))
    if (!all(ls %in% union(landmarks(b1), landmarks(b2))))
        stop("Landmarks must be the same in the 4 objects.")
    ls <- intersect(landmarks(a1), landmarks(b1))
    if (!all(ls %in% union(landmarks(a1), landmarks(b1))))
        stop("Landmarks must be the same in the 4 objects.")

    ds <- intersect(dimnames(a1)[[2L]], dimnames(a2)[[2L]])
    if (!all(ds %in% union(dimnames(a1)[[2L]], dimnames(a2)[[2L]])))
        stop("Dimensions must be the same in the 4 objects.")
    ds <- intersect(dimnames(b1)[[2L]], dimnames(b2)[[2L]])
    if (!all(ds %in% union(dimnames(b1)[[2L]], dimnames(b2)[[2L]])))
        stop("Dimensions must be the same in the 4 objects.")
    ds <- intersect(dimnames(a1)[[2L]], dimnames(b1)[[2L]])
    if (!all(ds %in% union(dimnames(a1)[[2L]], dimnames(b1)[[2L]])))
        stop("Dimensions must be the same in the 4 objects.")

    names(a1$data) <- paste0(ga1, "_", names(a1$data))
    names(a2$data) <- paste0(ga2, "_", names(a2$data))
    names(b1$data) <- paste0(gb1, "_", names(b1$data))
    names(b2$data) <- paste0(gb2, "_", names(b2$data))
    a1 <- a1[ls,ds,]
    a2 <- a2[ls,ds,]
    b1 <- b1[ls,ds,]
    b2 <- b2[ls,ds,]
    ab <- a1
    na1 <- length(a1$data)
    na2 <- length(a2$data)
    nb1 <- length(b1$data)
    nb2 <- length(b2$data)
    ab$name <- "data with 4 groups"
    ab$data <- c(a1$data, a2$data, b1$data, b2$data)
    ab$groups <- rep(1:4, c(na1, na2, nb1, nb2))
    ab
}

## dissimilarity matrix: pairwise distance is defined
## based on T-statistic (averaged on the log scale)
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

as.edma_data <- function(x, ...) UseMethod("as.edma_data")

as.edma_data.edma_data <- function(x, ...) x

as.edma_data.array <- function(x, ...) {
    d <- dim(x)
    K <- d[1L]
    D <- d[2L]
    n <- d[3L]
    dn <- list(dimnames(x)[[1L]], dimnames(x)[[2L]], dimnames(x)[[3L]])
    if (is.null(dn[[1L]]))
        dn[[1L]] <- paste0("L", seq_len(K))
    if (is.null(dn[[2L]]))
        dn[[2L]] <- c("X", "Y", "Z")[seq_len(D)]
    if (is.null(dn[[3L]]))
        dn[[3L]] <- paste0("S", seq_len(n))
    dimnames(x) <- dn
    DATA <- list()
    LM <- paste0("L", seq_len(K))
    for (i in seq_len(n)) {
        DATA[[i]] <- x[,,i]
    }
    names(DATA) <- dn[[3L]]
    out <- list(
        name="Landmark data",
        data=DATA
    )
    class(out) <- c("edma_data")
    out
}
