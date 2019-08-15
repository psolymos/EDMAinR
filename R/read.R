## read in xyz landmark data
## file: file name
## ...: passed to read.table for things like dec
## consider: what to do for 2D samples? Have only XY columns, or add Z=0?
read_xyz <- function(file, ...) {
    h <- readLines(file, n=4L)
    NAME <- h[1L]
    DIMS <- character(nchar(h[2L]))
    for (i in seq_along(DIMS))
        DIMS[i] <- substr(h[2L], i, i)
    tmp <- gsub("\\D+","", strsplit(h[3L], " ")[[1L]])
    K <- as.integer(tmp[1L])
    D <- as.integer(tmp[2L])
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
    NOTES <- try(read.table(file, header=FALSE, sep="\t", skip=4L+n*K,
        stringsAsFactors=FALSE)[[1L]], silent=TRUE)
    if (inherits(NOTES, "try-error"))
        NOTES <- NULL
    DATA <- list()
    for (i in seq_len(n)) {
        DATA[[i]] <- as.matrix(X[((i-1)*K+1):(i*K),,drop=FALSE])
        dimnames(DATA[[i]]) <- list(LABELS, DIMS)
    }
    out <- list(
        name=NAME,
        data=DATA,
        notes=NOTES
    )
    class(out) <- c("edma_data")
    out
}

## this turns the data into X expected by fitting functions
## is centering useful here?
stack.edma_data <- function(x, ...) {
    d <- x$data
    out <- do.call(rbind, d)
    rownames(out) <- paste0("rep", rep(seq_along(x$data),
        each=nrow(x$data[[1L]])), "_",
        rownames(x$data[[1L]]))
    out
}

## print function
print.edma_data <- function(x, ...) {
    cat("EDMA data: ", x$name, "\n",
        ncol(x$data[[1L]]), " dimensions, ",
        nrow(x$data[[1L]]), " landmarks, ",
        length(x$data), " replicates", sep="")
    invisible(x)
}

## subset replicates in the the data list
subset.edma_data <- function(x, subset, ...) {
    if (missing(subset))
        subset <- seq_along(x$data)
    x$data <- x$data[subset]
    x$notes <- x$notes[subset]
    x
}

## subset landmarks, dimensions, replicates
`[.edma_data` <- function (x, i, j, k) {
    if (missing(i))
        i <- seq_len(nrow(x$data[[1L]]))
    if (missing(j))
        j <- seq_len(ncol(x$data[[1L]]))
    if (missing(k))
        k <- seq_along(x$data)
    x <- subset(x, k)
    for (h in seq_along(x$data))
        x$data[[h]] <- x$data[[h]][i,j,drop=FALSE]
    x
}

## retrieve dimensions and dimnames
`dim.edma_data` <- function(x) {
    c(nrow(x$data[[1L]]), ncol(x$data[[1L]]), length(x$data))
}
`dimnames.edma_data` <- function(x) {
    c(dimnames(x$data[[1L]]), NULL)
}
landmark_names <- function (x, ...) UseMethod("landmark_names")
landmark_names.edma_data <- function(x, ...) dimnames(x)[[1L]]


## coercion methods
as.matrix.edma_data <- function(x, ...) stack(x, ...)
as.data.frame.edma_data <- function(x, ...) stack(x, ...)
as.array.edma_data <- function (x, ...) {
    out <- array(0, dim(x), dimnames(x))
    for (i in seq_along(x$data))
        out[,,i] <- x$data[[i]]
    out
}

edma_simulate_data <- function(n, M, SigmaK, H=NULL) {
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

    DATA <- list()
    for (i in seq_len(n)) {
        DATA[[i]] <- as.matrix(X[((i-1)*K+1):(i*K),,drop=FALSE])
        dimnames(DATA[[i]]) <- list(
            paste0("L", seq_len(K)),
            c("X", "Y", "Z")[seq_len(D)])
    }
    out <- list(
        name="Simulated landmark data",
        data=DATA,
        notes=NULL
    )
    class(out) <- c("edma_data")
    attr(out, "simulation_settings") <- list(M=M, SigmaK=SigmaK,
        SigmaKstar=SigmaKstar, H=H)
    out
}

