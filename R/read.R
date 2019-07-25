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
    tmp <- strsplit(h[3L], " ")[[1L]]
    K <- as.integer(substr(tmp[1L], 1L, nchar(tmp[1L])-1L))
    D <- as.integer(tmp[2L])
    n <- as.integer(tmp[3L])
    LABELS <- make.names(strsplit(h[4L], " ")[[1L]])
    X <- read.table(file, header=FALSE, sep=" ", skip=4L, nrows=n*K, ...)
    NOTES <- read.table(file, header=FALSE, sep="\t", skip=4L+n*K,
        stringsAsFactors=FALSE)[[1L]]
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

