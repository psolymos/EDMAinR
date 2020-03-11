## take design matrix and returns a factor for parameter matching
## NA --> 0
## other unique values will be factor levels
.mat2fac <- function(m) {
    K <- nrow(m)
    m[upper.tri(m)] <- m[lower.tri(m)]
    m <- factor(as.character(m))
    attr(m, "K") <- K
    m
}

## takes a named param vector and reconstructs matrix
.vec2mat <- function(parm, fac) {
    x <- matrix(0, attr(fac, "K"), attr(fac, "K"))
    fac <- droplevels(fac)
    if (!all(sort(names(parm)) == sort(levels(fac))))
        stop("names of parm and levels of fact must match")
    for (i in names(parm))
        x[!is.na(fac) & fac == i] <- parm[i]
    x
}

.SigmaK_fit <- function(SigmaKstar, H, pattern, init,
method = "Nelder-Mead", control = list()) {
    K <- nrow(SigmaKstar)
    pattern1 <- pattern
    pattern1[] <- NA
    diag(pattern1) <- diag(pattern)
    pattern0 <- pattern
    diag(pattern0) <- NA
    fac <- .mat2fac(pattern) # all parms
    lev1 <- levels(.mat2fac(pattern1)) # parms in diag
    lev0 <- levels(.mat2fac(pattern0)) # parms off-diag
    if (length(intersect(lev1, lev0)) > 0)
        stop("diagonal and off-diagonal parameters must not overlap")
    p <- nlevels(droplevels(fac))
    ## make sure diags are >0
    ## generate random starting values using runif()
    init0 <- structure(numeric(p), names=levels(fac))
    if (missing(init)) {
        init0[lev1] <- runif(length(lev1))
    } else {
        init <- init[names(init) %in% names(init0)]
        init0[names(init)] <- init
    }
    if (any(init0[lev1] <= 0))
        stop("inits for diagonal elements must be > 0")
    num_max <- .Machine$double.xmax^(1/3)
    fun <- function(parms){
        SigmaK <- .vec2mat(parms, fac)
        if (any(diag(SigmaK) <= 0))
            return(num_max)
        #10^4 * max((SigmaKstar - (H %*% SigmaK %*% H))^2)
        sum((SigmaKstar - (H %*% SigmaK %*% H))^2)
    }
    if (!is.null(control$fnscale) && control$fnscale < 0)
        stop("control$fnscale can not be negative")
    o <- suppressWarnings({
        optim(init0, fun, method=method, control=control, hessian=FALSE)
    })
    o$init <- init0
    o$SigmaK <- .vec2mat(o$par, fac)
    o$method <- method
    o$control <- control
    o
}

.check_pattern1 <- function(pattern) {
    utri <- pattern[upper.tri(pattern)]
    ltri <- t(pattern)[upper.tri(pattern)]
    if (!all(is.na(utri) == is.na(ltri)))
        stop("NA's (0's) are not symmertic in pattern matrix")
    if (!all(utri[!is.na(utri)] == ltri[!is.na(ltri)]))
        stop("non-NA's (non-0's) are not symmertic in pattern matrix")
    pattern[upper.tri(pattern)] <- NA
    pattern[upper.tri(pattern)] <- t(pattern)[upper.tri(pattern)]
    if (any(is.na(diag(pattern))))
        stop("pattern matrix must have parameters in the diagonal")
    if (dim(pattern)[1L] != dim(pattern)[1L])
        stop("pattern matrix must be square matrix")
    if (is.null(dimnames(pattern)))
        stop("pattern must have landmark names as dimnames")
    invisible(pattern)
}
.check_pattern <- function(object, pattern) {
    mode(pattern) <- "character"
    .check_pattern1(pattern)
    if (!all(landmarks(object) %in% rownames(pattern)) ||
        !all(landmarks(object) %in% colnames(pattern)))
        stop("dimnames of patterm must match landmark names")
    pattern <- pattern[rownames(object$SigmaKstar),
        colnames(object$SigmaKstar)]
    ## check spareness
    #UNK <- nlevels(factor(pattern))
    UNK <- sum(!is.na(pattern))
    MAX <- nrow(object$SigmaKstar)*(nrow(object$SigmaKstar)-1)/2
    if (UNK > MAX)
        stop(sprintf(
            "number of nonzero cells (%s) in pattern matrix must be <= %s",
            UNK, MAX))
    invisible(pattern)
}

read_pattern <- function(file, ...) {
    isCSV <- endsWith(file, ".csv")
    x <- if (isCSV) {
        read.csv(file, ...)
    } else {
        suppressMessages(read_excel(file, ...))
    }
    x <- as.data.frame(x)
    rownames(x) <- as.character(x[,1L])
    x[,1L] <- NULL
    x <- as.matrix(x)
    colnames(x) <- rownames(x)
    mode(x) <- "character"
    .check_pattern1(x)
    x
}

SigmaK_fit <- function(object, pattern, check_pattern=TRUE, ...) {
    if (check_pattern)
        pattern <- .check_pattern(object, pattern)
    o <- .SigmaK_fit(object$SigmaKstar, object$H, pattern, ...)
    object$SigmaK <- o$SigmaK
    o$SigmaK <- NULL
    object$pattern <- pattern
    object$results <- o
    if (!is.null(object$boot)) {
        for (i in seq_along(object$boot)) {
            object$boot[[i]][["SigmaK"]] <- .SigmaK_fit(
                object$boot[[i]][["SigmaKstar"]],
                object$boot[[i]][["H"]], pattern)$SigmaK#, ...)
        }
    }
    dimnames(object$SigmaK) <- dimnames(object$SigmaKstar)
    object$call <- match.call()
    class(object) <- c("edma_fit_p", "edma_fit", "edma_data")
    object
}

## print parametric fit object
print.edma_fit_p <- function(x, truncate=40, ...) {
    cat("EDMA parametric fit: ",
        .shorten_name(x$name, truncate), "\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        ncol(x$data[[1L]]), " dimensions, ",
        nrow(x$data[[1L]]), " landmarks, ",
        length(x$data), " replicates, ",
        if (length(x$boot))
            paste(length(x$boot), "bootstrap runs") else "no bootstrap",
        sep="")
    invisible(x)
}

## extract SigmaK estimate
SigmaK <- function (object, ...) UseMethod("SigmaK")
SigmaK.edma_fit_p <- function (object, ...) object[["SigmaK"]]

## evaluates sensitivity:
## par_* are parameters according to pattern matrix
## value is the loss function value evaluated at par_* from optim
sensitivity <- function (object, ...) UseMethod("sensitivity")
## this works because initial values are random, so replicate is fine
sensitivity.edma_fit_p <- function (object, m=10, ...) {
    if (m < 1)
        stop("m must be > 1")
    f <- function() {
        o <- .SigmaK_fit(object$SigmaKstar, object$H, object$pattern,
            method=object$results$method, control=object$results$control)
        unname(c(o$par, o$value))
    }
    out <- rbind(unname(c(object$results$par, object$results$value)),
        t(replicate(m, f())))
    colnames(out) <- c(paste0("par_", names(object$results$par)), "value")
    out
}

