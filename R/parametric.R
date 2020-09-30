## take design matrix and returns a factor for parameter matching
## NA --> 0
## other unique values will be factor levels
.mat2fac <- function(m) {
    K <- nrow(m)
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
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

.SigmaK_fit_old <- function(SigmaKstar, pattern, init,
method = "Nelder-Mead", control = list()) {
    K <- nrow(SigmaKstar)
    I <- diag(1, K)
    ones <- array(rep(1, K), c(1, K))
    H <- I - (1/K) * crossprod(ones, ones)
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
        stop("inits for diagonal elements must be > 0.")
    num_max <- .Machine$double.xmax^(1/3)
    fun <- function(parms){
        SigmaK <- .vec2mat(parms, fac)
        if (any(diag(SigmaK) <= 0))
            return(num_max)
        sum((SigmaKstar - (H %*% SigmaK %*% H))^2)
    }
    if (!is.null(control$fnscale) && control$fnscale < 0)
        stop("control$fnscale can not be negative.")
    o <- suppressWarnings({
        optim(init0, fun, method=method, control=control, hessian=FALSE)
    })
    o$init <- init0
    o$SigmaK <- .vec2mat(o$par, fac)
    dimnames(o$SigmaK) <- dimnames(SigmaKstar)
    o$method <- method
    o$control <- control
    o
}

make_Sigma <- function(params, pattern) {
    out <- .vec2mat(params, .mat2fac(pattern))
    K <- nrow(out)
    dimnames(out) <- list(paste0("L", seq_len(K)), paste0("L", seq_len(K)))
    out
}

.make_A <- function(pattern) {
    K <- as.integer(nrow(pattern))
    L <- cbind(rep(-1, K-1), diag(1, K-1, K-1))
    A.0 <- kronecker(L, L)
    Index.matrix0 <- matrix(seq(1:(K-1)^2),c(K-1,K-1))
    Index.vec0 <- 0
    for (i in 1:K-1){
        for (j in 1:K-1){
            if (i < j) {
                Index.vec0 <- c(Index.vec0, Index.matrix0[i,j])
            }
        }
    }
    Index.vec0 <- Index.vec0[-1]
    A.1 <- A.0[-Index.vec0,]
    Index.matrix1 <- matrix(seq(1:K^2),c(K,K))
    Index.mat1 <- c(0, 0)
    for (i in 1:K){
        for (j in 1:K){
            if (i < j) {
                Index.mat1 <- rbind(
                    Index.mat1,
                    c(Index.matrix1[i,j], Index.matrix1[j,i]))
            }
        }
    }
    Index.mat1 <- Index.mat1[-1,]
    A.2 <- A.1
    A.2[,Index.mat1[,2]] <- A.1[,Index.mat1[,2]] + A.1[,Index.mat1[,1]]
    A.2 <- A.2[,-Index.mat1[,1]]
    if (ncol(A.2) != K*(K+1)/2)
        stop("Something went wrong with constructing A.")
    Index.nonzero <- which(!is.na(pattern[lower.tri(pattern,diag=TRUE)]))
    A <- A.2[,Index.nonzero]
    A
}

.SigmaK_fit_full <- function(SigmaKstar, pattern) {
    K <- nrow(SigmaKstar)
    L <-  cbind(rep(-1, K-1), diag(1, K-1, K-1))
    A <- .make_A(pattern)

    IDs <- ifelse(!is.na(pattern), matrix(seq_len(length(pattern)), K, K), NA)
    IDs <- IDs[lower.tri(IDs, diag=TRUE)]
    IDs <- IDs[!is.na(IDs)]

    SigmaKc <- L %*% SigmaKstar %*% t(L)
    SigmaKcvec <- SigmaKc[lower.tri(SigmaKc, diag=TRUE)]
    AA <- try(solve(t(A) %*% A), silent=TRUE)
    if (inherits(AA, "try-error"))
        stop("Pattern matrix leads to non-invertible A matrix.")
    SigmaKhat1 <- AA %*% t(A) %*% SigmaKcvec

    SigmaKhat <- matrix(0, K, K)
    SigmaKhat[IDs] <- SigmaKhat1
    SigmaKhat <- t(SigmaKhat)
    SigmaKhat[IDs] <- SigmaKhat1
    dimnames(SigmaKhat) <- dimnames(SigmaKstar)
    SigmaKhat
}

.SigmaK_fit <- function(SigmaKstar, pattern, init,
method = "Nelder-Mead", control = list(), twostep=FALSE, check=TRUE) {
    K <- nrow(SigmaKstar)
    if (K < 3)
        stop("K must be at least 3")
    DF <- sum(!is.na(pattern[lower.tri(pattern,diag=TRUE)]))
    if (check && DF > K*(K-1)/2)
        stop(sprintf("Too many unknowns (%s > %s) in pattern matrix.",
            DF, K*(K-1)/2))
    I <- diag(1, K)
    ones <- array(rep(1, K), c(1, K))
    H <- I - (1/K) * crossprod(ones, ones)
    IDs <- ifelse(!is.na(pattern), matrix(seq_len(length(pattern)), K, K), NA)
#    IDs <- IDs[lower.tri(IDs, diag=TRUE)]
    IDs <- IDs[!is.na(IDs)]
    IDd <- which(diag(1, K, K) > 0)
    isDiag <- IDs %in% IDd
    names(IDs) <- pattern[IDs]
    uni <- !duplicated(names(IDs))
    attr(IDs, "diag") <- isDiag
    attr(IDs, "unique") <- uni
    isDiagUni <- isDiag[uni]

    SigmaKfull <- try(.SigmaK_fit_full(SigmaKstar, pattern), silent=TRUE)
    if (inherits(SigmaKfull, "try-error") && check)
        stop(as.character(SigmaKfull))
    if (inherits(SigmaKfull, "try-error")) {
        if (twostep)
            stop(as.character(SigmaKfull))
        SigmaKfull <- NULL
        parms_full <- rep(1, length(IDs))
    } else {
        parms_full <- SigmaKfull[IDs]
    }

    parms_uni <- parms_full[uni]
    names(parms_uni) <- names(IDs)[uni]
    for (i in names(parms_uni))
        parms_uni[i] <- mean(parms_full[names(IDs) == i])
    if (missing(init))
        init <- parms_uni
    if (length(init) != length(parms_uni))
        stop(sprintf("init length must be %s.", length(parms_uni)))
    init <- init[names(parms_uni)]
    if (any(is.na(init)))
        stop("Names of init must match pattern entries.")
    num_max <- .Machine$double.xmax^(1/3)
    if (!is.null(control$fnscale) && control$fnscale < 0)
        stop("control$fnscale can not be negative.")
    TMP <- matrix(0, K, K)
    re_match <- match(names(IDs), names(parms_uni))
    fun_c <- function(parms) {
        if (any(parms[isDiagUni] <= 0))
            return(num_max)
        sum((parms - parms_uni)^2)
    }
    fun_u <- function(parms) {
        if (any(parms[isDiagUni] <= 0))
            return(num_max)
        SigmaK <- TMP
        SigmaK[IDs] <- parms[re_match]
        sum((SigmaKstar - (H %*% SigmaK %*% H))^2)
    }
    fun <- if (twostep)
        fun_c else fun_u
    o <- suppressWarnings({
        optim(init, fun, method=method, control=control, hessian=FALSE)
    })

    names(o$par) <- names(parms_uni)
    parms_back <- o$par[re_match]

    SigmaKhat <- TMP
    SigmaKhat[IDs] <- o$par[re_match]

    o$init <- init
    o$SigmaK <- SigmaKhat
    o$SigmaKfull <- SigmaKfull
    o$method <- method
    o$control <- control
    o$id <- IDs
    o$twostep <- twostep
    o$check <- check
    o
}

.check_pattern1 <- function(pattern) {
    pattern1 <- pattern
    pattern1[] <- NA
    diag(pattern1) <- diag(pattern)
    pattern0 <- pattern
    diag(pattern0) <- NA
    fac <- .mat2fac(pattern) # all parms
    lev1 <- levels(.mat2fac(pattern1)) # parms in diag
    lev0 <- levels(.mat2fac(pattern0)) # parms off-diag
    if (length(intersect(lev1, lev0)) > 0)
        stop("Diagonal and off-diagonal parameters must not overlap.")
    utri <- pattern[upper.tri(pattern)]
    ltri <- t(pattern)[upper.tri(pattern)]
    if (!all(is.na(utri) == is.na(ltri)))
        stop("NA's (0's) are not symmertic in pattern matrix.")
    if (!all(utri[!is.na(utri)] == ltri[!is.na(ltri)]))
        stop("Non-NA's (non-0's) are not symmertic in pattern matrix.")
    pattern[upper.tri(pattern)] <- NA
    pattern[upper.tri(pattern)] <- t(pattern)[upper.tri(pattern)]
    if (any(is.na(diag(pattern))))
        stop("Pattern matrix must have parameters in the diagonal.")
    if (dim(pattern)[1L] != dim(pattern)[1L])
        stop("Pattern matrix must be square matrix.")
    if (is.null(dimnames(pattern)))
        stop("Pattern must have landmark names as dimnames.")
    invisible(pattern)
}
.check_pattern <- function(object, pattern) {
    mode(pattern) <- "character"
    .check_pattern1(pattern)
    if (!all(landmarks(object) %in% rownames(pattern)) ||
        !all(landmarks(object) %in% colnames(pattern)))
        stop("Dimnames of patterm must match landmark names.")
    pattern <- pattern[rownames(object$SigmaKstar),
        colnames(object$SigmaKstar)]
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

SigmaK_fit <- function(object, pattern, twostep=FALSE, check=TRUE, ...) {
    pattern <- .check_pattern(object, pattern)
    o <- .SigmaK_fit(object$SigmaKstar, pattern, twostep=twostep, check=check, ...)
    id <- o$id
    idu <- id[attr(id, "unique")]
    if (!is.null(object$boot)) {
        for (i in seq_along(object$boot)) {
            z <- .SigmaK_fit(
                object$boot[[i]][["SigmaKstar"]],
                pattern,
                twostep=twostep, check=check, ...)
            object$boot[[i]][["SigmaK"]] <- z$SigmaK[idu]
            object$boot[[i]][["SigmaKfull"]] <- z$SigmaKfull[id]
        }
    }
    object$SigmaK <- o$SigmaK
    o$SigmaK <- NULL
    object$SigmaKfull <- o$SigmaK
    o$SigmaKfull <- NULL
    object$twostep <- o$twostep
    o$SigmaK <- NULL
    object$pattern <- pattern
    object$results <- o
    object$call <- match.call()
    class(object) <- c("edma_fit_p", "edma_fit", "edma_data")
    object
}

## print parametric fit object
print.edma_fit_p <- function(x, truncate=40, ...) {
    cat("EDMA parametric fit (", if (x$twostep) "2-step" else "1-step", "): ",
        .shorten_name(x$name, truncate), "\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        nrow(x$data[[1L]]), " landmarks, ",
        ncol(x$data[[1L]]), " dimensions, ",
        length(x$data), " specimens, ",
        if (length(x$boot))
            paste(length(x$boot), "bootstrap runs") else "no bootstrap",
        sep="")
    invisible(x)
}

## extract SigmaK estimate
SigmaK <- function (object, ...) UseMethod("SigmaK")
SigmaK.edma_fit_p <- function (object, ...) object[["SigmaK"]]

## extract SigmaK estimate
SigmaKfull <- function (object, ...) UseMethod("SigmaKfull")
SigmaKfull.edma_fit_p <- function (object, ...) object[["SigmaKfull"]]

## evaluates sensitivity:
## par_* are parameters according to pattern matrix
## value is the loss function value evaluated at par_* from optim
sensitivity <- function (object, ...) UseMethod("sensitivity")

sensitivity.edma_fit_p <- function (object, m=10, scale=10, ...) {
    if (m < 1)
        stop("m must be > 1")
    init0 <- object$results$init
    f <- function() {
        o <- .SigmaK_fit(object$SigmaKstar, object$pattern,
            method=object$results$method, control=object$results$control,
            twostep=object$twostep, check=object$results$check,
            init=init0*runif(length(init0), 0.001, scale))
        unname(c(o$par, o$value))
    }
    out <- rbind(unname(c(object$results$par, object$results$value)),
        t(replicate(m, f())))
    colnames(out) <- c(paste0("par_", names(object$results$par)), "value")
    out
}

