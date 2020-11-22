## Total Least Squares
## use TLS to get scaling factor C
## Y=vec(FM1), X=vec(FM2)
.tlsXY <- function(X, Y) {
    if (missing(Y)) {
        XY <- X
    } else {
        XY <- cbind(X, Y)
    }
    lambda <- tail(eigen(t(XY) %*% XY)$values, 1L)
    as.numeric(solve(t(X) %*% X - lambda) %*% t(X) %*% Y)
}

## d_{ij,1}=c*d_{ij,2} for some c > 0 and for all {ij}
## Now FM is S: S2=C1*FM2, S1=C2*FM1
## C1=1, C2=d2_ij/d1_ij or based on TLS
## Shape difference matrix: S1-S2
.get_sdm <- function(M1, M2, log=TRUE, size=TRUE, C2=NULL) {
    S1 <- as.numeric(dist(M1))
    S2 <- as.numeric(dist(M2))
    ## C1 = 1
    ## C2 is supplied or estimate by TLS
    ## size=FALSE --> C2=1
    if (size) {
        C2 <- if (is.null(C2))
            .tlsXY(S1, S2) else C2
    } else {
        C2 <- 1
    }
    S1 <- C2 * S1
    if (log) {
        S1 <- log(S1)
        S2 <- log(S2)
    }
    SDM <- S1 - S2
    Range <- range(SDM)
    Zval <- Range[which.max(abs(Range))]
    list(sdm=SDM, Zval=Zval, Cval=C2)
}

## This function implements Lele & Cole's SDM test
.edma_sdm <- function(f1, f2, log=TRUE, size=TRUE, C2=NULL) {
    if (is.null(f1$boot) || is.null(f2$boot))
        stop("SDM requires bootstrapped EDMA fit objects")
    B <- min(length(f1$boot), length(f2$boot))
    res <- c(list(.get_sdm(Meanform(f1), Meanform(f2),
                           log=log, size=size, C2=C2)),
        lapply(seq_len(B), function(i) {
            .get_sdm(f1$boot[[i]]$M, f2$boot[[i]]$M,
                     log=log, size=size, C2=C2)
        }))
    SDM <- sapply(res, "[[","sdm")
    Zval <- sapply(res, "[[", "Zval")
    Cval <- sapply(res, "[[", "Cval")
    list(sdm=SDM, Zval=Zval, Cval=Cval)
}

## shape difference matrix with bootstrap
edma_sdm <- function(a, b, log=TRUE, size=TRUE, edge=NULL) {
    .compare_objects(a, b)
    if (is.null(edge)) {
        C2 <- NULL
    } else {
        if (is.numeric(edge))
            edge <- landmarks(a)[edge]
        edge <- as.character(edge)
        if (length(edge) != 2L)
            stop("Provide 2 landmarks for edge argument")
        da <- stack(dist(Meanform(a)))
        db <- stack(dist(Meanform(b)))
        i <- which(da$row == edge[1L] & da$col == edge[2L])
        if (!length(i))
            i <- which(da$row == edge[2L] & da$col == edge[1L])
        ## S1 = C2*FM1=(d2_ij/d1_ij)*FM1
        C2 <- db$dist[i]/da$dist[i]
    }
    d <- stack(dist(Meanform(a)))[,1:2]
    res <- .edma_sdm(a, b, log=log, size=size, C2)
    d$sdm <- res$sdm[,1L]
    out <- list(
        call=match.call(),
        a=a,
        b=b,
        dm=d,
        log=log,
        size=size,
        edge=edge,
        C2=C2,
        B=length(res$Zval)-1L,
        boot=res)
    class(out) <- c("edma_sdm", "edma_dm", class(out))
    out
}

print.edma_sdm <- function(x, level = 0.95, ...) {
    a <- c((1-level)/2, 1-(1-level)/2)
    Zci <- quantile(x$boot$Zval, a)
    Cci <- quantile(x$boot$Cval, a)
    cat("EDMA shape difference matrix\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", x$B, " bootstrap runs",
        if (x$log) " (difference of logarithms)" else "",
        "\n\n", sep="")
    print(rbind("Z (shape)"=Zci, "C (scale)"=Cci),
        digits=getOption("digits")-2L, ...)
    invisible(x)
}

## CI based on the 2x input object boot sample
confint.edma_sdm <- function (object, parm, level=0.95, ...) {
    d <- object$dm
    if (missing(parm))
        parm <- seq_len(nrow(d))
    a <- c((1-level)/2, 1-(1-level)/2)
    out <- t(apply(object$boot$sdm, 1, quantile, a))
    if (object$B < 1)
        out[] <- NA
    rownames(out) <- paste0(as.character(d$row), "-", as.character(d$col))
    out[parm,,drop=FALSE]
}

## this pulls out the stacked form difference matrix
get_sdm <- function (object, ...) UseMethod("get_sdm")
get_sdm.edma_sdm <- function (object, sort=FALSE,
level = 0.95, ...) {
    out <- object$dm
    ci <- confint(object, level=level)
    out$lower <- ci[,1L]
    out$upper <- ci[,2L]
    if (sort)
        out <- out[order(out$sdm, ...),]
    class(out) <- c("fdm", class(out))
    attr(out, "level") <- level
    out
}

Z_test <- function (object, ...) UseMethod("Z_test")
Z_test.edma_sdm <- function(object, level = 0.95, ...) {
    a <- c((1-level)/2, 1-(1-level)/2)
    Zci <- quantile(object$boot$Zval, a)
    Cci <- quantile(object$boot$Cval, a)
    cat("Bootstrap based EDMA Z-test\n", object$B,
        " bootstrap runs\n\n", sep="")
    print(rbind("Z (shape)"=Zci, "C (scale)"=Cci),
         digits=getOption("digits")-3, ...)
    invisible(object)
}

landmarks.edma_sdm <- function(x, ...)
    landmarks(x$a)
dimensions.edma_sdm <- function(x, ...)
    dimensions(x$a)

## influential landmarks
.influence2 <- function(i, object, statistic=c("Z", "C")) {
    ls <- landmarks(object)
    names(ls) <- ls
    i <- ls[i]
    lsd <- ls[!(ls %in% i)]
    a <- object$a
    b <- object$b
    a$M <- a$M[lsd,]
    b$M <- b$M[lsd,]
    for (j in seq_len(object$B)) {
        a$boot[[j]]$M <- a$boot[[j]]$M[lsd,]
        b$boot[[j]]$M <- b$boot[[j]]$M[lsd,]
    }
    val <- .edma_sdm(a, b, log=object$log)
    val[[paste0(statistic, "val")]]
}

## this is the quick version with CIs (no refitting)
get_influence.edma_sdm <- function (object, statistic=c("Z", "C"), level=0.95, ...) {
    statistic <- match.arg(statistic)
    ls <- landmarks(object)
    val <- if (statistic == "Z")
        object$boot$Zval[1L] else object$boot$Cval[1L]
    vals <- lapply(ls, .influence2, object=object, statistic=statistic)
    a <- c((1-level)/2, 1-(1-level)/2)
    vals <- t(sapply(vals, function(z) c(z[1L], quantile(z, a))))
    colnames(vals) <- c(paste0(statistic, "drop"), "lower", "upper")
    out <- data.frame(landmark=ls, vals)
    rownames(out) <- NULL
    attr(out, paste0(statistic, "val")) <- val
    attr(out, "level") <- level
    attr(out, "quick") <- FALSE
    attr(out, "statistic") <- statistic
    attr(out, "null") <- if (statistic == "Z") 0 else 1
    class(out) <- c("edma_influence", class(out))
    out
}

