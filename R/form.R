## compare data sets
## a and b are edma_fit objects
.compare_data <- function (a, b, text="", ...) {
    if (nrow(a$data[[1L]]) != nrow(b$data[[1L]]))
        stop(paste(text, "number of landmarks must be identical"))
    if (ncol(a$data[[1L]]) != ncol(b$data[[1L]]))
        stop(paste(text, "number of dimensions must be identical"))
    if (!all(rownames(a$data[[1L]]) == rownames(b$data[[1L]])))
        stop(paste(text, "landmark names and ordering must be identical"))
    if (!all(colnames(a$data[[1L]]) == colnames(b$data[[1L]])))
        stop(paste(text, "dimension names and ordering must be identical"))
    invisible(NULL)
}
## compare fit objects (imitating multiple dispatch for S3)
.compare_objects <- function (a, b, text="", ...) {
    if (!inherits(b, "edma_fit") || !inherits(a, "edma_fit"))
        stop("input must be edma_fit object")
    .compare_data(a, b, text, ...)
}

## inputs are meanform matrices
.formdiff <- function(a, b) {
    r <- dist(a) / dist(b)
    attr(r, "method") <- "euclidean_distance_ratio"
    attr(r, "call") <- NULL
    attr(r, "Tval") <- max(r) / min(r)
    r
}

## inputs are edma_fit objects
## this returns the distance matrix, not stacked
formdiff <- function (numerator, denominator, ...) {
    .compare_objects(numerator, denominator)
    .formdiff(Meanform(numerator), Meanform(denominator))
}

## calculate T-statistic
## x: object w bootstrap
## ref: reference object (w or w/o bootstrap)
.Ttest_fit <- function(x, ref, boot=TRUE) {
    if (!boot)
        x$boot <- NULL
    FMref <- as.numeric(as.dist(ref))
    FM <- if (is.null(x$boot)) {
        data.matrix(as.numeric(as.dist(x)))
    } else {
        cbind(as.numeric(as.dist(x)),
            sapply(x$boot, function(z) as.numeric(dist(z$M))))
    }
    FDM <- FM/FMref
    Tval <- apply(FDM, 2, function(z) max(z)/min(z))
    attr(FDM, "Tval") <- Tval
    FDM
}

## form difference matrix with bootstrap
edma_fdm <- function(numerator, denominator, ref_denom=TRUE) {
    .compare_objects(numerator, denominator)
    fd <- stack(formdiff(numerator, denominator))
    b <- if (ref_denom) {
        .Ttest_fit(numerator, denominator)
    } else {
        .Ttest_fit(denominator, numerator)
    }
    out <- list(
        call=match.call(),
        numerator=numerator,
        denominator=denominator,
        ref_denom=ref_denom,
        dm=fd,
        B=ncol(b)-1,
        boot=b)
    class(out) <- c("edma_fdm", "edma_dm", class(out))
    out
}

## this is an internal function for printing similar objects
## FDM, GM
.print_edma_fdm <- function(x, title="EDMA", truncate=20, ...) {
    H <- T_test(x)
    cat(title, "\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        if (x$B)
            paste(x$B, "bootstrap runs") else "no bootstrap",
        " (ref: ", if (x$ref_denom) "denominator" else "numerator", ")",
        "\nT=", round(H$statistic, 4),
        if (x$B)
            paste0(", p=", round(H$p.value, 4)) else "",
        sep="")
    invisible(x)
}
## print fmd objects
print.edma_fdm <- function(x, ...) {
    .print_edma_fdm(x, "EDMA form difference matrix", ...)
}

## this pulls out the stacked form difference matrix
get_fdm <- function (object, ...) UseMethod("get_fdm")
get_fdm.edma_fdm <- function (object, sort=FALSE,
level = 0.95, ...) {
    out <- object$dm
    ci <- confint(object, level=level)
    out$lower <- ci[,1L]
    out$upper <- ci[,2L]
    if (sort)
        out <- out[order(out$dist, ...),]
    class(out) <- c("fdm", class(out))
    attr(out, "level") <- level
    out
}

## T-test for 2 edma_fit objects
T_test <- function (object, ...) UseMethod("T_test")
.T_test <- function (object, DNAME="", ...) {
    METHOD <- "Bootstrap based EDMA T-test"
    Tval <- attr(object$boot, "Tval")
    PVAL <- sum(Tval[1L] <= Tval) / length(Tval)
    PARAMETER <- length(Tval) - 1L
    names(Tval) <- "T-value"
    names(PARAMETER) <- "B"
    out <- list(statistic = Tval[1L], parameter = PARAMETER,
        p.value = PVAL, method = METHOD, data.name = DNAME,
        Tvals=Tval)
    class(out) <- c("edma_test", "htest")
    out
}
T_test.edma_fdm <- function (object, ...)
    .T_test(object, DNAME="form difference matrix", ...)

## CI based on the 2x input object boot sample
confint.edma_dm <- function (object, parm, level=0.95, ...) {
    d <- object$dm
    if (missing(parm))
        parm <- seq_len(nrow(d))
    a <- c((1-level)/2, 1-(1-level)/2)
    out <- t(apply(object$boot, 1, quantile, a))
    if (object$B < 1)
        out[] <- NA
    rownames(out) <- rownames(d)
    out[parm,,drop=FALSE]
}

landmarks.edma_fdm <- function(x, ...)
    landmarks(x$numerator)
dimensions.edma_fdm <- function(x, ...)
    dimensions(x$numerator)

## influential landmarks
.influence <- function(i, object, quick=TRUE, level=0.95) {
    ls <- landmarks(object)
    names(ls) <- ls
    i <- ls[i]
    lsd <- ls[!(ls %in% i)]
    ci <- c(NA, NA)
    if (quick) {
        a <- c((1-level)/2, 1-(1-level)/2)
        keep <- !(object$dm$row %in% i) & !(object$dm$col %in% i)
        FDMkeep <- object$boot[keep,,drop=FALSE]
        Tvals <- apply(FDMkeep, 2, function(z) max(z)/min(z))
        Tval <- Tvals[1L]
        if (object$B > 0)
            ci <- unname(quantile(Tvals, a))
    } else {
        if (object$ref_denom) {
            x <- edma_fit(.get_data(object$numerator)[lsd,,])
            ref <- edma_fit(.get_data(object$denominator)[lsd,,])
        } else {
            x <- edma_fit(.get_data(object$denominator)[lsd,,])
            ref <- edma_fit(.get_data(object$numerator)[lsd,,])
        }
        Tval <- attr(.Ttest_fit(x, ref, boot=FALSE), "Tval")
    }
    c(Tval, ci)
}

## this is the quick version with CIs (no refitting)
get_influence <- function (object, ...) UseMethod("get_influence")
get_influence.edma_dm <- function (object, level=0.95, ...) {
    quick <- TRUE
    ls <- landmarks(object)
    Tval <- T_test(object)$statistic
    Tvals <- if (quick) {
        t(sapply(ls, .influence,
            object=object, quick=quick, level=level))
    } else {
        t(pbapply::pbsapply(ls, .influence,
            object=object, quick=quick, level=level))
    }
    colnames(Tvals) <- c("Tdrop", "lower", "upper")
    out <- data.frame(landmark=ls, Tvals)
    rownames(out) <- NULL
    attr(out, "Tval") <- Tval
    attr(out, "level") <- level
    attr(out, "quick") <- quick
    class(out) <- c("edma_influence", class(out))
    out
}
