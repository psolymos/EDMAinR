## compare data sets
## a and b are edma_fit objects
.compare_data <- function (a, b, text="", ...) {
    if (nrow(a$data[[1L]]) != nrow(b$data[[1L]]))
        stop(paste(text, "Number of landmarks must be identical."))
    if (ncol(a$data[[1L]]) != ncol(b$data[[1L]]))
        stop(paste(text, "Number of dimensions must be identical."))
    if (!all(rownames(a$data[[1L]]) == rownames(b$data[[1L]])))
        stop(paste(text, "Landmark names and ordering must be identical."))
    if (!all(colnames(a$data[[1L]]) == colnames(b$data[[1L]])))
        stop(paste(text, "Dimension names and ordering must be identical."))
    invisible(NULL)
}
## compare fit objects (imitating multiple dispatch for S3)
.compare_objects <- function (a, b, text="", ...) {
    if (!inherits(b, "edma_fit") || !inherits(a, "edma_fit"))
        stop("Input must be edma_fit object.")
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

## this does the fixed or mixed test once,
## d1: numerator, d2: denominator
## d1 and d2 are data objects
.Ttest_data <- function(d1, d2, ref_denom=TRUE, mix=FALSE) {
    n1 <- dim(d1)[3L]
    n2 <- dim(d2)[3L]
    s1 <- specimens(d1)
    s2 <- specimens(d2)
    if (mix) {
        ii <- c(-seq_len(n1), seq_len(n2))
        ii <- sample(ii, n1+n2, replace=TRUE)
        i1 <- ii[seq_len(n1)]
        i2 <- ii[seq_len(n2)+n1]
        d1$data <- c(d1$data[-i1[i1 < 0]], d2$data[i1[i1 > 0]])
        d2$data <- c(d1$data[-i2[i2 < 0]], d2$data[i2[i2 > 0]])
    } else {
        if (ref_denom) {
            d2$data <- d1$data[sample(seq_len(n1), n2, replace=TRUE)]
            d1$data <- d1$data[sample(seq_len(n1), n1, replace=TRUE)]
        } else {
            d1$data <- d2$data[sample(seq_len(n2), n1, replace=TRUE)]
            d2$data <- d2$data[sample(seq_len(n2), n2, replace=TRUE)]
        }
    }
    names(d1$data) <- s1
    names(d2$data) <- s2
    fit1 <- edma_fit(d1)
    fit2 <- edma_fit(d2)
    .formdiff(Meanform(fit1), Meanform(fit2))
}

.Ttest_fit <- function (d1, d2, B=0, ref_denom=TRUE, mix=FALSE,
ncores=getOption("Ncpus", 1L), ...) {
    fd <- .formdiff(Meanform(d1), Meanform(d2))
    d1 <- as.edma_data(d1)
    d2 <- as.edma_data(d2)
    .compare_data(d1, d2)
    bfd <- NULL
    if (B > 0) {
        ncores <- as.integer(ncores)
        if (ncores > 1L) {
            ncores <- min(ncores, parallel::detectCores(TRUE), na.rm=TRUE)
            cl <- parallel::makeCluster(ncores)
            on.exit(parallel::stopCluster(cl))
        } else {
            cl <- NULL
        }
        bfd <- pbapply::pbsapply(seq_len(B),
            function(i, d1, d2, mix, ref_denom) {
            as.numeric(
                EDMAinR::.Ttest_data(
                    d1=d1, d2=d2, mix=mix, ref_denom=ref_denom))
        }, d1=d1, d2=d2, mix=mix, ref_denom=ref_denom, cl=cl)
    }
    fds <- cbind(fd, bfd)
    Tval <- apply(fds, 2, function(z) max(z)/min(z))
    attr(fds, "Tval") <- unname(Tval)
    attr(fds, "mix") <- mix
    fds
}

## form difference matrix with bootstrap
edma_fdm <- function(numerator, denominator, B=0, ref_denom=TRUE, mix=FALSE) {
    .compare_objects(numerator, denominator)

    if (B > 0 && (is.null(numerator$boot) || is.null(denominator$boot)))
    stop("B > 0 requires bootstrapped EDMA fit objects")
    f1 <- if (ref_denom)
        numerator else denominator
    f2 <- if (ref_denom)
        denominator else numerator
    Bx <- min(length(f1$boot), length(f2$boot))
    #if (B > min(length(f1$boot), length(f2$boot)))
    #    stop(sprintf("B=%s is too large, inputs have %s and %x runs.",
    #                 B, length(f1$boot), length(f2$boot)))

    res <- cbind(.formdiff(Meanform(numerator), Meanform(denominator)))
    if (Bx > 0) {
        tmp <- sapply(seq_len(Bx), function(i) {
            .formdiff(numerator$boot[[i]]$M, denominator$boot[[i]]$M)
        })
        res <- cbind(res, tmp)
    }
    fd <- stack(.formdiff(Meanform(numerator), Meanform(denominator)))
    b <- .Ttest_fit(numerator, denominator, B=B, ref_denom=ref_denom, mix=mix)

    out <- list(
        call=match.call(),
        numerator=numerator,
        denominator=denominator,
        ref_denom=ref_denom,
        mix=mix,
        dm=fd,
        B=B,
        boot=b,    # mixed or reference bootstrap for global testing
        boot2=res) # 2-sample bootstrap for local testing (CI)
    class(out) <- c("edma_fdm", "edma_dm", class(out))
    out
}

## this is an internal function for printing similar objects
## FDM, GM
.print_edma_fdm <- function(x, title="EDMA", truncate=20, ...) {
    H <- Tobs_test(x)
    if (x$B) {
        fp <- format.pval(H$p.value)
        fp <- if (substr(fp, 1L, 1L) == "<")
            paste0(", p ", fp) else paste0(", p = ", fp)
    } else {
        fp <- ""
    }
    cat(title, "\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n",
        if (x$B)
            paste(x$B, "bootstrap runs") else "no bootstrap",
        " (ref: ", if (x$ref_denom) "denominator" else "numerator", ")",
        "\nTobs = ", format(H$statistic, digits = max(1L, getOption("digits") - 2L)),
        fp,
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
level = 0.95, what = "all", ...) {
    what <- match.arg(what, c("all", "less", "greater", "signif", "nonsignif"))
    out <- object$dm
    ci <- confint(object, level=level)
    out$lower <- ci[,1L]
    out$upper <- ci[,2L]
    if (sort)
        out <- out[order(out$dist, ...),]
    if (what == "less")
        out <- out[out$lower < 1 & out$upper < 1,]
    if (what == "greater")
        out <- out[out$lower > 1 & out$upper > 1,]
    if (what == "nonsignif")
        out <- out[out$lower < 1 & out$upper > 1,]
    if (what == "signif")
        out <- out[!(out$lower < 1 & out$upper > 1),]
    class(out) <- c("fdm", class(out))
    attr(out, "level") <- level
    out
}

## T-test for 2 edma_fit objects
Tobs_test <- function (object, ...) UseMethod("Tobs_test")
.Tobs_test <- function (object, DNAME="", ...) {
    METHOD <- if (attr(object$boot, "mix"))
        "Mixed bootstrap based EDMA T-test" else "Bootstrap based EDMA T-test"
    Tval <- attr(object$boot, "Tval")
    PVAL <- sum(Tval[1L] <= Tval[-1L]) / (length(Tval)-1L)
    PARAMETER <- length(Tval) - 1L
    names(Tval) <- "Tobs-value"
    names(PARAMETER) <- "B"
    out <- list(statistic = Tval[1L], parameter = PARAMETER,
        p.value = if (PARAMETER > 0) PVAL else NA,
        method = METHOD, data.name = DNAME, Tvals=Tval[-1L])
    class(out) <- c("edma_Ttest", "edma_test", "htest")
    out
}
Tobs_test.edma_fdm <- function (object, ...)
    .Tobs_test(object, DNAME="form difference matrix", ...)

## CI based on the 2x input object boot sample
confint.edma_dm <- function (object, parm, level=0.95, ...) {
    if (missing(parm))
        parm <- seq_len(nrow(object$dm))
    a <- c((1-level)/2, 1-(1-level)/2)
    out <- t(apply(object$boot2, 1, quantile, a))
    if (object$B < 1)
        out[] <- NA
    rownames(out) <- paste0(as.character(object$dm$row), "-",
        as.character(object$dm$col))
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
            x <- edma_fit(as.edma_data(object$numerator)[lsd,,])
            ref <- edma_fit(as.edma_data(object$denominator)[lsd,,])
        } else {
            x <- edma_fit(as.edma_data(object$denominator)[lsd,,])
            ref <- edma_fit(as.edma_data(object$numerator)[lsd,,])
        }
        Tval <- attr(.formdiff(Meanform(x), Meanform(ref)), "Tval")
    }
    c(Tval, ci)
}

## this is the quick version with CIs (no refitting)
get_influence <- function (object, ...) UseMethod("get_influence")
get_influence.edma_dm <- function (object, level=0.95, ...) {
    quick <- TRUE
    ls <- landmarks(object)
    Tval <- Tobs_test(object)$statistic
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
    attr(out, "statistic") <- "T"
    attr(out, "null") <- 1
    class(out) <- c("edma_influence", class(out))
    out
}
