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

## d_{ij,A}=c*d_{ij,B} for some c > 0 and for all {ij}
## Now FM is S: S1=FM1, S2=Cval*FM2
## Shape difference matrix: S1-S2
.get_sdm <- function(M1, M2, log=TRUE) {
    S1 <- as.numeric(dist(M1))
    S2 <- as.numeric(dist(M2))
    C2 <- .tlsXY(S1, S2) # C1 = 1
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
.edma_sdm <- function(f1, f2, log=TRUE) {
    if (is.null(f1$boot) || is.null(f2$boot))
        stop("SDM requires bootstrapped EDMA fit objects")
    B <- min(length(f1$boot), length(f2$boot))
    res <- c(list(.get_sdm(Meanform(f1), Meanform(f2), log=log)),
        lapply(seq_len(B), function(i) {
            .get_sdm(f1$boot[[i]]$M, f2$boot[[i]]$M, log=log)
        }))
    SDM <- sapply(res, "[[","sdm")
    Zval <- sapply(res, "[[", "Zval")
    Cval <- sapply(res, "[[", "Cval")
    list(sdm=SDM, Zval=Zval, Cval=Cval)
}

## Approximate composite likelihood approach
## for testing shape difference
##
## log-Composite likelihood function:
## This is based on the squared Euclidean distances.
## But it could be based on centered inner product matrix as well.
## We use the result that non-central chi-square distribution can be
## approximated by a normal distribution when the non-centrality
## parameter is large as compared to the variance.
## The results are easy to apply to any dimension.
## We will implement it to 2 and 3 dimensional objects.
##
## Notation:
## delta_lm : squared Euclidean distance between landmarks l and m.
## phi_lm : sigma_ll + sigma_mm - 2*sigma_lm.
##   These are the entries in SigmaK matrix.
##   SigmaD is assumed to be identity.
## e_lm : Squared Euclidean distance between landmarks l and m
##
## Under the Gaussian perturbation model, e_lm is approximately
## Normally distributed.
##
## 2 dimensional object:
##   Mean = 2*phi_lm + delta_lm
##   Variance = 4 * phi_lm^2 + 4 * delta_lm * phi_lm
##
## 3 dimensional object:
##   Mean = 3*phi_lm + delta_lm = alpha_1
##   Variance = 6 * phi_lm^2 + 4 * delta_lm * phi_lm = alpha_2
##
## We can use the method of moments to estimate these parameters.
## We can also (potentially) improve them by one-step maximum
## composite likelihood estimator.
##
## Notation: Let alpha_1_hat and alpha_2_hat denote the empirical
## mean and empirical variance of the observed squared Euclidean distances.
##
## 2 dimensional object:
##   delta_lm_hat = sqrt{(alpha_1_hat^2) - alpha_2_hat}
##   Be aware that for small samples, this can be negative.
##   phi_lm_hat = 0.5*(alpha_1_hat - delta_lm_hat)
##
## 3 dimensional object:
##   delta_lm_hat = sqrt{(alpha_1_hat^2) - 1.5*alpha_2_hat}
##   Be aware that for small samples, this can be negative.
##   phi_lm_hat = 0.33*(alpha_1_hat - delta_lm_hat)
##
## EDMAinR spits out the relevant information, namely,
## squared distances for the objects, their means and variances.
.edma_sdm_clr <- function(data1,data2) {
    ## Fit the data
    #fit1_edma = edma_fit(data1,B=1)
    #fit2_edma = edma_fit(data2,B=1)
    ## This is needed to estimate 'C_hat'
    #fit_sdm = edma_sdm(fit1_edma,fit2_edma)
    #C_hat = fit_sdm$boot$Cval[1]
    z_1 <- .edma_fit_np(data1, less=FALSE)
    z_2 <- .edma_fit_np(data2, less=FALSE)
    C_hat <- .get_sdm(z_1$M, z_2$M)$Cval
    K <- dim(data1)[1L]
    D <- dim(data1)[2L]
    n1 <- dim(data1)[3L]
    n2 <- dim(data2)[3L]
    ## Need to put these into a matrix with rows
    ## as the squared distances and columns as observations.
    ## There might be way to avoid this do loop.
    ## These are the observations for writing the composite likelihood.
    Eu_1 <- matrix(0, length(dist(z_1$M)), n1)
    for (i in seq_len(n1)) {
        Eu_1[,i] <- z_1$EuX[,,i][lower.tri(z_1$EuX[,,i])]
    }
    Eu_2 <- matrix(0, length(dist(z_2$M)), n2)
    for (i in seq_len(n2)) {
        Eu_2[,i] <- z_2$EuX[,,i][lower.tri(z_2$EuX[,,i])]
    }
    ## Estimated parameters for writing the composite likelihood
    delta_hat_1 <- as.vector(dist(z_1$M)^2)
    delta_hat_2 <- as.vector(dist(z_2$M)^2)
    phi_hat_1 <- (D/2)*(z_1$EuMean[lower.tri(z_1$EuMean)] - delta_hat_1)
    phi_hat_2 <- (D/2)*(z_2$EuMean[lower.tri(z_2$EuMean)] - delta_hat_2)
    ## Approximate composite likelihood can be computed using
    ## the 'dnorm' function or exact composite likelihood
    ## using the 'dchisq' function.
    ##
    ## Normal approximation based composite likelihood
    mean_1 <- D*phi_hat_1 + delta_hat_1
    var_1 <- 2*D*phi_hat_1^2 + 4*delta_hat_1*phi_hat_1
    CL_1 <- sum(dnorm(Eu_1,mean=mean_1,sd=sqrt(var_1),log=T))
    mean_2 <- D*phi_hat_2 + delta_hat_2
    var_2 <- 2*D*phi_hat_2^2 + 4*delta_hat_2*phi_hat_2
    CL_2 <- sum(dnorm(Eu_2,mean = mean_2,sd = sqrt(var_2),log=T))
    ## Constrained CL with the mean forms are proportional to each other.
    mean_2_const <- D*phi_hat_2 + (C_hat^2)*delta_hat_1
    var_2_const <- 2*D*phi_hat_2^2 + 4*(C_hat^2)*delta_hat_1*phi_hat_2
    CL_2_const <- sum(dnorm(Eu_2,
                            mean = mean_2_const,
                            sd = sqrt(var_2_const),
                            log=TRUE))
    ## Testing for shape difference can be done using the fact
    ## that under the null hypothesis M_2 = c* M_1,
    ## the corresponding parameters of the non-central chisquare
    ## can be computed. It only affects delta_lm.
    ## The method of moments estimators are available under
    ## the full model and restricted model.
    ## We can use parametric bootstrap to obtain the null distribution.
    CLR <- 2*(CL_2 - CL_2_const)
    list(CLR=CLR,
        C_hat=C_hat,
        M1=z_1$M,
        n1=n1,
        n2=n2,
        SigmaK1star=z_1$SigmaKstar,
        SigmaK2star=z_2$SigmaKstar)
}

## Now put a wrapper to do parametric bootstrap.
## This function only does 1 run
## x is list returned by .edma_sdm_clr
.edma_sdm_clr_null1 <- function(x) {
    M1star <- x$M1
    C_hat <- x$C_hat
    SigmaK1star <- x$SigmaK1star
    SigmaK2star <- x$SigmaK2star
    n1 <- x$n1
    n2 <- x$n2
    ## Generate the data.
    ## The input is the estimates from fitting the original data.
    ## SigmaKstar's are singular matrices.
    ##
    ## Generate centered matrix normal variates.
    K <- nrow(M1star)
    D <- ncol(M1star)
    x1.c <- MixMatrix::rmatrixnorm(n1,
        mean=M1star[-K,],
        U=SigmaK1star[1:K-1,1:K-1],
        V=diag(1, D))
    tmp <- array(rep(0,K*D*n1),c(K,D,n1))
    for (i in seq_len(n1)) {
        tmp[,,i] <- rbind(x1.c[,,i], -apply(x1.c[,,i],2,sum))
    }
    x1.c <- tmp
    ## This is the form matrix under the null.
    M2star = C_hat * M1star
    x2.c <- MixMatrix::rmatrixnorm(n1,
        mean=M2star[-K,],
        U=SigmaK2star[1:K-1,1:K-1],
        V=diag(1,D))
    tmp <- array(rep(0,K*D*n2), c(K,D,n2))
    for (i in seq_len(n2)) {
        tmp[,,i] <- rbind(x2.c[,,i],-apply(x2.c[,,i],2,sum))
    }
    x2.c <- tmp
    # Now input these to the CLR function.
    .edma_sdm_clr(as.edma_data(x1.c), as.edma_data(x2.c))$CLR
}

## CLR shape difference with parametric bootstrap
edma_sdm_clr <- function(a, b, B=0) {
    .compare_data(a, b)
    res <- .edma_sdm_clr(a, b)
    boot <- NULL
    if (B > 0) {
        boot <- pbapply::pbreplicate(B,
            .edma_sdm_clr_null1(res))
    }
    out <- list(
        call=match.call(),
        a=a,
        b=b,
        results=res,
        pvalue=if (B > 0) sum(boot > res$CLR) / B else NA,
        B=B,
        boot=boot)
    class(out) <- c("edma_sdm_clr", "edma_dm", class(out))
    out
}

CLR_test <- function (object, ...) UseMethod("CLR_test")
CLR_test.edma_sdm_clr <- function (object, ...) {
    METHOD <- "Parametric bootstrap based CLR shape test"
    CLR <- object$results$CLR
    PVAL <- object$pvalue
    PARAMETER <- length(object$boot)
    names(CLR) <- "CLR"
    names(PARAMETER) <- "B"
    out <- list(
        statistic = CLR,
        parameter = PARAMETER,
        p.value = PVAL,
        method = METHOD,
        data.name = "shape difference matrix", Null=object$boot)
    class(out) <- c("edma_CLRtest", "edma_test", "htest")
    out
}


# think about using htest
print.edma_sdm_clr <- function(x, ...) {
    d <- getOption("digits")-2L
    pval <- if (!is.null(x$boot)) {
        paste(", p-value =", format(x$pvalue, digits=d))
    } else ""
    cat("EDMA shape difference matrix\n",
        "Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", x$B, " parametric bootstrap runs\n",
        "CLR = ", format(x$results$CLR, digits=d),
        pval,
        "\n", sep="")
    invisible(x)
}


## shape difference matrix with bootstrap
edma_sdm <- function(a, b, log=TRUE) {
    .compare_objects(a, b)
    d <- stack(as.dist(a))[,1:2]
    res <- .edma_sdm(a, b, log=log)
    d$sdm <- res$sdm[,1L]
    out <- list(
        call=match.call(),
        a=a,
        b=b,
        dm=d,
        log=log,
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

