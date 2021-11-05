## FDM report
edma_fdm_report <- function(numerator, denominator, output="edma_output.txt",
    landmarks=NULL, B=0, level=0.95, ref_denom=TRUE, mix=FALSE,
    digits=4) {

    o0 <- options(max.print=99999)
    on.exit(options(o0), add=TRUE)

    read_fn <- function(x) {
        if (is.character(x)) {
            x <- read_xyz(x)
        }
        if (!inherits(x, "edma_data"))
            stop("Please check your inputs, something isn't right...")
        x
    }
    Cat <- function(...) {
        cat(..., file=output, append=TRUE, sep="\n")
    }
    Mat <- function(mat, ...) {
        mat[upper.tri(mat)] <- 0
        class(mat) <- "table"
        print(mat, digits=digits, zero.print = "", ...)
    }

    NUMERATOR_FILE <- is.character(numerator)
    DENOMINATOR_FILE <- is.character(denominator)
    NUMERATOR_LAB <- if (is.character(numerator))
            numerator else deparse(substitute(numerator))
    DENOMINATOR_LAB <- if (is.character(denominator))
            denominator else deparse(substitute(denominator))
    NUMERATOR_XYZ <- read_fn(numerator)
    DENOMINATOR_XYZ <- read_fn(denominator)
    .compare_data(NUMERATOR_XYZ, DENOMINATOR_XYZ)
    if (!is.null(landmarks)) {
        NUMERATOR_XYZ <- NUMERATOR_XYZ[landmarks,,]
        DENOMINATOR_XYZ <- DENOMINATOR_XYZ[landmarks,,]
    }
    NUMERATOR_FIT <- edma_fit(NUMERATOR_XYZ, B=B)
    DENOMINATOR_FIT <- edma_fit(DENOMINATOR_XYZ, B=B)
    FDM <- edma_fdm(
        numerator = NUMERATOR_FIT,
        denominator = DENOMINATOR_FIT,
        B=B, ref_denom=ref_denom, mix=mix)
    Dis <- .formdiff(
        Meanform(NUMERATOR_FIT),
        Meanform(DENOMINATOR_FIT))
    TV <- global_test(FDM)
    CI1 <- get_fdm(FDM, level=level, sort=FALSE)
    CI2 <- get_fdm(FDM, level=level, sort=TRUE)
    colnames(CI1) <- c("", "", "Estimate", "Low", "High")
    colnames(CI2) <- c("", "", "Estimate", "Low", "High")

    h <- hist(c(TV$statistic, TV$Tvals), plot=FALSE)
    MX <- 50
    NN <- round(MX * h$counts / max(h$counts))
    ST <- sapply(NN, function(z) paste0(rep("*", z), collapse=""))
    ii <- which.min(abs(TV$statistic-h$mids))
    #ST[ii] <- paste0(substr(ST[ii], 1, nchar(ST[ii])-1), "T")
    STEM <- paste0("  ", format(h$mids), "   |", ST)
    STEM[ii] <- gsub("   |", "  T|", STEM[ii], fixed=TRUE)
    ver <- read.dcf(file=system.file("DESCRIPTION", package="EDMAinR"),
        fields=c("Version", "Date"))

    sink(output, append=TRUE)
    on.exit(sink(), add=TRUE)
    op <- options(width = 10000, scipen=999)
    on.exit(options(op), add=TRUE)

    cat("# EDMAinR FDM output (v ", ver[1], " - ", ver[2], ")",
        "\n\nDate created: ", as.character(Sys.time()), "\n\n",
        file=output, sep="")
    Cat(
        "## Settings",
        "",
        paste0("NUMERATOR: ", NUMERATOR_XYZ$name),
        paste0("Number of specimens: ", dim(NUMERATOR_FIT)[3]),
        paste0(if (NUMERATOR_FILE) "Input file: " else "Input object: ",
            NUMERATOR_LAB),
        "",
        paste0("DENOMINATOR: ", DENOMINATOR_FIT$name),
        paste0("Number of specimens: ", dim(DENOMINATOR_FIT)[3]),
        paste0(if (DENOMINATOR_FILE) "Input file: " else "Input object: ",
            DENOMINATOR_LAB),
        "",
        paste0("Number of landmarks: ", dim(NUMERATOR_FIT)[1]),
        paste0("Number of dimensions:  ", dim(NUMERATOR_FIT)[2]),
        "",
        "## Form matrices",
        "",
        "XY(Z) coordinates of the NUMERATOR mean:",
        ""
    )
    print(Meanform(NUMERATOR_FIT), digits=digits)
    Cat(
        "",
        "XY(Z) coordinates of the DENOMINATOR mean:",
        ""
    )
    print(Meanform(DENOMINATOR_FIT), digits=digits)
    Cat(
        "",
        "Further-improved estimate of Sigma K* for the NUMERATOR:",
        ""
    )
    Mat(SigmaKstar(NUMERATOR_FIT))
    Cat(
        "",
        "Further-improved estimate of Sigma K* for the DENOMINATOR:",
        ""
    )
    Mat(SigmaKstar(DENOMINATOR_FIT))
    Cat(
        "",
        "Form matrix for the NUMERATOR:",
        ""
    )
    print(as.dist(NUMERATOR_FIT, diag=TRUE), digits=digits)
    Cat(
        "",
        "Form matrix for the DENOMINATOR:",
        ""
    )
    print(as.dist(DENOMINATOR_FIT, diag=TRUE), digits=digits)
    Cat(
        "",
        "## Form difference matrix",
        ""
    )
    print(as.dist(Dis, diag=TRUE), digits=digits)
    Cat(
        "",
        paste0("Tobs statistic (max/min): ", round(TV$statistic, digits)),
        "",
        paste0("The bootstrap reference sample is the ",
            if (ref_denom) "DENOMINATOR" else "NUMERATOR"),
        paste0("Number of bootstrap samples: ", B),
        paste0("Bootstrap type: ",
            if (mix) "mixed" else "reference is fixed"),
        paste0("Probability: ",
            format(max(TV$p.value, 10^-(digits+1)))),
        "",
        STEM,
        "",
        "Unsorted form-difference matrix with marginal confidence intervals:",
        ""
    )
    print(CI1, digits=4, row.names=FALSE, right=FALSE)
    Cat(
        "",
        paste0("Confidence intervals are based on ", B, " resamples."),
        paste0("Nonparametric resampling was used."),
        paste0("Alpha level: ", 1-level),
        "",
        "Sorted form-difference matrix:",
        ""
    )
    print(CI2, digits=4, row.names=FALSE, right=FALSE)
    Cat(
        "",
        "Significant form differences < 1:",
        ""
    )
    print(CI2[CI2$Low < 1 & CI2$High < 1,], digits=4,
        row.names=FALSE, right=FALSE)
    Cat(
        "",
        "Non-significant form differences:",
        ""
    )
    print(CI2[CI2$Low <= 1 & CI2$High >= 1,], digits=4,
        row.names=FALSE, right=FALSE)
    Cat(
        "",
        "Significant form differences > 1:",
        ""
    )
    print(CI2[CI2$Low > 1 & CI2$High > 1,], digits=4,
        row.names=FALSE, right=FALSE)
    Cat("")

    out <- list(
        NUMERATOR_LAB=NUMERATOR_LAB,
        DENOMINATOR_LAB=DENOMINATOR_LAB,
        NUMERATOR_FIT=NUMERATOR_FIT,
        DENOMINATOR_FIT=DENOMINATOR_FIT,
        FDM=FDM,
        TV=TV,
        CI1=CI1,
        CI2=CI2,
        STEM=STEM
    )
    invisible(out)
}

## GDM report
edma_gdm_report <- function(numerator_yng, numerator_old,
    denominator_yng, denominator_old,
    output="edma_output.txt",
    landmarks=NULL, B=0, level=0.95, ref_denom=TRUE, mix=FALSE,
    digits=4) {

    o0 <- options(max.print=99999)
    on.exit(options(o0), add=TRUE)

    read_fn <- function(x) {
        if (is.character(x)) {
            x <- read_xyz(x)
        }
        if (!inherits(x, "edma_data"))
            stop("Please check your inputs, something isn't right...")
        x
    }
    Cat <- function(...) {
        cat(..., file=output, append=TRUE, sep="\n")
    }
    Mat <- function(mat, ...) {
        mat[upper.tri(mat)] <- 0
        class(mat) <- "table"
        print(mat, digits=digits, zero.print = "", ...)
    }

    NUMERATOR_YNG_FILE <- is.character(numerator_yng)
    NUMERATOR_OLD_FILE <- is.character(numerator_old)
    DENOMINATOR_YNG_FILE <- is.character(denominator_yng)
    DENOMINATOR_OLD_FILE <- is.character(denominator_old)

    NUMERATOR_YNG_LAB <- if (is.character(numerator_yng))
            numerator_yng else deparse(substitute(numerator_yng))
    NUMERATOR_OLD_LAB <- if (is.character(numerator_old))
            numerator_old else deparse(substitute(numerator_old))
    DENOMINATOR_YNG_LAB <- if (is.character(denominator_yng))
            denominator_yng else deparse(substitute(denominator_yng))
    DENOMINATOR_OLD_LAB <- if (is.character(denominator_old))
            denominator_old else deparse(substitute(denominator_old))

    NUMERATOR_YNG_XYZ <- read_fn(numerator_yng)
    NUMERATOR_OLD_XYZ <- read_fn(numerator_old)
    DENOMINATOR_YNG_XYZ <- read_fn(denominator_yng)
    DENOMINATOR_OLD_XYZ <- read_fn(denominator_old)

    .compare_data(NUMERATOR_YNG_XYZ, NUMERATOR_OLD_XYZ)
    .compare_data(NUMERATOR_OLD_XYZ, NUMERATOR_YNG_XYZ)
    .compare_data(NUMERATOR_YNG_XYZ, DENOMINATOR_YNG_XYZ)
    .compare_data(NUMERATOR_OLD_XYZ, DENOMINATOR_OLD_XYZ)

    if (!is.null(landmarks)) {
        NUMERATOR_YNG_XYZ <- NUMERATOR_YNG_XYZ[landmarks,,]
        NUMERATOR_OLD_XYZ <- NUMERATOR_OLD_XYZ[landmarks,,]
        DENOMINATOR_YNG_XYZ <- DENOMINATOR_YNG_XYZ[landmarks,,]
        DENOMINATOR_OLD_XYZ <- DENOMINATOR_OLD_XYZ[landmarks,,]
    }
    NUMERATOR_YNG_FIT <- edma_fit(NUMERATOR_YNG_XYZ, B=B)
    NUMERATOR_OLD_FIT <- edma_fit(NUMERATOR_OLD_XYZ, B=B)
    DENOMINATOR_YNG_FIT <- edma_fit(DENOMINATOR_YNG_XYZ, B=B)
    DENOMINATOR_OLD_FIT <- edma_fit(DENOMINATOR_OLD_XYZ, B=B)

    GDM <- edma_gdm(
        a1=NUMERATOR_YNG_FIT,#yng_den
        a2=NUMERATOR_OLD_FIT,#yng_num
        b1=DENOMINATOR_YNG_FIT,#old_den
        b2=DENOMINATOR_OLD_FIT,#old_num
        B=B, ref_denom=ref_denom, mix=mix)

    ## GM numerator
    Dis_NUM <- .formdiff(
        Meanform(NUMERATOR_OLD_FIT),
        Meanform(NUMERATOR_YNG_FIT))
    ## GM denominator
    Dis_DEN <- .formdiff(
        Meanform(DENOMINATOR_OLD_FIT),
        Meanform(DENOMINATOR_YNG_FIT))
    ## GDM = (old_num/yng_num) / (old_den/yng_den)
    Dis <- Dis_NUM / Dis_DEN
    TV <- global_test(GDM)
    CI1 <- get_gdm(GDM, level=level, sort=FALSE)
    CI2 <- get_gdm(GDM, level=level, sort=TRUE)
    colnames(CI1) <- c("", "", "Estimate", "Low", "High")
    colnames(CI2) <- c("", "", "Estimate", "Low", "High")

    h <- hist(c(TV$statistic, TV$Tvals), plot=FALSE)
    MX <- 50
    NN <- round(MX * h$counts / max(h$counts))
    ST <- sapply(NN, function(z) paste0(rep("*", z), collapse=""))
    ii <- which.min(abs(TV$statistic-h$mids))
    #ST[ii] <- paste0(substr(ST[ii], 1, nchar(ST[ii])-1), "T")
    STEM <- paste0("  ", format(h$mids), "   |", ST)
    STEM[ii] <- gsub("   |", "  G|", STEM[ii], fixed=TRUE)
    ver <- read.dcf(file=system.file("DESCRIPTION", package="EDMAinR"),
        fields=c("Version", "Date"))

    sink(output, append=TRUE)
    on.exit(sink(), add=TRUE)
    op <- options(width = 10000, scipen=999)
    on.exit(options(op), add=TRUE)

    cat("# EDMAinR GDM output (v ", ver[1], " - ", ver[2], ")",
        "\n\nDate created: ", as.character(Sys.time()), "\n\n",
        file=output, sep="")
    Cat(
        "## Settings",
        "",
        paste0("OLDER NUMERATOR: ", NUMERATOR_OLD_XYZ$name),
        paste0("Number of specimens: ", dim(NUMERATOR_OLD_FIT)[3]),
        paste0(if (NUMERATOR_OLD_FILE) "Input file: " else "Input object: ",
            NUMERATOR_OLD_LAB),
        "",
        paste0("YOUNGER NUMERATOR: ", NUMERATOR_YNG_XYZ$name),
        paste0("Number of specimens: ", dim(NUMERATOR_YNG_FIT)[3]),
        paste0(if (NUMERATOR_YNG_FILE) "Input file: " else "Input object: ",
            NUMERATOR_YNG_LAB),
        "",
        paste0("OLDER DENOMINATOR: ", DENOMINATOR_OLD_FIT$name),
        paste0("Number of specimens: ", dim(DENOMINATOR_OLD_FIT)[3]),
        paste0(if (DENOMINATOR_OLD_FILE) "Input file: " else "Input object: ",
            DENOMINATOR_OLD_LAB),
        "",
        paste0("YOUNGER DENOMINATOR: ", DENOMINATOR_YNG_FIT$name),
        paste0("Number of specimens: ", dim(DENOMINATOR_YNG_FIT)[3]),
        paste0(if (DENOMINATOR_YNG_FILE) "Input file: " else "Input object: ",
            DENOMINATOR_YNG_LAB),
        "",
        paste0("Number of landmarks: ", dim(NUMERATOR_YNG_FIT)[1]),
        paste0("Number of dimensions:  ", dim(NUMERATOR_YNG_FIT)[2]),
        "",
        "## Form matrices",
        "",
        "XY(Z) coordinates of the OLDER NUMERATOR mean:",
        ""
    )
    print(Meanform(NUMERATOR_OLD_FIT), digits=digits)
    Cat(
        "",
        "XY(Z) coordinates of the YOUNGER NUMERATOR mean:",
        ""
    )
    print(Meanform(NUMERATOR_YNG_FIT), digits=digits)
    Cat(
        "",
        "XY(Z) coordinates of the OLDER DENOMINATOR mean:",
        ""
    )
    print(Meanform(DENOMINATOR_OLD_FIT), digits=digits)
    Cat(
        "",
        "XY(Z) coordinates of the YOUNGER DENOMINATOR mean:",
        ""
    )
    print(Meanform(DENOMINATOR_YNG_FIT), digits=digits)
    Cat(
        "",
        "Further-improved estimate of Sigma K* for the OLDER NUMERATOR:",
        ""
    )
    Mat(SigmaKstar(NUMERATOR_OLD_FIT))
    Cat(
        "",
        "Further-improved estimate of Sigma K* for the YOUNGER NUMERATOR:",
        ""
    )
    Mat(SigmaKstar(NUMERATOR_YNG_FIT))
    Cat(
        "",
        "Further-improved estimate of Sigma K* for the OLDER DENOMINATOR:",
        ""
    )
    Mat(SigmaKstar(DENOMINATOR_OLD_FIT))
    Cat(
        "",
        "Further-improved estimate of Sigma K* for the YOUNGER DENOMINATOR:",
        ""
    )
    Mat(SigmaKstar(DENOMINATOR_YNG_FIT))
    Cat(
        "",
        "Form matrix for the OLDER NUMERATOR:",
        ""
    )
    print(as.dist(NUMERATOR_OLD_FIT, diag=TRUE), digits=digits)
    Cat(
        "",
        "Form matrix for the YOUNGER NUMERATOR:",
        ""
    )
    print(as.dist(NUMERATOR_YNG_FIT, diag=TRUE), digits=digits)
    Cat(
        "",
        "Form matrix for the OLDER DENOMINATOR:",
        ""
    )
    print(as.dist(DENOMINATOR_OLD_FIT, diag=TRUE), digits=digits)
    Cat(
        "",
        "Form matrix for the YOUNGER DENOMINATOR:",
        ""
    )
    print(as.dist(DENOMINATOR_YNG_FIT, diag=TRUE), digits=digits)
    Cat(
        "",
        "## Growth matrix for the NUMERATOR",
        ""
    )
    print(as.dist(Dis_NUM, diag=TRUE), digits=digits)
    Cat(
        "",
        "## Growth matrix for the DENOMINATOR",
        ""
    )
    print(as.dist(Dis_DEN, diag=TRUE), digits=digits)
    Cat(
        "",
        "## Growth difference matrix",
        ""
    )
    print(as.dist(Dis, diag=TRUE), digits=digits)
    Cat(
        "",
        paste0("G statistic (max/min): ", round(TV$statistic, digits)),
        "",
        paste0("The bootstrap reference sample is the ",
            if (ref_denom) "DENOMINATOR" else "NUMERATOR"),
        paste0("Number of bootstrap samples: ", B),
        paste0("Bootstrap type: ",
            if (mix) "mixed" else "reference is fixed"),
        paste0("Probability: ",
            format(max(TV$p.value, 10^-(digits+1)))),
        "",
        STEM,
        "",
        "Unsorted form-difference matrix with marginal confidence intervals:",
        ""
    )
    print(CI1, digits=4, row.names=FALSE, right=FALSE)
    Cat(
        "",
        paste0("Confidence intervals are based on ", B, " resamples."),
        paste0("Nonparametric resampling was used."),
        paste0("Alpha level: ", 1-level),
        "",
        "Sorted form-difference matrix:",
        ""
    )
    print(CI2, digits=4, row.names=FALSE, right=FALSE)
    Cat(
        "",
        "Significant form differences < 1:",
        ""
    )
    ss <- CI2$Low < 1 & CI2$High < 1
    if (sum(ss) < 1) {
        cat("NONE\n")
    } else {
        print(CI2[ss,], digits=4,
            row.names=FALSE, right=FALSE)
    }
    Cat(
        "",
        "Non-significant form differences:",
        ""
    )
    ss <- CI2$Low <= 1 & CI2$High >= 1
    if (sum(ss) < 1) {
        cat("NONE\n")
    } else {
        print(CI2[ss,], digits=4,
            row.names=FALSE, right=FALSE)
    }
    Cat(
        "",
        "Significant form differences > 1:",
        ""
    )
    ss <- CI2$Low > 1 & CI2$High > 1
    if (sum(ss) < 1) {
        cat("NONE\n")
    } else {
        print(CI2[ss,], digits=4,
            row.names=FALSE, right=FALSE)
    }
    Cat("")

    out <- list(
        NUMERATOR_YNG_LAB=NUMERATOR_YNG_LAB,
        NUMERATOR_OLD_LAB=NUMERATOR_OLD_LAB,
        DENOMINATOR_YNG_LAB=DENOMINATOR_YNG_LAB,
        DENOMINATOR_OLD_LAB=DENOMINATOR_OLD_LAB,

        NUMERATOR_YNG_FIT=NUMERATOR_YNG_FIT,
        NUMERATOR_OLD_FIT=NUMERATOR_OLD_FIT,
        DENOMINATOR_YNG_FIT=DENOMINATOR_YNG_FIT,
        DENOMINATOR_OLD_FIT=DENOMINATOR_OLD_FIT,
        GDM=GDM,
        GM_NUM=Dis_NUM,
        GM_DEN=Dis_DEN,
        TV=TV,
        CI1=CI1,
        CI2=CI2,
        STEM=STEM
    )
    invisible(out)
}
