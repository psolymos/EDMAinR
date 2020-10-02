predict_landmarks <- function(A) {
    ## get dimensions from object A
    A0 <- A
    DIM <- dim(A)
    K <- DIM[1L]
    D <- DIM[2L]
    n <- DIM[3L]

    ## pairwise distances for specimens
    DM <- sapply(seq_len(n), function(i) as.vector(dist(A[,,i])))

    ## rows/cols
    tmp <- matrix(0, K, K)
    RC <- data.frame(
        row=row(tmp)[lower.tri(tmp)],
        col=col(tmp)[lower.tri(tmp)])
    RC$id <- seq_len(nrow(RC))

    ## where NAs are: row says which pair, col sais specimen
    NAS <- which(is.na(DM))
    NAS <- data.frame(id=NAS, row=row(DM)[NAS], col=col(DM)[NAS])

    ## we store this to keep predictions after SVM
    DM1 <- DM
    ## pick rows from NAS
    for (i in seq_len(nrow(NAS))) {
        ## complete obs based on dropping row/cols with NAs
        cdrop <- colSums(is.na(DM)) > 0
        #rdrop <- rowSums(is.na(DM[,-NAS$col[i],drop=FALSE])) > 0
        rdrop <- rowSums(is.na(DM[,!cdrop,drop=FALSE])) > 0
        ii <- RC$id == NAS$row[i]
        Y <- DM[ii, !cdrop, drop=TRUE]
        X <- DM[!ii & !rdrop, !cdrop, drop=FALSE]
        Xnew <- DM[!ii & !rdrop, NAS$col[i], drop=FALSE]
        ok <- !is.na(Xnew)
        df <- data.frame(Y=Y, X=t(X[ok,,drop=FALSE]))
        dfnew <- data.frame(X=t(Xnew[ok]))
        ## SVM fit
        m <- svm(Y ~ ., data=df)
        ## SVM predict
        pr <- unname(predict(m, newdata=dfnew))
        ## save output
        DM1[NAS$row[i], NAS$col[i]] <- pr
    }
    ## we store this to keep predictions after MDS
    DM2 <- DM1
    d0 <- dist(matrix(0, K, D))
    ## pick rows from NAS
    for (j in unique(NAS$col)) {
        z <- is.na(A[,,j])
        d <- d0
        d[] <- DM1[,j]
        xyz1 <- cmdscale(d, k=D)
        dd <- dist(xyz1)
        k <- NAS$row[NAS$col == j]
        DM2[k, j] <- dd[k]
    }
    ## use optim to find coordinates
    BX <- apply(A, 2, range, na.rm=TRUE)
    init <- colMeans(BX)
    ## parm: coordinates for missing landmark
    ## coord: corrdinates for non missing landmarks
    ## d: predicted distance between missing and
    ##    non missing landmarks (matching coordd)
    dfun <- function(parm, coord, pred) {
      sum((sqrt(colSums((t(coord) - parm)^2)) - pred)^2)
    }
    for (j in unique(NAS$col)) {
        Aj <- A[,,j]
        miss <- rowSums(is.na(Aj)) > 0
        coord <- Aj[!miss,,drop=FALSE]
        for (i in which(miss)) {
            a <- RC[RC$row == i | RC$col == i,]
            a$match <- ifelse(a$row == i, a$col, a$row)
            a$pred <- DM2[a$id,j]
            a <- a[a$match %in% which(!miss),]
            coord <- Aj[a$match,,drop=FALSE]
            o <- optim(init, dfun, coord=coord, pred=a$pred)
            A[i,,j] <- o$par
        }
    }
    list(A=A, dm=DM2, missing_A=which(is.na(A0)), missing_dm=which(is.na(DM)))
}
