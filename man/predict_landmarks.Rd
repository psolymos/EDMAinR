\name{predict_landmarks}
\alias{predict_landmarks}
\title{
Missing landmarks
}
\description{
Missing landmarks
}
\usage{
predict_landmarks(A)
}
\arguments{
  \item{A}{
K x D x n array of landmark coordinates.
}
}
\details{
To be added...
}
\author{
Peter Solymos,
Subhash R. Lele
}
\seealso{
Nonparametric fit: \code{\link{edma_fit}}

Parametric fit: \code{\link{SigmaK_fit}}

Form difference: \code{\link{edma_fdm}}

Growth difference: \code{\link{edma_gdm}}
}
\examples{
file1 <- system.file("extdata/growth/CZP0_wt_global.xyz",
    package="EDMAinR")
l <- c("amsph", "bas", "loci", "lpto", "lsqu",
        "lsyn", "roci", "rpto", "rsqu", "rsyn")
x1 <- read_xyz(file1)[l,,]

A <- as.array(x1)
A[1,,1] <- NA
A[1:2,,1] <- NA

pr <- predict_landmarks(A)

DM0 <- sapply(1:dim(x1)[3], function(i)
    as.vector(dist(as.array(x1)[,,i])))

Dist <- cbind(original=DM0[pr$missing_dm],
    estimate=pr$dm[pr$missing_dm])
Coord <- cbind(original=as.array(x1)[is.na(A)],
    estimate=pr$A[is.na(A)])

op <- par(mfrow=c(1,2))
plot(Dist, main="dist")
abline(0,1)

plot(Coord, main="coord")
abline(0,1)
par(op)
}
\keyword{ classif }
