\name{edma_class}
\alias{edma_class}
\alias{loo}
\title{
Classification
}
\description{
Classification
}
\usage{
edma_class(x, fit1, fit2, boot=FALSE,
    method=c("cip", "chisq", "bnorm", "norm"))

loo(fit1, fit2, B=0, level=0.95,
    method=c("cip", "chisq", "bnorm", "norm"))
}
\arguments{
  \item{x}{
matrix, new specimen landmark coordinates (K landmarks x D dimensions).
Can be an edma_data object with a single specimen.
}
  \item{fit1, fit2}{
edma_fit objects to be used to classify the new specimen \code{x}.
}
  \item{boot}{
logical, if bootstrap info is to be used for classification.
}
  \item{B}{
nonnegative integer, the number of bootstrap replicates.
}
  \item{level}{
numeric, between 0 and 1, alpha level for confidence interval.
}
  \item{method}{
character, classification method to use: centered inner product,
non-central Chi squared, bivariate normal, normal.
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
file2 <- system.file("extdata/growth/CZP0_mut_global.xyz",
    package="EDMAinR")
l <- c("amsph", "bas", "loci", "lpto", "lsqu",
        "lsyn", "roci", "rpto", "rsqu", "rsyn")

x1 <- read_xyz(file1)[l,,]
x2 <- read_xyz(file2)[l,,]

New <- x2[,,1]
x2 <- x2[,,-1]
x1
x2

fit1 <- edma_fit(x1, B=10)
fit2 <- edma_fit(x2, B=10)

## CLR > 0 => x belongs to fit2
edma_class(New, fit1, fit2, boot=TRUE)$complikr

L <- loo(fit1, fit2)
L
cm <- with(L, table(group, class))
sum(diag(cm)) / sum(cm)

}
\keyword{ classif }
