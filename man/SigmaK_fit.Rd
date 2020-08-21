\name{SigmaK_fit}
\alias{SigmaK_fit}
\alias{SigmaK}
\alias{print.edma_fit_p}
\alias{SigmaK.edma_fit_p}
\alias{read_pattern}
\alias{print_tb}
\alias{plot_tb}
\title{
Parametric EDMA fit
}
\description{
Estimate parameters of a structured SigmaK,
partly based on Hu (2007).
}
\usage{
## read pattern file
read_pattern(file, ...)

## estimation
SigmaK_fit(object, pattern, ...)

## methods
SigmaK(object, ...)
\method{print}{edma_fit_p}(x, truncate=40, ...)
\method{SigmaK}{edma_fit_p}(object, ...)

print_tb(x, ...)
plot_tb(x, mar=c(1, 1, 1, 4), ...)
}
\arguments{
  \item{object, x}{
EDMA fit object. For \code{print_tb} and \code{plot_tb} it
is a symmetric square matrix.
}
  \item{file}{
the name of the file which the data are to be read from
using \code{\link{read.csv}} or \code{\link[readxl]{read_excel}}.
}
  \item{truncate}{
numeric, number of characters to print for the object title.
}
  \item{pattern}{
a pattern matrix, see Details.
}
  \item{mar}{
plot margins, see \code{\link{par}}.
}
  \item{\dots}{
other arguments passed to methods.
}
}
\details{
The pattern matrix describes the unknowns in the SigmaK matrix.
Diagonal must be filled. Upper and lower diagonal values must
be symmetric or \code{NA}. See Examples.
}
\value{
\code{read_pattern} reads in a csv or xls/xlsx file with the
pattern matrix. See Examples.

\code{SigmaK_fit} estimates the SigmaK matrix and attaches
that to the EDMA fit object. The estimated SigmaK matrix
can be accessed by the \code{SigmaK} method.
Estimation is done by \code{\link{optim}}.
}
\references{
Hu, L., 2007. Euclidean Distance Matrix Analysis of Landmarks Data:
Estimation of Variance. Thesis, Master of Science in Statistics,
Department of Mathematical and Statistical Sciences,
University of Alberta, Edmonton, Alberta, Canada. Pp. 49.
}
\author{
Peter Solymos,
Subhash R. Lele
}
\seealso{
EDMA data: \code{\link{read_xyz}}

Nonparametric fit: \code{\link{edma_fit}}

Form difference: \code{\link{edma_fdm}}

Growth difference: \code{\link{edma_gdm}}
}
\examples{
## read pattern matrix in cs or xls/xlsx format
read_pattern(
    system.file("extdata/example.csv", package="EDMAinR"))
read_pattern(
    system.file("extdata/example.xlsx", package="EDMAinR"))

## here is a pattern matrix
m <- matrix(c(
    "a", NA, NA, NA,
    NA, "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, NA, "b"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.35)
print_tb(m)
plot_tb(m)

## simulate data, so we know what the true values were
M <- structure(c(-2.5, 7.5, -2.5, -2.5, -7.5, 2.5, 2.5, 4.5),
    .Dim=c(4L, 2L))
SigmaK <- EDMAinR:::.vec2mat(parm, EDMAinR:::.mat2fac(m))
sim <- edma_simulate_data(n=500, M, SigmaK)
dimnames(M) <- dimnames(sim$data[[1L]])
rownames(SigmaK) <- rownames(m) <- rownames(sim$data[[1L]])
colnames(SigmaK) <- colnames(m) <- rownames(sim$data[[1L]])

## nonparametric fit
fit <- edma_fit(sim)
## parametric fit
o <- SigmaK_fit(fit, m)
o

## estimates
cbind(true=parm, est=o$results$par)
SigmaK(o)
}
\keyword{ IO }
\keyword{ models }
