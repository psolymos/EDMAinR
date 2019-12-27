\name{SigmaK_fit}
\alias{SigmaK_fit}
\alias{SigmaK}
\alias{sensitivity}
\alias{print.edma_fit_p}
\alias{SigmaK.edma_fit_p}
\alias{sensitivity.edma_fit_p}
\alias{read_pattern}
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
SigmaK_fit(object, pattern, check_pattern = TRUE, ...)

## methods
SigmaK(object, ...)
sensitivity(object, ...)
\method{print}{edma_fit_p}(x, truncate = 40, ...)
\method{SigmaK}{edma_fit_p}(object, ...)
\method{sensitivity}{edma_fit_p}(object, m=10, ...)
}
\arguments{
  \item{object, x}{
EDMA fit object.
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
  \item{check_pattern}{
logical, if pattern should be checked for identifiability.
}
  \item{m}{
positive integer, number of times SigmaK estimation
is to be replicated to assess sensitivity.
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

\code{sensitivity} is used to re-estimate SigmaK matrix
multiple times with random starting values for numerical
optimization. When estimates vary while the loss function
value is the same, it is usually a sign of non-identifiability.
The return object has columns starting with
\code{par_*} that are the parameters according to pattern matrix,
\code{value} is the loss function value evaluated at
the estimates.
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
m

## simulate data, so we know what the true values were
M <- structure(c(-2.5, 7.5, -2.5, -2.5, -7.5, 2.5, 2.5, 4.5),
    .Dim = c(4L, 2L))
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

## sensitivity
s <- sensitivity(o)
summary(s)
}
\keyword{ IO }
\keyword{ models }
