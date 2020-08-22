\name{SigmaK_fit}
\alias{SigmaK_fit}
\alias{print.edma_fit_p}
\alias{SigmaK}
\alias{SigmaK.edma_fit_p}
\alias{SigmaKfull}
\alias{SigmaKfull.edma_fit_p}
\alias{read_pattern}
\alias{print_tb}
\alias{plot_tb}
\alias{make_Sigma}
\alias{sensitivity}
\alias{sensitivity.edma_fit_p}
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
SigmaK_fit(object, pattern,
    twostep=TRUE, check=TRUE, ...)

sensitivity(object, ...)

## methods
\method{print}{edma_fit_p}(x, truncate=40, ...)
SigmaK(object, ...)
\method{SigmaK}{edma_fit_p}(object, ...)
SigmaKfull(object, ...)
\method{SigmaKfull}{edma_fit_p}(object, ...)
\method{sensitivity}{edma_fit_p}(object, m=10, scale=10, ...)

print_tb(x, ...)
plot_tb(x, mar=c(1, 1, 1, 4), ...)

make_Sigma(params, pattern)
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
a pattern matrix, see Details and Examples.
}
  \item{twostep}{
logical, if the optimization is to be constrained by
the SigmaK estimate based on Hu (2007) approach,
see Details.
}
  \item{check}{
logical, if inputs should trigger errors when we know
the patterns are non identifiable, i.e. number of
non-zero elements > K(K-1)/2, or when the A matrix is non-invertible.
\code{check=FALSE} will still produce an error in such
cases when \code{twostep=TRUE}.
}

  \item{params}{
a named vector of parameters matching the non-\code{NA} values of the
pattern matrix, see Details and Examples.
}
  \item{m}{
positive integer, number of times SigmaK estimation
is to be replicated to assess sensitivity.
}
  \item{scale}{
positive numeric, scaling value for random numbers to be used in
sensitivity analysis. The original initial values are
multiplied by random Uniform(0.001, \code{scale}) values.
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

\code{SigmaKfull} calculates the 'full' parametric estimate according to
Hu (2007). \code{SigmaK} uses the \code{pattern} matrix to
constrain the parametric estimates in an additional optimization step
(thus the name 2-step). We minimize sum((SigmaKfull - SigmaK_hat)^2).

The 1-step approach minimizes
sum((SigmaKstar - (H \%*\% SigmaK_hat \%*\% H))^2)
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
p <- matrix(c(
    "a", NA, NA, NA,
    NA, "a", NA, NA,
    NA,  NA, "b", NA,
    NA,  NA, NA, "b"
), 4, 4, byrow=TRUE)
parm <- c(a=0.25, b=0.35)
print_tb(p)
plot_tb(p)
SigmaK <- make_Sigma(parm, p)
dimnames(p) <- dimnames(SigmaK)

## mean form
M <- matrix(c(-2.5, 7.5, -2.5, -2.5,
    -7.5, 2.5, 2.5, 4.5),
    4, 2)

## simulate data, so we know what the true values were
sim <- edma_simulate_data(n=500, M, SigmaK)

## nonparametric fit
fit <- edma_fit(sim)
## parametric fit
o <- SigmaK_fit(fit, p)
o

## estimates
cbind(true=parm, est=o$results$par[names(parm)])
SigmaK(o) # constrained parametric estimate
SigmaKfull(o) # full set of estimates

## sensitivity analysis
sensitivity(o)
}
\keyword{ IO }
\keyword{ models }
