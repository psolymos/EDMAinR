\name{edma_fit}
\alias{edma_fit}
\alias{Meanform}
\alias{SigmaKstar}
\alias{get_fm}
\alias{Meanform.edma_fit}
\alias{SigmaKstar.edma_fit}
\alias{get_fm.edma_fit}
\alias{confint.edma_fit}
\alias{as.dist.edma_fit}
\alias{print.edma_fit_np}
\alias{stack.dist}
\alias{plot.edma_fit}
\alias{plot_3d}
\alias{plot_2d.edma_fit}
\alias{plot_3d.edma_fit}
\alias{plot_ord.edma_fit}
\alias{plot_clust.edma_fit}
\alias{as.edma_data.edma_fit}
\title{
Nonparametric EDMA fit
}
\description{
Estimate mean form and SigmaKstar matrix based on
Lele (1991), Lele and Richtsmeier (1991) and Hu (2007).
}
\usage{
edma_fit(x, B=0, ncores=getOption("Ncpus", 1L))

## generics
Meanform(object, ...)
SigmaKstar(object, ...)
get_fm(object, ...)

## methods
\method{print}{edma_fit_np}(x, truncate=40, ...)
\method{Meanform}{edma_fit}(object, ...)
\method{SigmaKstar}{edma_fit}(object, ...)
\method{get_fm}{edma_fit}(object, sort=FALSE, level=0.95, ...)
\method{confint}{edma_fit}(object, parm, level=0.95, ...)
\method{as.edma_data}{edma_fit}(x, ...)

## plot methods
plot_3d(x, ...)
\method{plot}{edma_fit}(x, ...)
\method{plot_2d}{edma_fit}(x, ...)
\method{plot_3d}{edma_fit}(x, ...)
\method{plot_ord}{edma_fit}(x, ...)
\method{plot_clust}{edma_fit}(x, ...)

## distance manipulation
\method{as.dist}{edma_fit}(m, diag=FALSE, upper=FALSE)
\method{stack}{dist}(x, ...)
}
\arguments{
  \item{x, object, m}{
an EDMA data object of class 'edma_data'.
}
  \item{B}{
nonnegative integer, the number of bootstrap replicates.
}
  \item{ncores}{
positive integer, the number of cores to use when bootstrapping.
Use \code{options(Ncpus = 2)} to set it to 2 globally.
}
  \item{truncate}{
numeric, number of characters to print for the object title.
}
  \item{sort}{
logical, if stacked distances are to be sorted, see Examples.
}
  \item{level}{
numeric, between 0 and 1, alpha level for confidence interval.
}
  \item{parm}{
a specification of which parameters are to be given
confidence intervals, either a vector of numbers or a vector of names.
See \code{\link{confint}}.
}
  \item{diag, upper}{
logical, indicating whether the diagonal
and the upper triangle of the distance matrix should be printed.
See \code{\link{as.dist}}.
}
  \item{\dots}{
other arguments passed to methods. E.g.
for \code{plot_clust}, the \code{method}
describing the clustering agglomeration method to be used
by the \code{link{hclust}} function (default is \code{"ward.D2"}).
}
}
\details{
The function estimates mean form and
SigmaKstar matrix based on
Lele (1991), Lele and Richtsmeier (1991) and Hu (2007).
}
\value{
\code{edma_fit} returns and EDMA fit object of class 'edma_fit'.
\code{.edma_fit_np} is the workhorse function behind \code{edma_fit}.

\code{stack.dist} takes any distance matrix of class 'dist'
and turns that into a long form data frame with
columns \code{row} and \code{col} indicating the row and column
labels, \code{dist} giving the value in that cell.
Only returns the values from the lower triangle of the matrix.

\code{get_fm} is the intended user interface to extract
the form matrix (FM) from EDMA fit objects.
This has the stacked distances based on the mean form.
When the object has bootstrap replicates,
\code{get_fm} also returns confidence intervals
for the distances based on bootstrap and the \code{confint} method.

\code{Meanform} extracts the mean form (K x D) matrix,
\code{SigmaKstar} extracts the corresponding uncertainties
(K x K) based on the EDMA fit object.

\code{plot} and \code{plot_2d} produces a 2D plot of the mean form.
2D projection is used in case of 3D landmark data based on
metric multidimensional scaling.
\code{plot_3d} use the {rgl} package to make a 3D plot.
The sizes of the dots correspond to
square root of the SigmaKstar diagonal elements.

The \code{plot_ord} and \code{plot_clust}
produce plots based on dissimilarities among specimens,
see \code{\link{plot_ord.edma_data}} for details.
}
\references{
Lele, S. R., 1991.
Some comments on coordinate-free and scale-invariant
methods in morphometrics.
American Journal of Physical Anthropology 85:407--417.
<doi:10.1002/ajpa.1330850405>

Lele, S. R., and Richtsmeier, J. T., 1991.
Euclidean distance matrix analysis: A coordinate-free approach for
comparing biological shapes using landmark data.
American Journal of Physical Anthropology 86(3):415--27.
<doi:10.1002/ajpa.1330860307>

Hu, L., 2007. Euclidean Distance Matrix Analysis of Landmarks Data:
Estimation of Variance. Thesis, Master of Science in Statistics,
Department of Mathematical and Statistical Sciences,
University of Alberta, Edmonton, Alberta, Canada. Pp. 49.
}
\author{
Peter Solymos,
Subhash R. Lele,
Theodore M. Cole,
Liangyuan Hu,
Joan T. Richtsmeier
}
\seealso{
EDMA data: \code{\link{read_xyz}}

Parametric fit: \code{\link{SigmaK_fit}}

Form difference: \code{\link{edma_fdm}}

Growth difference: \code{\link{edma_gdm}}

Shape difference: \code{\link{edma_sdm}}
}
\examples{
file <- system.file(
    "extdata/crouzon/Crouzon_P0_Global_MUT.xyz",
    package="EDMAinR")
x <- read_xyz(file)
x <- x[,,1:10] # 10 specimens

## nonparametric fit
fit <- edma_fit(x, B=9)
fit
str(Meanform(fit))
str(SigmaKstar(fit))

## form matrix
str(as.dist(fit))
str(stack(as.dist(fit)))

head(get_fm(fit))
head(get_fm(fit, sort=TRUE, decreasing=TRUE))
head(get_fm(fit, sort=TRUE, decreasing=FALSE))

plot_ord(fit)
plot_clust(fit)
plot(fit)
plot_2d(fit)
if (interactive())
    plot_3d(fit)
}
\keyword{ models }
