\name{edma_fdm}
\alias{edma_fdm}
\alias{confint.edma_dm}
\alias{get_influence}
\alias{get_fdm}
\alias{get_influence.edma_dm}
\alias{print.edma_fdm}
\alias{T_test}
\alias{T_test.edma_fdm}
\alias{get_fdm.edma_fdm}
\alias{landmarks.edma_fdm}
\alias{dimensions.edma_fdm}
\alias{plot.edma_dm}
\alias{plot_2d.edma_dm}
\alias{plot_3d.edma_dm}
\alias{plot_ci}
\alias{plot_Ttest}
\alias{plot_Ttest.edma_dm}
\alias{plot_ord.edma_fdm}
\alias{plot_clust.edma_fdm}
\alias{plot_ci.edma_fdm}
\alias{plot.edma_influence}
\title{
Form difference
}
\description{
Form difference matrix based inference based on
Lele and Richtsmeier (1992, 1995).
}
\usage{
edma_fdm(numerator, denominator,
  B=0, ref_denom=TRUE, mix=FALSE)

get_influence(object, ...)
\method{get_influence}{edma_dm}(object, level=0.95, ...)
\method{plot}{edma_influence}(x, ...)

get_fdm(object, ...)
\method{get_fdm}{edma_fdm}(object, sort=FALSE, level=0.95, ...)
T_test(object, ...)
\method{T_test}{edma_fdm}(object, ...)
\method{confint}{edma_dm}(object, parm, level=0.95, ...)

\method{print}{edma_fdm}(x, ...)
\method{landmarks}{edma_fdm}(x, ...)
\method{dimensions}{edma_fdm}(x, ...)

plot_ci(x, ...)
plot_Ttest(x, ...)
\method{plot}{edma_dm}(x, ...)
\method{plot_2d}{edma_dm}(x, ...)
\method{plot_3d}{edma_dm}(x, ...)
\method{plot_Ttest}{edma_dm}(x, ...)
\method{plot_ci}{edma_fdm}(x, ...)
\method{plot_ord}{edma_fdm}(x, ...)
\method{plot_clust}{edma_fdm}(x, ...)
}
\arguments{
  \item{numerator, denominator}{
EDMA fit object to compare forms.
}
  \item{B}{
nonnegative integer, the number of bootstrap replicates.
}
  \item{ref_denom}{
logical, when \code{TRUE}, the \code{denominator} is used as
reference object (its form matrix is fixed when calculating
bootstrap comparing to the other object).
}
  \item{mix}{
logical, to use mixed bootstrap (numeriator and denominator
populations are mixed with replacement) or not
(only the non-reference population is resampled with replacement).
}
  \item{x, object}{
an EDMA FDM object of class 'edma_fdm'.
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
  \item{\dots}{
other arguments passed to methods.
}
}
\details{
Form difference (FDM) is calculated as the ratio of form matrices
(FM) from the numerator and denominator objects following
Lele and Richtsmeier (1992, 1995): FDM(A,B) = FM(B)/FM(A).
Form matrices are formed as pairwise Euclidean distances
between landmarks from EDMA fit objects using the estimated mean forms.

Bootstrap distribution is based on fixing the
reference FM and taking the ration between the reference FM
and the bootstrap FMs from the other object.

The T-statistic is based on the pairwise distanced in the FM,
taking the max/min of the distances. Confidence intervals for local testing
(via \code{confint}, \code{get_fdm}, and \code{plot_ci})
and T-test for global testing
(via \code{T_test}, and \code{plot_Ttest})
is based on the observed T-statistic and the bootstrap distribution.

Influential landmarks are identified by leaving one landmark out,
then comparing the T-statistic with the value based on all the
landmarks. The existing bootstrap distribution of
the mean form is used (i.e. no re-estimation of the mean form)
in \code{get_influence}.
}
\value{
\code{edma_fdm} compares two EDMA fit objects and calculates
form difference.

\code{confint} returns the confidence intervals for FDM,
the \code{get_fdm} extract the stacked FDM with confidence intervals,
the \code{plot_ci} visualizes the ordered form differences with
confidence intervals.

\code{get_influence} extracts landmark influence information,
the plot method visualizes this.

\code{T_test} presents the global T-test,
the bootstrap distribution and observed T-value is
visualized by \code{plot_Ttest}.

\code{plot} and \code{plot_2d} produces a 2D plot of the mean form
from the reference object ('prototype').
\code{plot_3d} use the {rgl} package to make a 3D plot using the same
mean form. Influential landmarks are colored red.
Lines represent distances between landmarks,
<1 differences are colored blue, >1 differences are colored red.

The \code{plot_ord} and \code{plot_clust}
produce plots based on dissimilarities among specimens
in the two objects.
}
\references{
Lele, S. R., and Richtsmeier, J. T., 1992.
On comparing biological shapes: detection of influential landmarks.
American Journal of Physical Anthropology 87:49--65.
<doi:10.1002/ajpa.1330870106>

Lele, S. R., and Richtsmeier, J. T., 1995.
Euclidean distance matrix analysis: confidence intervals for
form and growth differences.
American Journal of Physical Anthropology 98:73--86.
<doi:10.1002/ajpa.1330980107>
}
\author{
Peter Solymos,
Subhash R. Lele,
Theodore M. Cole,
Joan T. Richtsmeier
}
\seealso{
Nonparametric fit: \code{\link{edma_fit}}

Parametric fit: \code{\link{SigmaK_fit}}

Growth difference: \code{\link{edma_gdm}}
}
\examples{
file1 <- system.file("extdata/crouzon/Crouzon_P0_Global_MUT.xyz",
    package="EDMAinR")
x1 <- read_xyz(file1)

file2 <- system.file("extdata/crouzon/Crouzon_P0_Global_NON-MUT.xyz",
    package="EDMAinR")
x2 <- read_xyz(file2)

numerator <- edma_fit(x1, B=10)
denominator <- edma_fit(x2)

fdm <- edma_fdm(numerator, denominator, B=10)
fdm2 <- edma_fdm(numerator, denominator, B=10, ref_denom=FALSE)
fdm
fdm2

head(get_fdm(fdm))
head(get_fdm(fdm, sort=TRUE, decreasing=TRUE))
head(get_fdm(fdm, sort=TRUE, decreasing=FALSE))

T_test(fdm)
T_test(fdm2)

head(confint(fdm))

head(infl <- get_influence(fdm))
plot(infl)

plot_ord(fdm)
plot_clust(fdm)
plot_Ttest(fdm)
plot_ci(fdm)
plot_2d(fdm)
if (interactive())
    plot_3d(fdm)
}
\keyword{ models }
