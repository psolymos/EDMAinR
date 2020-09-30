\name{edma_gdm}
\alias{edma_gm}
\alias{edma_gdm}
\alias{print.edma_gm}
\alias{print.edma_gdm}
\alias{T_test.edma_gm}
\alias{T_test.edma_gdm}
\alias{plot_ci.edma_gm}
\alias{plot_ci.edma_gdm}
\alias{get_gm}
\alias{get_gm.edma_gm}
\alias{get_gdm}
\alias{get_gdm.edma_gdm}
\alias{landmarks.edma_gdm}
\alias{dimensions.edma_gdm}
\alias{plot_ord.edma_gdm}
\alias{plot_clust.edma_gdm}
\title{
Growth difference
}
\description{
Growth matrix and growth difference matrix based inference based on
Lele and Richtsmeier (1992, 1995).
}
\usage{
edma_gm(a1, a2, ...)
get_gm(object, ...)
\method{get_gm}{edma_gm}(object, sort=FALSE, level=0.95, ...)

edma_gdm(a1, a2, b1, b2, ...)
get_gdm(object, ...)
\method{get_gdm}{edma_gdm}(object, sort=FALSE, level=0.95, ...)

\method{print}{edma_gm}(x, ...)
\method{print}{edma_gdm}(x, ...)
\method{T_test}{edma_gm}(object, ...)
\method{T_test}{edma_gdm}(object, ...)
\method{landmarks}{edma_gdm}(x, ...)
\method{dimensions}{edma_gdm}(x, ...)

\method{plot_ord}{edma_gdm}(x, ...)
\method{plot_clust}{edma_gdm}(x, ...)
}
\arguments{
  \item{a1, a2, b1, b2}{
EDMA fit object to compare growths.
}
  \item{x, object}{
an EDMA GM or GDM objects.
}
  \item{sort}{
logical, if stacked distances are to be sorted, see Examples.
}
  \item{level}{
numeric, between 0 and 1, alpha level for confidence interval.
}
  \item{\dots}{
other arguments passed to \code{\link{edma_fdm}},
like \code{ref_denom}.
}
}
\details{
Growth matrix (GM) is calculated as the ratio of form matrices
(FM) from the numerator and denominator objects following
Lele and Richtsmeier (1992, 1995): GM(A1,A2) = FM(A2)/FM(A1).
Form matrices are formed as pairwise Euclidean distances
between landmarks from EDMA fit objects using the estimated mean forms.

Growth difference matrix (GDM) is calculated as
GDM(A1,A2,B1,B2) = GM(B1,B2) / GM(A1,A2).

Inference and visualization is similar to how it is done for FDMs.
}
\value{
\code{edma_gm} compares two EDMA fit objects and calculates GM.

\code{edma_gdm} compares 4 EDMA fit objects and calculates GDM.

The \code{plot_ord} and \code{plot_clust}
produce plots based on dissimilarities among specimens
in the 2 or 4 objects (for GM and GDM, respectively).
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

Form difference: \code{\link{edma_fdm}}

Shape difference: \code{\link{edma_sdm}}
}
\examples{
file_a1 <- system.file("extdata/growth/CZEM_wt_global.xyz",
    package="EDMAinR")
file_a2 <- system.file("extdata/growth/CZP0_wt_global.xyz",
    package="EDMAinR")

l <- c("amsph", "bas", "loci", "lpto", "lsqu",
        "lsyn", "roci", "rpto", "rsqu", "rsyn")

a1 <- read_xyz(file_a1)[l,,]
a2 <- read_xyz(file_a2)[l,,]
a1
a2

fit_a1 <- edma_fit(a1, B=10)
fit_a2 <- edma_fit(a2, B=10)

## --- growth matrix ---

gm <- edma_gm(a1=fit_a1, a2=fit_a2, B=10)
gm
T_test(gm)
head(confint(gm))
head(get_gm(gm))
head(get_gm(gm, sort=TRUE, decreasing=TRUE))
head(get_gm(gm, sort=TRUE, decreasing=FALSE))

plot_ord(gm)
plot_clust(gm)
plot_Ttest(gm)
plot_ci(gm)
plot_2d(gm)
if (interactive())
    plot_3d(gm)

## --- growth difference matrix ---

file_b1 <- system.file("extdata/growth/CZEM_mut_global.xyz",
    package="EDMAinR")
file_b2 <- system.file("extdata/growth/CZP0_mut_global.xyz",
    package="EDMAinR")

b1 <- read_xyz(file_b1)[l,,]
b2 <- read_xyz(file_b2)[l,,]
b1
b2

fit_b1 <- edma_fit(b1, B=10)
fit_b2 <- edma_fit(b2, B=10)

gdm <- edma_gdm(a1=fit_a1, a2=fit_a2, b1=fit_b1, b2=fit_b2, B=10)
gdm
T_test(gdm)
head(confint(gdm))
head(get_gdm(gdm))
head(get_gdm(gdm, sort=TRUE, decreasing=TRUE))
head(get_gdm(gdm, sort=TRUE, decreasing=FALSE))

plot_ord(gdm)
plot_clust(gdm)
plot_Ttest(gdm)
plot_ci(gdm)
plot_2d(gdm) # need real data
if (interactive())
    plot_3d(gdm)
}
\keyword{ models }
