\name{edma_sdm}
\alias{edma_sdm}
\alias{print.edma_sdm}
\alias{edma_sdm_clr}
\alias{print.edma_sdm_clr}
\alias{Z_test}
\alias{Z_test.edma_sdm}
\alias{plot_Ztest}
\alias{plot_Ztest.edma_sdm}
\alias{plot_ci.edma_sdm}
\alias{get_sdm}
\alias{get_sdm.edma_sdm}
\alias{confint.edma_sdm}
\alias{landmarks.edma_sdm}
\alias{dimensions.edma_sdm}
\alias{get_influence.edma_sdm}
\title{
Shape difference
}
\description{
Shape difference matrix based inference following
Lele and Cole (1996).
}
\usage{
edma_sdm_clr(a, b, B=0)
\method{print}{edma_sdm_clr}(x, level = 0.95, ...)


edma_sdm(a, b, log=TRUE)
get_sdm(object, ...)
\method{get_sdm}{edma_sdm}(object, sort=FALSE,
    level = 0.95, ...)

\method{print}{edma_sdm}(x, level = 0.95, ...)
Z_test(object, ...)
\method{Z_test}{edma_sdm}(object, level = 0.95, ...)
\method{landmarks}{edma_sdm}(x, ...)
\method{dimensions}{edma_sdm}(x, ...)

\method{confint}{edma_sdm}(object, parm, level=0.95, ...)
\method{get_influence}{edma_sdm}(object, statistic=c("Z", "C"),
    level=0.95, ...)

plot_Ztest(x, ...)
\method{plot_Ztest}{edma_sdm}(x, statistic=c("Z", "C"),
    level = 0.95, ...)
}
\arguments{
  \item{a, b}{
EDMA fit object to compare shapes.
}
  \item{B}{
number of parametric bootstrap runs.
}
  \item{x, object}{
a SDM object.
}
  \item{log}{
logical, if form matrix is to be log transformed before
calculating the differences.
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
  \item{statistic}{
character, the Z or C statistic to be plotted.
}
  \item{\dots}{
other arguments passed to other functions.
}
}
\details{
Shape difference matrix (SDM) is defined as
the difference between the scaled form matrices S(A) and S(B).
S(A) = C * FM(A), S(B) = FM(B),
where C is a scaling factor and is calculated using  TLS.
Shape difference matrix is S(A) - S(B) when \code{log=FALSE}
and log(S(A)) - log(S(B)) when \code{log=TRUE}.

Inference and visualization is similar to how it is done for FDMs.
}
\value{
\code{edma_sdm} compares 2 EDMA fit objects and calculates SDM.
}
\references{
Lele, S. R., and Cole, T. M. III., 1996.
A new test for shape differences when
variance-covariance matrices are unequal.
Journal of Human Evolution 31:193--212.
<doi:10.1006/jhev.1996.0057>
}
\author{
Peter Solymos,
Subhash R. Lele,
Theodore M. Cole
}
\seealso{
Nonparametric fit: \code{\link{edma_fit}}

Parametric fit: \code{\link{SigmaK_fit}}

Form difference: \code{\link{edma_fdm}}

Growth difference: \code{\link{edma_gdm}}
}
\examples{
file_a <- system.file("extdata/growth/CZP0_wt_global.xyz",
    package="EDMAinR")
file_b <- system.file("extdata/growth/CZP0_mut_global.xyz",
    package="EDMAinR")
l <- c("amsph", "bas", "loci", "lpto", "lsqu",
        "lsyn", "roci", "rpto", "rsqu", "rsyn")

a <- read_xyz(file_a)[l,,]
b <- read_xyz(file_b)[l,,]
a
b

fit_a <- edma_fit(a, B=10)
fit_b <- edma_fit(b, B=10)

sdm <- edma_sdm(a=fit_a, b=fit_b)
sdm
Z_test(sdm)
head(confint(sdm))
head(get_sdm(sdm))
head(get_sdm(sdm, sort=TRUE, decreasing=TRUE))
head(get_sdm(sdm, sort=TRUE, decreasing=FALSE))

get_influence(sdm)

plot_Ztest(sdm, "Z")
plot_Ztest(sdm, "C")
plot_ci(sdm)

plot(get_influence(sdm))
}
\keyword{ models }
