\name{edma_data}
\alias{read_xyz}
\alias{print.edma_data}
\alias{stack.edma_data}
\alias{subset.edma_data}
\alias{[.edma_data}
\alias{landmarks}
\alias{specimens}
\alias{dimensions}
\alias{as.edma_data}
\alias{as.edma_data.array}
\alias{landmarks.edma_data}
\alias{specimens.edma_data}
\alias{dimensions.edma_data}
\alias{dim.edma_data}
\alias{dimnames.edma_data}
\alias{as.array.edma_data}
\alias{as.data.frame.edma_data}
\alias{as.matrix.edma_data}
\alias{as.dist.edma_data}
\alias{edma_simulate_data}
\alias{plot.edma_data}
\alias{plot_2d}
\alias{plot_ord}
\alias{plot_clust}
\alias{plot_2d.edma_data}
\alias{plot_ord.edma_data}
\alias{plot_clust.edma_data}
\alias{combine_data}
\alias{combine_data4}
\title{
Functions for EDMA data objects
}
\description{
Functions for reading, simulating, and manipulating EDMA data.
}
\usage{
## read xyz files
read_xyz(file, ...)

## data generation
edma_simulate_data(n, M, SigmaK, H=NULL)

## print
\method{print}{edma_data}(x, truncate=40, ...)

## accessors
landmarks(x, ...)
dimensions(x, ...)
specimens(x, ...)
\method{landmarks}{edma_data}(x, ...)
\method{dimensions}{edma_data}(x, ...)
\method{specimens}{edma_data}(x, ...)
\method{dimnames}{edma_data}(x)
\method{dim}{edma_data}(x)

## subsetting
\method{subset}{edma_data}(x, subset, ...)
\method{[}{edma_data}(x, i, j, k)

## coercion
\method{stack}{edma_data}(x, ...)
\method{as.matrix}{edma_data}(x, ...)
\method{as.data.frame}{edma_data}(x, ...)
\method{as.array}{edma_data}(x, ...)
as.edma_data(x, ...)
\method{as.edma_data}{array}(x, ...)

combine_data(a, b,
    ga="G1", gb="G2")
combine_data4(a1, a2, b1, b2,
    ga1="A1", ga2="A2", gb1="B1", gb2="B2")

## plot methods
plot_2d(x, ...)
plot_ord(x, ...)
plot_clust(x, ...)
\method{plot}{edma_data}(x, which=NULL,
    ask=dev.interactive(), ...)
\method{plot_2d}{edma_data}(x, which=NULL, ...)
\method{plot_ord}{edma_data}(x, ...)
\method{plot_clust}{edma_data}(x, ...)

## dissimilarities
\method{as.dist}{edma_data}(m, diag=FALSE, upper=FALSE)
}
\arguments{
  \item{file}{
the name of the file which the data are to be read from,
see \code{\link{read.table}} for more details.
}
  \item{x, m}{
an EDMA data object of class 'edma_data'.
}
\item{which}{
if a subset of the specimens is required.
}
\item{ask}{
logical, if \code{TRUE}, the user is asked before each plot.
}
  \item{subset, i, j, k}{
subset is for subsetting specimens (e.g. for bootstrap).
[i, j, k] indices refer to [landmarks, dimensions, specimens].
}
  \item{n, M, SigmaK, H}{
number of specimens (n), mean form matrix (M, K x D),
variance-covariance matrix (K x K symmetric),
and H matrix (can be \code{NULL}) for simulations.
}
  \item{truncate}{
numeric, number of characters to print for the object title.
}
  \item{diag, upper}{
logical, indicating whether the diagonal
and the upper triangle of the distance matrix should be printed.
See \code{\link{as.dist}}.
}
  \item{a, b, a1, a2, b1, b2}{
EDMA data objects to be combined together. Landmarks must be
homologous (determined by dimension names).
}
  \item{ga, gb, ga1, ga2, gb1, gb2}{
character, group names that are prepended to the specimen names
to differentiate the groups.
}
  \item{\dots}{
other arguments passed to methods.
For \code{read_xyz}, arguments passed to \code{\link{read.table}}.
}
}
\details{
The xyz landmark data has the following structure, see Examples:

- Header: this is the description of the data.

- XYZ: indicates dimensions, XYZ means 3D landmark data.

- \code{42L 3 9}: dimensions, e.g. 42 landmarks (K), 3 dimensions (D), 9 specimens (n).

- Landmark names, separated by space.

- The stacked data of landmark coordinates, e.g. 3 columns,
space separated numeric values with K*n rows,
the K landmarks per individuals stacked n times.

- Blank line.

- Date on of scans for each specimen (n rows),
this part is also used to get specimen IDs.

After reading in or simulating and EDMA data object,
the methods help extracting info, or manipulate these objects.
See Values and Examples.

The EDMA data object (class 'edma_data') is a list with two
elements: \code{$name} is the data set name (header information from
the .xyz file), \code{$data} is a list of n matrices (the list can be
named if speciemen information is present),
each matrix is of dimension K x D, dimension names for the
matrices describing landmark names and coordinate names.
}
\value{
\code{edma_simulate_data} returns an EDMA data object of
class 'edma_data'.

\code{landmarks}, \code{dimensions}, and \code{specimens}
are dimensions names, \code{dimnames} returns these as a list.
Landmark names and dimensions are used to check
if landmarks are homogeneous among objects.

The \code{dim} returns the number of landmarks (K), dimensions (D),
and specimens (n) in a data object.

The \code{print} method prints info about the data object.

The methods \code{stack} and \code{as.matrix} return a stacked
2D array (K*n x D) with the landmark coordinates,
\code{as.data.frame} turns the same 2D stacked array into a data frame,
\code{as.array} returns a 3D array (K x D x n).
\code{as.edma_data} turns a 3D array to an EDMA data object,
this is useful to handle 3D array objects returned by many
functions of the geomorph package (i.e. after reding
Morphologika, NTS, TPS files).

\code{combine_data} and \code{combine_data4} combines
2 or 4 EDMA data sets together, landmarks must be homologous.

\code{as.dist} calculates the dissimilarity matrix (n x n, object
of class 'dist', see \code{\link{dist}}) containing
pairwise dissimilarities among the specimens.
Dissimilarity is based on the T-statistic (max/min distance)
averaged (so that it is symmetric) and on the log scale
(so that self dissimilarity is 0).

\code{subset} and \code{[i,j,k]} returns an EDMA data object
with the desired dimensions or permutations. See Examples.

\code{plot} and \code{plot_2d} produces a series of plots
as a side effect, returning the data object invisibly.
The function provides diagnostics for each specimen
or just the specimen selected by the \code{which} argument.
The 2D projection is used in case of 3D landmark data.
The convex hull of the specimens (excluding the one being
selected) is compared with the actual specimen's landmarks.
This allows easy spotting of erroneous data.

The \code{plot_ord} and \code{plot_clust} are based on the
dissimilarities among specimens and provide ordination
(metric multidimensional scaling using \code{\link{cmdscale}}
based on square rooted dissimilarities and Cailliez's correction).
and hierarchical cluster dendrogram (using the \code{\link{hclust}}
function with Ward's clustering method).
}
\author{
Peter Solymos
}
\seealso{
\code{\link{plot.edma_data}} for visualizing EDMA data objects.

\code{\link{edma_fit}} for EDMA analysis.

\code{\link{dist}} for dissimilarity matrices and
\code{\link{T_test}} for description of the T-statistic.
}
\examples{
## read xyz files
file <- system.file(
    "extdata/crouzon/Crouzon_P0_Global_MUT.xyz",
    package="EDMAinR")
x <- read_xyz(file)
x

## the orignal structure
l <- readLines(file)
cat(l[1:10], sep="\n")
cat(l[(length(l)-10):length(l)], sep="\n")

## plots
plot(x) # same as plot_2d()
plot_2d(x, which=2) # show specimen #2
plot_ord(x)
plot_clust(x)

## dimensions and names
dim(x)
dimnames(x)
landmarks(x)
specimens(x)
dimensions(x)

## subsets
x[1:10, 2:3, 1:5]
subset(x, 1:10)

## coercion
str(as.matrix(x))
str(as.data.frame(x))
str(stack(x))
str(as.array(x))
as.edma_data(as.array(x))

## simulate data
K <- 3 # number of landmarks
D <- 2 # dimension, 2 or 3
sig <- 0.75
rho <- 0
SigmaK <- sig^2*diag(1, K, K) + sig^2*rho*(1-diag(1, K, K))
M <- matrix(c(0,1,0,0,0,1), 3, 2)
M[,1] <- M[,1] - mean(M[,1])
M[,2] <- M[,2] - mean(M[,2])
M <- 10*M

edma_simulate_data(10, M, SigmaK)
}
\keyword{ IO }
\keyword{ manip }
\keyword{ methods }
\keyword{ datagen }
