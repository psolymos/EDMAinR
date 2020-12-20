\name{gpa_fit}
\alias{gpa_fit}
\alias{print.gpa_fit}
\title{
GPA
}
\description{
Fit GPA or WGPA to landmark data.
}
\usage{
gpa_fit(x, B = 0, ncores = getOption("Ncpus", 1L),
    weighted=FALSE, ...)
\method{print}{gpa_fit}(x, truncate=40, ...)

}
\arguments{
  \item{x}{
an EDMA data object of class 'edma_data'.
}
  \item{B}{
nonnegative integer, the number of bootstrap replicates.
}
  \item{weighted}{
logical, use \code{shapes::procWGPA} instead of \code{shapes::procGPA}.
}
  \item{ncores}{
positive integer, the number of cores to use when bootstrapping.
Use \code{options(Ncpus = 2)} to set it to 2 globally.
}
  \item{truncate}{
numeric, number of characters to print for the object title.
}
  \item{\dots}{
arguments passed to \code{shapes::procGPA} or \code{shapes::procWGPA}.
}
}
\value{
Returns only form matrix, SigmaKstar is \code{NA}.
}
\references{
Gower, J.C. (1975). Generalized Procrustes analysis, Psychometrika, 40, 33--50.
}
\author{
Peter Solymos wrote the wrapper for \code{shapes::procGPA}.
}
\examples{
file <- system.file(
    "extdata/crouzon/Crouzon_P0_Global_MUT.xyz",
    package="EDMAinR")
x <- read_xyz(file)
x <- x[,,1:10] # 10 specimens

## nonparametric fit
fit <- gpa_fit(x, B=9)
fit
str(Meanform(fit))
str(SigmaKstar(fit))
}
\keyword{ models }
