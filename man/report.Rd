\name{report}
\alias{edma_fdm_report}
\alias{edma_gdm_report}
\title{
EDMA Report
}
\description{
Writes output into a text file following close the WinEDMA implementation.
}
\usage{
edma_fdm_report(numerator, denominator, output="edma_output.txt",
  landmarks=NULL, B=0, level=0.95, ref_denom=TRUE, mix=FALSE,
  digits=4)

edma_gdm_report(numerator_yng, numerator_old,
    denominator_yng, denominator_old, output="edma_output.txt",
    landmarks=NULL, B=0, level=0.95, ref_denom=TRUE, mix=FALSE,
    digits=4)
}
\arguments{
  \item{numerator, denominator, numerator_yng, numerator_old, denominator_yng, denominator_old}{
input file names or EDMA data objects.
}
  \item{output}{
path and file name for the output file.
}
  \item{landmarks}{
a subset of landmarks to be specified, or \code{NULL} (use all landmarks).
}
  \item{B}{
number of bootstrap samples.
}
  \item{level}{
confidence level.
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
  \item{digits}{
significant digits to print.
}
}
\examples{
\dontrun{

edma_fdm_report(
    numerator = system.file(
        "extdata/crouzon/Crouzon_P0_Global_NON-MUT.xyz",
        package="EDMAinR"),
    denominator = system.file(
        "extdata/crouzon/Crouzon_P0_Global_MUT.xyz",
        package="EDMAinR"),
    output="edma_output.txt",
    landmarks=c("locc", "lpfl", "lpsq", "lpto", "lsqu",
        "rocc", "rpfl", "rpsq", "rpto", "rsqu"),
    B=1000,
    level=0.9,
    ref_denom=FALSE,
    mix=TRUE)

}
}
\keyword{ IO }
