\name{internals}
\alias{.edma_fit_np}
\alias{.Ttest_data}
\title{
Internal functions
}
\description{
Users should avoid directly using these functions. Not part of the public
user interface and can be changed any time without warning.
}
\usage{
.edma_fit_np(A, less=TRUE)
.Ttest_data(d1, d2, ref_denom=TRUE, mix=FALSE)
}
\arguments{
  \item{A}{
landmark data as an array (i.e. \code{as.array(x)}).
}
  \item{less}{
logical argument controlling what the function returns.
}
  \item{d1, d2}{
EDMA fit object to compare forms.
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
}
\keyword{ internal }
