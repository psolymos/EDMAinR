\name{edma_colors}
\alias{edma_colors}
\alias{plot_edma_colors}
\title{
Check and manipulate colors
}
\description{
Check and manipulate the default color values.
}
\usage{
edma_colors(n,
    type=c("diverging", "sequential", "qualitative"),
    alpha=1, rev=FALSE)

plot_edma_colors(n=9, maxq=9)
}
\arguments{
  \item{n}{
the number of colors (>0) to be in the palette.
}
  \item{type}{
the type of palette.
}
  \item{alpha}{
the alpha transparency, a number in [0,1].
}
  \item{rev}{
logical, should colors be reversed.
}
  \item{maxq}{
maximum number of qualitative colors to plot.
}
}
\details{
\code{edma_colors} create a vector of \code{n} colors
based on the settings in \code{getOption("edma_options")}.

The options can be set via \code{\link{options}} (see Examples).
The options can either be the name of a palette \code{\link{hcl.colors}}.
When the option is set to multiple values, those are treated
as colors to be interplolated with \code{\link{colorRampPalette}}.
For qualitative palettes, the color values are used directly
(recycled as needed).

Sequential palettes are produced as the higher half of the
diverging palette for consistency.

}
\value{
\code{edma_colors} returns a vector of hax color codes,
\code{plot_edma_colors} produces a plot with the
diverging, sequential, and qualitative default palettes
given settings in \code{getOption("edma_options")}.
}
\author{
Peter Solymos
}
\seealso{
\code{\link{hcl.colors}}, \code{\link{colorRampPalette}}, \code{\link{col2rgb}}
}
\examples{
## default palettes
plot_edma_colors(101)

## change default palettes
op <- options("edma_options" = list(
    diverging = "Green-Orange",
    qualitative = "Dark 2"))
plot_edma_colors(101)

## use color names
options("edma_options" = list(
    diverging = c("black", "grey", "pink"),
    qualitative = "Warm"))
plot_edma_colors(101)

## reset defaults
options(op)
plot_edma_colors(101)
}
\keyword{ hplot }
