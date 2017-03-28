\name{spplotLog}
\alias{spplotLog}

\title{Plot methods for spatial data with logarithmic scale}

\description{
This function extends \code{\link[sp]{spplot}}; it enabels logarithmic scaling.
}

\usage{
spplotLog(x, zcol, base, nTicks = 10, replace0 = FALSE, minNonzero = 0, ...)
}

\arguments{
  \item{x}{\code{\link{SpatialDataFrame-class}} object to be plotted}
  \item{zcol}{character; name(s) of attribute(s) to be plotted}
  \item{base}{integer; base for the logarithmic transformation, to be used in the ticks; if none is indicated it is selected automatically to fit the range of the data}
  \item{nTicks}{integer; number of ticks - approximate -}
  \item{replace0}{logical; if zeros are to be replaced by a value that can be displayed on logscale. The tick of the key at this value will show "0". If \code{replace0 = FALSE} 0 is plotted as missing value.}
  \item{minNonzero}{numeric >= 0; all values below \code{minNonzero} are replaced by 0. Can be used to avoid that very small values dominate the plot.}
  \item{...}{parameters to be forwarded to \code{\link[sp]{spplot}} for style, colour etc.}
)
}
\details{
The key refers to the original values although in the background log transformed data are used (except from \code{SpatialPointsDataFrame} where just the \code{cuts} are changed).
}

\value{
A plot.
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SPixelsDF)
data(SPointsDF)
# add arbitrary data to show the effect of different values (< 0, == 0, small, big)
SPixelsDF@data$a = c(-5, 0, 0.003, 10, 57, 320, 444, 1000, 10000)
# replace 0 and small values
spplotLog(SPixelsDF, replace0 = TRUE, minNonzero = 0.05)
# choose base an number of ticks manually; forward parameters to spplot
spplotLog(SPixelsDF, base = 5, nTicks = 20,
          col.regions = heat.colors,
          sp.layout = list("sp.points", SPointsDF[1,]))
}          

