\name{plot.Simulations}

\alias{plot.Simulations}
\alias{plot}

\title{
Plot overview of \code{simulations}.
}
\description{
Generate a common plot for the values associated with \code{locations} and with \code{plumes} and the \code{values}. 
}
\usage{
\method{plot}{Simulations}(x, ..., zcol = 1, main = "", 
                          col = terrain.colors(100), maxpixels = 1e+7)
}

\arguments{
  \item{x}{ \code{Simulations} object}
  \item{...}{ parameters to be forwarded }
  \item{zcol}{ integer, layer of \code{values} to be plotted }
  \item{main}{ character, header: should have short lines to fit into the upper right corner.}
  \item{col}{ color scheme }
  \item{maxpixels}{ if \code{values} (of one layer) fit into memory}
}
\details{
The plots show all factorial and numeric parameters, character parameters (e.g. names of \code{plumes}) are ignored. The colours cover min to max of each parameter.
If \code{values} (of one layer) fit into memory, they are plotted by \code{\link[graphics]{image}}. If they are bigger, we use \code{\link[raster]{plot}} and only plot a sample of size \code{maxpixels} (with warning).
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(radioactivePlumes)
plot(radioactivePlumes, zcol = 3, col = bpy.colors(), 
     main = "radioactive \n plumes \n area")
}

