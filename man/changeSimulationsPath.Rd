\name{changeSimulationsPath}

\alias{changeSimulationsPath}

\title{
Reset file paths in Simulations objects
}
\description{
Reset file paths in the \code{values} of Simulations objects and also in \code{\link[raster]{raster}} objects.
}
\usage{
changeSimulationsPath(simulations, path)
}

\arguments{
  \item{simulations}{\code{raster} or \code{Simulations}
  }
  \item{path}{\code{character} with new path(s) for each of the layers
  }
}

\value{
The function returns \code{simulations} with replaced paths to the files.
If \code{simulations} is in memory, it has no effect.
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de}}


\examples{
# if data are taken from files, the paths need to be updated
\dontrun{
library(sensors4plumesData)
data(radioactivePlumes_area)
radioactivePlumes_area = changeSimulationsPath(radioactivePlumes_area,
c(paste0(path.package("sensors4plumesData"),
          "/extdata/radioactivePlumes_area_finaldose.grd"),
  paste0(path.package("sensors4plumesData"),
          "/extdata/radioactivePlumes_area_maxdose.grd"),
  paste0(path.package("sensors4plumesData"),
          "/extdata/radioactivePlumes_area_time.grd"))
)
}

}
