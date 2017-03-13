\name{cbind.Simulations}

\alias{cbind.Simulations}% cannot be used directly as cbind, only as cbind.Simulations; why? (parameters differ)

\title{Combine plumes of Simulations objects with coinciding parameters}

\description{
\code{cbind}-like method for \code{\link{Simulations}}: if \code{locations} are the same and the names of \code{plumes} and the layer names of the \code{values} coincide, the plumes of two or more \code{Simulations} can be combined by combining the values of the \code{plumes} slot and of the \code{values} slot.
}
\usage{
cbind.Simulations(..., savePath = ".", overwrite = FALSE)
}
\arguments{
  \item{...}{\code{Simulations}}
  \item{savePath}{path to save resulting raster files in case they are too big to keep in memory and then contain the final data, may be a directory or a directory followed a name like \code{"/myDir/thisPbind"} as files are created by appending names to it, starting with \code{"_"}
  }
  \item{overwrite}{if files in \code{savePath} may be overwritten}

}

\value{
A \code{Simulations} object; if it is too big to keep in memory, it is saved in a file in \code{savePath}.
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SimulationsSmall) 
SimulationsSmall2 = cbind.Simulations(SimulationsSmall, SimulationsSmall)
}

