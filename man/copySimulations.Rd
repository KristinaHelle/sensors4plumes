\name{copySimulations}

\alias{copySimulations}

\title{
Copy Simulations (including raster files)
}
\description{
Copy or move the raster files associated with a \code{\link{Simulations}} object to another directory and update the paths in the \code{Simulations}. The resulting \code{Simulations} may be saved.
}
\usage{
copySimulations(simulations, newPath, newFile, 
  overwrite = FALSE, deleteOld = FALSE)
}

\arguments{
  \item{simulations}{ \code{Simulations}
  }
  \item{newPath}{path to the directory where to save the files; if it is missing, the raster files are not copied
  }
  \item{newFile}{filename where to save the resulting \code{simulations}; if missing, the \code{simulations} are not saved and just returned; the file is saved in the directory \code{newPath} if given and else in the current directory    
  }
  \item{overwrite}{\code{logical} if files may be overwritten when copying the raster files -- has no influence on saving of the \code{simulations} themselves as .Rdata file}
  \item{deleteOld}{\code{logical} if the old \code{raster} files are to be deleted (move files to new directory instead of copy)}
}

\value{
The \code{simulations} with paths in the \code{values} updated to the new \code{raster} files. In addition the \code{raster} files are copied or moved to \code{newPath} and the \code{simulations} may be saved in a .Rdata file. If the \code{values} of the \code{simulations} are in memory, there are no raster files to be copied. In this case the function will only save the \code{simulations} if a \code{newFile} is given.
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de}}

\examples{
\dontrun{
library(sensors4plumesData)
data(radioactivePlumes_area)
radioactivePlumes_area2 = 
  copySimulations(radioactivePlumes_area, 
  newPath = paste0(path.package("sensors4plumesData"), 
  newFile = "radioactivePlumes_area2")
}
}