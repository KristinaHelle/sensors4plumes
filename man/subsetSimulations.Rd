\name{subset.Simulations}

\alias{subset.Simulations}
\alias{subset}

\title{
Subsetting Simulations
}
\description{
Return subsets of \code{\link{Simulations}}, subsetting by \code{plumes}, \code{locations}, the kinds of \code{values} or by the columns of data associated with \code{plumes} or \code{locations} -- possibly combined.
}
\usage{
\method{subset}{Simulations}(x, ..., locations, plumes, kinds, dataLocations, dataPlumes, 
  saveName, saveDir = ".", overwrite = FALSE, valuesOnly = FALSE)
}

\arguments{
  \item{x}{\code{Simulations}}
  \item{...}{parameters to be forwarded}
  \item{locations}{indices of \code{locations} to keep}
  \item{plumes}{indices of \code{plumes} to keep}
  \item{kinds}{indices of kinds of \code{values} to keep}
  \item{dataLocations}{indices of the data associated with the \code{locations} to keep}
  \item{dataPlumes}{indices of the data associated with the \code{plumes} to keep}  
  \item{saveName}{name for saving the new \code{raster} files for the \code{values} if they cannot be kept in memory (without suffix)}
  \item{saveDir}{directory for saving if necessary}
  \item{overwrite}{\code{logical} if files may be overwritten by the raster files}
  \item{valuesOnly}{\code{logical} if only the subset of \code{values} is returned (if \code{FALSE} the result is the full new \code{Simulations})}
}

\details{
Multiple or invalid values in \code{locations} and \code{plumes} are ignored, only the order is taken into account. (Because in \code{SpatialIndexDataFrame} and \code{SpatialPolygridDataFrame} no locations can be selected multiply.)

However, if \code{valuesOnly = TRUE}, multiplity is taken into account and invalid indices result in \code{NA} values. In this case only the \code{values} are returned.
}

\value{
A \code{Simulations} object, subsetted by the given indices: subsetting \code{locations} or \code{plumes} results in subsetting of the \code{values}. If the resulting raster is too big to keep in memory, it is saved.

If \code{valuesOnly = TRUE}, only the new \code{values} are returned, i.e. a \code{\link[raster]{raster}} object.
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SimulationsSmall)

SimulationsSmall_a = subset(SimulationsSmall, locations = c(2,4,6,8))
SimulationsSmall_b = subset(SimulationsSmall, 
                            locations = 1:5, plumes = c(2,4), kinds = 1, 
                            dataPlumes = 1:2)
                            
# effect of 'valuesOnly'
data(radioactivePlumes_area)
locations = sample(nLocations(radioactivePlumes_area), 15)
plumes = sample(nPlumes(radioactivePlumes_area), 15)
subset_fm_RNA = subset(radioactivePlumes_area, 
                       plumes = c(plumes2, NA, 10000, NA, 0, plumes2), 
                       locations = c(locations2, NA, 10000, NA, 0,locations2))
                        
subset_mf_RNA_vO = subset(radioactivePlumes_area, 
                          plumes = c(plumes2, NA, 0, 10000, NA, plumes2), 
                          locations = c(locations2, NA, NA, 10000, 0, 0, locations2),
                          valuesOnly = TRUE)  
image(subset_fm_RNA@values, 1, col = rainbow(20))
image(subset_mf_RNA_vO, 1, col = rainbow(20))
}
